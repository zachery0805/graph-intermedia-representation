#include <stack>
#include <queue>
#include <algorithm>
#include <climits>
#include "gir.h"
#include "scc_common.h"
#include "graph.h"
#include "graph_undirected.h"
#include "wcc_cpu/wcc_common.h"
#include "block_cut_tree.h"
#include "dynamic_repair.h"
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <chrono>
#include <omp.h>

// -------------------- 构造/析构 --------------------

GIR::GIR(const std::string& path_directed,
         const std::string& path_undirected)
    : base_path_directed(path_directed),
      base_path_undirected(path_undirected) {
    std::string fw_bg_file = base_path_directed + "fw_begin.bin";
    std::string fw_aj_file = base_path_directed + "fw_adjacent.bin";
    std::string bw_bg_file = base_path_directed + "bw_begin.bin";
    std::string bw_aj_file = base_path_directed + "bw_adjacent.bin";
    g_directed = new graph(fw_bg_file.c_str(), fw_aj_file.c_str(),
                           bw_bg_file.c_str(), bw_aj_file.c_str());
}

GIR::~GIR() {
    if (g_directed) delete g_directed;
    if (g_undirected) delete g_undirected;
    for (auto* t : scc_bct) if (t) delete t;
}

// ==================== Dynamic edge overlay ====================
// See gir.h for the overall design note. Presence of a logical edge is
// always resolved through isDirectedEdgePresent/isUndirectedEdgePresent;
// every traversal helper below is built on top of those same rules so there
// is exactly one notion of "does this edge exist right now".

long long GIR::encodeDirected(int u, int v) {
    return (static_cast<long long>(u) << 32) | static_cast<unsigned int>(v);
}

long long GIR::encodeUndirected(int u, int v) {
    if (u > v) std::swap(u, v);
    return (static_cast<long long>(u) << 32) | static_cast<unsigned int>(v);
}

bool GIR::isDirectedEdgePresent(int u, int v) const {
    long long key = encodeDirected(u, v);
    if (removed_directed.count(key)) return false;
    auto it = extra_out_adj.find(u);
    if (it != extra_out_adj.end())
        for (int w : it->second) if (w == v) return true;
    if (g_directed && u >= 0 && (index_t)u < g_directed->vert_count) {
        for (index_t i = g_directed->fw_beg_pos[u]; i < g_directed->fw_beg_pos[u + 1]; ++i)
            if ((int)g_directed->fw_csr[i] == v) return true;
    }
    return false;
}

bool GIR::isUndirectedEdgePresent(int u, int v) const {
    long long key = encodeUndirected(u, v);
    if (removed_undirected.count(key)) return false;
    auto it = extra_undirected_adj.find(u);
    if (it != extra_undirected_adj.end())
        for (int w : it->second) if (w == v) return true;
    if (u >= 0 && (size_t)u < base_undirected_adj.size()) {
        const auto& adj = base_undirected_adj[u];
        if (std::binary_search(adj.begin(), adj.end(), v)) return true;
    }
    return false;
}

std::vector<int> GIR::directedOutNeighbors(int v) const {
    std::vector<int> result;
    if (g_directed && v >= 0 && (index_t)v < g_directed->vert_count) {
        for (index_t i = g_directed->fw_beg_pos[v]; i < g_directed->fw_beg_pos[v + 1]; ++i) {
            int w = (int)g_directed->fw_csr[i];
            if (!removed_directed.count(encodeDirected(v, w))) result.push_back(w);
        }
    }
    auto it = extra_out_adj.find(v);
    if (it != extra_out_adj.end()) result.insert(result.end(), it->second.begin(), it->second.end());
    return result;
}

std::vector<int> GIR::directedInNeighbors(int v) const {
    std::vector<int> result;
    if (g_directed && v >= 0 && (index_t)v < g_directed->vert_count) {
        for (index_t i = g_directed->bw_beg_pos[v]; i < g_directed->bw_beg_pos[v + 1]; ++i) {
            int w = (int)g_directed->bw_csr[i];
            if (!removed_directed.count(encodeDirected(w, v))) result.push_back(w);
        }
    }
    auto it = extra_in_adj.find(v);
    if (it != extra_in_adj.end()) result.insert(result.end(), it->second.begin(), it->second.end());
    return result;
}

std::vector<int> GIR::undirectedNeighbors(int v) const {
    std::vector<int> result;
    if (v >= 0 && (size_t)v < base_undirected_adj.size()) {
        for (int w : base_undirected_adj[v])
            if (!removed_undirected.count(encodeUndirected(v, w))) result.push_back(w);
    }
    auto it = extra_undirected_adj.find(v);
    if (it != extra_undirected_adj.end()) result.insert(result.end(), it->second.begin(), it->second.end());
    return result;
}

void GIR::addOverlayDirectedEdge(int u, int v) {
    long long key = encodeDirected(u, v);
    if (removed_directed.count(key)) { removed_directed.erase(key); return; }
    extra_out_adj[u].push_back(v);
    extra_in_adj[v].push_back(u);
}

void GIR::removeOverlayDirectedEdge(int u, int v) {
    long long key = encodeDirected(u, v);
    auto it = extra_out_adj.find(u);
    if (it != extra_out_adj.end()) {
        auto& vec = it->second;
        auto pos = std::find(vec.begin(), vec.end(), v);
        if (pos != vec.end()) {
            vec.erase(pos);
            auto it2 = extra_in_adj.find(v);
            if (it2 != extra_in_adj.end()) {
                auto& vec2 = it2->second;
                auto pos2 = std::find(vec2.begin(), vec2.end(), u);
                if (pos2 != vec2.end()) vec2.erase(pos2);
            }
            return;
        }
    }
    removed_directed.insert(key);
}

void GIR::addOverlayUndirectedEdge(int u, int v) {
    long long key = encodeUndirected(u, v);
    if (removed_undirected.count(key)) { removed_undirected.erase(key); return; }
    extra_undirected_adj[u].push_back(v);
    extra_undirected_adj[v].push_back(u);
}

void GIR::removeOverlayUndirectedEdge(int u, int v) {
    long long key = encodeUndirected(u, v);
    bool was_extra = false;
    auto it = extra_undirected_adj.find(u);
    if (it != extra_undirected_adj.end()) {
        auto& vec = it->second;
        auto pos = std::find(vec.begin(), vec.end(), v);
        if (pos != vec.end()) { vec.erase(pos); was_extra = true; }
    }
    auto it2 = extra_undirected_adj.find(v);
    if (it2 != extra_undirected_adj.end()) {
        auto& vec2 = it2->second;
        auto pos2 = std::find(vec2.begin(), vec2.end(), u);
        if (pos2 != vec2.end()) vec2.erase(pos2);
    }
    if (!was_extra) removed_undirected.insert(key);
}

// ==================== WCC layer: union-find ====================

int GIR::findWcc(int v) {
    while (wcc_parent[v] != v) {
        wcc_parent[v] = wcc_parent[wcc_parent[v]];
        v = wcc_parent[v];
    }
    return v;
}

void GIR::unionWcc(int a, int b) {
    a = findWcc(a);
    b = findWcc(b);
    if (a == b) return;
    if (wcc_size[a] < wcc_size[b]) std::swap(a, b);
    wcc_parent[b] = a;
    wcc_size[a] += wcc_size[b];
}

void GIR::rebuildWccComponent(const std::vector<int>& members) {
    if (members.empty()) return;
    int root = members[0];
    wcc_parent[root] = root;
    wcc_size[root] = (int)members.size();
    for (size_t i = 1; i < members.size(); ++i) wcc_parent[members[i]] = root;
}

bool GIR::probeUndirectedConnectivity(int u, int v, std::vector<int>& smaller_side_out) {
    if (u == v) { smaller_side_out = {u}; return true; }
    std::unordered_set<int> visited_u{u}, visited_v{v};
    std::vector<int> frontier_u{u}, frontier_v{v};
    while (!frontier_u.empty() && !frontier_v.empty()) {
        std::vector<int> next_u;
        for (int c : frontier_u) {
            for (int nb : undirectedNeighbors(c)) {
                if (visited_v.count(nb)) return true;
                if (visited_u.insert(nb).second) next_u.push_back(nb);
            }
        }
        frontier_u.swap(next_u);
        if (frontier_u.empty()) break;

        std::vector<int> next_v;
        for (int c : frontier_v) {
            for (int nb : undirectedNeighbors(c)) {
                if (visited_u.count(nb)) return true;
                if (visited_v.insert(nb).second) next_v.push_back(nb);
            }
        }
        frontier_v.swap(next_v);
        if (frontier_v.empty()) break;
    }
    if (frontier_u.empty()) smaller_side_out.assign(visited_u.begin(), visited_u.end());
    else smaller_side_out.assign(visited_v.begin(), visited_v.end());
    return false;
}

// -------------------- WCC 构建 --------------------
void GIR::buildWCC() {
    std::cout << "=== Build WCC ===" << std::endl;

    // NOTE: this used to call wcc_detection() (includes/wcc_cpu/wcc_core.cpp,
    // the parallel bfs_concurrent + color_propagation_fq pipeline). That
    // function was found to silently miscount WCCs even on the small
    // 18-vertex example_graph (reports 1 WCC when the graph actually has 2
    // disjoint components) while investigating the differential test added
    // alongside GIR's dynamic-update rewrite -- see docs/code-review.md.
    // Root-causing the parallel algorithm itself is a separate, larger
    // effort, so construction here uses a simple, unambiguously-correct
    // sequential union-find pass instead. It's asymptotically fine (union-
    // find is near-O(1) amortized per edge) even if not multi-threaded.
    vertex_t V = g_directed->vert_count;
    wcc_parent.assign(V, 0);
    wcc_size.assign(V, 1);
    for (vertex_t i = 0; i < V; ++i) wcc_parent[i] = (int)i;

    for (vertex_t u = 0; u < V; ++u) {
        for (index_t i = g_directed->fw_beg_pos[u]; i < g_directed->fw_beg_pos[u + 1]; ++i) {
            unionWcc((int)u, (int)g_directed->fw_csr[i]);
        }
    }

    std::unordered_set<int> roots;
    for (vertex_t i = 0; i < V; ++i) roots.insert(findWcc((int)i));
    std::cout << "Found " << roots.size() << " WCCs." << std::endl;
}

// ==================== SCC / condensation DAG layer ====================

void GIR::addDagEdge(int from_scc, int to_scc) {
    if (from_scc == to_scc) return;
    dag_adj_set[from_scc][to_scc]++;
    dag_radj_set[to_scc][from_scc]++;
}

// Decrements the reference count for this condensation edge, removing it
// only once no underlying vertex-level edge justifies it any more. Two
// SCCs can be connected by more than one real edge; deleting just one must
// not make the DAG edge disappear while another still exists.
void GIR::removeDagEdge(int from_scc, int to_scc) {
    auto it1 = dag_adj_set[from_scc].find(to_scc);
    if (it1 != dag_adj_set[from_scc].end() && --(it1->second) <= 0)
        dag_adj_set[from_scc].erase(it1);
    auto it2 = dag_radj_set[to_scc].find(from_scc);
    if (it2 != dag_radj_set[to_scc].end() && --(it2->second) <= 0)
        dag_radj_set[to_scc].erase(it2);
}

int GIR::newSccId() {
    for (size_t i = 0; i < scc_alive.size(); ++i) {
        if (!scc_alive[i]) {
            scc_alive[i] = 1;
            sccs[i].clear();
            dag_adj_set[i].clear();
            dag_radj_set[i].clear();
            return (int)i;
        }
    }
    int id = (int)sccs.size();
    sccs.push_back({});
    scc_alive.push_back(1);
    dag_adj_set.push_back({});
    dag_radj_set.push_back({});
    scc_rank.push_back(0);
    scc_bicc_blocks.push_back({});
    scc_ap.push_back({});
    scc_bct.push_back(nullptr);
    return id;
}

void GIR::renumberAllRanks() {
    std::vector<int> order;
    order.reserve(rank_to_scc.size());
    for (auto& kv : rank_to_scc) order.push_back(kv.second);
    rank_to_scc.clear();
    for (size_t i = 0; i < order.size(); ++i) {
        long long r = (long long)(i + 1) * RANK_GAP;
        scc_rank[order[i]] = r;
        rank_to_scc[r] = order[i];
    }
}

// Recomputes a valid topological order (and reassigns ranks) for exactly
// the SCCs whose rank currently falls strictly between lo_bound and
// hi_bound. Used after Case-2 merges/no-cycle backward inserts and Case-5b
// splits, where the touched region is exactly this interval -- SCCs outside
// it never need to move, since their existing rank already respects every
// real edge that isn't itself inside the interval.
void GIR::repairRankInterval(long long lo_bound, long long hi_bound) {
    for (int attempt = 0; attempt < 2; ++attempt) {
        auto lo_it = rank_to_scc.upper_bound(lo_bound);
        auto hi_it = rank_to_scc.lower_bound(hi_bound);
        std::vector<int> nodes;
        for (auto it = lo_it; it != hi_it; ++it) nodes.push_back(it->second);
        if (nodes.size() <= 1) return;

        long long span = hi_bound - lo_bound;
        if (span < (long long)(nodes.size() + 1)) {
            if (attempt == 1) return;  // already retried once; give up quietly
            renumberAllRanks();
            long long new_lo = scc_rank[nodes[0]], new_hi = scc_rank[nodes[0]];
            for (int n : nodes) {
                new_lo = std::min(new_lo, scc_rank[n]);
                new_hi = std::max(new_hi, scc_rank[n]);
            }
            lo_bound = new_lo - RANK_GAP / 2;
            hi_bound = new_hi + RANK_GAP / 2;
            continue;
        }

        std::unordered_set<int> node_set(nodes.begin(), nodes.end());
        std::unordered_map<int, int> indeg;
        for (int n : nodes) indeg[n] = 0;
        for (int n : nodes)
            for (auto& kv : dag_adj_set[n])
                if (node_set.count(kv.first)) indeg[kv.first]++;

        std::vector<int> queue_nodes;
        for (int n : nodes) if (indeg[n] == 0) queue_nodes.push_back(n);
        std::vector<int> order;
        order.reserve(nodes.size());
        size_t qi = 0;
        while (qi < queue_nodes.size()) {
            int n = queue_nodes[qi++];
            order.push_back(n);
            for (auto& kv : dag_adj_set[n])
                if (node_set.count(kv.first))
                    if (--indeg[kv.first] == 0) queue_nodes.push_back(kv.first);
        }
        if (order.size() != nodes.size()) return;  // touched region should be acyclic; guard only

        for (int n : nodes) rank_to_scc.erase(scc_rank[n]);
        size_t k = order.size();
        for (size_t i = 0; i < k; ++i) {
            long long r = lo_bound + (long long)((i + 1) * (span / (long long)(k + 1)));
            scc_rank[order[i]] = r;
            rank_to_scc[r] = order[i];
        }
        return;
    }
}

// Shifts an entire condensation-connected component (one former WCC's SCCs,
// used right after a Case-1 cross-WCC merge) below `must_be_below`,
// preserving its internal relative order. Bounded by that component's own
// SCC count, not the other side's.
void GIR::shiftComponentBelow(int start_scc, long long must_be_below) {
    std::vector<int> comp;
    std::unordered_set<int> visited;
    std::vector<int> stack = {start_scc};
    visited.insert(start_scc);
    long long max_rank = scc_rank[start_scc];
    while (!stack.empty()) {
        int c = stack.back();
        stack.pop_back();
        comp.push_back(c);
        max_rank = std::max(max_rank, scc_rank[c]);
        for (auto& kv : dag_adj_set[c]) if (!visited.count(kv.first)) { visited.insert(kv.first); stack.push_back(kv.first); }
        for (auto& kv : dag_radj_set[c]) if (!visited.count(kv.first)) { visited.insert(kv.first); stack.push_back(kv.first); }
    }
    if (max_rank < must_be_below) return;
    long long offset = must_be_below - max_rank - RANK_GAP;
    for (int c : comp) {
        rank_to_scc.erase(scc_rank[c]);
        scc_rank[c] += offset;
        rank_to_scc[scc_rank[c]] = c;
    }
}

// Mirror of shiftComponentBelow: shifts an entire condensation-connected
// component above `must_be_above` instead of below. Used when the *other*
// side of a Case-1 cross-WCC merge is the smaller one, so that side gets
// shifted instead (see handleInsert) -- avoids always shifting whichever
// component happens to contain the edge's source vertex, which is
// catastrophic on real graphs with one dominant giant WCC (shifting it on
// every insertion that touches it is effectively O(whole graph) instead of
// O(the small side)).
void GIR::shiftComponentAbove(int start_scc, long long must_be_above) {
    std::vector<int> comp;
    std::unordered_set<int> visited;
    std::vector<int> stack = {start_scc};
    visited.insert(start_scc);
    long long min_rank = scc_rank[start_scc];
    while (!stack.empty()) {
        int c = stack.back();
        stack.pop_back();
        comp.push_back(c);
        min_rank = std::min(min_rank, scc_rank[c]);
        for (auto& kv : dag_adj_set[c]) if (!visited.count(kv.first)) { visited.insert(kv.first); stack.push_back(kv.first); }
        for (auto& kv : dag_radj_set[c]) if (!visited.count(kv.first)) { visited.insert(kv.first); stack.push_back(kv.first); }
    }
    if (min_rank > must_be_above) return;
    long long offset = must_be_above - min_rank + RANK_GAP;
    for (int c : comp) {
        rank_to_scc.erase(scc_rank[c]);
        scc_rank[c] += offset;
        rank_to_scc[scc_rank[c]] = c;
    }
}

bool GIR::probeCondensationConnectivity(int su, int sv) {
    if (su == sv) return true;
    std::unordered_set<int> visited{su};
    std::queue<int> q;
    q.push(su);
    while (!q.empty()) {
        int c = q.front();
        q.pop();
        for (auto& kv : dag_adj_set[c]) {
            if (kv.first == sv) return true;
            if (visited.insert(kv.first).second) q.push(kv.first);
        }
        for (auto& kv : dag_radj_set[c]) {
            if (kv.first == sv) return true;
            if (visited.insert(kv.first).second) q.push(kv.first);
        }
    }
    return false;
}

std::vector<int> GIR::sccVertices(int scc_id) const {
    if (scc_id < 0 || (size_t)scc_id >= sccs.size()) return {};
    return sccs[scc_id];
}

void GIR::mergeSccs(const std::vector<int>& ids) {
    if (ids.size() < 2) return;
    int keep = ids[0];
    std::unordered_set<int> merge_set(ids.begin(), ids.end());

    for (size_t i = 1; i < ids.size(); ++i) {
        int other = ids[i];
        for (int v : sccs[other]) {
            nodeToScc[v] = keep;
            sccs[keep].push_back(v);
        }
        sccs[other].clear();
        scc_alive[other] = 0;
        rank_to_scc.erase(scc_rank[other]);
        if (scc_bct[other]) { delete scc_bct[other]; scc_bct[other] = nullptr; }
        scc_bicc_blocks[other].clear();
        scc_ap[other].clear();
    }

    // Accumulate (summed, not just union) external edge counts across all
    // merged SCCs -- two of them may have separately had edges to/from the
    // same external SCC, and the merged total must reflect both.
    std::unordered_map<int, int> new_out, new_in;
    for (int id : ids) {
        for (auto& kv : dag_adj_set[id]) if (!merge_set.count(kv.first)) new_out[kv.first] += kv.second;
        for (auto& kv : dag_radj_set[id]) if (!merge_set.count(kv.first)) new_in[kv.first] += kv.second;
        dag_adj_set[id].clear();
        dag_radj_set[id].clear();
    }
    dag_adj_set[keep] = new_out;
    dag_radj_set[keep] = new_in;
    for (auto& kv : new_out) {
        int s = kv.first;
        for (int id : ids) if (id != keep) dag_radj_set[s].erase(id);
        dag_radj_set[s][keep] = kv.second;
    }
    for (auto& kv : new_in) {
        int p = kv.first;
        for (int id : ids) if (id != keep) dag_adj_set[p].erase(id);
        dag_adj_set[p][keep] = kv.second;
    }

    long long min_rank = scc_rank[keep];
    for (int id : ids) min_rank = std::min(min_rank, scc_rank[id]);
    rank_to_scc.erase(scc_rank[keep]);
    scc_rank[keep] = min_rank;
    rank_to_scc[min_rank] = keep;

    rebuildSccBicc(keep);
}

void GIR::splitScc(int scc_id) {
    std::vector<int> old_verts = sccs[scc_id];
    std::unordered_set<int> old_vert_set(old_verts.begin(), old_verts.end());

    auto out_fn = [this, &old_vert_set](int x) {
        std::vector<int> nb = directedOutNeighbors(x);
        std::vector<int> filtered;
        for (int w : nb) if (old_vert_set.count(w)) filtered.push_back(w);
        return filtered;
    };
    auto fragments = dynamic_repair::find_scc_scoped(old_verts, out_fn);
    if (fragments.size() <= 1) {
        rebuildSccBicc(scc_id);
        return;
    }

    long long old_rank = scc_rank[scc_id];
    long long lo_bound = old_rank - RANK_GAP, hi_bound = old_rank + RANK_GAP;
    bool has_pred = false, has_succ = false;
    for (auto& kv : dag_radj_set[scc_id]) { if (!has_pred || scc_rank[kv.first] > lo_bound) lo_bound = scc_rank[kv.first]; has_pred = true; }
    for (auto& kv : dag_adj_set[scc_id]) { if (!has_succ || scc_rank[kv.first] < hi_bound) hi_bound = scc_rank[kv.first]; has_succ = true; }
    if (!has_pred) lo_bound = old_rank - RANK_GAP;
    if (!has_succ) hi_bound = old_rank + RANK_GAP;

    // scc_id is being fully retired (its fragments get freshly rewired via
    // addDagEdge below), so drop its entries from neighbors entirely
    // regardless of remaining count.
    for (auto& kv : dag_adj_set[scc_id]) dag_radj_set[kv.first].erase(scc_id);
    for (auto& kv : dag_radj_set[scc_id]) dag_adj_set[kv.first].erase(scc_id);
    dag_adj_set[scc_id].clear();
    dag_radj_set[scc_id].clear();
    rank_to_scc.erase(scc_rank[scc_id]);
    scc_alive[scc_id] = 0;
    sccs[scc_id].clear();
    if (scc_bct[scc_id]) { delete scc_bct[scc_id]; scc_bct[scc_id] = nullptr; }
    scc_bicc_blocks[scc_id].clear();
    scc_ap[scc_id].clear();

    std::vector<int> new_ids;
    new_ids.reserve(fragments.size());
    for (auto& frag : fragments) {
        int id = newSccId();
        sccs[id] = frag;
        for (int v : frag) nodeToScc[v] = id;
        new_ids.push_back(id);
    }

    {
        long long span = hi_bound - lo_bound;
        if (span < 1) span = RANK_GAP;
        size_t k = new_ids.size();
        for (size_t i = 0; i < k; ++i) {
            long long r = lo_bound + (long long)((i + 1) * (span / (long long)(k + 1)));
            scc_rank[new_ids[i]] = r;
            rank_to_scc[r] = new_ids[i];
        }
    }

    // Rewire cross-fragment edges. A real edge between two old_verts that
    // now land in *different* new fragments must get its own condensation
    // edge too, not just edges to genuinely external SCCs -- e.g. if the
    // split was triggered by deleting 9->2 but 2->9 still exists, that
    // surviving edge is now inter-fragment and must be registered. The
    // out-neighbor scan covers both cases (internal-to-different-fragment
    // and external) in one pass; the in-neighbor scan is restricted to
    // truly-external predecessors only, so an edge with both endpoints in
    // old_verts isn't double-counted (it's already captured once via its
    // source vertex's out-neighbor scan).
    for (int x : old_verts) {
        int sx = nodeToScc[x];
        for (int y : directedOutNeighbors(x)) {
            auto it = nodeToScc.find(y);
            if (it == nodeToScc.end()) continue;
            if (it->second != sx) addDagEdge(sx, it->second);
        }
        for (int y : directedInNeighbors(x)) {
            if (old_vert_set.count(y)) continue;
            auto it = nodeToScc.find(y);
            if (it != nodeToScc.end()) addDagEdge(it->second, sx);
        }
    }

    for (int id : new_ids) rebuildSccBicc(id);

    repairRankInterval(lo_bound, hi_bound);
}

// -------------------- BCT: per-SCC BiCC/AP + BC tree --------------------

void GIR::rebuildSccBicc(int scc_id) {
    const std::vector<int>& verts = sccs[scc_id];
    std::unordered_set<int> vert_set(verts.begin(), verts.end());
    auto neighbor_fn = [this, &vert_set](int x) {
        std::vector<int> nb = undirectedNeighbors(x);
        std::vector<int> filtered;
        for (int w : nb) if (vert_set.count(w)) filtered.push_back(w);
        return filtered;
    };
    auto result = dynamic_repair::find_bicc_scoped(verts, neighbor_fn);

    std::unordered_set<int> covered;
    for (auto& blk : result.blocks) for (int v : blk) covered.insert(v);
    for (int v : verts) if (!covered.count(v)) result.blocks.push_back({v});

    scc_bicc_blocks[scc_id] = std::move(result.blocks);
    scc_ap[scc_id] = std::move(result.articulation_points);

    if (scc_bct[scc_id]) delete scc_bct[scc_id];
    scc_bct[scc_id] = new BlockCutTree();
    scc_bct[scc_id]->build(scc_bicc_blocks[scc_id], scc_ap[scc_id]);
}

bool GIR::inSameBiccBlock(int scc_id, int u, int v) const {
    if (scc_id < 0 || (size_t)scc_id >= scc_bicc_blocks.size()) return false;
    for (auto& blk : scc_bicc_blocks[scc_id]) {
        bool hu = false, hv = false;
        for (int x : blk) { if (x == u) hu = true; if (x == v) hv = true; }
        if (hu && hv) return true;
    }
    return false;
}

bool GIR::isBridgeWithinScc(int scc_id, int u, int v) const {
    if (scc_id < 0 || (size_t)scc_id >= scc_bicc_blocks.size()) return true;
    for (auto& blk : scc_bicc_blocks[scc_id]) {
        bool hu = false, hv = false;
        for (int x : blk) { if (x == u) hu = true; if (x == v) hv = true; }
        if (hu && hv) return blk.size() == 2;
    }
    return true;
}

void GIR::buildBCT() {
    std::cout << "=== Build BCT (per SCC, per paper/task_parallel.tex model) ===" << std::endl;

    std::string undirected_bg_file = base_path_undirected + "undirected_begin.bin";
    std::string undirected_aj_file = base_path_undirected + "undirected_adjacent.bin";
    g_undirected = new graph_undirected(undirected_bg_file.c_str(), undirected_aj_file.c_str());

    // Symmetrize the base CSR once here (see gir.h's comment on
    // base_undirected_adj): the on-disk undirected_*.bin files are not
    // guaranteed to list both directions of every edge.
    base_undirected_adj.assign(g_undirected->vert_count, {});
    for (index_t u = 0; u < g_undirected->vert_count; ++u) {
        for (index_t i = g_undirected->beg_pos[u]; i < g_undirected->beg_pos[u + 1]; ++i) {
            int w = (int)g_undirected->csr[i];
            base_undirected_adj[u].push_back(w);
            base_undirected_adj[w].push_back((int)u);
        }
    }
    for (auto& adj : base_undirected_adj) {
        std::sort(adj.begin(), adj.end());
        adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
    }

    if (sccs.empty()) {
        std::cerr << "GIR::buildBCT: buildSCC_DAG() must be called first "
                     "(BiCC blocks are computed per SCC)." << std::endl;
        return;
    }

    scc_bicc_blocks.assign(sccs.size(), {});
    scc_ap.assign(sccs.size(), {});
    for (auto* t : scc_bct) if (t) delete t;
    scc_bct.assign(sccs.size(), nullptr);

    auto t0 = std::chrono::high_resolution_clock::now();
    for (size_t id = 0; id < sccs.size(); ++id) {
        if (!scc_alive[id] || sccs[id].empty()) continue;
        rebuildSccBicc((int)id);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_time = t1 - t0;
    std::cout << "Built per-SCC BC trees for " << sccs.size()
              << " SCCs in " << total_time.count() << " seconds\n";
}

// -------------------- SCC + DAG 构建 --------------------
void GIR::buildSCC_DAG() {
    std::cout << "=== Build SCC + DAG ===" << std::endl;

    // NOTE: this used to call scc_detection() (includes/scc_cpu/scc_core.cpp,
    // the parallel trim + elephant-fw/bw + mice-fw/bw pipeline). While
    // validating GIR's dynamic-update rewrite against a from-scratch
    // rebuild, that function was found to also under-merge SCCs on the
    // small 18-vertex example_graph -- e.g. vertices 2 and 9 have edges in
    // both directions (2->9 and 9->2), which should merge {0,1,2} and
    // {9,10,11} into one 6-vertex SCC, but scc_detection() reports them as
    // two separate SCCs (8 groups total instead of the correct 7). This is
    // on top of the SCC-label collision bug already fixed in fw_bw.h/
    // scc_core.cpp -- see docs/code-review.md. Given the risk of further
    // undiscovered bugs in that multi-stage pipeline, construction here
    // uses the same scoped Tarjan (dynamic_repair::find_scc_scoped) that
    // incremental repair already relies on, run once over the whole
    // vertex set, rather than trusting the parallel implementation.
    auto t0 = std::chrono::high_resolution_clock::now();
    {
        std::vector<int> all_verts(g_directed->vert_count);
        for (index_t i = 0; i < g_directed->vert_count; ++i) all_verts[i] = (int)i;
        auto out_fn = [this](int x) { return directedOutNeighbors(x); };
        sccs = dynamic_repair::find_scc_scoped(all_verts, out_fn);
        nodeToScc.clear();
        for (size_t i = 0; i < sccs.size(); ++i)
            for (int v : sccs[i]) nodeToScc[v] = (int)i;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_time1 = t1 - t0;
    std::cout << "Time to find SCC: " << total_time1.count() << " seconds\n";

    scc_alive.assign(sccs.size(), 1);

    dag_adj_set.assign(sccs.size(), {});
    dag_radj_set.assign(sccs.size(), {});
    scc_bicc_blocks.assign(sccs.size(), {});
    scc_ap.assign(sccs.size(), {});
    scc_bct.assign(sccs.size(), nullptr);

    for (index_t u = 0; u < g_directed->vert_count; ++u) {
        auto uit = nodeToScc.find(u);
        if (uit == nodeToScc.end()) continue;
        int su = uit->second;
        for (index_t i = g_directed->fw_beg_pos[u]; i < g_directed->fw_beg_pos[u + 1]; ++i) {
            vertex_t w = g_directed->fw_csr[i];
            auto wit = nodeToScc.find(w);
            if (wit == nodeToScc.end()) continue;
            int sw = wit->second;
            if (sw != su) { dag_adj_set[su][sw]++; dag_radj_set[sw][su]++; }
        }
    }

    std::vector<int> indeg(sccs.size(), 0);
    for (size_t s = 0; s < sccs.size(); ++s)
        for (auto& kv : dag_adj_set[s]) indeg[kv.first]++;
    std::vector<int> queue_nodes;
    for (size_t s = 0; s < sccs.size(); ++s) if (indeg[s] == 0) queue_nodes.push_back((int)s);
    std::vector<int> topo_order;
    topo_order.reserve(sccs.size());
    size_t qi = 0;
    while (qi < queue_nodes.size()) {
        int s = queue_nodes[qi++];
        topo_order.push_back(s);
        for (auto& kv : dag_adj_set[s]) if (--indeg[kv.first] == 0) queue_nodes.push_back(kv.first);
    }
    scc_rank.assign(sccs.size(), 0);
    rank_to_scc.clear();
    for (size_t i = 0; i < topo_order.size(); ++i) {
        long long r = (long long)(i + 1) * RANK_GAP;
        scc_rank[topo_order[i]] = r;
        rank_to_scc[r] = topo_order[i];
    }

    std::map<int, int> unused_map;
    std::unordered_map<int, std::vector<int>> unused_adj;
    generate_and_save_dag(unused_map, unused_adj, "dag.bin");

    std::cout << "Found " << sccs.size() << " SCCs, "
              << topo_order.size() << " ranked in initial condensation DAG." << std::endl;
}

void GIR::generate_and_save_dag(const std::map<int, int>& /*nodeToScc_unused*/,
                                const std::unordered_map<int, std::vector<int>>& /*adj_list_unused*/,
                                const std::string& output_filename) {
    std::ofstream outfile(output_filename, std::ios::out | std::ios::binary);
    if (!outfile) {
        std::cerr << "Error: cannot open file " << output_filename << std::endl;
        return;
    }
    int scc_count = (int)sccs.size();
    outfile.write(reinterpret_cast<const char*>(&scc_count), sizeof(int));
    for (int i = 0; i < scc_count; ++i) {
        std::vector<int> neighbors;
        neighbors.reserve(dag_adj_set[i].size());
        for (auto& kv : dag_adj_set[i]) neighbors.push_back(kv.first);
        int neighbor_count = (int)neighbors.size();
        outfile.write(reinterpret_cast<const char*>(&neighbor_count), sizeof(int));
        if (neighbor_count > 0)
            outfile.write(reinterpret_cast<const char*>(neighbors.data()), neighbor_count * sizeof(int));
    }
    outfile.close();
    std::cout << "DAG saved to " << output_filename << " (" << scc_count << " SCCs)" << std::endl;
}

// ==================== SPU dispatch (five-case classifier) ====================

void GIR::insertEdge(int u, int v) {
    bool has_d = isDirectedEdgePresent(u, v);
    bool has_u = isUndirectedEdgePresent(u, v);
    if (has_d && has_u) return;  // already fully present, no-op
    if (!has_d) addOverlayDirectedEdge(u, v);
    if (!has_u) addOverlayUndirectedEdge(u, v);
    handleInsert(u, v);
}

void GIR::deleteEdge(int u, int v) {
    bool had_d = isDirectedEdgePresent(u, v);
    bool had_u = isUndirectedEdgePresent(u, v);
    if (!had_d && !had_u) return;  // already fully absent, no-op
    if (had_d) removeOverlayDirectedEdge(u, v);
    if (had_u) removeOverlayUndirectedEdge(u, v);
    handleDelete(u, v);
}

void GIR::updateEdgeDirected(int u, int v, bool insert) {
    if (insert) insertEdge(u, v);
    else deleteEdge(u, v);
}

void GIR::updateEdgeUndirected(int u, int v, bool insert) {
    if (insert) insertEdge(u, v);
    else deleteEdge(u, v);
}

void GIR::handleInsert(int u, int v) {
    int wcc_u = findWcc(u), wcc_v = findWcc(v);
    int su = nodeToScc.at(u), sv = nodeToScc.at(v);

    if (wcc_u != wcc_v) {
        // Case 1: different WCCs. Shift whichever side has fewer vertices
        // (a cheap proxy for "smaller condensation graph", available in
        // O(1) via the union-find size *before* merging) -- shifting
        // whichever component happens to contain the edge's source vertex
        // unconditionally is fine for small graphs, but on any real graph
        // with one dominant giant WCC, about half of random cross-WCC
        // insertions would then shift that giant component, turning an
        // O(small side) operation into O(whole graph) each time.
        bool shift_u_side = wcc_size[wcc_u] <= wcc_size[wcc_v];
        unionWcc(u, v);
        if (shift_u_side) shiftComponentBelow(su, scc_rank[sv]);
        else shiftComponentAbove(sv, scc_rank[su]);
        addDagEdge(su, sv);
        return;
    }

    if (su != sv) {
        // Case 2: same WCC, different SCCs.
        if (scc_rank[su] <= scc_rank[sv]) {
            addDagEdge(su, sv);
            return;
        }
        long long lo_bound = scc_rank[sv], hi_bound = scc_rank[su];
        std::unordered_set<int> fwd;
        {
            std::queue<int> q;
            q.push(sv);
            fwd.insert(sv);
            while (!q.empty()) {
                int c = q.front();
                q.pop();
                for (auto& kv : dag_adj_set[c])
                    if (scc_rank[kv.first] <= hi_bound && fwd.insert(kv.first).second) q.push(kv.first);
            }
        }
        if (fwd.count(su)) {
            std::unordered_set<int> bwd;
            std::queue<int> q;
            q.push(su);
            bwd.insert(su);
            while (!q.empty()) {
                int c = q.front();
                q.pop();
                for (auto& kv : dag_radj_set[c])
                    if (scc_rank[kv.first] >= lo_bound && bwd.insert(kv.first).second) q.push(kv.first);
            }
            std::vector<int> merge_ids;
            for (int id : fwd) if (bwd.count(id)) merge_ids.push_back(id);
            if (merge_ids.size() < 2) merge_ids = {su, sv};
            mergeSccs(merge_ids);
        } else {
            addDagEdge(su, sv);
        }
        repairRankInterval(lo_bound, hi_bound);
        return;
    }

    // Same SCC: Case 3 or Case 4.
    if (inSameBiccBlock(su, u, v)) {
        return;  // Case 4: true no-op.
    }
    rebuildSccBicc(su);  // Case 3.
}

void GIR::handleDelete(int u, int v) {
    int su = nodeToScc.at(u), sv = nodeToScc.at(v);

    if (su == sv) {
        bool was_bridge = isBridgeWithinScc(su, u, v);
        if (!was_bridge) {
            rebuildSccBicc(su);  // Case 5a.
            return;
        }
        splitScc(su);  // Case 5b (splitScc itself no-ops the BiCC layer if it turns out not to actually split).
        return;
    }

    // Inter-SCC edge: it was itself a condensation-DAG edge.
    removeDagEdge(su, sv);
    if (probeCondensationConnectivity(su, sv)) {
        return;  // still connected via some other SCC -- nothing else to repair.
    }

    // Case 5c: the WCC actually splits.
    std::vector<int> side_a;
    probeUndirectedConnectivity(u, v, side_a);
    int old_root = findWcc(u);
    std::unordered_set<int> side_a_set(side_a.begin(), side_a.end());
    std::vector<int> side_b;
    for (size_t i = 0; i < wcc_parent.size(); ++i) {
        if (findWcc((int)i) == old_root && !side_a_set.count((int)i)) side_b.push_back((int)i);
    }
    rebuildWccComponent(side_a);
    rebuildWccComponent(side_b);
}

// -------------------- Getter --------------------
std::map<int, int> GIR::getNodeToScc() const { return nodeToScc; }

std::map<int, int> GIR::getNodeToWcc() const {
    std::map<int, int> result;
    for (size_t i = 0; i < wcc_parent.size(); ++i) {
        int v = (int)i;
        while (wcc_parent[v] != v) v = wcc_parent[v];
        result[(int)i] = v;
    }
    return result;
}

const std::vector<std::vector<int>>& GIR::getSCCs() const { return sccs; }

std::vector<std::vector<int>> GIR::getDAG() const {
    std::vector<std::vector<int>> result(dag_adj_set.size());
    for (size_t i = 0; i < dag_adj_set.size(); ++i) {
        result[i].reserve(dag_adj_set[i].size());
        for (auto& kv : dag_adj_set[i]) result[i].push_back(kv.first);
    }
    return result;
}

std::vector<std::vector<int>> GIR::getWCCs() const {
    std::unordered_map<int, std::vector<int>> groups;
    for (size_t i = 0; i < wcc_parent.size(); ++i) {
        int v = (int)i;
        while (wcc_parent[v] != v) v = wcc_parent[v];
        groups[v].push_back((int)i);
    }
    std::vector<std::vector<int>> result;
    result.reserve(groups.size());
    for (auto& kv : groups) result.push_back(std::move(kv.second));
    return result;
}

BlockCutTree* GIR::getBCT(int scc_id) const {
    if (scc_id < 0 || (size_t)scc_id >= scc_bct.size()) return nullptr;
    return scc_bct[scc_id];
}

const std::vector<std::vector<int>>& GIR::getBCTComponents(int scc_id) const {
    static const std::vector<std::vector<int>> empty;
    if (scc_id < 0 || (size_t)scc_id >= scc_bicc_blocks.size()) return empty;
    return scc_bicc_blocks[scc_id];
}

graph* GIR::getDirectedGraph() const { return g_directed; }
graph_undirected* GIR::getUndirectedGraph() const { return g_undirected; }

std::vector<int> GIR::getDirectedOutNeighbors(int v) const { return directedOutNeighbors(v); }
std::vector<int> GIR::getDirectedInNeighbors(int v) const { return directedInNeighbors(v); }
std::vector<int> GIR::getUndirectedNeighborsPublic(int v) const { return undirectedNeighbors(v); }
bool GIR::hasDirectedEdge(int u, int v) const { return isDirectedEdgePresent(u, v); }
bool GIR::hasUndirectedEdge(int u, int v) const { return isUndirectedEdgePresent(u, v); }
int GIR::vertexCount() const { return g_directed ? (int)g_directed->vert_count : 0; }

// ==================== Two-level parallel BFS (paper/task_parallel.tex) ====================

std::vector<int> GIR::bfsTaskParallel(int source) {
    const int V = vertexCount();
    std::vector<int> levels(V, -1);
    if (source < 0 || source >= V) return levels;
    auto src_it = nodeToScc.find(source);
    if (src_it == nodeToScc.end()) return levels;
    int src_scc = src_it->second;

    // Stage 1: prune to SCCs reachable from src_scc via the condensation DAG.
    // A hash-set membership test per edge is far more expensive than an O(1)
    // bitmap indexed by scc id (ids are dense indices into
    // `sccs`/`scc_alive`, see gir.h) -- but the real reason this runs FIRST,
    // before building any per-vertex lookup table below, is that its cost is
    // bounded by the size of the reachable region, not by V. Graphs with
    // many small/singleton SCCs (e.g. web-Google: 7.77M WCCs out of 8.64M
    // vertices -- a random BFS source typically reaches almost nothing) need
    // that: building any O(V) structure unconditionally before this would
    // dominate the whole call for such queries. See the note on scc_of
    // below for the concrete regression this caused.
    std::vector<char> reachable_sccs(sccs.size(), 0);
    {
        std::queue<int> q;
        q.push(src_scc);
        reachable_sccs[src_scc] = 1;
        while (!q.empty()) {
            int c = q.front();
            q.pop();
            for (auto& kv : dag_adj_set[c]) {
                if (!reachable_sccs[kv.first]) {
                    reachable_sccs[kv.first] = 1;
                    q.push(kv.first);
                }
            }
        }
    }

    // vertex -> scc id, restricted to vertices in reachable SCCs (via the
    // `sccs` vertex-lists, so this is bounded by the reachable region's
    // vertex count, not V). Needed to avoid an O(log V) nodeToScc lookup per
    // edge examined below -- but it must NOT be built by iterating all of
    // nodeToScc unconditionally: that was tried first and, while it fixed
    // task-parallel BFS being ~16x slower than sequential on twitter/
    // livejournal/RoadUSA (all dominated by one giant WCC, so an O(V) pass
    // is roughly the same order as the real work anyway), it made web-Google
    // ~220x slower (1.28ms -> 283ms) since almost every BFS there visits
    // only a handful of vertices and was paying a full 8.64M-entry
    // std::map iteration regardless. Scoping to reachable SCCs' vertex lists
    // fixes both: bounded by whichever is actually small for a given query.
    std::vector<int> scc_of(V, -1);
    for (int scc_id = 0; scc_id < (int)sccs.size(); ++scc_id) {
        if (!reachable_sccs[scc_id]) continue;
        for (int v : sccs[scc_id]) scc_of[v] = scc_id;
    }

    // Earlier versions of this function partitioned each level's frontier by
    // (scc, block) and dispatched one OpenMP task per group, to mirror the
    // paper's "task parallelism across blocks" model. Measured on real
    // datasets (see CHANGELOG) that grouping step -- and the block-index
    // precompute it depended on -- was actively harmful: most SCCs here have
    // no multi-vertex BiCC block (near-singleton SCCs dominate on e.g.
    // livejournal/twitter), so grouping degenerated into "one task per
    // vertex", and OpenMP's dynamic-scheduling overhead per task swamped the
    // ~14-edges-per-vertex of actual work each task did. Diagnostic timing
    // showed the "parallel" region taking about as long as the entire
    // sequential BFS baseline despite 128 cores being available -- i.e. the
    // task dispatch overhead alone matched a full single-threaded pass.
    // Replaced with plain data parallelism directly over the frontier array,
    // chunked to amortize scheduling overhead across a batch of vertices
    // instead of one at a time. This is no longer the paper's block-task
    // model -- see CHANGELOG for why that trade was made and what it means
    // for the task-parallel/data-parallel distinction the paper describes.
    const bool overlay_empty = extra_out_adj.empty() && removed_directed.empty();
    const int nthreads = omp_get_max_threads();
    levels[source] = 0;
    std::vector<int> frontier = {source};
    int round = 0;

    while (!frontier.empty()) {
        const size_t n = frontier.size();
        std::vector<std::vector<int>> per_thread_discovered(nthreads);
        const size_t chunk = 256;
        // Below this, OpenMP's own parallel-region overhead (thread-team
        // wake-up, per-level even when the team ends up doing almost no
        // work) costs more than just doing the level in the calling thread.
        // Real effect, not a hypothetical: on web-Google (7.77M WCCs, almost
        // all singleton/tiny -- a typical BFS from a random source only
        // reaches a handful of vertices per level) unconditional
        // parallelization made task-parallel BFS ~19x *slower* than
        // sequential even after the fixes above, purely from this fixed
        // per-level overhead. `if()` on the pragma serializes the loop in
        // the calling thread with no thread-team spin-up when below
        // threshold, rather than duplicating the loop body.
        static const size_t PAR_THRESHOLD = 2000;
        #pragma omp parallel for schedule(dynamic, chunk) if(n >= PAR_THRESHOLD)
        for (size_t fi = 0; fi < n; ++fi) {
            int u = frontier[fi];
            std::vector<int>& local_discovered = per_thread_discovered[omp_get_thread_num()];
            if (overlay_empty) {
                for (index_t i = g_directed->fw_beg_pos[u]; i < g_directed->fw_beg_pos[u + 1]; ++i) {
                    int w = (int)g_directed->fw_csr[i];
                    if (levels[w] != -1) continue;
                    int wscc = scc_of[w];
                    if (wscc < 0 || !reachable_sccs[wscc]) continue;
                    local_discovered.push_back(w);
                }
            } else {
                for (int w : directedOutNeighbors(u)) {
                    if (levels[w] != -1) continue;
                    int wscc = scc_of[w];
                    if (wscc < 0 || !reachable_sccs[wscc]) continue;
                    local_discovered.push_back(w);
                }
            }
        }

        std::vector<int> next_frontier;
        for (auto& dl : per_thread_discovered) {
            for (int w : dl) {
                if (levels[w] == -1) {  // sequential merge: no race, first writer per level wins
                    levels[w] = round + 1;
                    next_frontier.push_back(w);
                }
            }
        }
        frontier = std::move(next_frontier);
        ++round;
    }

    return levels;
}
