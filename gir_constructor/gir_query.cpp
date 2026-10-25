// Downstream analytics on GIR (paper Sec. 6 and 7.3): BFS, SSSP, BC.
//
//  bfs            direction-optimizing BFS. GIR's condensation DAG gives the
//                 exact set of vertices a source can reach, which is the
//                 candidate set of every bottom-up step (instead of all V).
//  bfsTwoLevel    the paper's two-level model: SCCs are processed as a
//                 wavefront over the condensation DAG (independent SCCs are
//                 concurrent tasks); inside the source SCC, BC-tree blocks are
//                 concurrent tasks started from their entry articulation point;
//                 large units use data parallelism inside.
//  sssp           Delta-stepping over the reachable region.
//  ssspTwoLevel   same DAG wavefront; singleton SCCs are settled directly,
//                 Delta-stepping runs only inside non-trivial SCCs.
//  bc             Brandes betweenness from a set of sources, sources processed
//                 as parallel tasks.
#include "gir.h"

#include <algorithm>
#include <parallel/algorithm>
#include <atomic>
#include <climits>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include <omp.h>

#include "graph.h"
#include "gir_iter.h"

namespace {

template <class T>
inline bool atomicMin(T* p, T v) {
    T old = __atomic_load_n(p, __ATOMIC_RELAXED);
    while (v < old)
        if (__atomic_compare_exchange_n(p, &old, v, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return true;
    return false;
}

inline void atomicAddDouble(double* p, double v) {
    #pragma omp atomic
    *p += v;
}


}  // namespace

uint32_t GIR::edgeWeight(int u, int v) {
    uint64_t x = ((uint64_t)(uint32_t)u << 32) | (uint32_t)v;
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return 1 + (uint32_t)(x % kMaxWeight);
}

std::vector<char> GIR::reachableSccs(int src_scc) const {
    std::vector<char> reach(sccs.size(), 0);
    std::vector<int> q{src_scc};
    reach[src_scc] = 1;
    for (size_t i = 0; i < q.size(); ++i)
        for (auto& kv : dag_out[q[i]])
            if (!reach[kv.first]) { reach[kv.first] = 1; q.push_back(kv.first); }
    return reach;
}

const GIR::ReachCache& GIR::reachInfo(int src_scc) const {
    if (reach_cache.version == struct_version && reach_cache.scc == src_scc) return reach_cache;
    reach_cache.sccs = reachableSccs(src_scc);
    reach_cache.verts.clear();
    for (size_t c = 0; c < sccs.size(); ++c)
        if (reach_cache.sccs[c]) reach_cache.verts.insert(reach_cache.verts.end(), sccs[c].begin(), sccs[c].end());
    long long e = 0;
    const auto& vs = reach_cache.verts;
    const index_t* fb = g_directed->fw_beg_pos;
    #pragma omp parallel for reduction(+ : e) schedule(static)
    for (size_t i = 0; i < vs.size(); ++i) e += (fb[vs[i] + 1] - fb[vs[i]]) + (long long)extra_out_adj[vs[i]].size();
    reach_cache.edges = e;
    reach_cache.pending.assign(sccs.size(), 0);
    for (size_t c = 0; c < sccs.size(); ++c) {
        if (!reach_cache.sccs[c]) continue;
        int cnt = 0;
        for (auto& kv : dag_in[c]) if (reach_cache.sccs[kv.first]) ++cnt;
        reach_cache.pending[c] = cnt;
    }
    reach_cache.scc = src_scc;
    reach_cache.version = struct_version;
    return reach_cache;
}

// ============================================================================
// Direction-optimizing BFS
// ============================================================================

namespace {
// Concatenate per-thread buffers into `out` in parallel.
void gatherParts(std::vector<std::vector<int>>& parts, std::vector<int>& out) {
    std::vector<size_t> off(parts.size() + 1, 0);
    for (size_t t = 0; t < parts.size(); ++t) off[t + 1] = off[t] + parts[t].size();
    out.resize(off.back());
    #pragma omp parallel for schedule(static, 1)
    for (size_t t = 0; t < parts.size(); ++t) {
        std::copy(parts[t].begin(), parts[t].end(), out.begin() + off[t]);
        parts[t].clear();
    }
}
}  // namespace

QueryVec<int> GIR::bfs(int source) {
    QueryVec<int> level(V);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < V; ++i) level[i] = -1;
    if (source < 0 || source >= V) return level;
    level[source] = 0;

    // GIR's condensation gives the exact size of the reachable region (for
    // the direction switch). It is used when already cached for this source
    // SCC; computing it costs about one BFS, so a query never builds it.
    const bool use_reach = pending.empty() && reach_cache.version == struct_version && reach_cache.scc == node_scc[source];
    const long long n_reach = use_reach ? (long long)reach_cache.verts.size() : (long long)V;
    const long long m_total = use_reach ? reach_cache.edges
                                        : (long long)g_directed->edge_count + extra_count - n_dead_dir;

    const index_t* fb = g_directed->fw_beg_pos;
    const vertex_t* fc = g_directed->fw_csr;
    const index_t* bb = g_directed->bw_beg_pos;
    const vertex_t* bc = g_directed->bw_csr;
    const bool clean = extra_count == 0 && n_dead_dir == 0;
    const bool chk_dead = n_dead_dir != 0;
    const bool no_extra = extra_count == 0;  // (then skip the per-vertex extra lists: a random read each)
    auto outdeg = [&](int v) { return (long long)(fb[v + 1] - fb[v]) + (no_extra ? 0LL : (long long)extra_out_adj[v].size()); };

    // Scratch reused across queries of this thread (already faulted in).
    struct Scratch {
        std::vector<int> frontier, next;
        std::vector<long long> pre;
        std::vector<uint64_t> fbits, nbits, vbits;
        std::vector<std::vector<int>> local;
    };
    static thread_local Scratch sc;
    auto& frontier = sc.frontier; auto& next = sc.next;
    frontier.assign(1, source);
    next.clear();
    const size_t W = ((size_t)V + 63) / 64;
    auto& fbits = sc.fbits; auto& nbits = sc.nbits;
    bool bits_valid = false;  // fbits holds the current frontier
    // Visited bitmap (V/8 bytes, cache-resident): a dense step skips a word
    // whose 64 vertices are all visited without touching `level`.
    auto& vbits = sc.vbits;
    vbits.resize(W);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < W; ++i) vbits[i] = 0;
    vbits[(size_t)source >> 6] |= 1ULL << (source & 63);
    const int nthreads = omp_get_max_threads();
    auto& local = sc.local;
    if ((int)local.size() < nthreads) local.resize(nthreads);
    for (auto& p : local) p.clear();
    std::vector<long long> m_part(nthreads);

    long long m_f = outdeg(source);
    long long m_unexplored = m_total - m_f;
    const long long alpha = 15, beta = 18;
    int L = 0;
    bool bottom_up = false;

    while (!frontier.empty()) {
        if (n_reach >= (long long)q_large_unit) {
            if (!bottom_up && m_f * alpha > m_unexplored) bottom_up = true;
            else if (bottom_up && (long long)frontier.size() * beta < n_reach) bottom_up = false;
        }
        std::fill(m_part.begin(), m_part.end(), 0);
        auto claim = [&](int w, std::vector<int>& out, long long& m) {
            // The bitmap answers "visited?" from cache; `level` is touched only for new vertices.
            if ((vbits[(size_t)w >> 6] >> (w & 63)) & 1) return;  // (visited before this level)
            if (__atomic_load_n(&level[w], __ATOMIC_RELAXED) != -1) return;
            int expect = -1;
            if (__atomic_compare_exchange_n(&level[w], &expect, L + 1, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
                out.push_back(w);
                m += outdeg(w);
            }
        };
        if (!bottom_up) {
            const bool par = frontier.size() >= q_par_frontier || m_f >= (1 << 16);
            if (par && m_f >= (1 << 16)) {
                // Edge-balanced (Ligra-style): the frontier's out-edges are cut
                // into equal chunks, so a few hubs don't serialize the level.
                // Position k of vertex u: base edge k, then u's extra edges.
                auto& pre = sc.pre;
                pre.resize(frontier.size() + 1);
                pre[0] = 0;
                for (size_t i = 0; i < frontier.size(); ++i) pre[i + 1] = pre[i] + outdeg(frontier[i]);
                const long long total = pre.back(), CH = 2048;
                const long long nch = (total + CH - 1) / CH;
                #pragma omp parallel
                {
                    const int t = omp_get_thread_num();
                    auto& out = local[t];
                    long long m = 0;
                    #pragma omp for schedule(dynamic, 4)
                    for (long long ch = 0; ch < nch; ++ch) {
                        long long lo = ch * CH, hi = std::min(total, lo + CH);
                        size_t i = std::upper_bound(pre.begin(), pre.end(), lo) - pre.begin() - 1;
                        for (long long p = lo; p < hi; ++i) {
                            const int u = frontier[i];
                            const long long base = fb[u + 1] - fb[u];
                            const long long end = std::min(hi, pre[i + 1]);
                            for (; p < end; ++p) {
                                const long long k = p - pre[i];
                                if (k < base) {
                                    const index_t e = fb[u] + k;
                                    if (chk_dead && bitGet(dead_fw, (long long)e)) continue;
                                    claim((int)fc[e], out, m);
                                } else {
                                    claim(extra_out_adj[u][k - base], out, m);
                                }
                            }
                        }
                    }
                    m_part[t] = m;
                }
            } else {
                #pragma omp parallel for schedule(dynamic, 64) if (par)
                for (size_t i = 0; i < frontier.size(); ++i) {
                    const int t = omp_get_thread_num();
                    auto& out = local[t];
                    long long& m = m_part[t];
                    int u = frontier[i];
                    if (clean) for (index_t k = fb[u]; k < fb[u + 1]; ++k) claim((int)fc[k], out, m);
                    else forOut(u, [&](int w) { claim(w, out, m); });
                }
            }
            bits_valid = false;
        } else {
            // Dense step: every unvisited vertex looks for a parent in the
            // frontier bitmap; the next frontier is written as a bitmap too,
            // one 64-vertex word per iteration (no atomics, no candidate lists).
            fbits.resize(W);
            nbits.resize(W);
            if (!bits_valid) {
                #pragma omp parallel for schedule(static)
                for (size_t i = 0; i < W; ++i) fbits[i] = 0;
                #pragma omp parallel for schedule(static)
                for (size_t i = 0; i < frontier.size(); ++i) {
                    int f = frontier[i];
                    __atomic_fetch_or(&fbits[(size_t)f >> 6], 1ULL << (f & 63), __ATOMIC_RELAXED);
                }
            }
            const uint64_t* fbp = fbits.data();
            auto inF = [fbp](int p) { return (fbp[(unsigned)p >> 6] >> (p & 63)) & 1; };
            #pragma omp parallel
            {
                const int t = omp_get_thread_num();
                auto& out = local[t];
                long long m = 0;
                #pragma omp for schedule(dynamic, 64) nowait
                for (size_t wd = 0; wd < W; ++wd) {
                    uint64_t bits = 0;
                    const uint64_t vis = vbits[wd];
                    const int v0 = (int)(wd << 6), v1 = std::min(V, v0 + 64);
                    if (vis == ~0ULL) { nbits[wd] = 0; continue; }
                    for (int v = v0; v < v1; ++v) {
                        if ((vis >> (v - v0)) & 1) continue;
                        bool hit = false;
                        if (clean) {
                            for (index_t k = bb[v]; k < bb[v + 1]; ++k) if (inF(bc[k])) { hit = true; break; }
                        } else {
                            for (index_t k = bb[v]; k < bb[v + 1]; ++k)
                                if (inF(bc[k]) && !(chk_dead && bitGet(dead_bw, (long long)k))) { hit = true; break; }
                            if (!hit) for (int p : extra_in_adj[v]) if (inF(p)) { hit = true; break; }
                        }
                        if (hit) { level[v] = L + 1; bits |= 1ULL << (v - v0); out.push_back(v); m += outdeg(v); }
                    }
                    nbits[wd] = bits;
                    vbits[wd] = vis | bits;
                }
                m_part[t] = m;
            }
            fbits.swap(nbits);
            bits_valid = true;
        }
        gatherParts(local, next);
        if (!bottom_up) {  // top-down claims set `level` only; mark them visited now
            #pragma omp parallel for schedule(static) if (next.size() >= q_par_frontier)
            for (size_t i = 0; i < next.size(); ++i)
                __atomic_fetch_or(&vbits[(size_t)next[i] >> 6], 1ULL << (next[i] & 63), __ATOMIC_RELAXED);
        }
        long long m_next = 0;
        for (long long m : m_part) m_next += m;
        m_unexplored -= m_next;
        m_f = m_next;
        frontier.swap(next);
        ++L;
    }
    return level;
}

QueryVec<int> GIR::bfsTaskParallel(int source) { return bfs(source); }

// ============================================================================
// Two-level model: DAG wavefront over SCCs, BC-tree wavefront inside the
// source SCC.
// ============================================================================

// Shared driver. `inner(c, is_source_scc)` finalizes dist[] inside SCC c given
// the tentative values its predecessors wrote; afterwards the driver relaxes
// c's outgoing condensation edges with `relax(x, w)` and releases successors
// whose predecessors are all done.
template <class Inner, class Relax>
void GIR::dagWavefront(int source, Inner inner, Relax relax) {
    const int s0 = node_scc[source];
    const ReachCache& ri = reachInfo(s0);
    // Remaining reachable-predecessor counts (cached with the reachable set).
    std::vector<std::atomic<int>> left(sccs.size());
    #pragma omp parallel for schedule(static)
    for (size_t c = 0; c < sccs.size(); ++c) left[c].store(ri.pending[c], std::memory_order_relaxed);

    std::vector<int> ready{s0};
    while (!ready.empty()) {
        std::vector<int> small, large;
        for (int c : ready) (sccs[c].size() >= q_large_unit ? large : small).push_back(c);
        // Large SCCs one at a time (data parallelism inside), small ones as
        // concurrent tasks.
        for (int c : large) inner(c, c == s0);
        #pragma omp parallel for schedule(dynamic, 16) if (small.size() > std::min<size_t>(64, q_par_frontier))
        for (size_t i = 0; i < small.size(); ++i) inner(small[i], small[i] == s0);

        const int nthreads = omp_get_max_threads();
        std::vector<std::vector<int>> nxt(nthreads);
        auto finish = [&](int c) {
            auto& out = nxt[omp_get_thread_num()];
            for (int x : sccs[c])
                forOut(x, [&](int w) { if (node_scc[w] != c) relax(x, w); });
            for (auto& kv : dag_out[c])
                if (left[kv.first].fetch_sub(1, std::memory_order_acq_rel) == 1) out.push_back(kv.first);
        };
        // Large SCCs: relax their (many) outgoing edges in parallel over
        // vertices; small ones: one task per SCC.
        for (int c : large) {
            const auto& vs = sccs[c];
            #pragma omp parallel for schedule(dynamic, 1024)
            for (size_t i = 0; i < vs.size(); ++i)
                forOut(vs[i], [&](int w) { if (node_scc[w] != c) relax(vs[i], w); });
            for (auto& kv : dag_out[c])
                if (left[kv.first].fetch_sub(1, std::memory_order_acq_rel) == 1) nxt[0].push_back(kv.first);
        }
        #pragma omp parallel for schedule(dynamic, 16) if (small.size() > std::min<size_t>(64, q_par_frontier))
        for (size_t i = 0; i < small.size(); ++i) finish(small[i]);
        ready.clear();
        for (auto& v : nxt) ready.insert(ready.end(), v.begin(), v.end());
    }
}

// Direction-optimizing BFS restricted to one unit (an SCC or a block):
// `verts` are the unit's vertices, in_scope(w) tests membership, and dist[]
// holds tentative levels for some of them (several seeds, possibly at
// different levels; INT_MAX = unknown). Seeds join the frontier when the
// sweep reaches their level. Data-parallel; used for the large units of the
// two-level model.
template <class InScope>
void GIR::scopedBfs(const std::vector<int>& verts, InScope in_scope, QueryVec<int>& dist) {
    std::map<int, std::vector<int>> seeds;
    for (int x : verts) if (dist[x] != INT_MAX) seeds[dist[x]].push_back(x);
    if (seeds.empty()) return;
    if (q_front.size() < (size_t)V) q_front.assign(V, 0);
    const index_t* fb = g_directed->fw_beg_pos;
    const index_t* bb = g_directed->bw_beg_pos;
    const vertex_t* bc = g_directed->bw_csr;
    const bool clean = extra_count == 0 && n_dead_dir == 0;
    auto outdeg = [&](int v) { return (long long)(fb[v + 1] - fb[v]) + (long long)extra_out_adj[v].size(); };
    long long m_unexplored = 0;
    #pragma omp parallel for reduction(+ : m_unexplored) schedule(static)
    for (size_t i = 0; i < verts.size(); ++i) m_unexplored += outdeg(verts[i]);
    const int nthreads = omp_get_max_threads();
    const long long alpha = 15, beta = 18, n = (long long)verts.size();
    std::vector<int> cand;
    bool bottom_up = false, have_cand = false;
    int L = seeds.begin()->first;
    std::vector<int> frontier = std::move(seeds.begin()->second);
    seeds.erase(seeds.begin());
    while (!frontier.empty() || !seeds.empty()) {
        if (frontier.empty()) {  // jump to the next seed level
            L = seeds.begin()->first;
            for (int x : seeds.begin()->second) if (dist[x] == L) frontier.push_back(x);
            seeds.erase(seeds.begin());
            continue;
        }
        long long m_f = 0;
        #pragma omp parallel for reduction(+ : m_f) schedule(static) if (frontier.size() >= q_par_frontier)
        for (size_t i = 0; i < frontier.size(); ++i) m_f += outdeg(frontier[i]);
        if (!bottom_up && m_f * alpha > m_unexplored) bottom_up = true;
        else if (bottom_up && (long long)frontier.size() * beta < n) bottom_up = false;
        if (bottom_up && !have_cand) {
            for (int x : verts) if (dist[x] == INT_MAX || dist[x] > L + 1) cand.push_back(x);
            have_cand = true;
        }
        std::vector<std::vector<int>> local(nthreads);
        if (!bottom_up) {
            #pragma omp parallel for schedule(dynamic, 256) if (frontier.size() >= q_par_frontier)
            for (size_t i = 0; i < frontier.size(); ++i) {
                auto& out = local[omp_get_thread_num()];
                forOut(frontier[i], [&](int w) {
                    if (dist[w] > L + 1 && in_scope(w) && atomicMin(&dist[w], L + 1)) out.push_back(w);
                });
            }
        } else {
            #pragma omp parallel for schedule(static) if (frontier.size() >= q_par_frontier)
            for (size_t i = 0; i < frontier.size(); ++i) q_front[frontier[i]] = 1;
            #pragma omp parallel for schedule(dynamic, 1024)
            for (size_t i = 0; i < cand.size(); ++i) {
                int v = cand[i];
                if (dist[v] <= L + 1) continue;
                bool hit = false;
                for (index_t k = bb[v]; k < bb[v + 1] && !hit; ++k) {
                    int p = bc[k];
                    if (q_front[p] && (clean || !(n_dead_dir && bitGet(dead_bw, (long long)k)))) hit = true;
                }
                if (!hit) for (int p : extra_in_adj[v]) if (q_front[p]) { hit = true; break; }
                if (hit) { dist[v] = L + 1; local[omp_get_thread_num()].push_back(v); }
            }
            #pragma omp parallel for schedule(static) if (frontier.size() >= q_par_frontier)
            for (size_t i = 0; i < frontier.size(); ++i) q_front[frontier[i]] = 0;
        }
        std::vector<int> next;
        for (auto& l : local) next.insert(next.end(), l.begin(), l.end());
        // Seeds sitting exactly at the next level (not reached more cheaply).
        auto it = seeds.find(L + 1);
        if (it != seeds.end()) {
            for (int x : it->second) if (dist[x] == L + 1) next.push_back(x);
            seeds.erase(it);
            std::sort(next.begin(), next.end());
            next.erase(std::unique(next.begin(), next.end()), next.end());
        }
        long long m_next = 0;
        #pragma omp parallel for reduction(+ : m_next) schedule(static) if (next.size() >= q_par_frontier)
        for (size_t i = 0; i < next.size(); ++i) m_next += outdeg(next[i]);
        m_unexplored -= m_next;
        if (have_cand && bottom_up) {
            size_t k = 0;
            for (int v : cand) if (dist[v] > L + 2 || dist[v] == INT_MAX) cand[k++] = v;
            cand.resize(k);
        }
        frontier.swap(next);
        ++L;
    }
}

// Bucketed BFS inside SCC c from whatever tentative levels its vertices hold
// (multi-source, sources at different levels). Data-parallel per bucket for
// large SCCs.
void GIR::bfsInsideScc(int c, QueryVec<int>& dist) {
    const auto& verts = sccs[c];
    if (verts.size() == 1) return;
    if (verts.size() >= q_large_unit) { scopedBfs(verts, [&](int w) { return node_scc[w] == c; }, dist); return; }
    int lo = INT_MAX;
    for (int x : verts) if (dist[x] < lo) lo = dist[x];
    if (lo == INT_MAX) return;
    std::vector<std::vector<int>> bucket(1);
    auto push = [&](int x) {
        size_t k = (size_t)(dist[x] - lo);
        if (bucket.size() <= k) bucket.resize(k + 1);
        bucket[k].push_back(x);
    };
    for (int x : verts) if (dist[x] != INT_MAX) push(x);
    const bool par = verts.size() >= q_large_unit;
    const int nthreads = par ? omp_get_max_threads() : 1;
    for (size_t k = 0; k < bucket.size(); ++k) {
        const int L = lo + (int)k;
        std::vector<int> fr;
        for (int x : bucket[k]) if (dist[x] == L) fr.push_back(x);
        std::sort(fr.begin(), fr.end());
        fr.erase(std::unique(fr.begin(), fr.end()), fr.end());
        std::vector<std::vector<int>> local(nthreads);
        #pragma omp parallel for schedule(dynamic, 256) if (par && fr.size() >= q_par_frontier)
        for (size_t i = 0; i < fr.size(); ++i) {
            auto& out = local[par ? omp_get_thread_num() : 0];
            forOut(fr[i], [&](int w) {
                if (node_scc[w] == c && atomicMin(&dist[w], L + 1)) out.push_back(w);
            });
        }
        for (auto& l : local)
            for (int w : l) if (dist[w] == L + 1) push(w);
    }
}

// Paper Sec. 6, Stages 2-4 inside the source SCC. Every path between two
// blocks passes through the articulation points between them, and a
// shortest path from an entry AP to a vertex of its block never leaves the
// block (it would have to re-enter through the same AP). So each block is
// one task: BFS from its entry AP restricted to the block, offset by the
// entry level. Blocks at the same BC-tree depth run concurrently; large
// blocks use data parallelism inside.
void GIR::bfsSourceSccByBlocks(int source, QueryVec<int>& dist) {
    const int c = node_scc[source];
    struct Task { int block, entry; };
    std::vector<char> seen_block(blk_uf.size(), 0);
    std::vector<Task> round;
    // Child index of this SCC's block-cut tree: (parent AP, block) pairs,
    // collected in parallel from the SCC's vertices (each block once).
    if (!(kids_cache.sv == struct_version && kids_cache.bv == block_version && kids_cache.scc == c)) {
        auto& kids = kids_cache.kids;
        kids.clear();
        const auto& vs = sccs[c];
        std::vector<char> mark(blk_uf.size(), 0);
        const int nt = omp_get_max_threads();
        std::vector<std::vector<std::pair<int, int>>> part(nt);
        #pragma omp parallel for schedule(static) if (vs.size() >= q_par_frontier)
        for (size_t i = 0; i < vs.size(); ++i) {
            int b = findBlkRO(v_block[vs[i]]);
            if (blk_parent_ap[b] < 0) continue;
            if (__atomic_exchange_n(&mark[b], (char)1, __ATOMIC_RELAXED)) continue;
            part[omp_get_thread_num()].push_back({blk_parent_ap[b], b});
        }
        for (auto& p : part) kids.insert(kids.end(), p.begin(), p.end());
        std::sort(kids.begin(), kids.end());
        kids_cache.sv = struct_version;
        kids_cache.bv = block_version;
        kids_cache.scc = c;
    }
    const auto& kids = kids_cache.kids;
    auto blocks_of = [&](int a, std::vector<int>& out) {
        out.clear();
        out.push_back(findBlkRO(v_block[a]));
        if (v_bcnt[a] < 2) return;
        auto lo = std::lower_bound(kids.begin(), kids.end(), std::make_pair(a, INT_MIN));
        for (auto it = lo; it != kids.end() && it->first == a; ++it) out.push_back(it->second);
    };
    std::vector<int> tmp;
    blocks_of(source, tmp);
    for (int b : tmp) { seen_block[b] = 1; round.push_back({b, source}); }
    dist[source] = 0;

    auto run_block = [&](const Task& t, std::vector<int>& new_aps, bool par) {
        const int B = t.block;
        auto in_b = [&](int w) { return findBlkRO(v_block[w]) == B || blk_parent_ap[B] == w; };
        if (par) {  // large block: data-parallel direction-optimizing BFS from the entry AP
            if (blk_stamp.size() < (size_t)V) blk_stamp.assign(V, 0);
            const int ep = ++blk_epoch;
            const auto& bv = blk_verts[B];
            #pragma omp parallel for schedule(static)
            for (size_t i = 0; i < bv.size(); ++i) blk_stamp[bv[i]] = ep;
            scopedBfs(bv, [&](int w) { return blk_stamp[w] == ep; }, dist);
            for (int w : blk_verts[B]) if (v_bcnt[w] >= 2 && w != t.entry) new_aps.push_back(w);
            return;
        }
        std::vector<int> fr{t.entry};
        const int base = dist[t.entry];
        int L = base;
        const int nthreads = par ? omp_get_max_threads() : 1;
        while (!fr.empty()) {
            std::vector<std::vector<int>> local(nthreads);
            #pragma omp parallel for schedule(dynamic, 256) if (par && fr.size() >= q_par_frontier)
            for (size_t i = 0; i < fr.size(); ++i) {
                auto& out = local[par ? omp_get_thread_num() : 0];
                forOut(fr[i], [&](int w) {
                    if (node_scc[w] != c || !in_b(w)) return;
                    if (atomicMin(&dist[w], L + 1)) out.push_back(w);
                });
            }
            std::vector<int> nx;
            for (auto& l : local)
                for (int w : l) if (dist[w] == L + 1) nx.push_back(w);
            std::sort(nx.begin(), nx.end());
            nx.erase(std::unique(nx.begin(), nx.end()), nx.end());
            for (int w : nx) if (v_bcnt[w] >= 2 && w != t.entry) new_aps.push_back(w);
            fr.swap(nx);
            ++L;
        }
    };

    while (!round.empty()) {
        std::vector<std::vector<int>> found_aps(round.size());
        std::vector<size_t> big, little;
        for (size_t i = 0; i < round.size(); ++i)
            (blk_verts[round[i].block].size() >= q_large_unit ? big : little).push_back(i);
        for (size_t i : big) run_block(round[i], found_aps[i], true);
        #pragma omp parallel for schedule(dynamic, 8) if (little.size() > std::min<size_t>(16, q_par_frontier))
        for (size_t j = 0; j < little.size(); ++j) run_block(round[little[j]], found_aps[little[j]], false);
        // Stage 2 for the next depth: blocks hanging off the APs just reached.
        // A 2-vertex block {a, w} inside an SCC means a<->w both exist, so
        // w's level is just a's + 1: settle it inline (most blocks of real
        // graphs are such pendants) instead of spawning a task.
        std::vector<Task> next;
        std::vector<int> work;
        for (auto& aps : found_aps) work.insert(work.end(), aps.begin(), aps.end());
        for (size_t wi = 0; wi < work.size(); ++wi) {
            int a = work[wi];
            blocks_of(a, tmp);
            for (int b : tmp) {
                if (seen_block[b]) continue;
                seen_block[b] = 1;
                const auto& bv = blk_verts[b];
                if (bv.size() == 2) {
                    int w = bv[0] == a ? bv[1] : bv[0];
                    if (dist[a] + 1 < dist[w]) dist[w] = dist[a] + 1;
                    if (v_bcnt[w] >= 2) work.push_back(w);
                } else {
                    next.push_back({b, a});
                }
            }
        }
        round.swap(next);
    }
}

QueryVec<int> GIR::bfsTwoLevel(int source) {
    sync();
    ensureBlocks();
    QueryVec<int> dist(V);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < V; ++i) dist[i] = INT_MAX;
    if (source < 0 || source >= V) { std::fill(dist.begin(), dist.end(), -1); return dist; }
    dist[source] = 0;
    dagWavefront(source,
                 [&](int c, bool is_src) {
                     if (is_src) bfsSourceSccByBlocks(source, dist);
                     else bfsInsideScc(c, dist);
                 },
                 [&](int x, int w) { atomicMin(&dist[w], dist[x] + 1); });
    for (auto& d : dist) if (d == INT_MAX) d = -1;
    return dist;
}

// ============================================================================
// SSSP
// ============================================================================

// Delta-stepping restricted to `scope` (all vertices if null), starting from
// the tentative distances already in dist. Buckets of width delta; a bucket
// is re-processed until nothing new lands in it.
template <class InScope>
void GIR::deltaStepping(QueryVec<uint64_t>& dist, std::vector<int> seeds, uint64_t delta, InScope in_scope,
                        bool par) {
    // Buckets of width delta; a bucket is re-processed until nothing new
    // lands in it. Each thread keeps its own buckets (no serial merge step);
    // a round gathers bucket k from all threads. Duplicate frontier entries
    // are skipped with a per-round marker instead of sorting the frontier.
    const int nthreads = par ? omp_get_max_threads() : 1;
    std::vector<std::vector<std::vector<int>>> tb(nthreads);
    std::vector<size_t> tmax(nthreads, 0);
    auto push = [&](int t, int x, uint64_t d) {
        size_t k = (size_t)(d / delta);
        auto& B = tb[t];
        if (B.size() <= k) B.resize(k + 1);
        B[k].push_back(x);
        if (k > tmax[t]) tmax[t] = k;
    };
    for (int x : seeds) push(0, x, dist[x]);
    if (par && round_mark.size() < (size_t)V) round_mark.assign(V, 0);
    std::vector<int> fr;
    std::vector<std::vector<int>> cur(nthreads);  // bucket k's lists, swapped out of tb
    std::vector<size_t> off(nthreads + 1);
    for (size_t k = 0;; ++k) {
        size_t kmax = 0;
        for (int t = 0; t < nthreads; ++t) kmax = std::max(kmax, tmax[t]);
        if (k > kmax) break;
        while (true) {
            // Bucket k's per-thread lists are swapped out (O(1)) and read in
            // place; pushes during the round go to fresh lists in tb. (A
            // gathering copy cost ~12 ms of a ~95 ms livejournal SSSP.)
            off[0] = 0;
            for (int t = 0; t < nthreads; ++t) {
                cur[t].clear();
                if (tb[t].size() > k) cur[t].swap(tb[t][k]);
                off[t + 1] = off[t] + cur[t].size();
            }
            if (off[nthreads] == 0) break;
            if (!par) fr.swap(cur[0]);
            // The shared marker is only used in data-parallel mode, where a
            // single call runs at a time; concurrent task-mode calls (small
            // SCCs in the two-level driver) de-duplicate locally instead.
            int round = 0;
            if (par) round = ++round_epoch;
            else { std::sort(fr.begin(), fr.end()); fr.erase(std::unique(fr.begin(), fr.end()), fr.end()); }
            const size_t m = par ? off[nthreads] : fr.size();
            #pragma omp parallel for schedule(dynamic, 256) if (par && m >= q_par_frontier)
            for (size_t i = 0; i < m; ++i) {
                int u;
                if (par) { const size_t s = std::upper_bound(off.begin(), off.end(), i) - off.begin() - 1; u = cur[s][i - off[s]]; }
                else u = fr[i];
                if (par && __atomic_exchange_n(&round_mark[u], round, __ATOMIC_RELAXED) == round) continue;
                uint64_t du = __atomic_load_n(&dist[u], __ATOMIC_RELAXED);
                if (du / delta != k) continue;  // stale entry
                const int t = par ? omp_get_thread_num() : 0;
                forOut(u, [&](int w) {
                    if (!in_scope(w)) return;
                    uint64_t nd = du + edgeWeight(u, w);
                    if (atomicMin(&dist[w], nd)) push(t, w, nd);
                });
            }
        }
        // (No freeing per bucket: the swap above hands each round's emptied
        // list back to tb, so list memory is reused rather than re-faulted.)
    }
}

QueryVec<uint64_t> GIR::sssp(int source, uint64_t delta) {
    QueryVec<uint64_t> dist(V);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < V; ++i) dist[i] = UINT64_MAX;
    if (source < 0 || source >= V) return dist;
    dist[source] = 0;
    deltaStepping(dist, {source}, delta, [](int) { return true; }, true);
    return dist;
}

QueryVec<uint64_t> GIR::ssspTwoLevel(int source, uint64_t delta) {
    sync();
    QueryVec<uint64_t> dist(V);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < V; ++i) dist[i] = UINT64_MAX;
    if (source < 0 || source >= V) return dist;
    dist[source] = 0;
    dagWavefront(source,
                 [&](int c, bool) {
                     if (sccs[c].size() == 1) return;  // settled by its predecessors already
                     std::vector<int> seeds;
                     for (int x : sccs[c]) if (dist[x] != UINT64_MAX) seeds.push_back(x);
                     deltaStepping(dist, seeds, delta, [&](int w) { return node_scc[w] == c; },
                                   sccs[c].size() >= q_large_unit);
                 },
                 [&](int x, int w) { atomicMin(&dist[w], dist[x] + edgeWeight(x, w)); });
    return dist;
}

// ============================================================================
// Betweenness centrality (Brandes, unweighted), sources as parallel tasks
// ============================================================================

QueryVec<double> GIR::bc(const std::vector<int>& sources) {
    QueryVec<double> score(V);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < V; ++i) score[i] = 0.0;
    // Per-vertex state packed in one record: the forward and backward passes
    // touch all three fields of each neighbour, so one cache line instead of
    // three per visit.
    struct St { int dist; double sigma, delta; };
    #pragma omp parallel
    {
        // Per-thread state kept across calls (every source resets what it
        // touched): a fresh V-sized array per thread and call cost 64 x 116 MB
        // of zero filling on livejournal, about a third of a 64-source run.
        static thread_local std::vector<St> st;
        if (st.size() != (size_t)V) st.assign(V, St{-1, 0.0, 0.0});
        static thread_local std::vector<int> order;
        order.reserve(1024);
        #pragma omp for schedule(dynamic, 1)
        for (size_t si = 0; si < sources.size(); ++si) {
            int s = sources[si];
            if (s < 0 || s >= V) continue;
            order.clear();
            st[s].dist = 0;
            st[s].sigma = 1;
            order.push_back(s);
            for (size_t i = 0; i < order.size(); ++i) {
                int u = order[i];
                const int du = st[u].dist;
                const double su = st[u].sigma;
                forOut(u, [&](int w) {
                    St& x = st[w];
                    if (x.dist < 0) { x.dist = du + 1; order.push_back(w); }
                    if (x.dist == du + 1) x.sigma += su;
                });
            }
            for (size_t i = order.size(); i-- > 0;) {
                int w = order[i];
                const St& y = st[w];
                const double coef = (1.0 + y.delta) / y.sigma;
                const int dw = y.dist;
                forIn(w, [&](int u) {
                    St& x = st[u];
                    if (x.dist >= 0 && x.dist == dw - 1) x.delta += x.sigma * coef;
                });
                if (w != s) atomicAddDouble(&score[w], y.delta);
            }
            for (int x : order) st[x] = St{-1, 0.0, 0.0};
        }
    }
    return score;
}

// ============================================================================
// Structural queries
// ============================================================================

// Every condensation edge goes from a lower to a higher rank, so a path
// SCC(u) ~> SCC(v) only visits SCCs with ranks in [rank(SCC(u)), rank(SCC(v))].
// Bidirectional search in that interval; each step expands the side whose
// next node has the smaller condensation degree (the giant SCC's adjacency
// has hundreds of thousands of entries).
bool GIR::reachable(int u, int v) {
    sync();
    if (u == v) return true;
    const int a = node_scc[u], b = node_scc[v];
    if (a == b) return true;
    const long long ra = scc_rank[a], rb = scc_rank[b];
    if (ra > rb) return false;
    if (pk_sf.size() < sccs.size()) { pk_sf.resize(sccs.size(), 0); pk_sb.resize(sccs.size(), 0); }
    if (++pk_ef == INT_MAX) { std::fill(pk_sf.begin(), pk_sf.end(), 0); pk_ef = 1; }
    if (++pk_eb == INT_MAX) { std::fill(pk_sb.begin(), pk_sb.end(), 0); pk_eb = 1; }
    const int ef = pk_ef, eb = pk_eb;
    std::vector<int> F{a}, B{b};
    pk_sf[a] = ef;
    pk_sb[b] = eb;
    size_t fi = 0, bi = 0;
    long long wf = 0, wb = 0;
    while (fi < F.size() && bi < B.size()) {
        const long long cf = (long long)dag_out[F[fi]].size(), cb = (long long)dag_in[B[bi]].size();
        if (wf + cf <= wb + cb) {
            const int x = F[fi++];
            wf += cf + 1;
            for (auto& kv : dag_out[x]) {
                const int y = kv.first;
                if (scc_rank[y] > rb) continue;
                if (pk_sb[y] == eb) return true;
                if (pk_sf[y] != ef) { pk_sf[y] = ef; F.push_back(y); }
            }
        } else {
            const int x = B[bi++];
            wb += cb + 1;
            for (auto& kv : dag_in[x]) {
                const int y = kv.first;
                if (scc_rank[y] < ra) continue;
                if (pk_sf[y] == ef) return true;
                if (pk_sb[y] != eb) { pk_sb[y] = eb; B.push_back(y); }
            }
        }
    }
    return false;
}

bool GIR::reachableNoIndex(int u, int v) {
    if (u == v) return true;
    const int ea = nextEpochA(), eb = nextEpochB();
    stamp_a[u] = ea;
    stamp_b[v] = eb;
    std::vector<int> qa{u}, qb{v};
    size_t ia = 0, ib = 0;
    long long wa = 0, wb = 0;
    bool met = false;
    auto odeg = [&](int x) { return (long long)(g_directed->fw_beg_pos[x + 1] - g_directed->fw_beg_pos[x]) + (long long)extra_out_adj[x].size(); };
    auto ideg = [&](int x) { return (long long)(g_directed->bw_beg_pos[x + 1] - g_directed->bw_beg_pos[x]) + (long long)extra_in_adj[x].size(); };
    while (!met && ia < qa.size() && ib < qb.size()) {
        const long long ca = odeg(qa[ia]), cb = ideg(qb[ib]);
        if (wa + ca <= wb + cb) {
            wa += ca + 1;
            anyOut(qa[ia++], [&](int w) {
                if (stamp_b[w] == eb) { met = true; return true; }
                if (stamp_a[w] != ea) { stamp_a[w] = ea; qa.push_back(w); }
                return false;
            });
        } else {
            wb += cb + 1;
            anyIn(qb[ib++], [&](int w) {
                if (stamp_a[w] == ea) { met = true; return true; }
                if (stamp_b[w] != eb) { stamp_b[w] = eb; qb.push_back(w); }
                return false;
            });
        }
    }
    return met;
}


// ============================================================================
// Incremental single-source BFS / SSSP (KickStarter-style dependence trees)
// ============================================================================

void GIR::trackIncremental(int source, bool weighted) {
    auto& S = inc[weighted];
    S.pending.clear();
    S.source = source;
    S.st.resize(V);
    auto w = [&](int a, int b) -> uint64_t { return weighted ? edgeWeight(a, b) : 1; };
    QueryVec<uint64_t> d(V);
    if (weighted) d = sssp(source);
    else {
        QueryVec<int> lv = bfs(source);
        #pragma omp parallel for schedule(static)
        for (int v = 0; v < V; ++v) d[v] = lv[v] < 0 ? UINT64_MAX : (uint64_t)lv[v];
    }
    // Parent: any in-neighbour on a shortest path.
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int v = 0; v < V; ++v) {
        if (d[v] == UINT64_MAX) { S.st[v] = ~0ULL; continue; }
        if (v == source) { S.st[v] = (uint64_t)(unsigned)source; continue; }
        int par = -1;
        anyIn(v, [&](int p) {
            if (d[p] != UINT64_MAX && d[p] + w(p, v) == d[v]) { par = p; return true; }
            return false;
        });
        S.st[v] = (d[v] << 32) | (unsigned)par;
    }
}

void GIR::updateIncremental(bool weighted) {
    auto& S = inc[weighted];
    if (S.source < 0 || S.pending.empty()) return;
    std::vector<std::pair<int, int>> P;
    P.swap(S.pending);
    if (P.size() >= 4096) __gnu_parallel::sort(P.begin(), P.end()); else std::sort(P.begin(), P.end());
    P.erase(std::unique(P.begin(), P.end()), P.end());
    auto& st = S.st;
    const int src = S.source;
    const uint64_t NONE = ~0ULL;
    auto w = [&](int a, int b) -> uint64_t { return weighted ? edgeWeight(a, b) : 1; };
    auto parentOf = [&](uint64_t k) { return (int)(unsigned)(k & 0xffffffffULL); };
    auto writeMin = [&](int v, uint64_t k) {
        uint64_t cur = __atomic_load_n(&st[v], __ATOMIC_RELAXED);
        while (k < cur)
            if (__atomic_compare_exchange_n(&st[v], &cur, k, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return true;
        return false;
    };
    if (inc_flag.size() < (size_t)V) inc_flag.assign(V, 0);
    char* flag = inc_flag.data();  // 1 = invalid, 2 = queued
    auto setFlag = [&](int v, char f) {
        char exp = 0;
        return __atomic_compare_exchange_n(&flag[v], &exp, f, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED);
    };
    const int T = omp_get_max_threads();
    std::vector<std::vector<int>> parts(T);
    auto gather = [&](std::vector<int>& out) {
        out.clear();
        for (auto& p : parts) { out.insert(out.end(), p.begin(), p.end()); p.clear(); }
    };
    std::vector<char> present(P.size());
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < P.size(); ++i) present[i] = isDirectedEdgePresent(P[i].first, P[i].second);

    // 1-2. Invalidate the dependence subtrees below deleted tree edges.
    std::vector<int> inv, fr, next;
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < P.size(); ++i) {
        const int u = P[i].first, v = P[i].second;
        if (present[i] || v == src || st[v] == NONE || parentOf(st[v]) != u) continue;
        if (setFlag(v, 1)) parts[omp_get_thread_num()].push_back(v);
    }
    gather(fr);
    while (!fr.empty()) {
        inv.insert(inv.end(), fr.begin(), fr.end());
        #pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < fr.size(); ++i) {
            const int x = fr[i];
            forOut(x, [&](int y) {
                const uint64_t k = st[y];
                if (y != src && k != NONE && parentOf(k) == x && setFlag(y, 1)) parts[omp_get_thread_num()].push_back(y);
            });
        }
        gather(next);
        fr.swap(next);
    }
    // 3-4. Invalid vertices restart from their valid in-neighbours.
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < inv.size(); ++i) st[inv[i]] = NONE;
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < inv.size(); ++i) {
        const int x = inv[i];
        uint64_t best = NONE;
        forIn(x, [&](int p) {
            if (flag[p] == 1) return;
            const uint64_t k = st[p];
            if (k == NONE) return;
            best = std::min(best, (((k >> 32) + w(p, x)) << 32) | (unsigned)p);
        });
        st[x] = best;
    }
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < inv.size(); ++i) flag[inv[i]] = 0;
    // Seeds: re-seeded invalid vertices, then targets improved by insertions.
    for (int x : inv) if (st[x] != NONE && setFlag(x, 2)) parts[0].push_back(x);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < P.size(); ++i) {
        if (!present[i]) continue;
        const int u = P[i].first, v = P[i].second;
        const uint64_t k = st[u];
        if (k == NONE) continue;
        if (writeMin(v, (((k >> 32) + w(u, v)) << 32) | (unsigned)u) && setFlag(v, 2)) parts[omp_get_thread_num()].push_back(v);
    }
    gather(fr);
    // 5. Propagate (frontier Bellman-Ford; a vertex improved again is requeued).
    while (!fr.empty()) {
        #pragma omp parallel for schedule(dynamic, 64)
        for (size_t i = 0; i < fr.size(); ++i) {
            const int x = fr[i];
            // Clear the queued mark, then read x's distance. A full fence
            // (exchange) between the two: a relaxed store could be passed by
            // the load, so a concurrent improvement of x would neither be read
            // here nor requeue x (its CAS on the flag still sees "queued").
            __atomic_exchange_n(&flag[x], (char)0, __ATOMIC_SEQ_CST);
            const uint64_t d = __atomic_load_n(&st[x], __ATOMIC_SEQ_CST) >> 32;
            forOut(x, [&](int y) {
                if (writeMin(y, ((d + w(x, y)) << 32) | (unsigned)x) && setFlag(y, 2)) parts[omp_get_thread_num()].push_back(y);
            });
        }
        gather(next);
        fr.swap(next);
    }
}
