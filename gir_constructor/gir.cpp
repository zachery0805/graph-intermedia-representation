// GIR: three-layer structural representation (WCC / SCC + condensation DAG /
// per-SCC block-cut forest) with incremental maintenance. See gir.h for the
// invariants every routine below relies on.
#include "gir.h"
#include <functional>

#include <algorithm>
#include <parallel/algorithm>
#include <chrono>
#include <climits>
#include <iostream>
#include <memory>
#include <queue>
#include <unordered_map>
#include <unordered_set>

#include <omp.h>

// jemalloc (linked when installed, see Makefile) reads this. Its default
// serves allocations over 8 MB from an arena that returns their pages at once,
// so every query's per-vertex arrays page-faulted again (BFS 31 -> 47 ms on
// livejournal); with huge pages for everything, the large per-vertex arrays
// also stay TLB-friendly under the random accesses of batched updates.
// Measured best of the settings tried. Ignored when jemalloc is not linked.
extern "C" {
const char* malloc_conf = "oversize_threshold:4294967296,thp:always,metadata_thp:always";
}

#include "block_cut_tree.h"
#include "graph.h"
#include "graph_undirected.h"
#include "gir_iter.h"

// ============================================================================
// Construction / destruction
// ============================================================================

GIR::GIR(const std::string& path_directed, const std::string& path_undirected)
    : base_path_directed(path_directed), base_path_undirected(path_undirected) {
    std::string d = base_path_directed;
    g_directed = new graph((d + "fw_begin.bin").c_str(), (d + "fw_adjacent.bin").c_str(),
                           (d + "bw_begin.bin").c_str(), (d + "bw_adjacent.bin").c_str());
    initAfterLoad();
}

GIR::GIR(graph* g) : g_directed(g) { initAfterLoad(); }

void GIR::initAfterLoad() {
    V = (int)g_directed->vert_count;
    normalizeBaseCsr();
    scc_pos.assign(V, 0);
    blk_pos.assign(V, 0);
    extra_out_adj.assign(V, {});
    extra_in_adj.assign(V, {});
    extra_und_adj.assign(V, {});
    stamp_a.assign(V, 0);
    stamp_b.assign(V, 0);
    scr_idx.assign(V, -1);
    scr_low.assign(V, 0);
    scr_flag.assign(V, 0);
    par_a.assign(V, -1);
    par_b.assign(V, -1);
}

GIR::~GIR() {
    delete g_directed;
    delete g_undirected;
}

graph_undirected* GIR::getUndirectedGraph() const {
    if (!g_undirected) {
        std::string p = base_path_undirected;
        g_undirected = new graph_undirected((p + "undirected_begin.bin").c_str(),
                                            (p + "undirected_adjacent.bin").c_str());
    }
    return g_undirected;
}

int GIR::nextEpochA() const {
    if (++epoch_a == INT_MAX) { std::fill(stamp_a.begin(), stamp_a.end(), 0); epoch_a = 1; }
    return epoch_a;
}
int GIR::nextEpochB() const {
    if (++epoch_b == INT_MAX) { std::fill(stamp_b.begin(), stamp_b.end(), 0); epoch_b = 1; }
    return epoch_b;
}

// ============================================================================
// Edge store
// ============================================================================

// Sort every adjacency list and drop duplicate edges (fw and bw alike). Duplicate
// base edges would otherwise be counted twice in the condensation-edge
// reference counts but removed once by a deletion.
static void normalizeOneCsr(index_t*& beg, vertex_t*& csr, index_t V, index_t& E) {
    std::vector<index_t> deg(V);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (index_t v = 0; v < V; ++v) {
        std::sort(csr + beg[v], csr + beg[v + 1]);
        deg[v] = (index_t)(std::unique(csr + beg[v], csr + beg[v + 1]) - (csr + beg[v]));
    }
    long long total = 0;
    for (index_t v = 0; v < V; ++v) total += deg[v];
    if (total == (long long)E) return;  // no duplicates: sorted in place
    index_t* nb = new index_t[V + 1];
    nb[0] = 0;
    for (index_t v = 0; v < V; ++v) nb[v + 1] = nb[v] + deg[v];
    vertex_t* nc = new vertex_t[total > 0 ? total : 1];
    #pragma omp parallel for schedule(dynamic, 1024)
    for (index_t v = 0; v < V; ++v) std::copy(csr + beg[v], csr + beg[v] + deg[v], nc + nb[v]);
    delete[] beg;
    delete[] csr;
    beg = nb;
    csr = nc;
    E = (index_t)total;
}

void GIR::normalizeBaseCsr() {
    index_t e_fw = g_directed->edge_count, e_bw = g_directed->edge_count;
    normalizeOneCsr(g_directed->fw_beg_pos, g_directed->fw_csr, g_directed->vert_count, e_fw);
    normalizeOneCsr(g_directed->bw_beg_pos, g_directed->bw_csr, g_directed->vert_count, e_bw);
    if (e_fw != e_bw) {
        std::cerr << "GIR: forward and backward CSR disagree after deduplication (" << e_fw << " vs " << e_bw
                  << " edges); the bw files are not the transpose of fw" << std::endl;
        exit(1);
    }
    if (e_fw != g_directed->edge_count)
        std::cout << "GIR: removed " << (g_directed->edge_count - e_fw) << " duplicate edges" << std::endl;
    g_directed->edge_count = e_fw;
}

long long GIR::fwPos(int u, int v) const {
    const vertex_t *b = g_directed->fw_csr + g_directed->fw_beg_pos[u], *e = g_directed->fw_csr + g_directed->fw_beg_pos[u + 1];
    const vertex_t* p = std::lower_bound(b, e, (vertex_t)v);
    return (p != e && *p == (vertex_t)v) ? (long long)(p - g_directed->fw_csr) : -1;
}
long long GIR::bwPos(int u, int v) const {
    const vertex_t *b = g_directed->bw_csr + g_directed->bw_beg_pos[v], *e = g_directed->bw_csr + g_directed->bw_beg_pos[v + 1];
    const vertex_t* p = std::lower_bound(b, e, (vertex_t)u);
    return (p != e && *p == (vertex_t)u) ? (long long)(p - g_directed->bw_csr) : -1;
}
long long GIR::symPos(int u, int v) const {
    if (sym_beg.empty()) return -1;
    auto b = sym_adj.begin() + sym_beg[u], e = sym_adj.begin() + sym_beg[u + 1];
    auto p = std::lower_bound(b, e, v);
    return (p != e && *p == v) ? (long long)(p - sym_adj.begin()) : -1;
}

bool GIR::baseHasDirected(int u, int v) const { return fwPos(u, v) >= 0; }
bool GIR::baseHasUndirected(int u, int v) const { return symPos(u, v) >= 0; }

bool GIR::isDirectedEdgePresent(int u, int v) const {
    if (u < 0 || v < 0 || u >= V || v >= V) return false;
    long long p = fwPos(u, v);
    if (p >= 0) return !n_dead_dir || !bitGet(dead_fw, p);
    const auto& L = extra_out_adj[u];
    return std::binary_search(L.begin(), L.end(), v);
}

bool GIR::isUndirectedEdgePresent(int u, int v) const {
    if (u < 0 || v < 0 || u >= V || v >= V) return false;
    long long p = symPos(u, v);
    if (p >= 0) return !n_dead_und || !bitGet(dead_sym, p);
    const auto& L = extra_und_adj[u];
    return std::binary_search(L.begin(), L.end(), v);
}

static void sortedInsert(std::vector<int>& vec, int val) {
    auto pos = std::lower_bound(vec.begin(), vec.end(), val);
    if (pos == vec.end() || *pos != val) vec.insert(pos, val);
}
static void sortedErase(std::vector<int>& vec, int val) {
    auto pos = std::lower_bound(vec.begin(), vec.end(), val);
    if (pos != vec.end() && *pos == val) vec.erase(pos);
}
// Sort a batch's pairs (parallel when large).
static void psort(std::vector<std::pair<int, int>>& l) {
    if (std::is_sorted(l.begin(), l.end())) return;
    if (l.size() >= 4096) __gnu_parallel::sort(l.begin(), l.end());
    else std::sort(l.begin(), l.end());
}
// Append per-thread buffers to `out` in thread order (keeps a static
// schedule's order).
template <class T>
static void concatParts(std::vector<std::vector<T>>& parts, std::vector<T>& out) {
    size_t tot = out.size();
    std::vector<size_t> off(parts.size());
    for (size_t t = 0; t < parts.size(); ++t) { off[t] = tot; tot += parts[t].size(); }
    out.resize(tot);
    #pragma omp parallel for schedule(static, 1)
    for (size_t t = 0; t < parts.size(); ++t) std::copy(parts[t].begin(), parts[t].end(), out.begin() + off[t]);
}

// The delta-list entries of a batch (edges that are not base edges, and
// undirected edges that are not base edges), built in parallel in order.
// Parallel stable compaction: for i in [0, n), emit(i, put) may call put(x)
// up to K times; outputs keep the order of i. Two passes over fixed chunks
// (count, then write at prefix offsets), so no per-thread buffers or merging.
template <int K, class T, class Emit>
static void parCompact(size_t n, std::vector<T>& out, Emit emit) {
    const int T_ = omp_get_max_threads();
    const size_t CH = std::max<size_t>(1024, (n + 4 * T_ - 1) / (4 * T_));
    const size_t nch = (n + CH - 1) / CH;
    std::vector<size_t> cnt(nch + 1, 0);
    #pragma omp parallel for schedule(static) if (n >= 8192)
    for (size_t ch = 0; ch < nch; ++ch) {
        size_t k = 0;
        for (size_t i = ch * CH, e = std::min(n, i + CH); i < e; ++i) emit(i, [&](const T&) { ++k; });
        cnt[ch + 1] = k;
    }
    for (size_t ch = 0; ch < nch; ++ch) cnt[ch + 1] += cnt[ch];
    const size_t base = out.size();
    out.resize(base + cnt[nch]);
    #pragma omp parallel for schedule(static) if (n >= 8192)
    for (size_t ch = 0; ch < nch; ++ch) {
        size_t k = base + cnt[ch];
        for (size_t i = ch * CH, e = std::min(n, i + CH); i < e; ++i) emit(i, [&](const T& x) { out[k++] = x; });
    }
}

static void buildDeltaLists(const std::vector<std::pair<int, int>>& es, const std::vector<long long>& pf,
                            const std::vector<char>& und, const std::vector<long long>& pu,
                            std::vector<std::pair<int, int>>& out_l, std::vector<std::pair<int, int>>& in_l,
                            std::vector<std::pair<int, int>>& und_l) {
    using E = std::pair<int, int>;
    out_l.clear(); in_l.clear(); und_l.clear();
    parCompact<1, E>(es.size(), out_l, [&](size_t i, auto put) { if (pf[i] < 0) put(es[i]); });
    parCompact<1, E>(es.size(), in_l, [&](size_t i, auto put) { if (pf[i] < 0) put(E{es[i].second, es[i].first}); });
    parCompact<2, E>(es.size(), und_l, [&](size_t i, auto put) {
        if (und[i] && pu[i] < 0) { put(es[i]); put(E{es[i].second, es[i].first}); }
    });
}

// Apply a sorted list of (vertex, neighbour) pairs to the sorted per-vertex
// delta lists: merged in (add) or subtracted (remove), one thread per vertex.
static void applyGrouped(std::vector<std::pair<int, int>>& l, std::vector<std::vector<int>>& adj, bool add) {
    psort(l);
    std::vector<size_t> starts;
    parCompact<1, size_t>(l.size(), starts, [&](size_t i, auto put) { if (i == 0 || l[i].first != l[i - 1].first) put(i); });
    starts.push_back(l.size());
    const size_t ng = starts.size() - 1;
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t g = 0; g < ng; ++g) {
        auto& a = adj[l[starts[g]].first];
        const size_t b = starts[g], e = starts[g + 1];
        if (add) {
            if (e - b <= 8) {
                // Few new entries (the usual case): insert in place.
                // (inplace_merge would allocate a buffer per vertex, and
                // 64 threads then contend on the allocator.)
                for (size_t i = b; i < e; ++i) sortedInsert(a, l[i].second);
            } else {
                const size_t m = a.size();
                for (size_t i = b; i < e; ++i) a.push_back(l[i].second);
                std::inplace_merge(a.begin(), a.begin() + m, a.end());
            }
        } else {
            size_t k = 0, i = b;
            for (size_t j = 0; j < a.size(); ++j) {
                while (i < e && l[i].second < a[j]) ++i;
                if (i < e && l[i].second == a[j]) continue;
                a[k++] = a[j];
            }
            a.resize(k);
        }
    }
}

void GIR::overlayAddDirected(int u, int v) {
    long long p = fwPos(u, v);
    if (p >= 0) {
        if (!dead_fw.empty() && bitGet(dead_fw, p)) {
            bitPut(dead_fw, p, false);
            bitPut(dead_bw, bwPos(u, v), false);
            --n_dead_dir;
        }
        return;
    }
    sortedInsert(extra_out_adj[u], v);
    sortedInsert(extra_in_adj[v], u);
    ++extra_count;
}
void GIR::overlayRemoveDirected(int u, int v) {
    long long p = fwPos(u, v);
    if (p >= 0) {
        if (dead_fw.empty()) {
            dead_fw.assign(((size_t)g_directed->fw_beg_pos[V] + 63) / 64, 0);
            dead_bw.assign(((size_t)g_directed->bw_beg_pos[V] + 63) / 64, 0);
        }
        if (!bitGet(dead_fw, p)) {
            bitPut(dead_fw, p, true);
            bitPut(dead_bw, bwPos(u, v), true);
            ++n_dead_dir;
        }
        return;
    }
    sortedErase(extra_out_adj[u], v);
    sortedErase(extra_in_adj[v], u);
    --extra_count;
}
void GIR::overlayAddUndirected(int u, int v) {
    long long p = symPos(u, v);
    if (p >= 0) {
        if (!dead_sym.empty() && bitGet(dead_sym, p)) {
            bitPut(dead_sym, p, false);
            bitPut(dead_sym, symPos(v, u), false);
            --n_dead_und;
        }
        return;
    }
    sortedInsert(extra_und_adj[u], v);
    sortedInsert(extra_und_adj[v], u);
}
void GIR::overlayRemoveUndirected(int u, int v) {
    long long p = symPos(u, v);
    if (p >= 0) {
        if (dead_sym.empty()) dead_sym.assign((sym_adj.size() + 63) / 64, 0);
        if (!bitGet(dead_sym, p)) {
            bitPut(dead_sym, p, true);
            bitPut(dead_sym, symPos(v, u), true);
            ++n_dead_und;
        }
        return;
    }
    sortedErase(extra_und_adj[u], v);
    sortedErase(extra_und_adj[v], u);
}

// The undirected edge {u,v} exists iff u->v or v->u exists. These two are the
// only places that change the edge set, so that invariant can't drift.
bool GIR::storeInsert(int u, int v, bool& und_changed) {
    und_changed = false;
    if (u == v || isDirectedEdgePresent(u, v)) return false;
    bool und_present = isUndirectedEdgePresent(u, v);
    overlayAddDirected(u, v);
    if (!und_present) { overlayAddUndirected(u, v); und_changed = true; }
    return true;
}
bool GIR::storeDelete(int u, int v, bool& und_changed) {
    und_changed = false;
    if (u == v || !isDirectedEdgePresent(u, v)) return false;
    overlayRemoveDirected(u, v);
    if (!isDirectedEdgePresent(v, u)) { overlayRemoveUndirected(u, v); und_changed = true; }
    return true;
}

void GIR::storeDeleteBulk(const std::vector<std::pair<int, int>>& es) {
    const size_t n = es.size();
    if (n == 0) return;
    if (dead_fw.empty()) {
        dead_fw.assign(((size_t)g_directed->fw_beg_pos[V] + 63) / 64, 0);
        dead_bw.assign(((size_t)g_directed->bw_beg_pos[V] + 63) / 64, 0);
    }
    if (dead_sym.empty() && !sym_adj.empty()) dead_sym.assign((sym_adj.size() + 63) / 64, 0);
    // (A batch from normalizeBatch is already sorted: search it directly.)
    const bool sorted = std::is_sorted(es.begin(), es.end());
    std::vector<long long> dk;
    if (!sorted) {
        dk.resize(n);
        for (size_t i = 0; i < n; ++i) dk[i] = keyD(es[i].first, es[i].second);
        if (dk.size() >= 4096) __gnu_parallel::sort(dk.begin(), dk.end()); else std::sort(dk.begin(), dk.end());
    }
    auto in_batch = [&](int a, int b) {
        return sorted ? std::binary_search(es.begin(), es.end(), std::make_pair(a, b))
                      : std::binary_search(dk.begin(), dk.end(), keyD(a, b));
    };
    // {u,v} disappears iff v->u is absent afterwards; if both directions go,
    // the smaller directed key removes it.
    std::vector<long long> pf(n), pu(n);
    std::vector<char> und_rm(n);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; ++i) {
        int u = es[i].first, v = es[i].second;
        pf[i] = fwPos(u, v);
        bool rev_in = in_batch(v, u);
        und_rm[i] = rev_in ? (u < v) : !isDirectedEdgePresent(v, u);
        pu[i] = und_rm[i] ? symPos(u, v) : -1;
    }
    auto set = [](std::vector<uint64_t>& m, long long i) {
        __atomic_fetch_or(&m[(size_t)i >> 6], 1ULL << (i & 63), __ATOMIC_RELAXED);
    };
    long long dd = 0, du = 0;
    #pragma omp parallel for schedule(static) reduction(+ : dd, du)
    for (size_t i = 0; i < n; ++i) {
        int u = es[i].first, v = es[i].second;
        if (pf[i] >= 0) { set(dead_fw, pf[i]); set(dead_bw, bwPos(u, v)); ++dd; }
        if (und_rm[i] && pu[i] >= 0) { set(dead_sym, pu[i]); set(dead_sym, symPos(v, u)); ++du; }
    }
    n_dead_dir += dd;
    n_dead_und += du;
    // Delta edges: erase from per-vertex lists, one thread per vertex group.
    std::vector<std::pair<int, int>> out_l, in_l, und_l;
    buildDeltaLists(es, pf, und_rm, pu, out_l, in_l, und_l);
    extra_count -= (long long)out_l.size();
    applyGrouped(out_l, extra_out_adj, false);
    applyGrouped(in_l, extra_in_adj, false);
    applyGrouped(und_l, extra_und_adj, false);
}

void GIR::storeInsertBulk(const std::vector<std::pair<int, int>>& es) {
    const size_t n = es.size();
    if (n == 0) return;
    // Classify against the current store (read-only, parallel).
    std::vector<long long> pf(n), pu(n);  // fw position (base) / sym position (base)
    std::vector<char> und_add(n);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; ++i) {
        int u = es[i].first, v = es[i].second;
        pf[i] = fwPos(u, v);
        und_add[i] = !isUndirectedEdgePresent(u, v);
        pu[i] = und_add[i] ? symPos(u, v) : -1;
    }
    // Both u->v and v->u may add {u,v}: keep one (the smaller directed key).
    {
        const bool sorted = std::is_sorted(es.begin(), es.end());  // (normalizeBatch output is)
        std::vector<long long> dk;
        if (!sorted) {
            dk.resize(n);
            for (size_t i = 0; i < n; ++i) dk[i] = keyD(es[i].first, es[i].second);
            if (dk.size() >= 4096) __gnu_parallel::sort(dk.begin(), dk.end()); else std::sort(dk.begin(), dk.end());
        }
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < n; ++i) {
            int u = es[i].first, v = es[i].second;
            if (!und_add[i] || u <= v) continue;
            bool rev = sorted ? std::binary_search(es.begin(), es.end(), std::make_pair(v, u))
                              : std::binary_search(dk.begin(), dk.end(), keyD(v, u));
            if (rev) und_add[i] = 0;
        }
    }
    // Base edges: clear tombstone bits (atomic: neighbouring bits share words).
    long long revived_d = 0, revived_u = 0;
    auto clr = [](std::vector<uint64_t>& m, long long i) {
        __atomic_fetch_and(&m[(size_t)i >> 6], ~(1ULL << (i & 63)), __ATOMIC_RELAXED);
    };
    #pragma omp parallel for schedule(static) reduction(+ : revived_d, revived_u)
    for (size_t i = 0; i < n; ++i) {
        int u = es[i].first, v = es[i].second;
        if (pf[i] >= 0) { clr(dead_fw, pf[i]); clr(dead_bw, bwPos(u, v)); ++revived_d; }
        if (und_add[i] && pu[i] >= 0) { clr(dead_sym, pu[i]); clr(dead_sym, symPos(v, u)); ++revived_u; }
    }
    n_dead_dir -= revived_d;
    n_dead_und -= revived_u;
    // Delta edges: append to per-vertex lists, one thread per vertex group.
    std::vector<std::pair<int, int>> out_l, in_l, und_l;
    buildDeltaLists(es, pf, und_add, pu, out_l, in_l, und_l);
    extra_count += (long long)out_l.size();
    applyGrouped(out_l, extra_out_adj, true);
    applyGrouped(in_l, extra_in_adj, true);
    applyGrouped(und_l, extra_und_adj, true);
}

GIR::Cursor GIR::outBegin(int v) const {
    return {(long long)g_directed->fw_beg_pos[v], (long long)g_directed->fw_beg_pos[v + 1], 0};
}
bool GIR::outNext(int v, Cursor& c, int& w) const {
    while (c.i < c.end) {
        w = (int)g_directed->fw_csr[c.i++];
        if (n_dead_dir && bitGet(dead_fw, c.i - 1)) continue;
        return true;
    }
    const auto& ex = extra_out_adj[v];
    if (c.j >= ex.size()) return false;
    w = ex[c.j++];
    return true;
}
GIR::Cursor GIR::undBegin(int v) const { return {sym_beg[v], sym_beg[v + 1], 0}; }
bool GIR::undNext(int v, Cursor& c, int& w) const {
    while (c.i < c.end) {
        w = sym_adj[c.i++];
        if (n_dead_und && bitGet(dead_sym, c.i - 1)) continue;
        return true;
    }
    const auto& ex = extra_und_adj[v];
    if (c.j >= ex.size()) return false;
    w = ex[c.j++];
    return true;
}

std::vector<int> GIR::getDirectedOutNeighbors(int v) const {
    std::vector<int> r; forOut(v, [&](int w) { r.push_back(w); }); return r;
}
std::vector<int> GIR::getDirectedInNeighbors(int v) const {
    std::vector<int> r; forIn(v, [&](int w) { r.push_back(w); }); return r;
}
std::vector<int> GIR::getUndirectedNeighborsPublic(int v) const {
    std::vector<int> r; forUnd(v, [&](int w) { r.push_back(w); }); return r;
}

// ============================================================================
// Generic iterative traversals
// ============================================================================

template <class InScope, class Emit>
void GIR::tarjanScc(const std::vector<int>& roots, InScope in_scope, Emit emit) {
    struct Frame { int v; Cursor c; };
    std::vector<int> touched, stk;
    std::vector<Frame> cs;
    int idx = 0;
    for (int r : roots) {
        if (scr_idx[r] != -1) continue;
        scr_idx[r] = scr_low[r] = idx++;
        touched.push_back(r);
        stk.push_back(r);
        scr_flag[r] = 1;
        cs.push_back({r, outBegin(r)});
        while (!cs.empty()) {
            int v = cs.back().v, w;
            if (outNext(v, cs.back().c, w)) {
                if (!in_scope(w)) continue;
                if (scr_idx[w] == -1) {
                    scr_idx[w] = scr_low[w] = idx++;
                    touched.push_back(w);
                    stk.push_back(w);
                    scr_flag[w] = 1;
                    cs.push_back({w, outBegin(w)});
                } else if (scr_flag[w]) {
                    scr_low[v] = std::min(scr_low[v], scr_idx[w]);
                }
            } else {
                cs.pop_back();
                if (!cs.empty()) scr_low[cs.back().v] = std::min(scr_low[cs.back().v], scr_low[v]);
                if (scr_low[v] == scr_idx[v]) {
                    std::vector<int> comp;
                    while (true) {
                        int x = stk.back(); stk.pop_back();
                        scr_flag[x] = 0;
                        comp.push_back(x);
                        if (x == v) break;
                    }
                    emit(comp);
                }
            }
        }
    }
    for (int x : touched) { scr_idx[x] = -1; scr_flag[x] = 0; }
}

// Vertex-stack variant of Hopcroft-Tarjan: when low[w] >= disc[v] for a tree
// edge v->w, the block is v plus everything above w on the stack.
template <class KeepEdge, class Emit>
void GIR::hopcroftTarjan(const std::vector<int>& roots, KeepEdge keep_edge, Emit emit) {
    struct Frame { int v, p; Cursor c; };
    std::vector<int> touched, stk;
    std::vector<Frame> cs;
    int t = 0;
    for (int r : roots) {
        if (scr_idx[r] != -1) continue;
        scr_idx[r] = scr_low[r] = t++;
        touched.push_back(r);
        stk.push_back(r);
        cs.push_back({r, -1, undBegin(r)});
        bool root_has_child = false;
        while (!cs.empty()) {
            int v = cs.back().v, p = cs.back().p, w;
            if (undNext(v, cs.back().c, w)) {
                if (w == v || w == p || !keep_edge(v, w)) continue;
                if (scr_idx[w] == -1) {
                    scr_idx[w] = scr_low[w] = t++;
                    touched.push_back(w);
                    stk.push_back(w);
                    if (v == r) root_has_child = true;
                    cs.push_back({w, v, undBegin(w)});
                } else {
                    scr_low[v] = std::min(scr_low[v], scr_idx[w]);
                }
            } else {
                cs.pop_back();
                if (p == -1) continue;
                scr_low[p] = std::min(scr_low[p], scr_low[v]);
                if (scr_low[v] >= scr_idx[p]) {
                    std::vector<int> block;
                    while (true) {
                        int x = stk.back(); stk.pop_back();
                        block.push_back(x);
                        if (x == v) break;
                    }
                    block.push_back(p);
                    emit(block);
                }
            }
        }
        stk.pop_back();  // r
        if (!root_has_child) { std::vector<int> single{r}; emit(single); }
    }
    for (int x : touched) scr_idx[x] = -1;
}

// ============================================================================
// WCC layer
// ============================================================================

void GIR::buildWCC() {
    std::cout << "=== Build WCC ===" << std::endl;

    // Symmetrized base adjacency (sorted, deduplicated, no self-loops). Each
    // vertex's range is just its fw list followed by its bw list, so this
    // is embarrassingly parallel.
    {
        std::vector<long long> deg(V + 1, 0);
        #pragma omp parallel for schedule(static)
        for (int u = 0; u < V; ++u)
            deg[u] = (long long)(g_directed->fw_beg_pos[u + 1] - g_directed->fw_beg_pos[u]) +
                     (long long)(g_directed->bw_beg_pos[u + 1] - g_directed->bw_beg_pos[u]);
        std::vector<long long> tmp_beg(V + 1, 0);
        for (int u = 0; u < V; ++u) tmp_beg[u + 1] = tmp_beg[u] + deg[u];
        std::vector<int> tmp(tmp_beg[V]);
        #pragma omp parallel for schedule(dynamic, 1024)
        for (int u = 0; u < V; ++u) {
            long long k = tmp_beg[u];
            for (index_t i = g_directed->fw_beg_pos[u]; i < g_directed->fw_beg_pos[u + 1]; ++i) tmp[k++] = g_directed->fw_csr[i];
            for (index_t i = g_directed->bw_beg_pos[u]; i < g_directed->bw_beg_pos[u + 1]; ++i) tmp[k++] = g_directed->bw_csr[i];
            auto b = tmp.begin() + tmp_beg[u], e = tmp.begin() + tmp_beg[u + 1];
            std::sort(b, e);
            e = std::unique(b, e);
            e = std::remove(b, e, u);
            deg[u] = e - b;
        }
        sym_beg.assign(V + 1, 0);
        for (int u = 0; u < V; ++u) sym_beg[u + 1] = sym_beg[u] + deg[u];
        sym_adj.resize(sym_beg[V]);
        #pragma omp parallel for schedule(dynamic, 1024)
        for (int u = 0; u < V; ++u)
            std::copy(tmp.begin() + tmp_beg[u], tmp.begin() + tmp_beg[u] + deg[u], sym_adj.begin() + sym_beg[u]);
    }
    // Tombstone bitmaps for the base CSRs (part of the edge store), allocated
    // here rather than on the first deletion, which would otherwise pay
    // ~30 MB of zero filling on livejournal.
    if (dead_fw.empty()) {
        dead_fw.assign(((size_t)g_directed->fw_beg_pos[V] + 63) / 64, 0);
        dead_bw.assign(((size_t)g_directed->bw_beg_pos[V] + 63) / 64, 0);
    }
    if (dead_sym.empty() && !sym_adj.empty()) dead_sym.assign((sym_adj.size() + 63) / 64, 0);

    std::vector<int> parent(V);
    for (int i = 0; i < V; ++i) parent[i] = i;
    auto find = [&](int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    for (int u = 0; u < V; ++u)
        for (index_t i = g_directed->fw_beg_pos[u]; i < g_directed->fw_beg_pos[u + 1]; ++i) {
            int a = find(u), b = find((int)g_directed->fw_csr[i]);
            if (a != b) parent[a] = b;
        }
    wcc_label.assign(V, 0);
    wcc_size.assign(V, 0);
    int count = 0;
    for (int i = 0; i < V; ++i) {
        int r = find(i);
        wcc_label[i] = r;
        if (wcc_size[r]++ == 0) ++count;
    }
    std::cout << "Found " << count << " WCCs." << std::endl;
}

void GIR::relabelComponent(int start, int from_label, int to_label) {
    std::vector<int> q{start};
    wcc_label[start] = to_label;
    int moved = 0;
    for (size_t i = 0; i < q.size(); ++i) {
        ++moved;
        forUnd(q[i], [&](int w) {
            if (wcc_label[w] == from_label) { wcc_label[w] = to_label; q.push_back(w); }
        });
    }
    wcc_size[from_label] -= moved;
    wcc_size[to_label] += moved;
}

bool GIR::undirectedConnected(int u, int v, std::vector<int>& side) {
    side.clear();
    if (u == v) return true;
    int ea = nextEpochA(), eb = nextEpochB();
    stamp_a[u] = ea;
    stamp_b[v] = eb;
    std::vector<int> qa{u}, qb{v};
    size_t ia = 0, ib = 0;
    long long work_a = 0, work_b = 0;
    bool met = false;
    // Expand the side whose next vertex keeps its total work lower (so an
    // isolated endpoint is seen before a hub's list is scanned).
    auto deg = [&](int x) {
        return (sym_beg.empty() ? 0LL : (long long)(sym_beg[x + 1] - sym_beg[x])) + (long long)extra_und_adj[x].size();
    };
    while (ia < qa.size() && ib < qb.size() && !met) {
        if (work_a + deg(qa[ia]) <= work_b + deg(qb[ib])) {
            forUnd(qa[ia++], [&](int w) {
                ++work_a;
                if (stamp_b[w] == eb) met = true;
                else if (stamp_a[w] != ea) { stamp_a[w] = ea; qa.push_back(w); }
            });
        } else {
            forUnd(qb[ib++], [&](int w) {
                ++work_b;
                if (stamp_a[w] == ea) met = true;
                else if (stamp_b[w] != eb) { stamp_b[w] = eb; qb.push_back(w); }
            });
        }
    }
    if (met) return true;
    side = (ia >= qa.size()) ? std::move(qa) : std::move(qb);
    return false;
}

// ============================================================================
// SCC layer + condensation DAG
// ============================================================================

int GIR::newSccId() {
    ++struct_version;
    if (!free_scc_ids.empty()) {
        int id = free_scc_ids.back(); free_scc_ids.pop_back();
        scc_alive[id] = 1;  // its list and maps were released by retireScc
        return id;
    }
    int id = (int)sccs.size();
    sccs.emplace_back();
    scc_alive.push_back(1);
    dag_out.emplace_back();
    dag_in.emplace_back();
    scc_rank.push_back(0);
    return id;
}

void GIR::retireScc(int id) {
    ++struct_version;
    auto it = rank_to_scc.find(scc_rank[id]);
    if (it != rank_to_scc.end() && it->second == id) rank_to_scc.erase(it);
    scc_alive[id] = 0;
    // Release, not clear(): clear() costs O(bucket count) on every reuse, and
    // an id that once held a large SCC keeps a huge bucket array.
    std::vector<int>().swap(sccs[id]);
    std::unordered_map<int, int>().swap(dag_out[id]);
    std::unordered_map<int, int>().swap(dag_in[id]);
    free_scc_ids.push_back(id);
}

void GIR::addDagEdge(int a, int b) {
    if (a == b) return;
    if (dag_out[a][b]++ == 0) ++struct_version;
    dag_in[b][a]++;
}
void GIR::removeDagEdge(int a, int b) {
    auto it = dag_out[a].find(b);
    if (it != dag_out[a].end() && --it->second <= 0) { dag_out[a].erase(it); ++struct_version; }
    auto it2 = dag_in[b].find(a);
    if (it2 != dag_in[b].end() && --it2->second <= 0) dag_in[b].erase(it2);
}
bool GIR::dagEdgeExists(int a, int b) const { return dag_out[a].count(b) > 0; }

void GIR::setRank(int scc, long long r) {
    scc_rank[scc] = r;
    rank_to_scc[r] = scc;
}

void GIR::respreadAround(long long r, long long need) {
    auto it = rank_to_scc.find(r);
    if (it == rank_to_scc.end()) { renumberAllRanks(); return; }
    for (size_t w = 32;; w *= 4) {
        auto lo_it = it;
        size_t below = 0;
        while (below < w && lo_it != rank_to_scc.begin()) { --lo_it; ++below; }
        auto hi_it = it;
        size_t above = 0;
        while (above < w && std::next(hi_it) != rank_to_scc.end()) { ++hi_it; ++above; }
        const bool at_begin = lo_it == rank_to_scc.begin(), at_end = std::next(hi_it) == rank_to_scc.end();
        const long long count = (long long)(below + above + 1);
        const long long pad = RANK_GAP * (count + need + 1);
        const long long Lb = at_begin ? lo_it->first - pad : std::prev(lo_it)->first;
        const long long Rb = at_end ? hi_it->first + pad : std::next(hi_it)->first;
        const long long step = (Rb - Lb) / (count + 2 * need + 1);
        if (step < 2) {
            if (at_begin && at_end) { renumberAllRanks(); return; }
            continue;
        }
        const int anchor = it->second;
        std::vector<int> ids;
        for (auto p = lo_it;; ++p) { ids.push_back(p->second); if (p == hi_it) break; }
        for (int id : ids) rank_to_scc.erase(scc_rank[id]);
        long long pos = 1;
        for (int id : ids) {
            if (id == anchor) { pos += need; setRank(id, Lb + step * pos); pos += need + 1; }
            else { setRank(id, Lb + step * pos); ++pos; }
        }
        return;
    }
}

void GIR::renumberAllRanks() {
    std::vector<int> order;
    order.reserve(rank_to_scc.size());
    for (auto& kv : rank_to_scc) order.push_back(kv.second);
    rank_to_scc.clear();
    for (size_t i = 0; i < order.size(); ++i) setRank(order[i], (long long)(i + 1) * RANK_GAP);
}

void GIR::buildSCC_DAG() {
    std::cout << "=== Build SCC + DAG ===" << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    sccs.clear();
    // SCC ids never exceed V plus a few temporaries (retired ids are
    // reused). Reserving up front keeps newSccId from reallocating these
    // per-SCC vectors (millions of hash maps) mid-update.
    const size_t cap = (size_t)V + 1024;
    sccs.reserve(cap);
    std::vector<int> all(V);
    for (int i = 0; i < V; ++i) all[i] = i;
    tarjanScc(all, [](int) { return true; }, [&](std::vector<int>& c) { sccs.push_back(std::move(c)); });
    std::chrono::duration<double> dt = std::chrono::high_resolution_clock::now() - t0;
    std::cout << "Time to find SCC: " << dt.count() << " seconds\n";

    const int S = (int)sccs.size();
    node_scc.assign(V, -1);
    for (int s = 0; s < S; ++s)
        for (size_t i = 0; i < sccs[s].size(); ++i) { node_scc[sccs[s][i]] = s; scc_pos[sccs[s][i]] = (int)i; }
    scc_alive.reserve(cap);
    dag_out.reserve(cap);
    dag_in.reserve(cap);
    scc_rank.reserve(cap);
    scc_alive.assign(S, 1);
    free_scc_ids.clear();
    dag_out.assign(S, {});
    dag_in.assign(S, {});
    for (int u = 0; u < V; ++u) {
        int su = node_scc[u];
        for (index_t i = g_directed->fw_beg_pos[u]; i < g_directed->fw_beg_pos[u + 1]; ++i) {
            int sw = node_scc[g_directed->fw_csr[i]];
            if (sw != su) { dag_out[su][sw]++; dag_in[sw][su]++; }
        }
    }
    // Tarjan emits SCCs in reverse topological order, so reversing the
    // emission order is already a valid global topological order.
    scc_rank.assign(S, 0);
    rank_to_scc.clear();
    for (int s = 0; s < S; ++s) setRank(s, (long long)(S - s) * RANK_GAP);
    std::cout << "Found " << S << " SCCs, " << S << " ranked in initial condensation DAG." << std::endl;
    pickHubCandidates();  // (for the batch certificates; one parallel pass over V)
}

// Dynamic topological order with SCC merging.
//
// Search: F = descendants of b with rank <= rank(a), B = ancestors of a with
// rank >= rank(b), in lockstep. a is never expanded forward and b never
// backward (a's descendants rank above it, b's ancestors below). Nodes with a
// huge condensation adjacency ("hubs", typically the giant SCC) are not
// expanded on the first attempt.
//
//  - If one side completes without meeting a hub, only that side is
//    re-ranked, into the free gap above rank(a) (F) or below rank(b) (B).
//    Every edge leaving F ends above rank(a) and every edge entering B starts
//    below rank(b), so nothing outside the moved side is violated. A cycle
//    exists iff a is in F (equivalently b in B); the side nodes on it merge.
//  - If both sides stall on the same single hub h, the cycle b ~> h ~> a
//    closes through h, and the merge set is computable without expanding h:
//    a node x != h on a cycle through a->b has a path b ~> x or x ~> a that
//    avoids h (otherwise h ~> x ~> h and x would already be inside h). The
//    merged SCC keeps rank(h); only nodes that would now violate the order
//    move: descendants (below rank h) of members ranked below h go into the
//    gap just above rank(h), ancestors (above rank h) of members ranked above
//    h go into the gap just below.
//  - Otherwise hubs are expanded and the search is repeated.
// Returns 1 if SCCs merged, 0 if not, -1 if `cross_wcc` and both sides
// stalled on hubs (the caller then uses the smaller-WCC move instead).
int GIR::insertDagEdgePK(int a, int b, bool cross_wcc) {
    if (scc_rank[a] < scc_rank[b]) { addDagEdge(a, b); return 0; }
    // An SCC without condensation edges (a vertex joining a component, the
    // common cross-WCC case) has no order constraint but this edge: it goes
    // to the very front of the order (a isolated, before b) or the very end
    // (b isolated, after a). Not into the gap next to the other endpoint:
    // that is often the giant SCC, and thousands of halvings there used up
    // the gaps the remaining Pearce-Kelly moves need (their merges went from
    // ~6.6 to ~15 us). Ranks may go negative.
    {
        const bool a_iso = dag_in[a].empty() && dag_out[a].empty();
        const bool b_iso = !a_iso && dag_in[b].empty() && dag_out[b].empty();
        if (a_iso || b_iso) {
            const int x = a_iso ? a : b;
            rank_to_scc.erase(scc_rank[x]);
            const long long nr = a_iso ? rank_to_scc.begin()->first - RANK_GAP : std::prev(rank_to_scc.end())->first + RANK_GAP;
            scc_rank[x] = nr;
            rank_to_scc.emplace_hint(a_iso ? rank_to_scc.begin() : rank_to_scc.end(), nr, x);
            addDagEdge(a, b);
            order_stats[0]++;
            return 0;
        }
    }
    // Two-node cycles without a search. If b is a's only predecessor, every
    // path b ~> a is the edge b -> a, so with that edge the new edge a -> b
    // closes exactly {a, b}, and the merged node can take b's rank (a's other
    // neighbours are successors, after rank(a) > rank(b)). Symmetrically if a
    // is b's only successor (the merged node takes a's rank). On livejournal
    // most merges are a small SCC joining the giant one this way.
    if (!cross_wcc) {
        long long r = 0;
        bool two = true;
        if (dag_in[a].size() == 1 && dag_in[a].begin()->first == b) r = scc_rank[b];
        else if (dag_out[b].size() == 1 && dag_out[b].begin()->first == a) r = scc_rank[a];
        else two = false;
        if (two) {
            // (rank_to_scc is a tree with an entry per SCC: retiring the
            // absorbed SCC drops its entry; if the survivor keeps its own
            // rank, which is the giant-SCC case, nothing else changes.)
            std::vector<int> M{a, b};
            const int keep = chooseKeep(M);
            mergeBlockForests(M, keep);
            const int s = mergeSccs(M, keep);
            if (scc_rank[s] != r) { rank_to_scc.erase(scc_rank[s]); setRank(s, r); }
            order_stats[0]++;
            return 1;
        }
    }
    const long long ra = scc_rank[a], rb = scc_rank[b];
    if (pk_sf.size() < sccs.size()) { pk_sf.resize(sccs.size(), 0); pk_sb.resize(sccs.size(), 0); }

    for (int attempt = 0; attempt < 2; ++attempt) {
        const size_t HUB = attempt == 0 ? hub_degree : SIZE_MAX;
        if (++pk_ef == INT_MAX) { std::fill(pk_sf.begin(), pk_sf.end(), 0); pk_ef = 1; }
        if (++pk_eb == INT_MAX) { std::fill(pk_sb.begin(), pk_sb.end(), 0); pk_eb = 1; }
        const int ef = pk_ef, eb = pk_eb;
        auto inF = [&](int x) { return pk_sf[x] == ef; };
        auto inB = [&](int x) { return pk_sb[x] == eb; };
        std::vector<int> F{b}, B{a}, hubF, hubB;
        pk_sf[b] = ef;
        pk_sb[a] = eb;
        size_t fi = 0, bi = 0;
        long long wf = 0, wb = 0;
        int complete = -1;  // 0 = F, 1 = B
        auto check_complete = [&]() {
            if (fi >= F.size() && hubF.empty()) complete = 0;
            else if (bi >= B.size() && hubB.empty()) complete = 1;
        };
        while (fi < F.size() || bi < B.size()) {
            bool pickF = fi < F.size() && (bi >= B.size() || wf <= wb);
            if (pickF) {
                int x = F[fi++];
                if (x != a) {
                    if (dag_out[x].size() > HUB) hubF.push_back(x);  // includes b itself
                    else {
                        wf += 1 + (long long)dag_out[x].size();
                        for (auto& kv : dag_out[x])
                            if (scc_rank[kv.first] <= ra && !inF(kv.first)) { pk_sf[kv.first] = ef; F.push_back(kv.first); }
                    }
                }
            } else {
                int x = B[bi++];
                if (x != b) {
                    if (dag_in[x].size() > HUB) hubB.push_back(x);  // includes a itself
                    else {
                        wb += 1 + (long long)dag_in[x].size();
                        for (auto& kv : dag_in[x])
                            if (scc_rank[kv.first] >= rb && !inB(kv.first)) { pk_sb[kv.first] = eb; B.push_back(kv.first); }
                    }
                }
            }
            check_complete();
            if (complete >= 0) break;
        }
        if (complete < 0) check_complete();

        if (complete >= 0) {
            const bool forward = complete == 0;
            std::vector<int>& side = forward ? F : B;
            const bool cycle = forward ? inF(a) : inB(b);
            auto gap = [&](long long& lo, long long& hi) {
                if (forward) {
                    lo = scc_rank[a];
                    auto nx = std::next(rank_to_scc.find(lo));
                    hi = (nx == rank_to_scc.end()) ? lo + RANK_GAP : nx->first;
                } else {
                    hi = scc_rank[b];
                    auto it = rank_to_scc.find(hi);
                    lo = (it == rank_to_scc.begin()) ? hi - RANK_GAP : std::prev(it)->first;
                }
            };
            long long lo, hi;
            gap(lo, hi);
            if (hi - lo < (long long)side.size() + 1) { respreadAround(forward ? scc_rank[a] : scc_rank[b], (long long)side.size() + 1); gap(lo, hi); }
            std::sort(side.begin(), side.end(), [&](int x, int y) { return scc_rank[x] < scc_rank[y]; });
            addDagEdge(a, b);
            std::vector<int> order;  // topological: old rank order, merged node in place of its members
            if (cycle) {
                // Side nodes on a cycle through a->b: those that reach a
                // (forward side) or are reachable from b (backward side).
                std::unordered_set<int> Ms;
                if (forward) {
                    for (size_t i = side.size(); i-- > 0;) {
                        int x = side[i];
                        if (x == a) { Ms.insert(x); continue; }
                        for (auto& kv : dag_out[x]) if (Ms.count(kv.first)) { Ms.insert(x); break; }
                    }
                } else {
                    for (int x : side) {
                        if (x == b) { Ms.insert(x); continue; }
                        for (auto& kv : dag_in[x]) if (Ms.count(kv.first)) { Ms.insert(x); break; }
                    }
                }
                std::vector<int> M(Ms.begin(), Ms.end());
                // The moved nodes are re-sorted topologically (Kahn) once the
                // members have collapsed into one node.
                for (int x : side) rank_to_scc.erase(scc_rank[x]);
                int keep = chooseKeep(M);
                mergeBlockForests(M, keep);
                int survivor = mergeSccs(M, keep);
                std::vector<int> nodes;
                for (int x : side) if (!Ms.count(x)) nodes.push_back(x);
                nodes.push_back(survivor);
                order = kahnOrder(nodes);
            } else {
                for (int x : side) rank_to_scc.erase(scc_rank[x]);
                order = side;
            }
            const long long step = (hi - lo) / ((long long)order.size() + 1);
            for (size_t i = 0; i < order.size(); ++i) setRank(order[i], lo + step * (long long)(i + 1));
            order_stats[attempt == 0 ? 0 : 2]++;
            order_moved[0] += (long long)order.size();
            order_moved[1] = std::max<long long>(order_moved[1], (long long)order.size());
            if (order.size() > 10000) order_moved[2]++;
            return cycle ? 1 : 0;
        }

        // Both sides stalled on hubs.
        std::sort(hubF.begin(), hubF.end()); hubF.erase(std::unique(hubF.begin(), hubF.end()), hubF.end());
        std::sort(hubB.begin(), hubB.end()); hubB.erase(std::unique(hubB.begin(), hubB.end()), hubB.end());
        if (hubF.size() == 1 && hubB.size() == 1 && hubF[0] == hubB[0] && hubF[0] != a && hubF[0] != b) {
            hubMerge(a, b, hubF[0], F, B, inF, inB);
            order_stats[1]++;
            return 1;
        }
        if (cross_wcc) { order_stats[3]++; return -1; }
    }
    return 0;  // unreachable: the second attempt has no hubs
}

template <class InF, class InB>
void GIR::hubMerge(int a, int b, int h, std::vector<int>& F, std::vector<int>& B, InF inF, InB inB) {
    const long long rh = scc_rank[h];
    auto by_rank = [&](int x, int y) { return scc_rank[x] < scc_rank[y]; };
    std::sort(F.begin(), F.end(), by_rank);
    std::sort(B.begin(), B.end(), by_rank);
    std::unordered_set<int> Ms{h};
    for (int x : F) if (inB(x)) Ms.insert(x);
    {  // F nodes reaching h (decreasing rank; successors first)
        std::unordered_set<int> r;
        for (size_t i = F.size(); i-- > 0;) {
            int x = F[i];
            if (x == h || x == a) continue;
            for (auto& kv : dag_out[x]) if (kv.first == h || r.count(kv.first)) { r.insert(x); break; }
        }
        Ms.insert(r.begin(), r.end());
    }
    {  // B nodes reachable from h (increasing rank; predecessors first)
        std::unordered_set<int> r;
        for (int x : B) {
            if (x == h || x == b) continue;
            for (auto& kv : dag_in[x]) if (kv.first == h || r.count(kv.first)) { r.insert(x); break; }
        }
        Ms.insert(r.begin(), r.end());
    }
    std::vector<int> M(Ms.begin(), Ms.end());
    // Nodes that would violate the order once the members collapse to rank(h).
    std::vector<int> Y, Z;
    std::unordered_set<int> Ys, Zs;
    for (int x : M) {
        if (x == h) continue;
        if (scc_rank[x] < rh) {
            for (auto& kv : dag_out[x])
                if (!Ms.count(kv.first) && scc_rank[kv.first] < rh && Ys.insert(kv.first).second) Y.push_back(kv.first);
        } else {
            for (auto& kv : dag_in[x])
                if (!Ms.count(kv.first) && scc_rank[kv.first] > rh && Zs.insert(kv.first).second) Z.push_back(kv.first);
        }
    }
    for (size_t i = 0; i < Y.size(); ++i)
        for (auto& kv : dag_out[Y[i]])
            if (scc_rank[kv.first] < rh && !Ms.count(kv.first) && Ys.insert(kv.first).second) Y.push_back(kv.first);
    for (size_t i = 0; i < Z.size(); ++i)
        for (auto& kv : dag_in[Z[i]])
            if (scc_rank[kv.first] > rh && !Ms.count(kv.first) && Zs.insert(kv.first).second) Z.push_back(kv.first);
    std::sort(Y.begin(), Y.end(), by_rank);
    std::sort(Z.begin(), Z.end(), by_rank);

    addDagEdge(a, b);
    for (int x : M) rank_to_scc.erase(scc_rank[x]);
    for (int x : Y) rank_to_scc.erase(scc_rank[x]);
    for (int x : Z) rank_to_scc.erase(scc_rank[x]);
    int keep = chooseKeep(M);
    mergeBlockForests(M, keep);
    int survivor = mergeSccs(M, keep);
    setRank(survivor, rh);
    auto place = [&]() -> bool {
        long long r = scc_rank[survivor];
        auto it = rank_to_scc.find(r);
        auto nx = std::next(it);
        long long hi_y = (nx == rank_to_scc.end()) ? r + RANK_GAP : nx->first;
        long long lo_z = (it == rank_to_scc.begin()) ? r - RANK_GAP : std::prev(it)->first;
        if (hi_y - r < (long long)Y.size() + 1 || r - lo_z < (long long)Z.size() + 1) return false;
        long long sy = (hi_y - r) / ((long long)Y.size() + 1), sz = (r - lo_z) / ((long long)Z.size() + 1);
        for (size_t i = 0; i < Y.size(); ++i) setRank(Y[i], r + sy * (long long)(i + 1));
        for (size_t i = 0; i < Z.size(); ++i) setRank(Z[i], lo_z + sz * (long long)(i + 1));
        return true;
    };
    if (!place()) { respreadAround(scc_rank[survivor], (long long)std::max(Y.size(), Z.size()) + 1); place(); }
}

int GIR::chooseKeep(const std::vector<int>& M) const {
    // Survivor: the member whose data is costliest to move (vertex list plus
    // condensation adjacency), so merging only ever touches the small ones.
    int keep = M[0];
    auto cost = [&](int x) { return sccs[x].size() + dag_out[x].size() + dag_in[x].size(); };
    for (int x : M) if (cost(x) > cost(keep)) keep = x;
    return keep;
}

std::vector<int> GIR::kahnOrder(const std::vector<int>& nodes) {
    std::unordered_map<int, int> indeg;
    for (int x : nodes) indeg[x] = 0;
    auto succ_in = [&](int x, auto&& f) {
        if (dag_out[x].size() <= nodes.size()) {
            for (auto& kv : dag_out[x]) if (indeg.count(kv.first)) f(kv.first);
        } else {
            for (int y : nodes) if (dag_out[x].count(y)) f(y);
        }
    };
    for (int x : nodes) succ_in(x, [&](int y) { indeg[y]++; });
    std::vector<int> order;
    for (int x : nodes) if (indeg[x] == 0) order.push_back(x);
    for (size_t i = 0; i < order.size(); ++i)
        succ_in(order[i], [&](int y) { if (--indeg[y] == 0) order.push_back(y); });
    return order;
}

void GIR::insertCrossWccEdge(int a, int b, bool move_a_side) {
    if (scc_rank[a] < scc_rank[b]) { addDagEdge(a, b); return; }
    std::vector<int> nodes{move_a_side ? a : b};
    std::unordered_set<int> seen(nodes.begin(), nodes.end());
    const long long ra = scc_rank[a], rb = scc_rank[b];
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (move_a_side) {
            for (auto& kv : dag_in[nodes[i]])
                if (scc_rank[kv.first] >= rb && seen.insert(kv.first).second) nodes.push_back(kv.first);
        } else {
            for (auto& kv : dag_out[nodes[i]])
                if (scc_rank[kv.first] <= ra && seen.insert(kv.first).second) nodes.push_back(kv.first);
        }
    }
    std::sort(nodes.begin(), nodes.end(), [&](int x, int y) { return scc_rank[x] < scc_rank[y]; });
    // Gap strictly below rank(b) (moving a's ancestors) or strictly above
    // rank(a) (moving b's descendants). No moved node owns the neighbouring
    // key: a's ancestors rank >= rank(b), b's descendants rank <= rank(a).
    auto gap = [&](long long& lo, long long& hi) {
        if (move_a_side) {
            hi = scc_rank[b];
            auto it = rank_to_scc.find(hi);
            lo = (it == rank_to_scc.begin()) ? hi - RANK_GAP : std::prev(it)->first;
        } else {
            lo = scc_rank[a];
            auto nx = std::next(rank_to_scc.find(lo));
            hi = (nx == rank_to_scc.end()) ? lo + RANK_GAP : nx->first;
        }
    };
    const long long need = (long long)nodes.size() + 1;
    long long lo, hi;
    gap(lo, hi);
    if (hi - lo < need) { respreadAround(move_a_side ? scc_rank[b] : scc_rank[a], need); gap(lo, hi); }
    for (int x : nodes) rank_to_scc.erase(scc_rank[x]);
    const long long step = (hi - lo) / need;
    for (size_t i = 0; i < nodes.size(); ++i) setRank(nodes[i], lo + step * (long long)(i + 1));
    addDagEdge(a, b);
}

int GIR::mergeSccs(const std::vector<int>& ids, int keep) {
    if (keep < 0) keep = chooseKeep(ids);
    // (Most merges have two members: a linear test instead of building a set.)
    std::unordered_set<int> Mset;
    const bool small = ids.size() <= 8;
    if (!small) Mset.insert(ids.begin(), ids.end());
    struct { const std::vector<int>& ids; const std::unordered_set<int>& s; bool small;
             bool count(int x) const { return small ? std::find(ids.begin(), ids.end(), x) != ids.end() : s.count(x) > 0; } }
        M{ids, Mset, small};
    for (int o : ids) {
        if (o == keep) continue;
        for (auto& kv : dag_out[o]) {
            if (M.count(kv.first)) continue;
            dag_in[kv.first].erase(o);
            dag_out[keep][kv.first] += kv.second;
            dag_in[kv.first][keep] += kv.second;
        }
        for (auto& kv : dag_in[o]) {
            if (M.count(kv.first)) continue;
            dag_out[kv.first].erase(o);
            dag_in[keep][kv.first] += kv.second;
            dag_out[kv.first][keep] += kv.second;
        }
        dag_out[keep].erase(o);
        dag_in[keep].erase(o);
        for (int v : sccs[o]) { node_scc[v] = keep; scc_pos[v] = (int)sccs[keep].size(); sccs[keep].push_back(v); }
        retireScc(o);
    }
    return keep;
}

// A 2-hop detour u->w->v (w inside SCC scc if scc >= 0): very common on
// social graphs. Iterates the shorter sorted base list and binary-searches the
// longer one; read-only, so safe in parallel.
bool GIR::twoHop(int u, int v, int scc, int* mid) const {
    const index_t* fb = g_directed->fw_beg_pos; const vertex_t* fc = g_directed->fw_csr;
    const index_t* bb = g_directed->bw_beg_pos; const vertex_t* bc = g_directed->bw_csr;
    const vertex_t *a0 = fc + fb[u], *a1 = fc + fb[u + 1];   // out(u)
    const vertex_t *b0 = bc + bb[v], *b1 = bc + bb[v + 1];   // in(v)
    bool a_short = (a1 - a0) <= (b1 - b0);
    const vertex_t *s0 = a_short ? a0 : b0, *s1 = a_short ? a1 : b1;
    const vertex_t *l0 = a_short ? b0 : a0, *l1 = a_short ? b1 : a1;
    int budget = 4096;
    const bool chk = n_dead_dir != 0;
    for (const vertex_t* p = s0; p < s1 && budget-- > 0; ++p) {
        int w = *p;
        if (w == u || w == v || (scc >= 0 && node_scc[w] != scc)) continue;
        const vertex_t* q = std::lower_bound(l0, l1, (vertex_t)w);
        if (q == l1 || *q != (vertex_t)w) continue;
        // u->w sits in fw_csr, w->v in bw_csr, whichever list each pointer walks.
        if (chk) {
            const vertex_t* pa = a_short ? p : q;  // in out(u)
            const vertex_t* pb = a_short ? q : p;  // in in(v)
            if (bitGet(dead_fw, pa - fc) || bitGet(dead_bw, pb - bc)) continue;
        }
        if (mid) *mid = w;
        return true;
    }
    return false;
}

bool GIR::reachableWithinScc(int u, int v, int scc, std::vector<int>* side, bool* fwd) {
    if (u == v) return true;
    // Cheap certificate first: a 2-hop detour u->w->v inside the SCC.
    if (twoHop(u, v, scc)) return true;
    int ea = nextEpochA(), eb = nextEpochB();
    stamp_a[u] = ea;
    stamp_b[v] = eb;
    std::vector<int> qa{u}, qb{v};
    size_t ia = 0, ib = 0;
    long long work_a = 0, work_b = 0;
    bool met = false;
    auto odeg = [&](int x) { return (long long)(g_directed->fw_beg_pos[x + 1] - g_directed->fw_beg_pos[x]) + (long long)extra_out_adj[x].size(); };
    auto ideg = [&](int x) { return (long long)(g_directed->bw_beg_pos[x + 1] - g_directed->bw_beg_pos[x]) + (long long)extra_in_adj[x].size(); };
    // Hubs (degree above search_hub_deg) are not expanded when their turn
    // comes but deferred: every vertex the other side reaches is tested for
    // an edge to them (binary search). A side expands its deferred hubs only
    // once its queue is empty, so an exhausted side is still everything it
    // can reach, and `side` stays exact.
    std::vector<int> ha, hb;
    size_t ja = 0, jb = 0;
    auto visitA = [&](int w) {
        ++work_a;
        if (met || node_scc[w] != scc) return;
        if (stamp_b[w] == eb) { met = true; return; }
        if (stamp_a[w] == ea) return;
        stamp_a[w] = ea;
        qa.push_back(w);
        for (int h : hb) if (isDirectedEdgePresent(w, h)) { met = true; return; }
    };
    auto visitB = [&](int w) {
        ++work_b;
        if (met || node_scc[w] != scc) return;
        if (stamp_a[w] == ea) { met = true; return; }
        if (stamp_b[w] == eb) return;
        stamp_b[w] = eb;
        qb.push_back(w);
        for (int h : ha) if (isDirectedEdgePresent(h, w)) { met = true; return; }
    };
    while (!met && (ia < qa.size() || ja < ha.size()) && (ib < qb.size() || jb < hb.size())) {
        // The cheaper next step wins. A side left with only deferred hubs is
        // charged 8x the hub's degree: the other side usually meets the hub
        // through the adjacency tests long before that (measured faster on
        // twitter than charging the plain degree), and the hub is still
        // expanded eventually (no budget here, so exhaustion stays exact).
        const long long ca = ia < qa.size() ? odeg(qa[ia]) : 8 * odeg(ha[ja]);
        const long long cb = ib < qb.size() ? ideg(qb[ib]) : 8 * ideg(hb[jb]);
        const bool pick_a = work_a + ca <= work_b + cb;
        if (pick_a) {
            int x;
            if (ia < qa.size()) {
                x = qa[ia++];
                if (ca > search_hub_deg && ha.size() < 32) { ha.push_back(x); work_a += 32; continue; }
            } else {
                x = ha[ja++];
            }
            forOut(x, visitA);
        } else {
            int x;
            if (ib < qb.size()) {
                x = qb[ib++];
                if (cb > search_hub_deg && hb.size() < 32) { hb.push_back(x); work_b += 32; continue; }
            } else {
                x = hb[jb++];
            }
            forIn(x, visitB);
        }
    }
    if (met) return true;
    if (side) {
        bool f = ia >= qa.size() && ja >= ha.size();
        *fwd = f;
        *side = f ? std::move(qa) : std::move(qb);
    }
    return false;
}

bool GIR::removeVerticesFromBlocks(const std::vector<int>& X, const std::vector<std::pair<int, int>>& extra_att) {
    if (blocks_on_demand) {
        for (int x : X) markBlocksDirty(x);
        for (auto& e : extra_att) markBlocksDirty(e.second);
        return true;
    }
    ++block_version;
    std::unordered_set<int> Xs(X.begin(), X.end());
    std::unordered_map<int, std::vector<int>> by_block;
    std::vector<int> order;
    // An articulation point x in X: R stays connected without x, so R lies
    // entirely on one side of x: either x's parent side (x's child subtrees
    // then lie inside X and vanish) or a single child subtree of x (the
    // parent side vanishes and that child block becomes the new root; see
    // the parent-AP check below).
    for (int x : X) {
        int b = findBlk(v_block[x]);
        if (!by_block.count(b)) order.push_back(b);
        by_block[b].push_back(x);
    }
    // Child blocks of a departing articulation point contain it only as their
    // parent AP, so they are not reached through v_block: add them explicitly.
    for (int x : X) {
        if (v_bcnt[x] < 2) continue;
        forUnd(x, [&](int w) {
            if (node_scc[w] != node_scc[x] && !Xs.count(w)) return;
            int b = commonBlock(x, w);
            if (b >= 0 && blk_parent_ap[b] == x && !by_block.count(b)) { by_block[b]; order.push_back(b); }
        });
    }
    // A child block reached only through a just-deleted edge (still in the
    // blocks, no longer in the graph) is found through the deleted edges.
    for (auto& xy : extra_att) {
        int x = xy.first;
        if (!Xs.count(x) || v_bcnt[x] < 2) continue;
        int b = commonBlock(x, xy.second);
        if (b >= 0 && blk_parent_ap[b] == x && !by_block.count(b)) { by_block[b]; order.push_back(b); }
    }
    // Each touched block is handled on its own: removing X can only split
    // blocks (never merge them), and b \ X stays connected because the rest
    // of the SCC is still strongly connected.
    for (int b : order) {
        // Count against the block's own vertex list: its parent AP may be in
        // X without being grouped under it.
        // Members of X in b: those grouped under b (v_block) plus b's parent
        // AP if it is in X -- no scan of b's (possibly huge) vertex list.
        std::vector<int> xb = by_block[b];
        auto& bv = blk_verts[b];
        if (blk_parent_ap[b] >= 0 && Xs.count(blk_parent_ap[b])) xb.push_back(blk_parent_ap[b]);
        std::unordered_set<int> gone(xb.begin(), xb.end());
        size_t rest = bv.size() - xb.size();
        if (rest >= 2) {
            // The block hung off an articulation point that leaves with X: it
            // becomes the root of what remains of the tree.
            if (blk_parent_ap[b] >= 0 && Xs.count(blk_parent_ap[b])) blk_parent_ap[b] = -1;
            if (bv.size() >= fastpath_min_block) {
                std::vector<int> seeds, extra;
                for (int x : xb) {
                    forUnd(x, [&](int w) { if (!Xs.count(w) && inBlock(w, b)) seeds.push_back(w); });
                    for (auto& xy : extra_att)
                        if (xy.first == x) { seeds.push_back(xy.second); extra.push_back(xy.second); }
                }
                if (peelAndCertify(b, gone, seeds, extra)) continue;
            }
            recomputeBlockWithout(b, gone);
            continue;
        }
        if (rest == 0) {  // a block hanging off an articulation point of X: gone entirely
            std::vector<int>().swap(bv);
            continue;
        }
        // rest == 1: b was {x, w}; it disappears and w loses one block.
        int w = -1;
        for (int y : bv) if (!Xs.count(y)) w = y;
        std::vector<int>().swap(bv);
        v_bcnt[w]--;
        if (blk_parent_ap[b] == w) continue;  // b hung off w
        // b was w's parent block (b was a root). Promote one of w's child
        // blocks to root; if w has none left, it becomes a singleton.
        int child = -1;
        forUnd(w, [&](int z) {
            if (child >= 0 || Xs.count(z) || node_scc[z] != node_scc[w]) return;
            int cb = commonBlock(w, z);
            if (cb >= 0 && cb != b) child = cb;
        });
        if (child < 0) {
            int nb = newBlock();
            blk_verts[nb] = {w};
            v_block[w] = nb;
            v_bcnt[w] = 1;
        } else {
            blk_parent_ap[child] = -1;
            v_block[w] = child;
        }
    }
    return true;
}

// Single deletion u->v broke SCC s; the probe exhausted the small side X
// (= SCC(v) if it was the backward side, SCC(u) if forward). Peel small pieces
// off until the rest R is certified strongly connected, then detach them.
//
// Certificate (any anchor a in R): with U = the removed vertices, let Z be the
// R-vertices with an edge into U and E those with an edge out of U, counting
// the deleted edge u->v if it ran between R and U. Every r in R reaches some
// z in Z and is reached from some e in E inside R (follow an old path to or
// from U up to its first/last visit of U). So R is strongly connected iff every
// z reaches a and a reaches every e. When a probe fails, the side it exhausted
// is small and therefore disjoint from the giant SCC (the closure of a giant
// vertex contains the whole giant SCC): it joins U, and if it contains the
// anchor, the probe's other endpoint becomes the anchor. Pieces may end up on
// either side of R in the condensation; a topological sort of the small DAG
// (pieces + R) decides. Falls back to a scoped Tarjan over all of s when U
// would exceed half of s or the boundary is large.
void GIR::splitAfterDelete(int s, int u, int v, std::vector<int>& X, bool X_is_forward) {
    const size_t C = sccs[s].size();
    if (force_resplit || X.size() * 2 > C) { split_stats[1]++; resplitScc(s, true, deletedEdgeAttachment(u, v)); return; }
    const int tmp = newSccId();
    std::vector<int> U = X;
    for (int x : U) node_scc[x] = tmp;  // probes below stay inside R
    auto inR = [&](int w) { return node_scc[w] == s; };
    auto restore = [&]() { for (int x : U) node_scc[x] = s; retireScc(tmp); };
    int a = X_is_forward ? v : u;  // the deleted edge's endpoint outside X
    bool ok = false;
    for (int round = 0; round < 256; ++round) {
        std::vector<int> Z, E;
        {
            std::unordered_set<int> zs, es;
            for (int x : U) {
                forIn(x, [&](int w) { if (inR(w) && zs.insert(w).second) Z.push_back(w); });
                forOut(x, [&](int w) { if (inR(w) && es.insert(w).second) E.push_back(w); });
            }
            if (inR(u) && node_scc[v] == tmp && zs.insert(u).second) Z.push_back(u);
            if (inR(v) && node_scc[u] == tmp && es.insert(v).second) E.push_back(v);
        }
        if (Z.size() + E.size() > 512) { split_stats[2]++; break; }
        std::vector<int> grow;
        int new_anchor = -1;
        for (int z : Z) {
            if (z == a) continue;
            std::vector<int> side;
            bool fwd = false;
            if (reachRO(z, a, s) || reachableWithinScc(z, a, s, &side, &fwd)) continue;
            grow = std::move(side);           // fwd: what z reaches (cannot reach a)
            if (!fwd) new_anchor = z;         // backward side of a: a itself leaves
            break;
        }
        if (grow.empty())
            for (int e2 : E) {
                if (e2 == a) continue;
                std::vector<int> side;
                bool fwd = false;
                if (reachRO(a, e2, s) || reachableWithinScc(a, e2, s, &side, &fwd)) continue;
                grow = std::move(side);       // !fwd: what reaches e2 (unreachable from a)
                if (fwd) new_anchor = e2;     // forward side of a: a itself leaves
                break;
            }
        if (grow.empty()) { ok = true; break; }
        if ((U.size() + grow.size()) * 2 > C) { split_stats[1]++; break; }
        for (int x : grow) node_scc[x] = tmp;
        U.insert(U.end(), grow.begin(), grow.end());
        if (new_anchor >= 0) a = new_anchor;
    }
    if (!ok) { split_stats[3]++; restore(); resplitScc(s, true, deletedEdgeAttachment(u, v)); return; }

    std::vector<std::vector<int>> below, above;
    orderPieces(s, tmp, U, below, above);
    for (int x : U) node_scc[x] = s;
    retireScc(tmp);
    // If the deleted edge's undirected copy is gone and it ran between U and
    // the rest, it was an edge of the block u and v shared.
    detachPieces(s, below, above, deletedEdgeAttachment(u, v));
}

void GIR::orderPieces(int s, int tmp, const std::vector<int>& U, std::vector<std::vector<int>>& below,
                      std::vector<std::vector<int>>& above) {
    auto inR = [&](int w) { return node_scc[w] == s; };
    // SCCs of U, then their order relative to the rest R: Kahn on the small
    // DAG of pieces plus one node for R.
    std::vector<std::vector<int>> frags;
    tarjanScc(U, [&](int w) { return node_scc[w] == tmp; },
              [&](std::vector<int>& comp) { frags.push_back(std::move(comp)); });
    const int K = (int)frags.size(), RN = K;
    std::unordered_map<int, int> fid;
    for (int i = 0; i < K; ++i) for (int x : frags[i]) fid[x] = i;
    std::vector<std::unordered_set<int>> succ(K + 1);
    for (int i = 0; i < K; ++i)
        for (int x : frags[i]) {
            forOut(x, [&](int w) {
                auto it = fid.find(w);
                if (it != fid.end()) { if (it->second != i) succ[i].insert(it->second); }
                else if (inR(w)) succ[i].insert(RN);
            });
            forIn(x, [&](int w) { if (inR(w)) succ[RN].insert(i); });
        }
    std::vector<int> indeg(K + 1, 0), order;
    for (auto& sset : succ) for (int t : sset) indeg[t]++;
    for (int i = 0; i <= K; ++i) if (!indeg[i]) order.push_back(i);
    for (size_t i = 0; i < order.size(); ++i)
        for (int t : succ[order[i]]) if (--indeg[t] == 0) order.push_back(t);
    bool seen_r = false;
    for (int i : order) {
        if (i == RN) { seen_r = true; continue; }
        (seen_r ? above : below).push_back(std::move(frags[i]));
    }
}

// Detach the given pieces (each a new SCC, listed in topological order) from
// SCC s, whose remaining vertices keep id, rank and blocks. `below` pieces
// precede the rest of s, `above` pieces follow it. Only the pieces' vertices
// and edges are touched. `extra` = (x, y): a just-deleted undirected edge
// between piece vertex x and rest vertex y (block attachment), or (-1, -1).
// (x, y) for the just-deleted edge u->v if its undirected copy is gone and
// u, v shared a block (either orientation; detachPieces keeps the one whose
// first vertex leaves).
std::pair<int, int> GIR::deletedEdgeAttachment(int u, int v) {
    if (isUndirectedEdgePresent(u, v) || commonBlock(u, v) < 0) return {-1, -1};
    return {u, v};
}

void GIR::detachPieces(int s, std::vector<std::vector<int>>& below, std::vector<std::vector<int>>& above,
                       std::vector<std::pair<int, int>> extras) {
    std::vector<int> U;
    for (auto* L : {&below, &above}) for (auto& f : *L) U.insert(U.end(), f.begin(), f.end());
    std::unordered_set<int> Us(U.begin(), U.end());
    // Keep (piece vertex, rest vertex) pairs, oriented so the first one leaves.
    std::vector<std::pair<int, int>> ex;
    for (auto e : extras) {
        if (!Us.count(e.first) && Us.count(e.second)) std::swap(e.first, e.second);
        if (Us.count(e.first) && !Us.count(e.second)) ex.push_back(e);
    }
    bool blocks_ok = removeVerticesFromBlocks(U, ex);
    for (int x : U) sccListRemove(s, x);
    if (!sccs[s].empty()) markBlocksDirty(sccs[s][0]);  // (on demand) what stays of s lost vertices
    std::vector<int> ids_below, ids_above;
    auto make = [&](std::vector<std::vector<int>>& L, std::vector<int>& ids) {
        for (auto& f : L) {
            int id = newSccId();
            for (size_t j = 0; j < f.size(); ++j) { node_scc[f[j]] = id; scc_pos[f[j]] = (int)j; }
            sccs[id] = std::move(f);
            ids.push_back(id);
        }
    };
    make(below, ids_below);
    make(above, ids_above);
    // Condensation edges: only the pieces' edges change owner.
    for (int x : U) {
        int fx = node_scc[x];
        forOut(x, [&](int w) {
            if (Us.count(w)) { if (node_scc[w] != fx) addDagEdge(fx, node_scc[w]); return; }
            int t = node_scc[w];
            if (t != s) removeDagEdge(s, t);
            addDagEdge(fx, t);
        });
        forIn(x, [&](int w) {
            if (Us.count(w)) return;
            int t = node_scc[w];
            if (t != s) removeDagEdge(t, s);
            addDagEdge(t, fx);
        });
    }
    // Ranks: below-pieces in the gap just below s, above-pieces just above.
    // External predecessors of a piece were predecessors of s (rank below
    // s) and external successors were successors of s (rank above s).
    const long long need = (long long)std::max(ids_below.size(), ids_above.size()) + 1;
    auto place = [&]() -> bool {
        long long r = scc_rank[s];
        auto it = rank_to_scc.find(r);
        long long lo = (it == rank_to_scc.begin()) ? r - RANK_GAP : std::prev(it)->first;
        auto nx = std::next(it);
        long long hi = (nx == rank_to_scc.end()) ? r + RANK_GAP : nx->first;
        if (r - lo < need || hi - r < need) return false;
        long long sb = (r - lo) / ((long long)ids_below.size() + 1), sa = (hi - r) / ((long long)ids_above.size() + 1);
        for (size_t i = 0; i < ids_below.size(); ++i) setRank(ids_below[i], lo + sb * (long long)(i + 1));
        for (size_t i = 0; i < ids_above.size(); ++i) setRank(ids_above[i], r + sa * (long long)(i + 1));
        return true;
    };
    if (!place()) { respreadAround(scc_rank[s], need); place(); }
    for (int id : ids_below) rebuildSccBlocks(id);
    for (int id : ids_above) rebuildSccBlocks(id);
    if (!blocks_ok) rebuildSccBlocks(s);
    split_stats[blocks_ok ? 0 : 4]++;
}

bool GIR::resplitScc(int s, bool incremental, std::pair<int, int> extra) {
    std::vector<std::vector<int>> frags;
    tarjanScc(sccs[s], [&](int w) { return node_scc[w] == s; },
              [&](std::vector<int>& c) { frags.push_back(std::move(c)); });
    if (frags.size() <= 1) return false;
    std::reverse(frags.begin(), frags.end());  // now topological

    // One dominant fragment: it keeps the old id, rank and blocks; only the
    // small fragments are detached (O(their size) instead of O(|s|)).
    size_t g = 0;
    for (size_t i = 1; i < frags.size(); ++i) if (frags[i].size() > frags[g].size()) g = i;
    if (incremental && frags[g].size() * 2 > sccs[s].size()) {
        std::vector<std::vector<int>> below(std::make_move_iterator(frags.begin()), std::make_move_iterator(frags.begin() + g));
        std::vector<std::vector<int>> above(std::make_move_iterator(frags.begin() + g + 1), std::make_move_iterator(frags.end()));
        detachPieces(s, below, above, extra);
        return true;
    }

    // Fragments go into the rank gap around the old SCC: every external
    // predecessor ranks below it and every external successor above it.
    const long long need = (long long)frags.size() + 1;
    auto gap = [&](long long& lo, long long& hi) {
        long long r = scc_rank[s];
        auto it = rank_to_scc.find(r);
        lo = (it == rank_to_scc.begin()) ? r - RANK_GAP : std::prev(it)->first;
        auto nx = std::next(it);
        hi = (nx == rank_to_scc.end()) ? r + RANK_GAP : nx->first;
    };
    long long lo, hi;
    gap(lo, hi);
    if (hi - lo < need) { respreadAround(scc_rank[s], need); gap(lo, hi); }

    std::vector<int> old_verts = std::move(sccs[s]);
    for (auto& kv : dag_out[s]) dag_in[kv.first].erase(s);
    for (auto& kv : dag_in[s]) dag_out[kv.first].erase(s);
    retireScc(s);

    std::vector<int> ids;
    const long long step = (hi - lo) / need;
    for (size_t i = 0; i < frags.size(); ++i) {
        int id = newSccId();
        for (size_t j = 0; j < frags[i].size(); ++j) { node_scc[frags[i][j]] = id; scc_pos[frags[i][j]] = (int)j; }
        sccs[id] = std::move(frags[i]);
        setRank(id, lo + step * (long long)(i + 1));
        ids.push_back(id);
    }
    int ea = nextEpochA();
    for (int x : old_verts) stamp_a[x] = ea;
    for (int x : old_verts) {
        int sx = node_scc[x];
        forOut(x, [&](int w) { if (node_scc[w] != sx) addDagEdge(sx, node_scc[w]); });
        forIn(x, [&](int w) { if (stamp_a[w] != ea) addDagEdge(node_scc[w], sx); });
    }
    for (int id : ids) rebuildSccBlocks(id);
    return true;
}

// ============================================================================
// Block-cut forest
// ============================================================================

int GIR::findBlk(int b) const {
    int r = b;
    while (blk_uf[r] != r) r = blk_uf[r];
    while (blk_uf[b] != r) { int n = blk_uf[b]; blk_uf[b] = r; b = n; }
    return r;
}

int GIR::newBlock() {
    ++block_version;
    if (!free_blk_ids.empty()) {
        int id = free_blk_ids.back();
        free_blk_ids.pop_back();
        blk_uf[id] = id;
        blk_parent_ap[id] = -1;
        blk_verts[id].clear();
        return id;
    }
    int id = (int)blk_uf.size();
    blk_uf.push_back(id);
    blk_parent_ap.push_back(-1);
    blk_verts.emplace_back();
    return id;
}

bool GIR::inBlock(int x, int b) const { return findBlk(v_block[x]) == b || blk_parent_ap[b] == x; }

// Read-only variants (no path compression) for parallel filters.
void GIR::sccListRemove(int s, int x) {
    auto& L = sccs[s];
    size_t i = (size_t)scc_pos[x];
    if (i >= L.size() || L[i] != x) i = (size_t)(std::find(L.begin(), L.end(), x) - L.begin());
    if (i >= L.size()) return;
    int last = L.back();
    L[i] = last;
    scc_pos[last] = (int)i;
    L.pop_back();
}

void GIR::blkListRemove(int b, int x) {
    auto& L = blk_verts[b];
    size_t i = (size_t)blk_pos[x];
    if (i >= L.size() || L[i] != x) i = (size_t)(std::find(L.begin(), L.end(), x) - L.begin());
    if (i >= L.size()) return;
    int last = L.back();
    L[i] = last;
    blk_pos[last] = (int)i;  // only meaningful if b is last's parent block; else a harmless stale hint
    L.pop_back();
}

void GIR::setBlkPos(int b) {
    const auto& L = blk_verts[b];
    for (size_t i = 0; i < L.size(); ++i) blk_pos[L[i]] = (int)i;  // parent APs get re-hinted by their own parent later
}

int GIR::nextRmEpoch(std::vector<int>& st, int& ep) {
    if (st.size() < (size_t)V) st.assign(V, 0);
    if (++ep == INT_MAX) { std::fill(st.begin(), st.end(), 0); ep = 1; }
    return ep;
}

int GIR::findBlkRO(int b) const {
    while (blk_uf[b] != b) b = blk_uf[b];
    return b;
}
int GIR::commonBlockRO(int x, int y) const {
    int bx = findBlkRO(v_block[x]), by = findBlkRO(v_block[y]);
    if (bx == by) return bx;
    if (blk_parent_ap[bx] == y) return bx;
    if (blk_parent_ap[by] == x) return by;
    return -1;
}

int GIR::commonBlock(int x, int y) const {
    int bx = findBlk(v_block[x]), by = findBlk(v_block[y]);
    if (bx == by) return bx;
    if (blk_parent_ap[bx] == y) return bx;
    if (blk_parent_ap[by] == x) return by;
    return -1;
}

template <class BaseCnt>
void GIR::installBlocks(std::vector<std::vector<int>>& blocks, BaseCnt base_cnt_fn, int p, int reuse_id) {
    ++block_version;
    const size_t k = blocks.size();
    std::vector<int> ids(k);
    for (size_t i = 0; i < k; ++i) {
        ids[i] = (i == 0 && reuse_id >= 0) ? reuse_id : newBlock();
        blk_uf[ids[i]] = ids[i];
        blk_parent_ap[ids[i]] = -1;
        blk_verts[ids[i]] = std::move(blocks[i]);
    }
    // Occurrence counts inside the region (scr_low reused as a counter).
    int eb = nextEpochB();
    for (size_t i = 0; i < k; ++i)
        for (int x : blk_verts[ids[i]]) {
            if (stamp_b[x] != eb) { stamp_b[x] = eb; scr_low[x] = 0; }
            scr_low[x]++;
        }
    std::unordered_map<int, std::vector<int>> shared;
    std::vector<int> p_blocks;
    for (size_t i = 0; i < k; ++i)
        for (int x : blk_verts[ids[i]]) {
            if (scr_low[x] >= 2) shared[x].push_back((int)i);
            if (x == p) p_blocks.push_back((int)i);
        }
    int ea = nextEpochA();
    for (size_t i = 0; i < k; ++i)
        for (int x : blk_verts[ids[i]])
            if (stamp_a[x] != ea) { stamp_a[x] = ea; v_bcnt[x] = base_cnt_fn(x) + scr_low[x]; }

    // Root each tree and set parent pointers by BFS over shared vertices.
    ea = nextEpochA();  // now: "parent block already assigned"
    std::vector<char> visited(k, 0);
    std::vector<int> q;
    if (p >= 0) {
        stamp_a[p] = ea;  // p keeps its parent outside the region
        for (int i : p_blocks) { visited[i] = 1; blk_parent_ap[ids[i]] = p; q.push_back(i); }
    }
    size_t qi = 0;
    for (size_t start = 0; start <= k; ++start) {
        while (qi < q.size()) {
            int i = q[qi++];
            int id = ids[i];
            for (int x : blk_verts[id]) {
                if (x == blk_parent_ap[id] || stamp_a[x] == ea) continue;
                stamp_a[x] = ea;
                v_block[x] = id;
                auto it = shared.find(x);
                if (it == shared.end()) continue;
                for (int j : it->second)
                    if (!visited[j]) { visited[j] = 1; blk_parent_ap[ids[j]] = x; q.push_back(j); }
            }
        }
        if (start < k && !visited[start]) { visited[start] = 1; q.push_back((int)start); }
    }
    for (size_t i = 0; i < k; ++i) {
        const auto& L = blk_verts[ids[i]];
        for (size_t j = 0; j < L.size(); ++j) if (v_block[L[j]] == ids[i]) blk_pos[L[j]] = (int)j;
    }
}

// Tarjan-Vishkin. T = BFS tree of the SCC's undirected view rooted at r,
// pre = preorder numbers of T, nd = subtree sizes. low(v) / high(v) = the
// smallest / largest pre among v's subtree and the non-tree neighbours of its
// subtree. Tree edges (named by their child vertex) are joined when
//   (1) a non-tree edge {v,w} links two unrelated vertices: e(v) ~ e(w);
//   (2) for tree edge (u,v), u not the root, the subtree of v reaches outside
//       u's subtree: low(v) < pre(u) or high(v) >= pre(u) + nd(u): e(v) ~ e(u).
// The classes are the blocks (as edge sets); a block's vertices are its tree
// edges' endpoints.
void GIR::parallelSccBlocks(int s) {
    const std::vector<int>& S = sccs[s];
    const int n = (int)S.size();

    // Local index of each vertex (in s) via a V-sized scratch map.
    static std::vector<int> loc;
    if (loc.size() < (size_t)V) loc.assign(V, -1);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) loc[S[i]] = i;
    // loc is -1 outside s, so loc[w] >= 0 doubles as the membership test.
    std::vector<int> par(n, -1), lvl_of(n, -1);
    // Parallel BFS from S[0] (local ids).
    std::vector<std::vector<int>> levels{{0}};
    lvl_of[0] = 0;
    par[0] = 0;
    const int T = omp_get_max_threads();
    // Direction-optimizing: a level whose frontier is large is found bottom
    // up (every unvisited vertex looks for a neighbour in the frontier and
    // stops at the first), which on twitter skips most of the edge scans.
    // Either way the result is a BFS tree, which is all Tarjan-Vishkin needs.
    auto undDegree = [&](int x) { return (long long)(sym_beg[x + 1] - sym_beg[x]) + (long long)extra_und_adj[x].size(); };
    long long visited = 1, total_edges = 0;
    #pragma omp parallel for schedule(static) reduction(+ : total_edges)
    for (int i = 0; i < n; ++i) total_edges += undDegree(S[i]);
    while (!levels.back().empty()) {
        const auto& fr = levels.back();
        const int L = (int)levels.size();
        std::vector<std::vector<int>> parts(T);
        long long fr_edges = 0;
        #pragma omp parallel for schedule(static) reduction(+ : fr_edges)
        for (size_t i = 0; i < fr.size(); ++i) fr_edges += undDegree(S[fr[i]]);
        // (Beamer's rule: the frontier's edges against the unexplored ones.)
        const long long rest = (long long)((double)total_edges * (n - visited) / n);
        const bool bottom_up = fr_edges * 14 > rest || fr.size() * 20 > (size_t)(n - visited);
        if (bottom_up) {
            #pragma omp parallel for schedule(dynamic, 1024)
            for (int y = 0; y < n; ++y) {
                if (lvl_of[y] != -1) continue;
                anyUnd(S[y], [&](int w) {
                    int x = loc[w];
                    if (x < 0 || lvl_of[x] != L - 1) return false;
                    par[y] = x;
                    parts[omp_get_thread_num()].push_back(y);
                    return true;
                });
            }
            for (auto& p : parts) for (int y : p) lvl_of[y] = L;  // (after the scan: L-1 tests stay exact)
        } else {
            #pragma omp parallel for schedule(dynamic, 64)
            for (size_t i = 0; i < fr.size(); ++i) {
                const int x = fr[i];
                auto& p = parts[omp_get_thread_num()];
                forUnd(S[x], [&](int w) {
                    int y = loc[w];
                    if (y < 0) return;
                    int exp = -1;
                    if (__atomic_load_n(&lvl_of[y], __ATOMIC_RELAXED) == -1 &&
                        __atomic_compare_exchange_n(&lvl_of[y], &exp, L, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
                        par[y] = x;
                        p.push_back(y);
                    }
                });
            }
        }
        std::vector<int> next;
        for (auto& p : parts) next.insert(next.end(), p.begin(), p.end());
        visited += (long long)next.size();
        levels.push_back(std::move(next));
    }
    levels.pop_back();
    // Subtree sizes (bottom-up), preorder numbers (top-down).
    std::vector<int> nd(n, 1), pre(n, 0), acc(n, 0);
    for (size_t L = levels.size(); L-- > 1;) {
        const auto& lv = levels[L];
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < lv.size(); ++i) __atomic_fetch_add(&nd[par[lv[i]]], nd[lv[i]], __ATOMIC_RELAXED);
    }
    for (size_t L = 1; L < levels.size(); ++L) {
        const auto& lv = levels[L];
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < lv.size(); ++i) {
            int v = lv[i], u = par[v];
            pre[v] = pre[u] + 1 + __atomic_fetch_add(&acc[u], nd[v], __ATOMIC_RELAXED);
        }
    }
    // Union-find over tree edges (id = child vertex; the root has none).
    std::vector<int> uf(n);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) uf[i] = i;
    auto find = [&](int x) {
        while (true) {
            int p = __atomic_load_n(&uf[x], __ATOMIC_RELAXED);
            if (p == x) return x;
            int gp = __atomic_load_n(&uf[p], __ATOMIC_RELAXED);
            if (gp != p) __atomic_compare_exchange_n(&uf[x], &p, gp, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED);
            x = gp;
        }
    };
    auto unite = [&](int a, int b) {
        while (true) {
            a = find(a);
            b = find(b);
            if (a == b) return;
            if (a < b) std::swap(a, b);
            int exp = a;
            if (__atomic_compare_exchange_n(&uf[a], &exp, b, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return;
        }
    };
    // low / high: own non-tree neighbours, then aggregated bottom-up.
    // Rule (1) of the union-find (a non-tree edge between two unrelated
    // vertices joins their tree edges) needs only pre and nd, so it runs in
    // the same adjacency scan.
    // (pre, nd, par packed: one cache line per neighbour in this scan.)
    // (pre and subtree size packed: 8 bytes per neighbour; the parent is read
    // only for a neighbour inside v's subtree, which a tree child must be.)
    struct TreeNode { int pre, nd; };
    std::vector<TreeNode> tn(n);
    #pragma omp parallel for schedule(static)
    for (int v = 0; v < n; ++v) tn[v] = {pre[v], nd[v]};
    std::vector<int> low(n), high(n);
    #pragma omp parallel for schedule(dynamic, 256)
    for (int v = 0; v < n; ++v) {
        const int pv = par[v], pre_v = pre[v], end_v = pre[v] + nd[v];
        int lo = pre_v, hi = pre_v;
        forUnd(S[v], [&](int w) {
            const int y = loc[w];
            if (y < 0) return;
            const TreeNode t = tn[y];
            const bool below = pre_v < t.pre && t.pre < end_v;
            if (y == pv || (below && par[y] == v)) return;  // tree edge
            lo = std::min(lo, t.pre);
            hi = std::max(hi, t.pre);
            const bool related = below || y == v || (t.pre <= pre_v && pre_v < t.pre + t.nd);
            if (y > v && !related) unite(v, y);
        });
        low[v] = lo;
        high[v] = hi;
    }
    for (size_t L = levels.size(); L-- > 1;) {
        const auto& lv = levels[L];
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < lv.size(); ++i) {
            int v = lv[i], u = par[v];
            int lo = low[v], hi = high[v];
            int cur = __atomic_load_n(&low[u], __ATOMIC_RELAXED);
            while (lo < cur && !__atomic_compare_exchange_n(&low[u], &cur, lo, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {}
            cur = __atomic_load_n(&high[u], __ATOMIC_RELAXED);
            while (hi > cur && !__atomic_compare_exchange_n(&high[u], &cur, hi, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {}
        }
    }
    // Rule (2) for tree edge (par[v], v).
    #pragma omp parallel for schedule(static)
    for (int v = 1; v < n; ++v) {
        int u = par[v];
        if (u != 0 && (low[v] < pre[u] || high[v] >= pre[u] + nd[u])) unite(v, u);
    }
    // Blocks: tree edge e(v) contributes v and par[v] to its class.
    std::vector<std::pair<int, int>> mem;  // (class, local vertex)
    mem.resize(2 * (size_t)(n - 1));
    std::vector<int> cls(n, -1);
    #pragma omp parallel for schedule(static)
    for (int v = 1; v < n; ++v) {
        int c = find(v);
        cls[v] = c;
        mem[2 * (size_t)(v - 1)] = {c, v};
        mem[2 * (size_t)(v - 1) + 1] = {c, par[v]};
    }
    __gnu_parallel::sort(mem.begin(), mem.end());
    mem.erase(std::unique(mem.begin(), mem.end()), mem.end());
    // Rooted block-cut forest straight from T: every vertex x but the root
    // lies in the block of its own tree edge as a non-head member, so that
    // block is x's parent block; a block's parent articulation point is its
    // head, the tree parent of its earliest tree edge (smallest pre). Blocks
    // headed by the root hang off it, except one, which is the forest root.
    std::vector<size_t> start;
    for (size_t i = 0; i < mem.size(); ++i) if (i == 0 || mem[i].first != mem[i - 1].first) start.push_back(i);
    const size_t k = start.size();
    start.push_back(mem.size());
    std::vector<int> bid(n, -1);  // class root (local vertex) -> block index
    for (size_t b = 0; b < k; ++b) bid[mem[start[b]].first] = (int)b;
    std::vector<long long> headkey(k, LLONG_MAX);
    #pragma omp parallel for schedule(static)
    for (int v = 1; v < n; ++v) {
        long long key = ((long long)pre[v] << 32) | (unsigned)v;
        long long* p = &headkey[bid[cls[v]]];
        long long cur = __atomic_load_n(p, __ATOMIC_RELAXED);
        while (key < cur && !__atomic_compare_exchange_n(p, &cur, key, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {}
    }
    std::vector<int> ids(k);
    for (size_t b = 0; b < k; ++b) ids[b] = newBlock();
    int root_block = -1;
    for (size_t b = 0; b < k && root_block < 0; ++b)
        if (par[(int)(headkey[b] & 0xffffffffLL)] == 0) root_block = (int)b;
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t b = 0; b < k; ++b) {
        const int id = ids[b];
        auto& L = blk_verts[id];
        L.clear();
        for (size_t i = start[b]; i < start[b + 1]; ++i) L.push_back(S[mem[i].second]);
        const int head = par[(int)(headkey[b] & 0xffffffffLL)];
        blk_uf[id] = id;
        blk_parent_ap[id] = (int)b == root_block ? -1 : S[head];
    }
    #pragma omp parallel for schedule(static)
    for (int v = 0; v < n; ++v) {
        v_block[S[v]] = ids[v == 0 ? root_block : bid[cls[v]]];
        v_bcnt[S[v]] = 0;
    }
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < mem.size(); ++i) __atomic_fetch_add(&v_bcnt[S[mem[i].second]], 1, __ATOMIC_RELAXED);
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t b = 0; b < k; ++b) {
        const auto& L = blk_verts[ids[b]];
        for (size_t j = 0; j < L.size(); ++j) if (v_block[L[j]] == ids[b]) blk_pos[L[j]] = (int)j;
    }
    ++block_version;
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) loc[S[i]] = -1;
}

void GIR::rebuildSccBlocks(int s) {
    if (blocks_on_demand) { for (int x : sccs[s]) markBlocksDirty(x); return; }
    if (sccs[s].size() >= par_bcc_min && sccs[s].size() > 1) { parallelSccBlocks(s); return; }
    std::vector<std::vector<int>> blocks;
    hopcroftTarjan(sccs[s], [&](int, int w) { return node_scc[w] == s; },
                   [&](std::vector<int>& b) { blocks.push_back(std::move(b)); });
    installBlocks(blocks, [](int) { return 0; }, -1, -1);
}

// Blocks on demand: the SCCs holding a dirty vertex get their blocks rebuilt.
// Every SCC whose blocks changed since the last rebuild holds one (updates mark
// their endpoints, splits mark every vertex they move), and so does every
// SCC descended from it, so the old blocks of those SCCs touch no other SCC:
// retire them (reached through v_block of their members), reset the members'
// block fields, and rebuild each SCC.
void GIR::ensureBlocks() {
    if (blk_dirty_v.empty()) return;
    std::vector<int> D;
    {
        std::unordered_set<int> seen;
        for (int x : blk_dirty_v) if (seen.insert(node_scc[x]).second) D.push_back(node_scc[x]);
        std::vector<int>().swap(blk_dirty_v);
    }
    for (int s : D)
        for (int x : sccs[s]) {
            if (v_block[x] < 0) continue;
            int b = findBlk(v_block[x]);
            if (blk_verts[b].empty() && blk_parent_ap[b] == -2) continue;  // already retired
            std::vector<int>().swap(blk_verts[b]);
            blk_parent_ap[b] = -2;  // (marks it retired)
            free_blk_ids.push_back(b);
        }
    for (int s : D) for (int x : sccs[s]) { v_block[x] = -1; v_bcnt[x] = 0; }
    const bool od = blocks_on_demand;
    blocks_on_demand = false;  // (rebuildSccBlocks must not just mark)
    for (int s : D) rebuildSccBlocks(s);
    blocks_on_demand = od;
}

void GIR::buildBCT() {
    std::cout << "=== Build BCT (per SCC) ===" << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    v_block.assign(V, -1);
    v_bcnt.assign(V, 0);
    blk_uf.clear();
    blk_parent_ap.clear();
    blk_verts.clear();
    blk_pos.assign(V, 0);
    // Large SCCs in parallel (Tarjan-Vishkin), the rest with one sequential
    // Hopcroft-Tarjan pass.
    std::vector<int> small;
    std::vector<int> big;
    for (int s = 0; s < (int)sccs.size(); ++s) {
        if (sccs[s].size() >= par_bcc_min && sccs[s].size() > 1) big.push_back(s);
        else small.insert(small.end(), sccs[s].begin(), sccs[s].end());
    }
    std::vector<std::vector<int>> blocks;
    hopcroftTarjan(small, [&](int v, int w) { return node_scc[v] == node_scc[w]; },
                   [&](std::vector<int>& b) { blocks.push_back(std::move(b)); });
    size_t nb = blocks.size();
    installBlocks(blocks, [](int) { return 0; }, -1, -1);
    for (int s : big) parallelSccBlocks(s);
    std::chrono::duration<double> dt = std::chrono::high_resolution_clock::now() - t0;
    std::cout << "Built per-SCC BC trees (" << nb << " blocks over " << sccs.size()
              << " SCCs) in " << dt.count() << " seconds\n";
}

void GIR::repairBlock(int b) {
    if (blocks_on_demand) { if (!blk_verts[findBlk(b)].empty()) markBlocksDirty(blk_verts[findBlk(b)][0]); return; }
    ++block_version;
    b = findBlk(b);
    std::vector<int> verts = blk_verts[b];
    int ea = nextEpochA();
    for (int x : verts) stamp_a[x] = ea;
    std::vector<std::vector<int>> sub;
    // Two vertices share at most one block, so an edge with both endpoints
    // in b is an edge of b.
    hopcroftTarjan(verts, [&](int, int w) { return stamp_a[w] == ea; },
                   [&](std::vector<int>& blk) { sub.push_back(std::move(blk)); });
    if (sub.size() == 1 && sub[0].size() == verts.size()) {
        blk_verts[b] = std::move(sub[0]);
        for (size_t j = 0; j < blk_verts[b].size(); ++j) {
            int x = blk_verts[b][j];
            if (findBlk(v_block[x]) == b) blk_pos[x] = (int)j;
        }
        return;  // still biconnected
    }
    int p = blk_parent_ap[b];
    installBlocks(sub, [&](int x) { return v_bcnt[x] - 1; }, p, b);
}

bool GIR::blockPath(int b, int s, int t, const std::unordered_set<int>& banned, bool skip_direct,
                    long long budget, std::vector<int>& interior, std::vector<int>* ordered) {
    interior.clear();
    int ea = nextEpochA(), eb = nextEpochB();
    stamp_a[s] = ea; par_a[s] = -1;
    stamp_b[t] = eb; par_b[t] = -1;
    std::vector<int> qa{s}, qb{t};
    size_t ia = 0, ib = 0;
    long long work_a = 0, work_b = 0;
    int meet_a = -1, meet_b = -1;  // edge meet_a (a-side) -- meet_b (b-side)
    auto ok = [&](int w) { return !banned.count(w) && inBlock(w, b); };
    while (ia < qa.size() && ib < qb.size() && meet_a < 0 && work_a + work_b < budget) {
        if (work_a + undDeg(qa[ia]) <= work_b + undDeg(qb[ib])) {
            int x = qa[ia++];
            forUnd(x, [&](int w) {
                if (meet_a >= 0) return;
                ++work_a;
                if (skip_direct && x == s && w == t) return;
                if (stamp_b[w] == eb) { meet_a = x; meet_b = w; return; }
                if (stamp_a[w] == ea || !ok(w)) return;
                stamp_a[w] = ea; par_a[w] = x; qa.push_back(w);
            });
        } else {
            int x = qb[ib++];
            forUnd(x, [&](int w) {
                if (meet_a >= 0) return;
                ++work_b;
                if (skip_direct && x == t && w == s) return;
                if (stamp_a[w] == ea) { meet_a = w; meet_b = x; return; }
                if (stamp_b[w] == eb || !ok(w)) return;
                stamp_b[w] = eb; par_b[w] = x; qb.push_back(w);
            });
        }
    }
    if (meet_a < 0) return false;
    for (int x = meet_a; x != s; x = par_a[x]) interior.push_back(x);
    const size_t na = interior.size();
    for (int x = meet_b; x != t; x = par_b[x]) interior.push_back(x);
    if (ordered) {
        ordered->assign(1, s);
        for (size_t i = na; i-- > 0;) ordered->push_back(interior[i]);
        for (size_t i = na; i < interior.size(); ++i) ordered->push_back(interior[i]);
        ordered->push_back(t);
    }
    return true;
}

// Open-addressing int -> int map, cleared in O(used); per-thread scratch for
// the parallel certificates.
struct GIR::LocalSearch {
    struct Map {
        std::vector<int> key, val, used;
        explicit Map(int cap) : key(cap, -1), val(cap) {}
        int* find(int k) {
            size_t m = key.size() - 1, h = ((unsigned)k * 2654435761u) & m;
            while (key[h] != -1) { if (key[h] == k) return &val[h]; h = (h + 1) & m; }
            return nullptr;
        }
        bool put(int k, int v) {  // false if full
            if (used.size() * 2 >= key.size()) return false;
            size_t m = key.size() - 1, h = ((unsigned)k * 2654435761u) & m;
            while (key[h] != -1) { if (key[h] == k) { val[h] = v; return true; } h = (h + 1) & m; }
            key[h] = k; val[h] = v; used.push_back((int)h);
            return true;
        }
        void clear() { for (int h : used) key[h] = -1; used.clear(); }
    };
    // Small tables first (cache-resident); a search that hits their limit or
    // the work budget is retried once with the large ones (hub-heavy graphs).
    Map sa{1 << 15}, sb{1 << 15}, ban{1 << 9};
    std::unique_ptr<Map> la, lb;
    std::vector<int> qa, qb, path;
    static LocalSearch& get() { static thread_local LocalSearch ls; return ls; }

    bool ok(const GIR& g, int w, SearchMode mode, int filt, const std::unordered_set<int>* ext) {
        if (mode != SearchMode::UndBlock) { if (filt >= 0 && g.node_scc[w] != filt) return false; }
        else if (mode == SearchMode::UndBlock) {
            if (g.findBlkRO(g.v_block[w]) != filt && g.blk_parent_ap[filt] != w) return false;
        }
        return !ban.find(w) && !(ext && ext->count(w));
    }

    // Two common neighbours w1 != w2 of s and t (base edges, alive, allowed)
    // give the disjoint paths s-w1-t and s-w2-t.
    // Number of such common neighbours found (0, 1 or 2); the first one is
    // left in `common`.
    int common = -1;
    int twoCommon(const GIR& g, int s, int t, SearchMode mode, int filt, const std::unordered_set<int>* ext) {
        const int* base = g.sym_adj.data();
        const int *a0 = base + g.sym_beg[s], *a1 = base + g.sym_beg[s + 1];
        const int *b0 = base + g.sym_beg[t], *b1 = base + g.sym_beg[t + 1];
        if ((a1 - a0) > (b1 - b0)) { std::swap(a0, b0); std::swap(a1, b1); }
        const bool chk = g.n_dead_und != 0;
        int found = 0, budget = 4096;
        for (const int* p = a0; p < a1 && budget-- > 0; ++p) {
            int w = *p;
            if (w == s || w == t) continue;
            const int* q = std::lower_bound(b0, b1, w);
            if (q == b1 || *q != w) continue;
            if (chk && (bitGet(g.dead_sym, p - base) || bitGet(g.dead_sym, q - base))) continue;
            if (!ok(g, w, mode, filt, ext)) continue;
            if (found == 0) common = w;
            if (++found == 2) return 2;
        }
        return found;
    }

    // Bidirectional search u ~> v (directed: out-edges from u, in-edges into
    // v). Vertices in `ban` or *ext are skipped; skip_direct forbids the edge
    // u-v. On success the interior vertices of the path are left in `path`.
    // Vertices above GIR::search_hub_deg are not expanded (on twitter one
    // hub's list has millions of entries). They become targets instead:
    // every vertex the other side reaches is tested for an edge to them by
    // binary search in the sorted list.
    std::vector<int> big_a, big_b;

    bool search(const GIR& g, int u, int v, SearchMode mode, int filt, bool skip_direct,
                const std::unordered_set<int>* ext) {
        if (u == v) { path.clear(); return true; }
        int r = run(g, u, v, mode, filt, skip_direct, ext, sa, sb, 1 << 17);
        if (r >= 0) return r == 1;
        if (!la) { la.reset(new Map(1 << 18)); lb.reset(new Map(1 << 18)); }
        return run(g, u, v, mode, filt, skip_direct, ext, *la, *lb, 1 << 20) == 1;
    }

    // 1 = path found (interior in `path`), 0 = none, -1 = gave up (limits).
    int run(const GIR& g, int u, int v, SearchMode mode, int filt, bool skip_direct,
            const std::unordered_set<int>* ext, Map& ma, Map& mb, const long long budget) {
        ma.clear(); mb.clear();
        ma.put(u, -1); mb.put(v, -1);
        qa.assign(1, u); qb.assign(1, v);
        big_a.clear(); big_b.clear();
        size_t ia = 0, ib = 0, ja = 0, jb = 0;  // queue heads; next target to expand
        int meet_a = -1, meet_b = -1;
        bool full = false;
        long long work_a = 0, work_b = 0;
        const bool dir = mode == SearchMode::Directed;
        auto allowed = [&](int w) { return ok(g, w, mode, filt, ext); };
        auto deg = [&](int x, bool side_a) -> long long {
            if (!dir) return (g.sym_beg[x + 1] - g.sym_beg[x]) + (long long)g.extra_und_adj[x].size();
            if (side_a) return (long long)(g.g_directed->fw_beg_pos[x + 1] - g.g_directed->fw_beg_pos[x]) + (long long)g.extra_out_adj[x].size();
            return (long long)(g.g_directed->bw_beg_pos[x + 1] - g.g_directed->bw_beg_pos[x]) + (long long)g.extra_in_adj[x].size();
        };
        // a ~> b edge between an a-side vertex and a b-side vertex.
        auto edge = [&](int a, int b) {
            if (skip_direct && a == u && b == v) return false;
            return dir ? g.isDirectedEdgePresent(a, b) : g.isUndirectedEdgePresent(a, b);
        };
        // A newly reached vertex w on one side: does it touch a target of the other side?
        auto hits_big = [&](int w, bool side_a) -> int {
            const auto& L = side_a ? big_b : big_a;
            for (int t : L) if (side_a ? edge(w, t) : edge(t, w)) return t;
            return -1;
        };
        // A side is exhausted only when its queue is empty and all its
        // targets have been expanded too (targets are expanded last).
        auto alive_a = [&]() { return ia < qa.size() || ja < big_a.size(); };
        auto alive_b = [&]() { return ib < qb.size() || jb < big_b.size(); };
        while (alive_a() && alive_b() && meet_a < 0 && !full && work_a + work_b < budget) {
            const bool side_a = work_a <= work_b;
            int x;
            bool expand_big = false;
            if (side_a) { if (ia < qa.size()) x = qa[ia++]; else { x = big_a[ja++]; expand_big = true; } }
            else { if (ib < qb.size()) x = qb[ib++]; else { x = big_b[jb++]; expand_big = true; } }
            if (!expand_big && deg(x, side_a) > g.search_hub_deg && (side_a ? big_a : big_b).size() < 32) {
                // Target instead of expanding: check what the other side already reached.
                (side_a ? big_a : big_b).push_back(x);
                auto& other = side_a ? mb : ma;
                // (Vertices the other side reaches from now on are checked as they come.)
                if (other.used.size() <= 4096) for (int h : other.used) {
                    int y = other.key[h];
                    if (side_a ? edge(x, y) : edge(y, x)) {
                        if (side_a) { meet_a = x; meet_b = y; } else { meet_a = y; meet_b = x; }
                        break;
                    }
                }
                (side_a ? work_a : work_b) += 64;
                continue;
            }
            auto visit = [&](int w) {
                if (meet_a >= 0 || full) return;
                if (++(side_a ? work_a : work_b) + (side_a ? work_b : work_a) >= budget) { full = true; return; }
                if (side_a) {
                    if (skip_direct && x == u && w == v) return;
                    if (mb.find(w)) { meet_a = x; meet_b = w; return; }
                    if (ma.find(w) || !allowed(w)) return;
                    if (!ma.put(w, x)) { full = true; return; }
                    qa.push_back(w);
                    int t = hits_big(w, true);
                    if (t >= 0) { meet_a = w; meet_b = t; }
                } else {
                    if (skip_direct && x == v && w == u) return;
                    if (ma.find(w)) { meet_a = w; meet_b = x; return; }
                    if (mb.find(w) || !allowed(w)) return;
                    if (!mb.put(w, x)) { full = true; return; }
                    qb.push_back(w);
                    int t = hits_big(w, false);
                    if (t >= 0) { meet_a = t; meet_b = w; }
                }
            };
            // Early-exit scans: stop as soon as the sides meet or a limit hits
            // (a hub's list can have millions of entries).
            auto step = [&](int w) { visit(w); return meet_a >= 0 || full; };
            if (dir) { if (side_a) g.anyOut(x, step); else g.anyIn(x, step); }
            else g.anyUnd(x, step);
        }
        // Otherwise one side was exhausted: no path within the restrictions.
        if (meet_a < 0) return (full || work_a + work_b >= budget) ? -1 : 0;
        path.clear();
        for (int x = meet_a; x != u; x = *ma.find(x)) path.push_back(x);
        std::reverse(path.begin(), path.end());  // (interior in order u -> v: callers may re-check it edge by edge)
        for (int x = meet_b; x != v; x = *mb.find(x)) path.push_back(x);
        return 1;
    }
};

bool GIR::reachRO(int u, int v, int scc, std::vector<int>* path) const {
    int mid = -1;
    if (u != v && twoHop(u, v, scc, &mid)) { if (path) path->assign(1, mid); return true; }
    auto& ls = LocalSearch::get();
    ls.ban.clear();
    if (!ls.search(*this, u, v, SearchMode::Directed, scc, false, nullptr)) return false;
    if (path) *path = ls.path;
    return true;
}
void GIR::pickHubCandidates() {
    if (!hub_by_in.empty()) return;
    const int K = 64;
    // Top K by degree: per-thread top K, then merged (one pass over V).
    auto top = [&](const index_t* beg) {
        const int T = omp_get_max_threads();
        std::vector<std::vector<std::pair<long long, int>>> part(T);
        #pragma omp parallel
        {
            auto& h = part[omp_get_thread_num()];
            #pragma omp for schedule(static)
            for (int x = 0; x < V; ++x) {
                long long d = beg[x + 1] - beg[x];
                if ((int)h.size() < K) { h.push_back({d, x}); std::push_heap(h.begin(), h.end(), std::greater<>()); }
                else if (d > h.front().first) {
                    std::pop_heap(h.begin(), h.end(), std::greater<>());
                    h.back() = {d, x};
                    std::push_heap(h.begin(), h.end(), std::greater<>());
                }
            }
        }
        std::vector<std::pair<long long, int>> all;
        for (auto& h : part) all.insert(all.end(), h.begin(), h.end());
        std::sort(all.begin(), all.end(), std::greater<>());
        if ((int)all.size() > K) all.resize(K);
        std::vector<int> r;
        for (auto& p : all) r.push_back(p.second);
        return r;
    };
    hub_by_in = top(g_directed->bw_beg_pos);
    hub_by_out = top(g_directed->fw_beg_pos);
}

bool GIR::hubsCover(int a, int b, size_t C) const {
    const long long din = g_directed->bw_beg_pos[a + 1] - g_directed->bw_beg_pos[a];
    const long long dout = g_directed->fw_beg_pos[b + 1] - g_directed->fw_beg_pos[b];
    return din * 128 >= (long long)C && dout * 128 >= (long long)C;
}

int GIR::hubProbe(int x, int h, bool fwd, int scc, long long budget, std::vector<int>* path, bool big) const {
    auto hit = [&](int w) { return w == h || (fwd ? isDirectedEdgePresent(w, h) : isDirectedEdgePresent(h, w)); };
    if (hit(x)) { if (path) path->clear(); return 1; }
    auto& ls = LocalSearch::get();
    if (big && !ls.la) { ls.la.reset(new LocalSearch::Map(1 << 18)); ls.lb.reset(new LocalSearch::Map(1 << 18)); }
    auto& m = big ? *ls.la : ls.sa;
    m.clear();
    m.put(x, -1);  // (values: BFS parent)
    ls.qa.assign(1, x);
    long long work = 0;
    int r = 0;
    for (size_t i = 0; i < ls.qa.size() && r == 0; ++i) {
        const int p = ls.qa[i];
        auto visit = [&](int w) {
            if (++work >= budget) { r = -1; return true; }
            if (node_scc[w] != scc) return false;
            if (hit(w)) {
                r = 1;
                if (path) { path->assign(1, w); for (int y = p; y != -1; y = *m.find(y)) path->push_back(y); }
                return true;
            }
            if (!m.find(w)) {
                if (!m.put(w, p)) { r = -1; return true; }
                ls.qa.push_back(w);
            }
            return false;
        };
        if (fwd) anyOut(p, visit); else anyIn(p, visit);
    }
    return r;
}

bool GIR::reachSccUndRO(int x, int scc) const {
    if (node_scc[x] == scc) return true;
    auto& ls = LocalSearch::get();
    auto& m = ls.sa;
    m.clear();
    m.put(x, -1);
    ls.qa.assign(1, x);
    long long work = 0;
    bool found = false, stop = false;
    for (size_t i = 0; i < ls.qa.size() && !found && !stop; ++i) {
        anyUnd(ls.qa[i], [&](int w) {
            if (++work >= (1 << 16)) { stop = true; return true; }
            if (node_scc[w] == scc) { found = true; return true; }
            if (!m.find(w)) {
                if (!m.put(w, 0)) { stop = true; return true; }  // table full: give up
                ls.qa.push_back(w);
            }
            return false;
        });
    }
    return found;
}

bool GIR::reachUndRO(int u, int v) const {
    auto& ls = LocalSearch::get();
    ls.ban.clear();
    return ls.search(*this, u, v, SearchMode::UndScc, -1, false, nullptr);
}

// Two internally disjoint undirected s-t paths avoiding *banned, restricted
// to SCC filt (UndScc) or block filt (UndBlock). Greedy: one path, then a
// second one avoiding the first's interior (and the direct edge if the first
// path is that edge), so false is not a proof that no two paths exist.
bool GIR::twoPathsRO(int s, int t, SearchMode mode, int filt, const std::unordered_set<int>* banned) const {
    auto& ls = LocalSearch::get();
    ls.ban.clear();
    const int k = ls.twoCommon(*this, s, t, mode, filt, banned);
    if (k >= 2) return true;
    if (k == 1) {
        // s-w-t is one path; the edge s-t (if any) is a second, otherwise
        // search for one that avoids w.
        if (isUndirectedEdgePresent(s, t)) return true;
        ls.ban.put(ls.common, 0);
        return ls.search(*this, s, t, mode, filt, false, banned);
    }
    if (!ls.search(*this, s, t, mode, filt, false, banned)) return false;
    const bool direct = ls.path.empty();
    for (int x : ls.path) if (!ls.ban.put(x, 0)) return false;
    return ls.search(*this, s, t, mode, filt, direct, banned);
}

// The directed u ~> v path is also an undirected u-v path inside s, so the
// block certificate needs one more search: a second path avoiding its
// interior. Two internally disjoint u-v paths inside s put u and v on a
// common cycle of s, i.e. in one block.
int GIR::certifyKeep(int u, int v, int s, bool need_block) const {
    if (blocks_on_demand) need_block = false;
    auto& ls = LocalSearch::get();
    std::vector<int> p1;
    int w = -1;
    if (twoHop(u, v, s, &w)) {
        p1.assign(1, w);
    } else {
        ls.ban.clear();
        if (!ls.search(*this, u, v, SearchMode::Directed, s, false, nullptr)) return 0;
        p1 = ls.path;
    }
    if (!need_block) return blocks_on_demand ? 2 : 1;
    ls.ban.clear();
    if (ls.twoCommon(*this, u, v, SearchMode::UndScc, s, nullptr) >= 2) return 2;
    for (int x : p1) if (!ls.ban.put(x, 0)) return 1;
    if (ls.search(*this, u, v, SearchMode::UndScc, s, false, nullptr)) return 2;
    // The greedy choice of the first path can block the second one; retry
    // once with a shortest undirected path as the first.
    return twoPathsRO(u, v, SearchMode::UndScc, s, nullptr) ? 2 : 1;
}

bool GIR::twoDisjointPaths(int b, int s, int t, const std::unordered_set<int>& banned, long long budget) {
    // Cheap certificate first: two common neighbours w1 != w2 of s and t
    // inside the block give the disjoint paths s-w1-t and s-w2-t.
    if (!sym_beg.empty()) {
        const int *a0 = sym_adj.data() + sym_beg[s], *a1 = sym_adj.data() + sym_beg[s + 1];
        const int *b0 = sym_adj.data() + sym_beg[t], *b1 = sym_adj.data() + sym_beg[t + 1];
        if ((a1 - a0) > (b1 - b0)) { std::swap(a0, b0); std::swap(a1, b1); }
        const bool chk = n_dead_und != 0;
        const int* base = sym_adj.data();
        int found = 0, budget2 = 4096;
        for (const int* p = a0; p < a1 && budget2-- > 0 && found < 2; ++p) {
            int w = *p;
            if (w == s || w == t || banned.count(w)) continue;
            const int* q = std::lower_bound(b0, b1, w);
            if (q == b1 || *q != w) continue;
            if (chk && (bitGet(dead_sym, p - base) || bitGet(dead_sym, q - base))) continue;
            if (!inBlock(w, b)) continue;
            ++found;
        }
        if (found >= 2) return true;
    }
    std::vector<int> p1, p2;
    if (!blockPath(b, s, t, banned, false, budget, p1)) return false;
    std::unordered_set<int> banned2(banned);
    banned2.insert(p1.begin(), p1.end());
    // If the first path is the direct edge, the second one must not reuse it.
    return blockPath(b, s, t, banned2, p1.empty(), budget, p2);
}

bool GIR::carveChunk(int b, const std::vector<int>& chunk, int cut, const std::unordered_set<int>& gone) {
    const int p = blk_parent_ap[b];
    for (int x : chunk) if (x != cut && x == p) return false;
    ++block_version;
    ++carve_count;
    int ea = nextEpochA();
    for (int x : chunk) stamp_a[x] = ea;
    // The chunk's vertices other than the cut only have b-edges inside the
    // chunk or to departing vertices, so its own blocks come from its induced
    // subgraph.
    std::vector<std::vector<int>> sub;
    hopcroftTarjan(chunk, [&](int, int w) { return stamp_a[w] == ea && !gone.count(w); },
                   [&](std::vector<int>& blk) { sub.push_back(std::move(blk)); });
    for (int x : chunk) if (x != cut) blkListRemove(b, x);
    // Base counts: chunk vertices lose block b; the cut keeps it.
    installBlocks(sub, [&](int y) { return y == cut ? v_bcnt[y] : v_bcnt[y] - 1; }, cut, -1);
    return true;
}

bool GIR::peelAndCertify(int b, const std::unordered_set<int>& gone, const std::vector<int>& seeds,
                         const std::vector<int>& extra_att) {
    ++block_version;
    const long long BUDGET = 1LL << 22;
    std::unordered_set<int> out(gone);  // gone + peeled: not in the core
    std::vector<std::pair<int, int>> peeled;  // (y, z): y hangs off z
    auto in_core = [&](int w) { return !out.count(w) && inBlock(w, b); };
    std::vector<int> q(seeds);
    for (size_t i = 0; i < q.size(); ++i) {
        int y = q[i];
        if (out.count(y) || !inBlock(y, b)) continue;
        int deg = 0, z = -1;
        anyUnd(y, [&](int w) { if (w != z && in_core(w)) { ++deg; z = w; } return deg >= 2; });
        if (deg == 0) return false;
        if (deg == 1) {
            if (y == blk_parent_ap[b]) return false;  // would re-hang b itself; leave to the slow path
            out.insert(y);
            peeled.push_back({y, z});
            q.push_back(z);
        }
    }
    const size_t core = blk_verts[b].size() - out.size();
    if (core < 2) return false;
    if (core >= 3) {
        // Attachment vertices: core vertices that lost an edge. If c were a
        // cut vertex of the core, every component of core - c would contain
        // one (otherwise c already cut b). Two disjoint paths from att[0] to
        // every other attachment rule out any c != att[0]; c == att[0] would
        // separate two other attachments, which the pairs from att[1] catch.
        std::vector<int> att;
        std::unordered_set<int> seen;
        auto add = [&](int w) { if (in_core(w) && seen.insert(w).second) att.push_back(w); };
        for (int x : out) forUnd(x, add);
        for (int x : extra_att) add(x);
        if (att.size() > (1u << 16)) return false;
        // On a failing pair, the exact Menger search returns the small side of
        // the minimum vertex cut: carve it out as its own subtree hanging off
        // the cut vertex (which then counts as an attachment of the rest) and
        // certify again.
        // The greedy pair checks are read-only and run in parallel; a pair
        // they can't certify goes to the exact Menger search.
        // A pair certified once stays certified after a carve: the carved
        // chunk hangs off one cut vertex, so no simple path between two
        // vertices outside it passes through it.
        std::unordered_set<long long> done;
        for (int iter = 0; iter < 64; ++iter) {
            std::vector<std::pair<int, int>> pairs;
            for (size_t i = 1; i < att.size(); ++i) pairs.push_back({att[0], att[i]});
            for (size_t i = 2; i < att.size(); ++i) pairs.push_back({att[1], att[i]});
            if (!done.empty()) {
                size_t k = 0;
                for (auto& p : pairs) if (!done.count(keyU(p.first, p.second))) pairs[k++] = p;
                pairs.resize(k);
            }
            std::vector<char> okp(pairs.size(), 0);
            #pragma omp parallel for schedule(dynamic, 1) if (pairs.size() >= 4)
            for (size_t i = 0; i < pairs.size(); ++i)
                okp[i] = twoPathsRO(pairs[i].first, pairs[i].second, SearchMode::UndBlock, b, &out);
            for (size_t i = 0; i < pairs.size(); ++i) if (okp[i]) done.insert(keyU(pairs[i].first, pairs[i].second));
            int bad_s = -1, bad_t = -1, r = 1;
            std::vector<int> chunk;
            int cut = -1;
            bool from_s = false;
            for (size_t i = 0; i < pairs.size() && r == 1; ++i) {
                if (okp[i]) continue;
                bad_s = pairs[i].first; bad_t = pairs[i].second;
                r = mengerSplit(b, bad_s, bad_t, chunk, cut, from_s, BUDGET, &out);
                if (r == 1) done.insert(keyU(bad_s, bad_t));
            }
            if (r == 1) break;  // every pair has two disjoint paths
            if (r < 0 || chunk.size() * 2 > blk_verts[b].size()) return false;
            if (!carveChunk(b, chunk, cut, out)) return false;
            std::unordered_set<int> inchunk(chunk.begin(), chunk.end());
            std::vector<int> natt;
            for (int x : att) if (!inchunk.count(x) || x == cut) natt.push_back(x);
            if (std::find(natt.begin(), natt.end(), cut) == natt.end()) natt.push_back(cut);
            att.swap(natt);
            if (att.size() < 2) break;
            if (iter == 63) return false;
        }
    }
    auto& bv = blk_verts[b];
    for (int x : out) blkListRemove(b, x);  // O(|out|) via position hints
    (void)bv;
    for (auto& yz : peeled) {
        int nb = newBlock();
        blk_verts[nb] = {yz.first, yz.second};
        blk_parent_ap[nb] = yz.second;
        v_block[yz.first] = nb;
        blk_pos[yz.first] = 0;
        v_bcnt[yz.second]++;
    }
    return true;
}

void GIR::recomputeBlockWithout(int b, const std::unordered_set<int>& gone) {
    std::vector<int> keep;
    for (int y : blk_verts[b]) if (!gone.count(y)) keep.push_back(y);
    int ea = nextEpochA();
    for (int y : keep) stamp_a[y] = ea;
    std::vector<std::vector<int>> sub;
    hopcroftTarjan(keep, [&](int, int w) { return stamp_a[w] == ea; },
                   [&](std::vector<int>& blk) { sub.push_back(std::move(blk)); });
    installBlocks(sub, [&](int y) { return v_bcnt[y] - 1; }, blk_parent_ap[b], b);
}

// Ordered s..t path inside block b (undirected); empty if none within budget.
std::vector<int> GIR::blockPathOrdered(int b, int s, int t, long long budget, const std::unordered_set<int>* banned) {
    static const std::unordered_set<int> none;
    std::vector<int> interior, path;
    if (!blockPath(b, s, t, banned ? *banned : none, false, budget, interior, &path)) return {};
    return path;
}

// Exact Menger test with vertex capacities 1 inside block b, after one s-t
// path P1 is known: look for an augmenting path in the residual network,
// searching from s and from t in lockstep. Returns 1 if a second disjoint
// path exists. Returns 0 if not; then `chunk` is the side of the minimum
// vertex cut closest to the end whose search exhausted first (that end is
// reported in `from_s`), including the cut vertex `cut`. Returns -1 if the
// budget ran out.
long long GIR::undDeg(int x) const {
    return (sym_beg.empty() ? 0LL : (long long)(sym_beg[x + 1] - sym_beg[x])) + (long long)extra_und_adj[x].size();
}

int GIR::mengerSplit(int b, int s, int t, std::vector<int>& chunk, int& cut, bool& from_s, long long budget,
                     const std::unordered_set<int>* banned) {
    std::vector<int> P = blockPathOrdered(b, s, t, budget, banned);
    if (P.size() < 2) { menger_nopath++; return -1; }
    if (blk_verts[b].size() == 2) return 1;  // {s, t} joined by an edge is already a block
    std::unordered_map<int, int> pos;  // index on P (oriented s -> t)
    for (size_t i = 0; i < P.size(); ++i) pos[P[i]] = (int)i;
    const int last = (int)P.size() - 1;
    auto ip = [&](int v) { auto it = pos.find(v); return it == pos.end() ? -1 : it->second; };
    auto interior = [&](int v) { int i = ip(v); return i > 0 && i < last; };
    // Residual network of one unit of flow along P, vertex capacities 1
    // (state = 2*v + out). Forward search from (s,out) follows residual arcs;
    // backward search from (t,in) follows them in reverse. They stop as soon
    // as they share a state (augmenting path => second disjoint path).
    std::unordered_set<long long> F, Bk;
    std::vector<long long> fq{2LL * s + 1}, bq{2LL * t};
    F.insert(fq[0]);
    Bk.insert(bq[0]);
    size_t fi = 0, bi = 0;
    long long wf = 0, wb = 0;
    int exhausted = -1;  // 0 = forward, 1 = backward
    bool met = false;
    // Hubs (degree above search_hub_deg) are deferred instead of expanded:
    // forward (h,out) states in hubF, backward (h,in) states in hubB. A state
    // the other side adds is tested for the residual arc to them, (x,out) ->
    // (y,in) for an edge {x,y} not saturated along P. Deferred states are
    // expanded once their side has nothing else, so exhaustion (and the cut
    // read off it) is exact.
    std::vector<int> hubF, hubB;
    size_t jf = 0, jb = 0;
    auto arc = [&](int x, int y) {
        if (!isUndirectedEdgePresent(x, y)) return false;
        int ix = ip(x), iy = ip(y);
        return !(ix >= 0 && iy >= 0 && iy == ix + 1);
    };
    auto pushF = [&](long long st) {
        if (Bk.count(st)) { met = true; return; }
        if (!F.insert(st).second) return;
        fq.push_back(st);
        if (st & 1) for (int h : hubB) if (arc((int)(st >> 1), h)) { met = true; return; }
    };
    auto pushB = [&](long long st) {
        if (F.count(st)) { met = true; return; }
        if (!Bk.insert(st).second) return;
        bq.push_back(st);
        if (!(st & 1)) for (int h : hubF) if (arc(h, (int)(st >> 1))) { met = true; return; }
    };
    // Cost of expanding a state: the vertex's list for (v,out) forward and
    // (v,in) backward, nothing for the others. Expanding the side that keeps
    // its total lower lets a small side run out before a hub is scanned.
    auto costF = [&](long long st) { return (st & 1) ? undDeg((int)(st >> 1)) : 1LL; };
    auto costB = [&](long long st) { return (st & 1) ? 1LL : undDeg((int)(st >> 1)); };
    while (!met) {
        const bool fq_has = fi < fq.size(), bq_has = bi < bq.size();
        if (!fq_has && jf >= hubF.size()) { exhausted = 0; break; }
        if (!bq_has && jb >= hubB.size()) { exhausted = 1; break; }
        if (wf + wb > budget) { menger_budget++; return -1; }
        // The cheaper next step wins; a side left with only deferred hubs
        // pays the hub's degree, so deferring only postpones a hub (always
        // preferring the other side can drive a huge side into the budget).
        const long long cf = fq_has ? costF(fq[fi]) : undDeg(hubF[jf]);
        const long long cb = bq_has ? costB(bq[bi]) : undDeg(hubB[jb]);
        const bool fwd = wf + cf <= wb + cb;  // (both sides are alive here)
        if (fwd) {
            long long st;
            if (fq_has) {
                st = fq[fi++];
                if ((st & 1) && cf > search_hub_deg && hubF.size() < 32) { hubF.push_back((int)(st >> 1)); wf += 32; continue; }
            } else {
                st = 2LL * hubF[jf++] + 1;
            }
            int v = (int)(st >> 1);
            if (st & 1) {  // (v,out) -> (w,in) except along P; (v,out) -> (v,in) if v interior
                int iv = ip(v);
                forUnd(v, [&](int w) {
                    if (met) return;
                    ++wf;
                    if (!inBlock(w, b) || (banned && banned->count(w))) return;
                    if (iv >= 0 && ip(w) == iv + 1) return;
                    pushF(2LL * w);
                });
                if (interior(v)) pushF(2LL * v);
            } else {       // (v,in) -> (v,out) if not interior; (x_i,in) -> (x_{i-1},out) if on P
                if (!interior(v)) { if (v != t) pushF(2LL * v + 1); }
                else pushF(2LL * P[ip(v) - 1] + 1);
            }
        } else {
            long long st;
            if (bq_has) {
                st = bq[bi++];
                if (!(st & 1) && cb > search_hub_deg && hubB.size() < 32) { hubB.push_back((int)(st >> 1)); wb += 32; continue; }
            } else {
                st = 2LL * hubB[jb++];
            }
            int v = (int)(st >> 1);
            if (!(st & 1)) {  // predecessors of (v,in): (w,out) unless w->v is along P; (v,out) if v interior
                int iv = ip(v);
                forUnd(v, [&](int w) {
                    if (met) return;
                    ++wb;
                    if (!inBlock(w, b) || (banned && banned->count(w))) return;
                    if (iv >= 0 && ip(w) == iv - 1) return;
                    pushB(2LL * w + 1);
                });
                if (interior(v)) pushB(2LL * v + 1);
            } else {          // predecessors of (v,out): (v,in) if not interior; (x_{i+1},in) if v = x_i on P
                int iv = ip(v);
                if (!interior(v)) { if (v != s) pushB(2LL * v); }
                if (iv >= 0 && iv < last) pushB(2LL * P[iv + 1]);
            }
        }
    }
    if (met) return 1;
    // Max flow is 1. The exhausted side is the source side of the min cut
    // closest to its end: forward -> the P vertex with "in" reached but
    // "out" not; backward -> "out" reached (can reach t) but "in" not.
    const auto& Sx = exhausted == 0 ? F : Bk;
    std::unordered_set<int> verts;
    for (long long st : Sx) verts.insert((int)(st >> 1));
    cut = -1;
    for (int i = 1; i < last; ++i) {
        int v = P[i];
        bool in = Sx.count(2LL * v), out = Sx.count(2LL * v + 1);
        if (exhausted == 0 ? (in && !out) : (out && !in)) { cut = v; break; }
    }
    from_s = exhausted == 0;
    if (cut < 0) {
        // Edge cut at the exhausted end's first P arc: that end hangs off the
        // next P vertex by a bridge.
        if (verts.size() == 1 && P.size() > 2) {
            cut = exhausted == 0 ? P[1] : P[last - 1];
            chunk = {exhausted == 0 ? s : t, cut};
            return 0;
        }
        menger_nocut++;
        return -1;
    }
    chunk.assign(verts.begin(), verts.end());
    return 0;
}

// Block b loses everything of `chunk` except `cut`, which becomes the
// articulation point joining b and a new block for the chunk. b keeps its id
// (it is the big part), so only the chunk's vertices are touched.
void GIR::splitChunkOff(int b, const std::vector<int>& chunk, int cut) {
    ++block_version;
    const int p = blk_parent_ap[b];
    std::unordered_set<int> K(chunk.begin(), chunk.end());
    const bool root_in_chunk = p >= 0 && p != cut && K.count(p);
    int nk = newBlock();
    blk_verts[nk] = chunk;
    for (int x : chunk) {
        if (x == cut) continue;
        if (x != p && findBlk(v_block[x]) == b) v_block[x] = nk;
        blkListRemove(b, x);
    }
    for (size_t i = 0; i < chunk.size(); ++i) if (findBlk(v_block[chunk[i]]) == nk) blk_pos[chunk[i]] = (int)i;
    v_bcnt[cut]++;
    if (root_in_chunk) {
        // The chunk holds b's parent AP: it takes over b's place in the tree
        // and b now hangs off the cut vertex.
        blk_parent_ap[nk] = p;
        blk_parent_ap[b] = cut;
        v_block[cut] = nk;
    } else {
        blk_parent_ap[nk] = cut;
    }
}

void GIR::repairBlockAfterDelete(int b, int u, int v, bool skip_quick) {
    if (blocks_on_demand) { markBlocksDirty(u); markBlocksDirty(v); return; }
    b = findBlk(b);
    if (blk_verts[b].size() >= fastpath_min_block) {
        static const std::unordered_set<int> none;
        // Common case: two common neighbours / two disjoint paths at once.
        if (!force_exact_cert && !skip_quick && twoPathsRO(u, v, SearchMode::UndBlock, b, nullptr)) { block_stats[0]++; return; }
        // Otherwise the block minus {u,v} is a chain of blocks between u and
        // v. Peel chain links off whichever end is small, one min cut at a
        // time, until the rest is certified 2-connected.
        int s = u, t = v;
        for (int iter = 0; iter < 256; ++iter) {
            std::vector<int> chunk;
            int cut;
            bool from_s;
            // After a carve, the rest between s and t is usually one block
            // already: the bounded certificate saves the confirming Menger run.
            if (iter > 0 && !force_exact_cert && twoPathsRO(s, t, SearchMode::UndBlock, b, nullptr)) {
                block_stats[1]++;
                return;
            }
            int r = mengerSplit(b, s, t, chunk, cut, from_s, 1LL << 22);
            if (r == 1) { block_stats[1]++; return; }
            // A piece bigger than half a large block is cheaper to rebuild.
            if (r < 0 || (chunk.size() * 2 > blk_verts[b].size() && blk_verts[b].size() > 1024)) break;
            splitChunkOff(b, chunk, cut);
            block_stats[2]++;
            (from_s ? s : t) = cut;
            if (s == t) return;
        }
    }
    block_stats[3]++;
    repairBlock(b);
}

void GIR::mergeBlockPath(int u, int v) {
    if (blocks_on_demand) { markBlocksDirty(u); markBlocksDirty(v); return; }
    ++block_version;
    // Node encoding: block b -> 2b, articulation point a -> 2a+1.
    auto start = [&](int x) -> long long {
        return v_bcnt[x] >= 2 ? 2LL * x + 1 : 2LL * findBlk(v_block[x]);
    };
    auto up = [&](long long n) -> long long {
        if (n & 1) return 2LL * findBlk(v_block[(int)(n >> 1)]);
        int ap = blk_parent_ap[(int)(n >> 1)];
        return ap < 0 ? -1 : 2LL * ap + 1;
    };
    std::vector<long long> px{start(u)}, py{start(v)};
    std::unordered_map<long long, size_t> mx{{px[0], 0}}, my{{py[0], 0}};
    long long L = -1;
    if (px[0] == py[0]) L = px[0];
    while (L < 0) {
        bool moved = false;
        long long n = up(px.back());
        if (n >= 0) {
            moved = true;
            px.push_back(n);
            auto it = my.find(n);
            if (it != my.end()) { L = n; py.resize(it->second + 1); break; }
            mx[n] = px.size() - 1;
        }
        n = up(py.back());
        if (n >= 0) {
            moved = true;
            py.push_back(n);
            auto it = mx.find(n);
            if (it != mx.end()) { L = n; px.resize(it->second + 1); break; }
            my[n] = py.size() - 1;
        }
        if (!moved) {
            std::cerr << "GIR::mergeBlockPath: endpoints in different block-cut trees" << std::endl;
            return;
        }
    }
    std::vector<long long> path = px;  // u-side ... L
    for (size_t i = py.size() - 1; i-- > 0;) path.push_back(py[i]);

    std::vector<int> blocks;
    std::unordered_set<int> pc;
    for (long long n : path) {
        if (n & 1) pc.insert((int)(n >> 1));
        else blocks.push_back((int)(n >> 1));
    }
    int N = blocks[0];
    for (int b : blocks) if (blk_verts[b].size() > blk_verts[N].size()) N = b;
    const int new_parent = (L & 1) ? (int)(L >> 1) : blk_parent_ap[(int)(L >> 1)];

    for (size_t i = 0; i < path.size(); ++i) {
        if (!(path[i] & 1)) continue;
        int c = (int)(path[i] >> 1);
        bool adj_to_N = (i > 0 && path[i - 1] == 2LL * N) || (i + 1 < path.size() && path[i + 1] == 2LL * N);
        if (!adj_to_N) { blk_pos[c] = (int)blk_verts[N].size(); blk_verts[N].push_back(c); }
        if (i > 0 && i + 1 < path.size()) v_bcnt[c]--;  // two of its blocks became one
    }
    for (int b : blocks) {
        if (b == N) continue;
        for (int x : blk_verts[b]) if (!pc.count(x)) { blk_pos[x] = (int)blk_verts[N].size(); blk_verts[N].push_back(x); }
        std::vector<int>().swap(blk_verts[b]);
        blk_uf[b] = N;
    }
    blk_parent_ap[N] = new_parent;
}

void GIR::rerootAt(int b) {
    ++block_version;
    std::vector<int> chain_blocks{findBlk(v_block[b])}, chain_aps{b};
    while (blk_parent_ap[chain_blocks.back()] != -1) {
        int c = blk_parent_ap[chain_blocks.back()];
        chain_aps.push_back(c);
        chain_blocks.push_back(findBlk(v_block[c]));
    }
    for (size_t i = 0; i < chain_blocks.size(); ++i) {
        blk_parent_ap[chain_blocks[i]] = chain_aps[i];
        if (i > 0) v_block[chain_aps[i]] = chain_blocks[i - 1];
    }
}

void GIR::linkTrees(int a, int b) {
    if (blocks_on_demand) { markBlocksDirty(a); markBlocksDirty(b); return; }
    ++block_version;
    int Bb = findBlk(v_block[b]);
    if (blk_verts[Bb].size() == 1) {  // b was an isolated singleton block
        std::vector<int>().swap(blk_verts[Bb]);
        v_bcnt[b] = 0;
    } else {
        rerootAt(b);
    }
    int Ba = findBlk(v_block[a]);
    int N;
    if (blk_verts[Ba].size() == 1) {  // a was an isolated singleton block: grow it
        N = Ba;
        blk_pos[b] = (int)blk_verts[N].size();
        blk_verts[N].push_back(b);
        blk_parent_ap[N] = -1;
    } else {
        N = newBlock();
        blk_verts[N] = {a, b};
        blk_parent_ap[N] = a;
        blk_pos[b] = 1;
        v_bcnt[a]++;
    }
    v_block[b] = N;
    v_bcnt[b]++;
}

void GIR::mergeBlockForests(const std::vector<int>& old_ids, int keep) {
    if (blocks_on_demand) { if (!sccs[keep].empty()) markBlocksDirty(sccs[keep][0]); return; }
    std::unordered_map<int, int> tpar, tsize;
    for (int s : old_ids) { tpar[s] = s; tsize[s] = (int)sccs[s].size(); }
    auto tfind = [&](int s) {
        while (tpar[s] != s) { tpar[s] = tpar[tpar[s]]; s = tpar[s]; }
        return s;
    };
    for (int s : old_ids) {
        if (s == keep) continue;
        for (int x : sccs[s]) {
            forUnd(x, [&](int w) {
                int t = node_scc[w];
                if (t == s || !tpar.count(t)) return;
                int tx = tfind(s), tw = tfind(t);
                if (tx != tw) {
                    if (tsize[tx] >= tsize[tw]) linkTrees(x, w); else linkTrees(w, x);
                    if (tsize[tx] < tsize[tw]) std::swap(tx, tw);
                    tpar[tw] = tx;
                    tsize[tx] += tsize[tw];
                } else if (commonBlock(x, w) < 0) {
                    mergeBlockPath(x, w);
                }
            });
        }
    }
}

// ============================================================================
// Updates
// ============================================================================

void GIR::handleInsert(int u, int v, bool und_new) {
    int su = node_scc[u], sv = node_scc[v];
    if (wcc_label[u] != wcc_label[v]) {
        // Case 1. Relabel the smaller WCC and, if needed, move the smaller
        // side's part of the order. No cycle can close across WCCs and no
        // block changes.
        int lu = wcc_label[u], lv = wcc_label[v];
        bool u_smaller = wcc_size[lu] <= wcc_size[lv];
        if (u_smaller) relabelComponent(u, lu, lv);
        else relabelComponent(v, lv, lu);
        if (insertDagEdgePK(su, sv, /*cross_wcc=*/true) < 0) insertCrossWccEdge(su, sv, u_smaller);
        case_counts[C1_CROSS_WCC]++;
        return;
    }
    if (su != sv) {
        bool merged = insertDagEdgePK(su, sv) > 0;
        case_counts[merged ? C2_SCC_MERGE : C2_DAG_EDGE_NO_CYCLE]++;
        return;
    }
    // Same SCC. If the undirected edge already existed (v->u present) the
    // endpoints are adjacent in the SCC and therefore share a block.
    if (blocks_on_demand) {  // same SCC: only blocks could change
        if (und_new) { markBlocksDirty(u); markBlocksDirty(v); }
        case_counts[C4_NOOP]++;
        return;
    }
    if (!und_new || commonBlock(u, v) >= 0) { case_counts[C4_NOOP]++; return; }
    mergeBlockPath(u, v);
    case_counts[C3_BLOCK_MERGE]++;
}

void GIR::attributeTime(const std::array<long long, NUM_CASES>& before, double secs) {
    for (int c = 0; c < NUM_CASES; ++c)
        if (case_counts[c] != before[c]) { case_seconds[c] += secs; return; }
}

void GIR::insertEdge(int u, int v) {
    sync();
    auto before = case_counts;
    auto t0 = std::chrono::steady_clock::now();
    bool und_new;
    if (!storeInsert(u, v, und_new)) case_counts[C4_NOOP]++;
    else { logChanges({{u, v}}); handleInsert(u, v, und_new); }
    attributeTime(before, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
}

void GIR::deleteEdge(int u, int v) {
    sync();
    auto before = case_counts;
    auto t0 = std::chrono::steady_clock::now();
    if (isDirectedEdgePresent(u, v)) logChanges({{u, v}});
    deleteEdgeImpl(u, v);
    attributeTime(before, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
}

void GIR::deleteEdgeImpl(int u, int v, int hint) {
    bool und_gone;
    if (!storeDelete(u, v, und_gone)) return;
    int su = node_scc[u], sv = node_scc[v];
    if (su != sv) {
        // A condensation edge. SCCs and blocks are unaffected; only the WCC
        // can split, and only if no other u'->v' edge joins these SCCs and no
        // other undirected path joins u and v.
        removeDagEdge(su, sv);
        std::vector<int> side;
        if (dagEdgeExists(su, sv) || reachUndRO(u, v) || undirectedConnected(u, v, side)) {
            case_counts[D_INTER_SCC]++;
            return;
        }
        int old_label = wcc_label[u];
        int new_label = (int)wcc_size.size();
        wcc_size.push_back(0);
        for (int x : side) wcc_label[x] = new_label;
        wcc_size[new_label] = (int)side.size();
        wcc_size[old_label] -= (int)side.size();
        case_counts[D_WCC_SPLIT]++;
        return;
    }
    if (und_gone) { markBlocksDirty(u); markBlocksDirty(v); }  // (on demand only)
    // Inside one SCC. The SCC survives iff u still reaches v: every old
    // path that used u->v can then be rerouted. Bridge status in the
    // undirected view says nothing here (any edge of a directed cycle can
    // break strong connectivity). The WCC cannot split: if the undirected
    // edge disappeared, v->u was absent, so a v ~> u path that avoids u->v
    // still connects the endpoints.
    // Common case in a large block: the SCC survives and the block needs no
    // repair. The bounded searches certify both; anything they can't
    // certify takes the exact path. (Running the two concurrently in a
    // 2-thread OpenMP team was measured slower: the fork costs more than
    // the ~15 us it overlaps.)
    bool ro_failed = hint == 2;
    if (hint == 1 || hint == 3) {  // u ~> v holds: the SCC stays; only the block needs repair
        // (1: the bounded block certificate already failed; 3: not tried)
        if (und_gone) repairBlockAfterDelete(commonBlock(u, v), u, v, hint == 1);
        case_counts[D_INTRA_SCC_KEEP]++;
        return;
    }
    if (und_gone && !force_exact_cert && !ro_failed) {
        int b = commonBlock(u, v);  // (compresses the union-find before the read-only searches)
        if (b >= 0 && blk_verts[b].size() >= fastpath_min_block) {
            const int cert = certifyKeep(u, v, su, true);
            const bool r_ok = cert >= 1, b_ok = cert == 2;
            if (r_ok) {
                if (b_ok) block_stats[0]++;
                else repairBlockAfterDelete(b, u, v);
                case_counts[D_INTRA_SCC_KEEP]++;
                return;
            }
            ro_failed = true;
        }
    }
    std::vector<int> side;
    bool side_fwd = false;
    if ((!ro_failed && reachRO(u, v)) || reachableWithinScc(u, v, su, &side, &side_fwd)) {
        if (und_gone) repairBlockAfterDelete(commonBlock(u, v), u, v);
        case_counts[D_INTRA_SCC_KEEP]++;
        return;
    }
    auto before = split_stats;
    auto t0 = std::chrono::steady_clock::now();
    splitAfterDelete(su, u, v, side, side_fwd);
    double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    for (int i = 0; i < 5; ++i) if (split_stats[i] != before[i]) { split_secs[i] += dt; break; }
    case_counts[D_SCC_SPLIT]++;
}

void GIR::updateEdgeDirected(int u, int v, bool insert) { insert ? insertEdge(u, v) : deleteEdge(u, v); }
void GIR::updateEdgeUndirected(int u, int v, bool insert) { insert ? insertEdge(u, v) : deleteEdge(u, v); }

// Filter-Merge-Update.
//  Filter: per edge only the last operation matters; ops that leave the edge
//          in its current state (cancelling pairs, duplicates) are dropped.
//  Merge:  deletions are grouped by SCC and by block, so each touched SCC
//          gets one split check and each touched block one repair.
//  Update: deletions first (the final graph doesn't depend on order), then
//          insertions through the single-edge path.
template <class Hit>
int GIR::parBfs(int x, int dir, int scc, Hit hit, std::vector<int>* path) {
    if (hit(x)) { if (path) path->clear(); return 1; }
    if (pa_stamp.size() < (size_t)V) { pa_stamp.assign(V, 0); pa_parent.assign(V, -1); }
    if (++pa_epoch == INT_MAX) { std::fill(pa_stamp.begin(), pa_stamp.end(), 0); pa_epoch = 1; }
    const int ep = pa_epoch;
    pa_stamp[x] = ep;
    pa_parent[x] = -1;
    std::vector<int> fr{x}, next;
    std::vector<std::vector<int>> parts(omp_get_max_threads());
    int found = -1;
    while (!fr.empty() && found < 0) {
        #pragma omp parallel for schedule(dynamic, 16)
        for (size_t i = 0; i < fr.size(); ++i) {
            if (__atomic_load_n(&found, __ATOMIC_RELAXED) >= 0) continue;
            const int p = fr[i];
            auto& out = parts[omp_get_thread_num()];
            auto visit = [&](int w) {
                if (scc >= 0 && node_scc[w] != scc) return false;
                int exp = __atomic_load_n(&pa_stamp[w], __ATOMIC_RELAXED);
                if (exp == ep || !__atomic_compare_exchange_n(&pa_stamp[w], &exp, ep, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return false;
                pa_parent[w] = p;
                if (hit(w)) { int none = -1; __atomic_compare_exchange_n(&found, &none, w, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED); return true; }
                out.push_back(w);
                return false;
            };
            if (dir == 0) anyOut(p, visit); else if (dir == 1) anyIn(p, visit); else anyUnd(p, visit);
        }
        next.clear();
        for (auto& q : parts) { next.insert(next.end(), q.begin(), q.end()); q.clear(); }
        fr.swap(next);
    }
    if (found < 0) return 0;
    if (path) { path->clear(); for (int y = found; y != x && y != -1; y = pa_parent[y]) path->push_back(y); }
    return 1;
}

// FMU: the batch's uncertified deletions inside large SCC s, handled at once.
// Most of them split off single vertices that lost their last out-edge (or
// in-edge) inside s. Trim on the graph after the batch: repeatedly take out
// vertices of s with no remaining in-neighbour or no remaining out-neighbour
// in s; call the removed set T.
// Reason in G, the graph with every uncertified deletion still present: the
// structure is valid for G (certified deletions change nothing), so s is
// strongly connected in G. Put back the uncertified deletions inside R = s \ T
// and keep those touching T deleted: no edge between two R vertices is
// missing, so, as in splitAfterDelete, R is strongly connected iff every z in
// Z reaches an anchor a and a reaches every e in E (Z, E: R-vertices with a
// G-edge into / out of T). The probes run in parallel; when one fails, the
// side it exhausted joins T and the check repeats. T's own SCCs and their
// order come from orderPieces, and T is detached in one step. If anything
// fails, nothing has changed and the deletions take the sequential path.
std::vector<int> GIR::bulkTrimSplit(int s, const std::vector<std::pair<int, int>>& dels,
                                    const std::vector<char>& no_reach,
                                    const std::vector<std::pair<int, int>>& all_dels,
                                    std::vector<std::pair<int, int>>& wcc_probe, bool& certified) {
    certified = false;
    const size_t C = sccs[s].size();
    const int tmp = newSccId();
    auto inS = [&](int w) { return node_scc[w] == s; };
    std::vector<int> T, work;
    for (auto& e : dels) { work.push_back(e.first); work.push_back(e.second); }
    bool give_up = false;
    for (size_t i = 0; i < work.size() && !give_up; ++i) {
        int x = work[i];
        if (!inS(x)) continue;
        bool has_out = anyOut(x, [&](int w) { return w != x && inS(w); });
        if (has_out && anyIn(x, [&](int w) { return w != x && inS(w); })) continue;
        node_scc[x] = tmp;
        T.push_back(x);
        if (T.size() * 4 > C || T.size() > (1u << 16) || work.size() > (1u << 22)) { give_up = true; break; }
        forOut(x, [&](int w) { if (inS(w)) work.push_back(w); });
        forIn(x, [&](int w) { if (inS(w)) work.push_back(w); });
    }
    auto undo = [&]() { for (int x : T) node_scc[x] = s; retireScc(tmp); T.clear(); };
    bool any_split = false;
    for (char f : no_reach) any_split |= f != 0;
    if (give_up || (T.empty() && !any_split)) { undo(); return {}; }
    // Deletions inside R that only failed the block certificate go back into
    // the graph (their blocks are repaired one by one afterwards). Those that
    // failed reachability stay deleted and become probes u ~> v below.
    for (size_t i = 0; i < dels.size(); ++i)
        if (!no_reach[i] && inS(dels[i].first) && inS(dels[i].second)) {
            bool und_new;
            storeInsert(dels[i].first, dels[i].second, und_new);
        }
    // Certify R = s \ T, growing T by the side a failing probe exhausts
    // (as in splitAfterDelete; such a side is small and cannot be part of
    // the giant SCC, and if it holds the anchor, the anchor moves).
    // Anchors: the in-anchor a (every z must reach it) and the out-anchor b
    // (it must reach every e), plus a ~> b and b ~> a. That is the single
    // anchor argument with a ~> b ~> e in place of a ~> e. In a big SCC they
    // are its top hubs, and the probes stop at their neighbourhoods.
    int a = -1, b = -1;
    bool hubs = false;  // a and b are hub anchors (else a == b, a local anchor)
    const bool use_hubs = C >= hub_anchor_min;
    if (use_hubs) pickHubCandidates();
    // Anchor set (where no single hub covers enough, e.g. livejournal): the
    // SCC's top in- and out-degree vertices, shown mutually reachable inside
    // R each round (probes h ~> a0, a0 ~> h). Then "z reaches some anchor"
    // and "some anchor reaches e" give the certificate as with one anchor;
    // a probe stops at the first vertex with an edge into (out of) an
    // anchor, marked in anc_in (anc_out). ANY stands for the set in probes.
    const int ANY = -2;
    bool multi = false, multi_off = false;
    std::vector<int> anchors;  // the anchors of anc_list still in R this round
    auto clearAnchors = [&]() { anchors.clear(); multi = false; };
    auto buildAnchors = [&]() {
        clearAnchors();
        if (anc_list.empty()) {  // first use: the global top hubs and their marks
            anc_in.assign(V, 0);
            anc_out.assign(V, 0);
            std::unordered_set<int> seen;
            for (size_t i = 0; i < hub_by_in.size() && anc_list.size() < 32; ++i)
                for (int h : {hub_by_in[i], i < hub_by_out.size() ? hub_by_out[i] : -1})
                    if (h >= 0 && anc_list.size() < 32 && seen.insert(h).second) anc_list.push_back(h);
            const int K = (int)anc_list.size();
            #pragma omp parallel for schedule(dynamic, 1)
            for (int k = 0; k < K; ++k) {
                const int h = anc_list[k];
                auto put = [&](std::vector<int>& M, int w) {
                    int exp = 0;
                    __atomic_compare_exchange_n(&M[w], &exp, k + 1, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED);
                };
                put(anc_in, h);
                put(anc_out, h);
                forIn(h, [&](int w) { put(anc_in, w); });
                forOut(h, [&](int w) { put(anc_out, w); });
            }
        }
        long long cover = 0;
        for (int h : anc_list)
            if (inS(h)) { anchors.push_back(h); cover += undDeg(h); }
        if (anchors.size() < 2 || (!force_anchor_set && cover * 64 < (long long)C)) { anchors.clear(); return; }
        multi = true;
    };
    // z ~> ANY (fwd) or ANY ~> e: bounded BFS to the first marked vertex.
    auto anchorProbe = [&](int x, bool fwd, std::vector<int>& path, bool big) -> int {
        const std::vector<int>& M = fwd ? anc_in : anc_out;
        auto& ls = LocalSearch::get();
        if (big && !ls.la) { ls.la.reset(new LocalSearch::Map(1 << 18)); ls.lb.reset(new LocalSearch::Map(1 << 18)); }
        auto& m = big ? *ls.la : ls.sa;
        const long long budget = big ? (1 << 20) : seq_probe_budget;
        // w's mark names an anchor h; it counts if h is in R and the edge
        // w -> h (fwd) / h -> w is still there.
        auto valid = [&](int w) {
            const int k = M[w] - 1;
            if (k < 0) return false;
            const int h = anc_list[k];
            if (!inS(h)) return false;
            return h == w || (fwd ? isDirectedEdgePresent(w, h) : isDirectedEdgePresent(h, w));
        };
        auto finish = [&](int w, int p) {
            path.clear();
            if (w != x) { path.push_back(w); for (int y = p; y != -1 && y != x; y = *m.find(y)) path.push_back(y); }
            const int h = anc_list[M[w] - 1];
            if (h != w) path.push_back(h);
            return 1;
        };
        if (valid(x)) { m.clear(); return finish(x, -1); }
        m.clear();
        m.put(x, -1);
        ls.qa.assign(1, x);
        long long work = 0;
        int r = 0, hw = -1, hp = -1;
        for (size_t q = 0; q < ls.qa.size() && r == 0; ++q) {
            const int p = ls.qa[q];
            auto visit = [&](int w) {
                if (++work >= budget) { r = -1; return true; }
                if (!inS(w) || m.find(w)) return false;
                if (valid(w)) { r = 1; hw = w; hp = p; return true; }
                if (!m.put(w, p)) { r = -1; return true; }
                ls.qa.push_back(w);
                return false;
            };
            if (fwd) anyOut(p, visit); else anyIn(p, visit);
        }
        return r == 1 ? finish(hw, hp) : r;
    };
    // The same search with all threads (level-synchronous BFS), for the few
    // probes that outgrow anchorProbe's sequential budget: one such probe
    // otherwise decides a round's time (tens of ms on livejournal).
    auto parAnchorProbe = [&](int x, bool fwd, std::vector<int>& path) -> int {
        const std::vector<int>& M = fwd ? anc_in : anc_out;
        auto valid = [&](int w) {
            const int k = M[w] - 1;
            if (k < 0) return false;
            const int h = anc_list[k];
            if (!inS(h)) return false;
            return h == w || (fwd ? isDirectedEdgePresent(w, h) : isDirectedEdgePresent(h, w));
        };
        if (!parBfs(x, fwd ? 0 : 1, s, valid, &path)) return 0;  // exhausted: no anchor
        const int w = path.empty() ? x : path.front();
        const int h = anc_list[M[w] - 1];
        if (h != w) path.push_back(h);
        return 1;
    };
    auto pickAnchors = [&]() {
        a = b = -1;
        hubs = false;
        if (!use_hubs) return;
        for (int h : hub_by_in) if (inS(h)) { a = h; break; }
        for (int h : hub_by_out) if (inS(h)) { b = h; break; }
        if (a < 0 || b < 0) { a = b = -1; return; }
        if (!force_anchor_set && hubsCover(a, b, C)) { hubs = true; return; }
        a = b = -1;
        if (!multi_off) buildAnchors();
        if (multi) a = b = anchors[0];
    };
    std::unordered_map<unsigned long long, std::vector<int>> passed;  // probe -> its path
    int near = -1;  // a vertex of R next to T (the WCC probes below search from it)
    for (int round = 0;; ++round) {
        if (round == 64) { undo(); return {}; }
        // Boundary of T in R, counting the batch's deleted edges.
        std::vector<int> Z, E;
        {
            std::unordered_set<int> zs, es;
            for (int x : T) {
                forIn(x, [&](int w) { if (inS(w) && zs.insert(w).second) Z.push_back(w); });
                forOut(x, [&](int w) { if (inS(w) && es.insert(w).second) E.push_back(w); });
            }
            for (auto& e : all_dels) {
                if (inS(e.first) && node_scc[e.second] == tmp && zs.insert(e.first).second) Z.push_back(e.first);
                if (node_scc[e.first] == tmp && inS(e.second) && es.insert(e.second).second) E.push_back(e.second);
            }
        }
        // Probes: z ~> a, a ~> e, and u ~> v for each deleted edge left
        // inside R (the old paths through it must be rerouted inside R).
        // With hub anchors a deleted edge u -> v left inside R becomes
        // u ~> a and b ~> v (its old paths are rerouted through the anchors).
        std::vector<std::pair<int, int>> pr;
        if (multi && !std::all_of(anchors.begin(), anchors.end(), inS)) { clearAnchors(); multi_off = true; a = b = -1; }
        if (a < 0 || !inS(a) || b < 0 || !inS(b)) pickAnchors();
        for (size_t i = 0; i < dels.size(); ++i)
            if (no_reach[i] && inS(dels[i].first) && inS(dels[i].second)) {
                if (hubs) { pr.push_back({dels[i].first, a}); pr.push_back({b, dels[i].second}); }
                else if (multi) { pr.push_back({dels[i].first, ANY}); pr.push_back({ANY, dels[i].second}); }
                else pr.push_back(dels[i]);
            }
        if (Z.empty() && E.empty() && pr.empty()) { certified = true; break; }
        if (multi) {
            for (size_t k = 1; k < anchors.size(); ++k) {
                pr.push_back({anchors[k], anchors[0]});
                pr.push_back({anchors[0], anchors[k]});
            }
            for (int z : Z) pr.push_back({z, ANY});
            for (int e : E) pr.push_back({ANY, e});
        } else {
            if (!hubs) {
                if (a < 0 || !inS(a)) a = !Z.empty() ? Z[0] : !E.empty() ? E[0] : pr[0].first;
                b = a;
            } else if (a != b) {
                pr.push_back({a, b});
                pr.push_back({b, a});
            }
            for (int z : Z) pr.push_back({z, a});
            for (int e : E) pr.push_back({b, e});
        }
        anchor_set_calls += multi;
        hub_anchor_calls += hubs;
        for (int x : Z) if (inS(x)) { near = x; break; }
        if (near < 0 || !inS(near)) for (int x : E) if (inS(x)) { near = x; break; }
        auto tgt = [&](int x) { return x == ANY ? anchors[0] : x; };  // (exact fallbacks use anchor 0)
        auto inS2 = [&](int x) { return x == ANY || inS(x); };
        // A probe that passed in an earlier round stays certified while its
        // path avoids T (the edges do not change during the rounds).
        std::vector<char> ok(pr.size(), 0), fresh(pr.size(), 0), longp(pr.size(), 0);
        std::vector<std::vector<int>> paths(pr.size());
        auto key = [](int x, int y) { return (unsigned long long)(unsigned)x << 32 | (unsigned)y; };
        #pragma omp parallel for schedule(dynamic, 4)
        for (size_t i = 0; i < pr.size(); ++i) {
            const int x = pr[i].first, y = pr[i].second;
            auto it = passed.find(key(x, y));
            if (it != passed.end() && std::all_of(it->second.begin(), it->second.end(), inS)) { ok[i] = 1; continue; }
            if (multi && x != ANY && y != ANY) {  // anchor pair: a path from an earlier batch, if intact
                auto pt = anc_pair_path.find(key(x, y));
                if (pt != anc_pair_path.end()) {
                    const auto& P = pt->second;
                    bool good = std::all_of(P.begin(), P.end(), inS);
                    for (size_t q = 0; good && q <= P.size(); ++q) {
                        const int from = q == 0 ? x : P[q - 1], to = q == P.size() ? y : P[q];
                        good = isDirectedEdgePresent(from, to);
                    }
                    if (good) { paths[i] = P; ok[i] = fresh[i] = 1; continue; }
                }
            }
            int r = -1;
            if (hubs && y == a && x != b) { r = hubProbe(x, a, true, s, seq_probe_budget, &paths[i]); if (r < 0) { longp[i] = 2; continue; } }
            else if (hubs && x == b && y != a) { r = hubProbe(y, b, false, s, seq_probe_budget, &paths[i]); if (r < 0) { longp[i] = 2; continue; } }
            else if (y == ANY || x == ANY) {
                r = anchorProbe(y == ANY ? x : y, y == ANY, paths[i], false);
                if (r < 0) { longp[i] = 1; continue; }  // second pass, all threads
            }
            ok[i] = r == 1 || (r < 0 && reachRO(tgt(x), tgt(y), s, &paths[i]));
            fresh[i] = ok[i];
        }
        for (size_t i = 0; i < pr.size(); ++i) {
            if (!longp[i]) continue;
            const int x = pr[i].first, y = pr[i].second;
            if (longp[i] == 1) { ok[i] = fresh[i] = parAnchorProbe(y == ANY ? x : y, y == ANY, paths[i]) == 1; continue; }
            // Hub probe on all threads: x ~> a, or b ~> y.
            const bool f = y == a;
            const int src = f ? x : y, h = f ? a : b;
            auto hit = [&](int w) { return w == h || (f ? isDirectedEdgePresent(w, h) : isDirectedEdgePresent(h, w)); };
            ok[i] = fresh[i] = parBfs(src, f ? 0 : 1, s, hit, &paths[i]) == 1;
            if (ok[i]) paths[i].push_back(h);
        }
        for (size_t i = 0; i < pr.size(); ++i) {
            if (!fresh[i]) continue;
            const int x = pr[i].first, y = pr[i].second;
            if (multi && x != ANY && y != ANY && !longp[i]) anc_pair_path[key(x, y)] = paths[i];  // (reachRO paths are in order)
            passed[key(x, y)] = std::move(paths[i]);
        }
        // Failed probes, one by one against the current T (earlier ones may
        // already have moved what this one would find). Growing T only
        // shrinks R, so passed probes need no recheck until the next round.
        size_t grown = 0;
        bool moved = false;
        for (size_t i = 0; i < pr.size() && !moved; ++i) {
            if (ok[i] || !inS2(pr[i].first) || !inS2(pr[i].second)) continue;
            std::vector<int> side;
            bool fwd = false;
            if (reachableWithinScc(tgt(pr[i].first), tgt(pr[i].second), s, &side, &fwd)) continue;
            if ((T.size() + side.size()) * 2 > C) { undo(); return {}; }
            // The exhausted side: what the source reaches (fwd) or what
            // reaches the target. If it holds the anchor, the anchor moves to
            // the probe's other end, and the anchor probes are redone.
            bool has_anchor = false;
            if (multi) {
                std::unordered_set<int> as(anchors.begin(), anchors.end());
                for (int x : side) if (as.count(x)) { has_anchor = true; break; }
            }
            if (has_anchor) {  // back to a single local anchor for this call
                clearAnchors();
                multi_off = true;
                a = b = -1;
                moved = true;
            } else if (std::find(side.begin(), side.end(), a) != side.end() ||
                std::find(side.begin(), side.end(), b) != side.end()) {
                if (hubs) a = b = -1;  // re-picked next round
                else a = fwd ? pr[i].second : pr[i].first;
                moved = true;
            }
            for (int x : side) { node_scc[x] = tmp; T.push_back(x); }
            grown += side.size();
        }
        if (grown == 0) { certified = true; break; }
    }
    if (T.empty()) { undo(); return {}; }  // nothing splits
    // Deletions put back above whose endpoint has since joined T are taken
    // out again (removing edges that touch T leaves R's certificate intact).
    // Deletions still inside R did not split it: they go back too and get
    // their block repair on the normal path.
    for (auto& e : dels) {
        bool touches_t = node_scc[e.first] == tmp || node_scc[e.second] == tmp;
        bool present = isDirectedEdgePresent(e.first, e.second);
        bool flag;
        if (touches_t && present) storeDelete(e.first, e.second, flag);
        else if (!touches_t && !present) storeInsert(e.first, e.second, flag);
    }
    std::vector<std::vector<int>> below, above;
    orderPieces(s, tmp, T, below, above);
    for (int x : T) node_scc[x] = s;
    retireScc(tmp);
    // Deleted undirected edges between T and R were edges of their block.
    std::unordered_set<int> Ts(T.begin(), T.end());
    std::vector<std::pair<int, int>> extras;
    for (auto& e : all_dels) {
        int x = e.first, y = e.second;
        if (Ts.count(y)) std::swap(x, y);
        if (!Ts.count(x) || Ts.count(y) || !inS(y)) continue;
        if (!isUndirectedEdgePresent(x, y) && commonBlock(x, y) >= 0) extras.push_back({x, y});
    }
    detachPieces(s, below, above, extras);
    // A vertex that lost all its edges into R may have left the WCC.
    for (int x : T)
        if (!anyUnd(x, [&](int w) { return node_scc[w] == s; })) wcc_probe.push_back({x, near >= 0 && inS(near) ? near : a});
    return std::vector<int>(Ts.begin(), Ts.end());
}

void GIR::normalizeBatch(const std::vector<Update>& batch, std::vector<std::pair<int, int>>& ins,
                         std::vector<std::pair<int, int>>& del, long long& noops,
                         std::vector<std::pair<long long, char>>* before) {
    // Last operation per edge wins: sort (key, position), keep the last entry
    // of each key. Invalid operations get a key that sorts last.
    const size_t n = batch.size();
    const long long BAD = LLONG_MAX;
    std::vector<std::pair<long long, int>> kv(n);  // (key, position)
    #pragma omp parallel for schedule(static) if (n >= 4096)
    for (size_t i = 0; i < n; ++i) {
        const auto& op = batch[i];
        bool bad = op.u == op.v || op.u < 0 || op.v < 0 || op.u >= V || op.v >= V;
        kv[i] = {bad ? BAD : keyD(op.u, op.v), (int)i};
    }
    if (n >= 4096) __gnu_parallel::sort(kv.begin(), kv.end());
    else std::sort(kv.begin(), kv.end());
    size_t m = n;
    while (m > 0 && kv[m - 1].first == BAD) --m;
    // Each key's last entry, classified against the current store. Results
    // stay in key order (ins/del sorted by source vertex: cache-friendly for
    // the store updates).
    // Code per sorted entry (parallel lookups), then one sequential pass
    // collects them (cheaper than per-thread buffers for batch sizes).
    // 0 = skip, 1 = insert, 2 = delete, 3 = no-op insert.
    std::vector<char> code(m);
    #pragma omp parallel for schedule(static) if (m >= 4096)
    for (size_t i = 0; i < m; ++i) {
        if (i + 1 < m && kv[i + 1].first == kv[i].first) { code[i] = 0; continue; }  // not the last op on this edge
        const long long k = kv[i].first;
        const bool want = batch[kv[i].second].insert;
        const bool present = isDirectedEdgePresent((int)(k >> 32), (int)(unsigned)(k & 0xffffffffLL));
        code[i] = want == present ? (want ? 3 : 0) : want ? 1 : 2;
    }
    // One counting pass (per chunk: inserts, deletions, no-ops), then the
    // three outputs written at their prefix offsets.
    using E = std::pair<int, int>;
    const int T_ = omp_get_max_threads();
    const size_t CH = std::max<size_t>(1024, (m + 4 * T_ - 1) / (4 * T_)), nch = (m + CH - 1) / CH;
    std::vector<size_t> ci(nch + 1, 0), cd(nch + 1, 0);
    long long nn = 0;
    #pragma omp parallel for schedule(static) reduction(+ : nn) if (m >= 8192)
    for (size_t ch = 0; ch < nch; ++ch) {
        size_t a = 0, d = 0;
        for (size_t i = ch * CH, e = std::min(m, i + CH); i < e; ++i) { a += code[i] == 1; d += code[i] == 2; nn += code[i] == 3; }
        ci[ch + 1] = a; cd[ch + 1] = d;
    }
    for (size_t ch = 0; ch < nch; ++ch) { ci[ch + 1] += ci[ch]; cd[ch + 1] += cd[ch]; }
    ins.resize(ci[nch]);
    del.resize(cd[nch]);
    const size_t b0 = before ? before->size() : 0;
    if (before) before->resize(b0 + ci[nch] + cd[nch]);
    #pragma omp parallel for schedule(static) if (m >= 8192)
    for (size_t ch = 0; ch < nch; ++ch) {
        size_t a = ci[ch], d = cd[ch], b = b0 + ci[ch] + cd[ch];
        for (size_t i = ch * CH, e = std::min(m, i + CH); i < e; ++i) {
            if (code[i] != 1 && code[i] != 2) continue;
            const long long k = kv[i].first;
            const E ed{(int)(k >> 32), (int)(unsigned)(k & 0xffffffffLL)};
            if (code[i] == 1) ins[a++] = ed; else del[d++] = ed;
            if (before) (*before)[b++] = {k, (char)(code[i] == 2)};
        }
    }
    noops = nn;
}

void GIR::ingestBatch(const std::vector<Update>& batch) {
    std::vector<std::pair<int, int>> ins, del;
    long long noops;
    normalizeBatch(batch, ins, del, noops, &pending);
    logChanges(del);
    logChanges(ins);
    storeDeleteBulk(del);
    storeInsertBulk(ins);
}

// Net change since the last sync: edges whose presence now differs from
// their presence then. Put the store back to the synced state for exactly
// those edges, then run them through applyBatch (which updates the store
// again, plus every layer).
void GIR::sync() {
    if (pending.empty()) return;
    std::vector<std::pair<long long, char>> p;
    p.swap(pending);
    if (p.size() >= (1u << 15)) __gnu_parallel::stable_sort(p.begin(), p.end(), [](auto& a, auto& b) { return a.first < b.first; });
    else std::stable_sort(p.begin(), p.end(), [](auto& a, auto& b) { return a.first < b.first; });
    std::vector<std::pair<long long, char>> first;  // first record per edge = presence at the last sync
    for (size_t i = 0; i < p.size(); ++i) if (i == 0 || p[i].first != p[i - 1].first) first.push_back(p[i]);
    std::vector<char> now(first.size());
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < first.size(); ++i) {
        long long k = first[i].first;
        now[i] = isDirectedEdgePresent((int)(k >> 32), (int)(unsigned)(k & 0xffffffffLL));
    }
    std::vector<std::pair<int, int>> added, removed;
    for (size_t i = 0; i < first.size(); ++i) {
        if ((bool)now[i] == (bool)first[i].second) continue;
        long long k = first[i].first;
        std::pair<int, int> e{(int)(k >> 32), (int)(unsigned)(k & 0xffffffffLL)};
        (now[i] ? added : removed).push_back(e);
    }
    storeDeleteBulk(added);
    storeInsertBulk(removed);
    std::vector<Update> ops;
    ops.reserve(added.size() + removed.size());
    for (auto& e : removed) ops.push_back({e.first, e.second, false});
    for (auto& e : added) ops.push_back({e.first, e.second, true});
    in_sync = true;
    applyBatch(ops);
    in_sync = false;
}

void GIR::applyBatch(const std::vector<Update>& batch) {
    sync();
    std::vector<std::pair<int, int>> ins, del;
    {
        long long noops;
        normalizeBatch(batch, ins, del, noops);
        case_counts[C4_NOOP] += noops;
    }
    logChanges(del);
    logChanges(ins);

    // ---- Certify (parallel, read-only) against the graph after the whole batch ----
    // Apply every deletion to the edge store first, so certificates are
    // checked in the final graph F. A deleted edge u->v inside SCC s changes no
    // layer if u still reaches v in F (that path also exists in every graph
    // between F and the old one, and its vertices are in s), and, if {u,v}
    // disappears too, u and v keep two internally disjoint paths inside s.
    // Let C be the certified deletions. Every C edge is replaced by such paths
    // in old-minus-C, so SCCs, blocks and the block-cut forest of old-minus-C
    // equal the old ones (a cycle through two edges maps to a closed walk in
    // one block). An inter-SCC deletion is certified if {u,v} stays or u and
    // v stay connected. The uncertified deletions are put back and run through
    // the normal path, now against old-minus-C, which the structure matches.
    std::vector<char> del_noop(del.size(), 0), del_und(del.size(), 0), del_reach(del.size(), 0), del_long(del.size(), 0);
    storeDeleteBulk(del);
    // Whether {u,v} is gone in F (both directions may be in the batch).
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < del.size(); ++i) del_und[i] = !isUndirectedEdgePresent(del[i].first, del[i].second);
    // Hub shortcut for the largest SCC: with hubs a, b in it and a ~> b in
    // F, u ~> a and b ~> v (short searches, see hubProbe) give u ~> v.
    int hub_s = -1, hub_a = -1, hub_b = -1;
    if (!del.empty() && V >= (int)hub_anchor_min) {
        pickHubCandidates();
        for (int h : hub_by_in)
            if (sccs[node_scc[h]].size() >= hub_anchor_min) { hub_a = h; hub_s = node_scc[h]; break; }
        if (hub_s >= 0)
            for (int h : hub_by_out) if (node_scc[h] == hub_s) { hub_b = h; break; }
        if (hub_b < 0 || !hubsCover(hub_a, hub_b, sccs[hub_s].size()) ||
            (hub_a != hub_b && !reachRO(hub_a, hub_b, hub_s))) hub_s = -1;
    }
    // Bounded undirected BFS from x until the hub or one of its neighbours.
    auto undToHub = [&](int x) {
        const int h = hub_a;
        if (x == h || isUndirectedEdgePresent(x, h)) return true;
        auto& ls = LocalSearch::get();
        auto& m = ls.sa;
        m.clear();
        m.put(x, 0);
        ls.qa.assign(1, x);
        long long work = 0;
        int r = 0;
        for (size_t q = 0; q < ls.qa.size() && r == 0; ++q)
            anyUnd(ls.qa[q], [&](int w) {
                if (++work >= std::min<long long>(4096, seq_probe_budget)) { r = -1; return true; }
                if (m.find(w)) return false;
                if (w == h || isUndirectedEdgePresent(w, h)) { r = 1; return true; }
                if (!m.put(w, 0)) { r = -1; return true; }
                ls.qa.push_back(w);
                return false;
            });
        return r == 1;
    };
    #pragma omp parallel for schedule(dynamic, 16)
    for (size_t i = 0; i < del.size(); ++i) {
        int u = del[i].first, v = del[i].second;
        int s = node_scc[u];
        if (s != node_scc[v]) {
            // {u,v} gone: u and v stay connected if both are still connected
            // to the hub vertex in the graph after the batch (twitter: most
            // vertices are a hop or two away). (Not "to the hub's SCC": this
            // batch's deletions may split it.)
            if (!del_und[i]) { del_noop[i] = 1; continue; }
            if (hub_s >= 0 && undToHub(u) && undToHub(v)) { del_noop[i] = 1; continue; }
            del_noop[i] = reachUndRO(u, v);
            continue;
        }
        int cert = -1;
        if (s == hub_s && (blocks_on_demand || !del_und[i])) {
            const int h1 = hubProbe(u, hub_a, true, s, seq_probe_budget);
            const int h2 = h1 == 1 ? hubProbe(v, hub_b, false, s, seq_probe_budget) : 0;
            if (h1 == 1 && h2 == 1) cert = blocks_on_demand ? 2 : 1;  // (what certifyKeep returns for a path without block check)
            else if (h1 < 0 || h2 < 0) { del_long[i] = 1; continue; }  // second pass, all threads
        }
        if (cert < 0) cert = certifyKeep(u, v, s, del_und[i]);
        del_reach[i] = cert >= 1;
        del_noop[i] = cert == 2 || (cert == 1 && !del_und[i]);
    }
    // Hub certificates that outgrew the sequential budget: again in
    // parallel with the large tables and a larger budget, and what is still
    // left one at a time on all threads (one such search otherwise sets the
    // batch's time; doing all of them that way costs ~0.1 ms each).
    {
        std::vector<size_t> L;
        for (size_t i = 0; i < del.size(); ++i) if (del_long[i]) L.push_back(i);
        #pragma omp parallel for schedule(dynamic, 1)
        for (size_t j = 0; j < L.size(); ++j) {
            const size_t i = L[j];
            const int u = del[i].first, v = del[i].second, s = node_scc[u];
            const int h1 = hubProbe(u, hub_a, true, s, seq_probe_budget << 4, nullptr, true);
            const int h2 = h1 == 1 ? hubProbe(v, hub_b, false, s, seq_probe_budget << 4, nullptr, true) : 0;
            if (h1 == 1 && h2 == 1) { del_long[i] = 0; del_reach[i] = 1; del_noop[i] = 1; }
            else if (h1 == 0 || (h1 == 1 && h2 == 0)) del_long[i] = 0;  // exhausted: not certified
        }
    }
    for (size_t i = 0; i < del.size(); ++i) {
        if (!del_long[i]) continue;
        const int u = del[i].first, v = del[i].second, s = node_scc[u];
        auto to_a = [&](int w) { return w == hub_a || isDirectedEdgePresent(w, hub_a); };
        auto from_b = [&](int w) { return w == hub_b || isDirectedEdgePresent(hub_b, w); };
        const bool r = parBfs(u, 0, s, to_a) == 1 && parBfs(v, 1, s, from_b) == 1;
        del_reach[i] = r;
        del_noop[i] = r;  // (hub pass only runs with blocks on demand or when {u,v} stays)
    }
    if (blocks_on_demand)  // an undirected edge leaving an SCC may change its blocks
        for (size_t i = 0; i < del.size(); ++i)
            if (del_und[i] && node_scc[del[i].first] == node_scc[del[i].second]) {
                markBlocksDirty(del[i].first);
                markBlocksDirty(del[i].second);
            }
    // Certified inter-SCC deletions drop their condensation edge now, before
    // any split re-derives DAG counts from the edge store.
    for (size_t i = 0; i < del.size(); ++i) {
        if (!del_noop[i]) continue;
        int u = del[i].first, v = del[i].second;
        if (node_scc[u] != node_scc[v]) { removeDagEdge(node_scc[u], node_scc[v]); case_counts[D_INTER_SCC]++; }
        else case_counts[D_INTRA_SCC_KEEP]++;
    }
    // Uncertified deletions inside a large SCC are first tried in bulk
    // (bulkTrimSplit); every other uncertified deletion is put back and runs
    // through the normal path below.
    std::vector<std::pair<int, int>> wcc_probe;
    std::vector<std::pair<int, int>> wcc_trimmed;  // (trimmed vertex, rest vertex) from bulkTrimSplit
    std::vector<char> del_bulk(del.size(), 0);
    {
        std::unordered_map<int, std::vector<size_t>> big;
        for (size_t i = 0; i < del.size(); ++i) {
            if (del_noop[i]) continue;
            int u = del[i].first, v = del[i].second, s = node_scc[u];
            if (!force_resplit && s == node_scc[v] && sccs[s].size() >= fastpath_min_block) big[s].push_back(i);
            else { bool und_new; storeInsert(u, v, und_new); }
        }
        for (auto& kv : big) {
            std::vector<std::pair<int, int>> d;
            std::vector<char> nr;
            for (size_t i : kv.second) { d.push_back(del[i]); nr.push_back(!del_reach[i]); }
            bool certified;
            std::vector<int> T = bulkTrimSplit(kv.first, d, nr, del, wcc_trimmed, certified);
            std::unordered_set<int> Ts(T.begin(), T.end());
            for (size_t i : kv.second) {
                if (Ts.count(del[i].first) || Ts.count(del[i].second)) { del_bulk[i] = 1; case_counts[D_SCC_SPLIT]++; }
                // R = s \ T is strongly connected without all its deletions
                // (bulkTrimSplit's certificate), so each of them keeps the SCC.
                else if (certified && !del_reach[i]) del_reach[i] = 2;  // (block certificate not tried yet)
            }
        }
        for (size_t i = 0; i < del.size(); ++i)
            if (!del_noop[i] && !del_bulk[i]) { bool und_new; storeInsert(del[i].first, del[i].second, und_new); }
    }

    // ---- deletions: edge store + condensation edges, repairs deferred ----
    std::unordered_map<int, std::vector<std::pair<int, int>>> intra;  // scc -> deleted edges
    std::vector<std::pair<int, int>> dirty_blocks;  // (scc, block)
    std::vector<std::pair<int, int>> dirty_edges;   // parallel to dirty_blocks
    for (size_t di = 0; di < del.size(); ++di) {
        auto& e = del[di];
        int u = e.first, v = e.second;
        if (del_noop[di] || del_bulk[di]) continue;
        // Inside a large SCC a bulk re-split would cost O(|SCC|) per batch;
        // the single-edge path with its local shortcuts is far cheaper. It is
        // safe to run it here: it only inspects paths inside this SCC and the
        // edges of the vertices it moves, and every deferred change so far
        // (other SCCs, condensation edges) is already in the edge store and
        // the DAG counts.
        if (node_scc[u] == node_scc[v] && sccs[node_scc[u]].size() >= fastpath_min_block) {
            deleteEdgeImpl(u, v, del_reach[di] == 1 ? 1 : del_reach[di] == 2 ? 3 : 2);
            continue;
        }
        bool und_gone;
        storeDelete(u, v, und_gone);
        int su = node_scc[u], sv = node_scc[v];
        if (su != sv) {
            removeDagEdge(su, sv);
            if (und_gone && !dagEdgeExists(su, sv)) wcc_probe.push_back(e);
            else case_counts[D_INTER_SCC]++;
        } else {
            intra[su].push_back(e);
            if (und_gone) {
                dirty_blocks.push_back({su, commonBlock(u, v)});
                dirty_edges.push_back(e);
                wcc_probe.push_back(e);  // only probed if the SCC ends up split
            }
        }
    }
    // ---- one split check per touched SCC ----
    std::unordered_set<int> split;
    for (auto& kv : intra) {
        int s = kv.first;
        bool intact = true;
        for (auto& e : kv.second)
            if (!reachableWithinScc(e.first, e.second, s)) { intact = false; break; }
        if (!intact && resplitScc(s)) {
            split.insert(s);
            case_counts[D_SCC_SPLIT] += kv.second.size();
        } else {
            case_counts[D_INTRA_SCC_KEEP] += kv.second.size();
        }
    }
    // ---- one repair per touched block (split SCCs were rebuilt) ----
    {
        // Block -> deleted edges inside it. A block with exactly one deleted
        // edge can use the single-deletion shortcuts (changes elsewhere in
        // the batch can't affect paths inside it).
        std::unordered_map<int, std::vector<std::pair<int, int>>> per_block;
        std::vector<int> todo;
        for (size_t i = 0; i < dirty_blocks.size(); ++i) {
            if (split.count(dirty_blocks[i].first)) continue;
            int b = findBlk(dirty_blocks[i].second);
            auto& lst = per_block[b];
            if (lst.empty()) todo.push_back(b);
            lst.push_back(dirty_edges[i]);
        }
        for (int b : todo) {
            auto& lst = per_block[b];
            if (lst.size() == 1) repairBlockAfterDelete(b, lst[0].first, lst[0].second);
            else repairBlock(b);
        }
    }
    // ---- WCC splits ----
    // A batch can split one WCC into several pieces. Every piece contains an
    // endpoint of some deleted edge, so it suffices to make each such
    // endpoint's final component a label class of its own: per old label,
    // compare every endpoint against a representative and split off
    // whichever side the lockstep probe exhausts.
    {
        std::unordered_map<int, std::vector<int>> by_label;
        std::unordered_set<int> seen;
        for (auto* L : {&wcc_probe, &wcc_trimmed})
            for (auto& e : *L) {
                if (node_scc[e.first] == node_scc[e.second]) continue;  // SCC intact => connected
                for (int x : {e.first, e.second})
                    if (seen.insert(x).second) by_label[wcc_label[x]].push_back(x);
            }
        // Cheap cases first, in parallel: an endpoint with no edges left is a
        // component of its own; one that a bounded search connects to its
        // label's representative shares the representative's fate below.
        {
            std::vector<std::pair<int, int>> cand;  // (label, endpoint)
            for (auto& kv : by_label) {
                auto& lst = kv.second;
                int r = lst[0];  // representative: endpoint in the largest SCC
                for (int x : lst) if (sccs[node_scc[x]].size() > sccs[node_scc[r]].size()) r = x;
                std::swap(*std::find(lst.begin(), lst.end(), r), lst[0]);
                for (size_t i = 1; i < lst.size(); ++i) cand.push_back({kv.first, lst[i]});
            }
            // 1 = isolated, 2 = joined to the representative, 3 = the bounded
            // search ran out of vertices: what it visited is x's whole new
            // component (it does not reach the representative's SCC), which
            // is split off without the exact sequential search. On sparse,
            // nearly acyclic graphs (cit-Patents) almost every deletion splits
            // a small piece off this way.
            // A label with many endpoints (the giant WCC) first grows a ball
            // around its representative, a parallel undirected BFS of up to
            // 64K vertices (on cit-Patents 2^15-2^17 measured best; 2M cost
            // more to build than it saved): all in the representative's component, so an
            // endpoint whose search meets the ball stays. Without it, on
            // graphs whose SCCs are single vertices (cit-Patents) "meets the
            // representative's SCC" means meeting one vertex, which half of
            // the searches never did, and each went to the exact sequential
            // search (~190 us).
            std::unordered_map<int, int> ball_of;  // label -> epoch
            for (auto& kv : by_label) {
                // (only where the representative's SCC is a poor target: in
                // a giant SCC, as on livejournal and twitter, the searches
                // reach it anyway and building the ball cost ~20% of a batch)
                if (kv.second.size() < 64 || sccs[node_scc[kv.second[0]]].size() >= 1024) continue;
                if (ball_mark.size() < (size_t)V) ball_mark.assign(V, 0);
                if (++ball_epoch == INT_MAX) { std::fill(ball_mark.begin(), ball_mark.end(), 0); ball_epoch = 1; }
                const int ep = ball_epoch, r = kv.second[0];
                ball_of[kv.first] = ep;
                ball_mark[r] = ep;
                // The cap is enforced inside a level too (a shared counter, and
                // scans stop at it): from a hub one level reaches millions of
                // vertices on twitter, and checking only between levels made
                // the "64K" ball cost tens of ms per batch there.
                const long long CAP = 1 << 16;
                std::vector<int> fr{r}, next;
                std::vector<std::vector<int>> parts(omp_get_max_threads());
                long long total = 1;
                while (!fr.empty() && total < CAP) {
                    #pragma omp parallel for schedule(dynamic, 64)
                    for (size_t i = 0; i < fr.size(); ++i) {
                        if (__atomic_load_n(&total, __ATOMIC_RELAXED) >= CAP) continue;
                        auto& out = parts[omp_get_thread_num()];
                        anyUnd(fr[i], [&](int w) {
                            int e = __atomic_load_n(&ball_mark[w], __ATOMIC_RELAXED);
                            if (e != ep && __atomic_compare_exchange_n(&ball_mark[w], &e, ep, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
                                out.push_back(w);
                                if (__atomic_add_fetch(&total, 1, __ATOMIC_RELAXED) >= CAP) return true;
                            }
                            return false;
                        });
                    }
                    next.clear();
                    for (auto& p : parts) { next.insert(next.end(), p.begin(), p.end()); p.clear(); }
                    fr.swap(next);
                }
            }
            std::vector<char> st(cand.size(), 0);
            std::vector<std::vector<int>> side(cand.size());
            #pragma omp parallel for schedule(dynamic, 16)
            for (size_t i = 0; i < cand.size(); ++i) {
                const int x = cand[i].second;
                if (!anyUnd(x, [](int) { return true; })) { st[i] = 1; continue; }
                const int target = node_scc[by_label.at(cand[i].first)[0]];
                auto bit = ball_of.find(cand[i].first);
                const int bep = bit == ball_of.end() ? -1 : bit->second;
                if (bep >= 0 && ball_mark[x] == bep) { st[i] = 2; continue; }
                auto& ls = LocalSearch::get();
                auto& m = ls.sa;
                m.clear();
                m.put(x, 0);
                ls.qa.assign(1, x);
                long long work = 0;
                int r = 0;  // 0 exhausted, 1 reached, -1 gave up
                for (size_t q = 0; q < ls.qa.size() && r == 0; ++q)
                    anyUnd(ls.qa[q], [&](int w) {
                        if (++work >= (1 << 16)) { r = -1; return true; }
                        if (node_scc[w] == target || (bep >= 0 && ball_mark[w] == bep)) { r = 1; return true; }
                        if (!m.find(w)) {
                            if (!m.put(w, 0)) { r = -1; return true; }
                            ls.qa.push_back(w);
                        }
                        return false;
                    });
                if (r == 1) st[i] = 2;
                else if (r == 0) { st[i] = 3; side[i] = ls.qa; }
            }
            std::unordered_set<int> drop;
            for (size_t i = 0; i < cand.size(); ++i) {
                if (!st[i]) continue;
                int x = cand[i].second;
                drop.insert(x);
                if (st[i] == 1 || (st[i] == 3 && wcc_label[x] == cand[i].first)) {
                    // (with 3, an endpoint of the same new component may have
                    // labelled it already: its search found the same set)
                    const int L = cand[i].first, nl = (int)wcc_size.size();
                    if (st[i] == 1) side[i].assign(1, x);
                    wcc_size.push_back((int)side[i].size());
                    for (int y : side[i]) wcc_label[y] = nl;
                    wcc_size[L] -= (int)side[i].size();
                    case_counts[D_WCC_SPLIT]++;
                }
            }
            if (!drop.empty())
                for (auto& kv : by_label) {
                    auto& lst = kv.second;
                    lst.erase(std::remove_if(lst.begin() + 1, lst.end(), [&](int x) { return drop.count(x) > 0; }), lst.end());
                }
        }
        for (auto& kv : by_label) {
            const int L = kv.first;
            std::vector<int> rest = kv.second;
            while (rest.size() >= 2) {
                int r = rest[0];
                std::vector<int> again;
                for (size_t i = 1; i < rest.size(); ++i) {
                    int x = rest[i];
                    if (!again.empty() || wcc_label[r] != L) {  // r was split off: re-queue
                        if (wcc_label[x] == L) again.push_back(x);
                        continue;
                    }
                    if (wcc_label[x] != L) continue;  // x's component already split off
                    std::vector<int> side;
                    if (undirectedConnected(r, x, side)) continue;
                    int nl = (int)wcc_size.size();
                    wcc_size.push_back((int)side.size());
                    for (int y : side) wcc_label[y] = nl;
                    wcc_size[L] -= (int)side.size();
                    case_counts[D_WCC_SPLIT]++;
                    if (wcc_label[r] != L && wcc_label[x] == L) again.push_back(x);  // r's side left; x still unresolved
                }
                rest.swap(again);
            }
        }
        for (auto& e : wcc_probe)
            if (node_scc[e.first] != node_scc[e.second] && wcc_label[e.first] == wcc_label[e.second])
                case_counts[D_INTER_SCC]++;
    }
    // ---- insertions ----
    // Filter (parallel, read-only): an insertion whose endpoints already
    // share a block (or are already adjacent inside one SCC) changes no
    // layer, and stays that way under the other insertions, which only merge.
    // 1 = no-op (inside an SCC, no block change), 2 = forward condensation
    // edge (same WCC, rank(su) < rank(sv): no cycle, the order stays valid).
    std::vector<char> ins_noop(ins.size(), 0);
    #pragma omp parallel for schedule(dynamic, 256)
    for (size_t i = 0; i < ins.size(); ++i) {
        int u = ins[i].first, v = ins[i].second;
        const int su = node_scc[u], sv = node_scc[v];
        if (su != sv) {
            if (wcc_label[u] == wcc_label[v] && scc_rank[su] < scc_rank[sv]) ins_noop[i] = 2;
            continue;
        }
        ins_noop[i] = blocks_on_demand || isUndirectedEdgePresent(u, v) || commonBlockRO(u, v) >= 0;
    }
    // Whether each insertion adds its undirected edge, decided before any of
    // them is stored (ins is in key order, as the sequential path processed
    // it: of u->v and v->u in one batch, the smaller key adds {u,v}).
    std::vector<char> und_new(ins.size());
    #pragma omp parallel for schedule(static) if (ins.size() >= 4096)
    for (size_t i = 0; i < ins.size(); ++i) {
        const int u = ins[i].first, v = ins[i].second;
        und_new[i] = !isUndirectedEdgePresent(u, v) && !(u > v && std::binary_search(ins.begin(), ins.end(), std::make_pair(v, u)));
    }
    // No-op insertions and forward condensation edges change no SCC, WCC or
    // block. Their condensation-edge counts are added in parallel, grouped by
    // SCC so that each map is changed by one thread (dag_out by source, then
    // dag_in by target).
    {
        using E = std::pair<int, int>;
        std::vector<E> fwd;  // (source SCC, target SCC)
        parCompact<1, E>(ins.size(), fwd, [&](size_t i, auto put) {
            if (ins_noop[i] == 2) put(E{node_scc[ins[i].first], node_scc[ins[i].second]});
        });
        long long noop = 0;
        #pragma omp parallel for schedule(static) reduction(+ : noop) if (ins.size() >= 4096)
        for (size_t i = 0; i < ins.size(); ++i) noop += ins_noop[i] == 1;
        auto addGrouped = [&](std::vector<E>& l, bool out) {
            psort(l);
            std::vector<size_t> st;
            parCompact<1, size_t>(l.size(), st, [&](size_t i, auto put) { if (i == 0 || l[i].first != l[i - 1].first) put(i); });
            st.push_back(l.size());
            int grew = 0;
            const size_t ng = st.size() - 1;
            #pragma omp parallel for schedule(dynamic, 16) reduction(| : grew) if (l.size() >= 1024)
            for (size_t g = 0; g < ng; ++g) {
                auto& m = (out ? dag_out : dag_in)[l[st[g]].first];
                for (size_t i = st[g]; i < st[g + 1]; ++i) grew |= m[l[i].second]++ == 0;
            }
            return grew;
        };
        std::vector<E> rev(fwd.size());
        #pragma omp parallel for schedule(static) if (fwd.size() >= 4096)
        for (size_t i = 0; i < fwd.size(); ++i) rev[i] = {fwd[i].second, fwd[i].first};
        if (addGrouped(fwd, true)) ++struct_version;
        addGrouped(rev, false);
        if (blocks_on_demand) {  // a new undirected edge inside an SCC may change its blocks
            std::vector<int> d;
            parCompact<2, int>(ins.size(), d, [&](size_t i, auto put) {
                if (ins_noop[i] == 1 && und_new[i]) { put(ins[i].first); put(ins[i].second); }
            });
            blk_dirty_v.insert(blk_dirty_v.end(), d.begin(), d.end());
        }
        case_counts[C4_NOOP] += noop;
        case_counts[C2_DAG_EDGE_NO_CYCLE] += (long long)rev.size();
    }
    // With blocks on demand every insertion goes into the store in one bulk
    // write, and only the structural steps stay sequential (none of them
    // reads edges another insertion of the batch could add across a label:
    // the WCC relabel walks one old label, the condensation steps use the
    // DAG maps). Eager block maintenance reads the store while it merges
    // blocks, so there the other insertions are stored one by one as before.
    if (blocks_on_demand) {
        storeInsertBulk(ins);
        for (size_t i = 0; i < ins.size(); ++i)
            if (!ins_noop[i]) handleInsert(ins[i].first, ins[i].second, und_new[i]);
    } else {
        std::vector<std::pair<int, int>> bulk;
        parCompact<1, std::pair<int, int>>(ins.size(), bulk, [&](size_t i, auto put) { if (ins_noop[i]) put(ins[i]); });
        storeInsertBulk(bulk);
        for (size_t i = 0; i < ins.size(); ++i) {
            if (ins_noop[i]) continue;
            bool und;
            if (!storeInsert(ins[i].first, ins[i].second, und)) continue;
            handleInsert(ins[i].first, ins[i].second, und);
        }
    }
}

// ============================================================================
// Getters and invariant checks
// ============================================================================

std::map<int, int> GIR::getNodeToScc() const {
    std::map<int, int> m;
    for (int v = 0; v < V; ++v) m[v] = node_scc[v];
    return m;
}
std::map<int, int> GIR::getNodeToWcc() const {
    std::map<int, int> m;
    for (int v = 0; v < V; ++v) m[v] = wcc_label[v];
    return m;
}
std::vector<std::vector<int>> GIR::getDAG() const {
    std::vector<std::vector<int>> r(dag_out.size());
    for (size_t i = 0; i < dag_out.size(); ++i)
        for (auto& kv : dag_out[i]) r[i].push_back(kv.first);
    return r;
}
std::vector<std::vector<int>> GIR::getWCCs() const {
    std::unordered_map<int, std::vector<int>> g;
    for (int v = 0; v < V; ++v) g[wcc_label[v]].push_back(v);
    std::vector<std::vector<int>> r;
    for (auto& kv : g) r.push_back(std::move(kv.second));
    return r;
}
std::vector<std::vector<int>> GIR::getBCTComponents(int s) const {
    std::vector<std::vector<int>> r;
    if (s < 0 || (size_t)s >= sccs.size() || !scc_alive[s]) return r;
    std::unordered_set<int> seen;
    for (int x : sccs[s]) {
        int b = findBlk(v_block[x]);
        if (seen.insert(b).second) r.push_back(blk_verts[b]);
    }
    return r;
}
BlockCutTree* GIR::getBCT(int s) const {
    auto comps = getBCTComponents(s);
    std::unordered_set<int> aps;
    for (auto& c : comps) for (int x : c) if (v_bcnt[x] >= 2) aps.insert(x);
    auto* t = new BlockCutTree();
    t->build(comps, aps);
    return t;
}

bool GIR::checkBlockForestInvariants(std::string* why) const {
    auto fail = [&](const std::string& m) { if (why) *why = m; return false; };
    // Global topological order: distinct ranks, every condensation edge
    // goes up in rank, rank_to_scc indexes exactly the live SCCs.
    size_t alive = 0;
    for (size_t s = 0; s < sccs.size(); ++s) {
        if (!scc_alive[s]) continue;
        ++alive;
        auto it = rank_to_scc.find(scc_rank[s]);
        if (it == rank_to_scc.end() || it->second != (int)s) return fail("rank_to_scc out of sync");
        for (auto& kv : dag_out[s]) {
            if (!scc_alive[kv.first]) return fail("DAG edge to dead SCC");
            if (scc_rank[s] >= scc_rank[kv.first]) return fail("DAG edge violates rank order");
        }
    }
    if (alive != rank_to_scc.size()) return fail("rank_to_scc has stale entries");
    for (size_t s = 0; s < sccs.size(); ++s) {
        if (!scc_alive[s]) continue;
        std::unordered_map<int, int> occ;
        std::unordered_set<int> blocks;
        for (int x : sccs[s]) {
            if (node_scc[x] != (int)s) return fail("node_scc mismatch");
            int b = findBlk(v_block[x]);
            blocks.insert(b);
        }
        int roots = 0;
        for (int b : blocks) {
            if (blk_parent_ap[b] == -1) ++roots;
            std::unordered_set<int> uniq(blk_verts[b].begin(), blk_verts[b].end());
            if (uniq.size() != blk_verts[b].size()) return fail("duplicate vertex in block");
            for (int x : blk_verts[b]) {
                if (node_scc[x] != (int)s) return fail("block vertex outside its SCC");
                if (!inBlock(x, b)) return fail("block lists a vertex whose pointers disagree");
                occ[x]++;
            }
            int p = blk_parent_ap[b];
            if (p >= 0) {
                if (!uniq.count(p)) return fail("parent AP not in block");
                if (findBlk(v_block[p]) == b) return fail("AP's parent is its own child");
            }
            // climb to root, must terminate
            int cur = b;
            for (size_t steps = 0; blk_parent_ap[cur] != -1; ++steps) {
                if (steps > blocks.size()) return fail("cycle in block-cut tree");
                cur = findBlk(v_block[blk_parent_ap[cur]]);
            }
        }
        if (roots != 1) return fail("SCC " + std::to_string(s) + " has " + std::to_string(roots) + " roots");
        for (int x : sccs[s])
            if (occ[x] != v_bcnt[x]) return fail("v_bcnt mismatch at vertex " + std::to_string(x));
    }
    return true;
}

