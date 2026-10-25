#ifndef GIR_H
#define GIR_H

#include <algorithm>
#include <array>
#include <memory>
#include <type_traits>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class graph;
class graph_undirected;
class BlockCutTree;

// Allocator that leaves elements uninitialized (default-init): query results
// are V-sized and filled in parallel right after allocation, so the zero fill
// of std::allocator (sequential, ~8 ms for 4.8M ints, mostly page faults on
// one thread) is pure overhead. Ligra/Terrace allocate their per-query arrays
// the same way.
template <class T>
struct DefaultInitAllocator : std::allocator<T> {
    template <class U> struct rebind { using other = DefaultInitAllocator<U>; };
    DefaultInitAllocator() = default;
    template <class U> DefaultInitAllocator(const DefaultInitAllocator<U>&) noexcept {}
    template <class U> void construct(U* p) noexcept(std::is_nothrow_default_constructible<U>::value) { ::new ((void*)p) U; }
    template <class U, class... A> void construct(U* p, A&&... a) { ::new ((void*)p) U(std::forward<A>(a)...); }
};
template <class T> using QueryVec = std::vector<T, DefaultInitAllocator<T>>;
// (found by ADL: the allocator lives in this namespace)
template <class T> bool operator==(const QueryVec<T>& a, const std::vector<T>& b) { return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin()); }
template <class T> bool operator==(const std::vector<T>& a, const QueryVec<T>& b) { return b == a; }
template <class T> bool operator!=(const QueryVec<T>& a, const std::vector<T>& b) { return !(a == b); }
template <class T> bool operator!=(const std::vector<T>& a, const QueryVec<T>& b) { return !(b == a); }

// GIR maintains three structural layers over a directed graph G and keeps
// them consistent under edge insertions and deletions:
//
//   WCC  : weakly connected components of G.
//   SCC  : strongly connected components plus the condensation DAG, with one
//          *global* topological rank per SCC (Pearce-Kelly maintenance).
//   BiCC : for every SCC C, the biconnected blocks of the undirected graph
//          induced by C, organized as a rooted block-cut forest.
//
// Invariants that the update code relies on (and the differential test in
// gir_test.cpp checks):
//  - The undirected view is *always* the symmetrization of the current
//    directed edge set: {u,v} exists iff u->v or v->u exists. It is never
//    loaded independently from disk (the on-disk undirected_*.bin files are
//    not trustworthy, see docs/code-review.md), so it cannot drift.
//  - Ranks are globally distinct and every condensation edge a->b satisfies
//    rank(a) < rank(b). Because the order is global, a cross-WCC insertion
//    needs no special rank rebasing: it is the no-cycle case of the same
//    Pearce-Kelly reorder used for intra-WCC insertions.
//  - Block-cut forest representation: every vertex x stores v_block[x], the
//    (union-find) id of its *parent* block, and v_bcnt[x], the number of
//    blocks containing it (x is an articulation point iff v_bcnt[x] >= 2).
//    Every block stores the articulation point it hangs from (-1 for a
//    root). A vertex belongs to block B iff find(v_block[x]) == B or
//    blk_parent_ap[B] == x. Blocks merged along a BC-tree path are unioned
//    in blk_uf, so merging never has to touch vertices outside the path.
//
// Two update entry points:
//  - insertEdge / deleteEdge: one edge at a time (Scope-Prune-Update).
//  - applyBatch: Filter-Merge-Update. Cancelling insert/delete pairs and
//    duplicates are dropped; all deletions are applied first with SCC split
//    checks and block repairs deferred to once per touched SCC/block; then
//    insertions run through the single-edge insert path.
class GIR {
public:
    struct Update {
        int u, v;
        bool insert;
    };

    // Indices into caseCounts().
    enum Case {
        C1_CROSS_WCC = 0,      // insert, different WCCs
        C2_DAG_EDGE_NO_CYCLE,  // insert, same WCC, different SCCs, no cycle
        C2_SCC_MERGE,          // insert, closes a cycle, SCCs merged
        C3_BLOCK_MERGE,        // insert, same SCC, different blocks
        C4_NOOP,               // insert, same block (or edge already present)
        D_INTER_SCC,           // delete of a condensation edge, WCC intact
        D_WCC_SPLIT,           // delete of a condensation edge, WCC split
        D_INTRA_SCC_KEEP,      // delete inside an SCC, SCC stays strongly connected
        D_SCC_SPLIT,           // delete inside an SCC, SCC split
        NUM_CASES
    };

    GIR(const std::string& path_directed, const std::string& path_undirected);
    // Takes ownership of an already-built graph (e.g. graph::fromEdges).
    explicit GIR(graph* g);
    ~GIR();

    // Construction. buildWCC, then buildSCC_DAG, then buildBCT.
    void buildWCC();
    void buildSCC_DAG();
    void buildBCT();

    // Updates.
    void insertEdge(int u, int v);
    void deleteEdge(int u, int v);
    void applyBatch(const std::vector<Update>& batch);
    // Deferred maintenance. ingestBatch applies a batch to the edge store
    // only (parallel, like a plain graph store); the WCC/SCC/block layers lag
    // behind until sync(), which applies the net change since the last sync
    // through applyBatch. While changes are pending, the structure getters
    // (sccOf, getSCCs, ...) describe the graph as of the last sync; the
    // updates above and the structure-based queries sync first. BFS, SSSP
    // and BC read only the edge store and see every ingested change.
    void ingestBatch(const std::vector<Update>& batch);
    void sync();
    bool synced() const { return pending.empty(); }
    // Kept for older call sites; both route through insertEdge/deleteEdge.
    void updateEdgeDirected(int u, int v, bool insert);
    void updateEdgeUndirected(int u, int v, bool insert);

    // Cheap per-vertex queries.
    int sccOf(int v) const { return node_scc[v]; }
    int wccOf(int v) const { return wcc_label[v]; }
    long long sccRank(int s) const { return scc_rank[s]; }
    bool isArticulationPoint(int v) const { return v_bcnt[v] >= 2; }

    // Structural queries answered from the index (sync first; the block ones
    // also ensureBlocks()). reachable(u, v): u ~> v via SCC ids, topological
    // ranks and a rank-pruned search on the condensation DAG.
    bool reachable(int u, int v);
    bool sameScc(int u, int v) { sync(); return node_scc[u] == node_scc[v]; }
    bool sameWcc(int u, int v) { sync(); return wcc_label[u] == wcc_label[v]; }
    bool articulationPoint(int v) { sync(); ensureBlocks(); return v_bcnt[v] >= 2; }
    // u and v on a common cycle of the undirected view of one SCC (same block).
    bool biconnected(int u, int v) { sync(); ensureBlocks(); return u == v || commonBlock(u, v) >= 0; }
    // The same reachability answer without the index: bidirectional BFS on
    // the edge store (what a plain graph store can do per query).
    bool reachableNoIndex(int u, int v);
    int vertexCount() const { return V; }

    // Materialized views (allocate; meant for tests and tools).
    std::map<int, int> getNodeToScc() const;
    std::map<int, int> getNodeToWcc() const;
    const std::vector<std::vector<int>>& getSCCs() const { return sccs; }
    std::vector<std::vector<int>> getDAG() const;
    std::vector<std::vector<int>> getWCCs() const;
    // BiCC blocks (vertex sets) of one SCC.
    std::vector<std::vector<int>> getBCTComponents(int scc_id) const;
    // Builds a BlockCutTree object for one SCC on demand (caller owns it).
    BlockCutTree* getBCT(int scc_id) const;
    // Checks the internal block-cut forest invariants of every live SCC
    // (parent pointers, AP counts, membership). Test-only; O(V + blocks).
    bool checkBlockForestInvariants(std::string* why = nullptr) const;

    graph* getDirectedGraph() const { return g_directed; }
    // The raw undirected CSR from disk, loaded lazily. Only the legacy
    // ParallelBFS tests use it; GIR itself never reads it.
    graph_undirected* getUndirectedGraph() const;

    std::vector<int> getDirectedOutNeighbors(int v) const;
    std::vector<int> getDirectedInNeighbors(int v) const;
    std::vector<int> getUndirectedNeighborsPublic(int v) const;
    bool hasDirectedEdge(int u, int v) const { return isDirectedEdgePresent(u, v); }
    bool hasUndirectedEdge(int u, int v) const { return isUndirectedEdgePresent(u, v); }

    const std::array<long long, NUM_CASES>& caseCounts() const { return case_counts; }
    // Wall time (seconds) spent in single-edge insertEdge/deleteEdge calls,
    // attributed to the case each call landed in.
    const std::array<double, NUM_CASES>& caseSeconds() const { return case_seconds; }
    // Blocks at least this large try the path-certificate shortcuts before a
    // full block recompute on deletion (tests set 0 to exercise them).
    void setFastPathMinBlock(size_t n) { fastpath_min_block = n; }
    // Condensation degree above which order-maintenance searches treat a
    // node as a hub (tests set it tiny to exercise the hub path).
    void setHubDegree(size_t d) { hub_degree = d; }
    // Bounded searches don't expand vertices above this degree (tests shrink it).
    void setSearchHubDegree(long long d) { search_hub_deg = d; }
    // SCCs from this size on use hub anchors in bulkTrimSplit (tests: 2).
    void setHubAnchorMin(size_t n) { hub_anchor_min = n; }
    // Tests: never use a single hub anchor (forces the anchor-set path).
    void setForceAnchorSet(bool on) { force_anchor_set = on; }
    // Tests: edge budget of the sequential hub/anchor probes (small values
    // send most probes to the parallel second pass).
    void setSeqProbeBudget(long long b) { seq_probe_budget = b; }
    // Test knob: skip the quick block certificate so every intra-block
    // deletion goes through the exact Menger search.
    void setForceExactBlockCert(bool f) { force_exact_cert = f; }
    // Test knob: skip the split fast path, always use the (incremental) resplit.
    void setForceResplit(bool f) { force_resplit = f; }
    // Block layer maintenance. Eager (default): every update keeps the
    // per-SCC block-cut forests exact. On demand: updates maintain WCC, SCC
    // and the condensation DAG only and record which SCCs' blocks went stale;
    // ensureBlocks() rebuilds exactly those (the block getters below describe
    // the last ensureBlocks() until then; bfsTwoLevel calls it itself).
    void setBlocksOnDemand(bool on) { if (!on) ensureBlocks(); blocks_on_demand = on; }
    void setParallelBlocksMin(size_t n) { par_bcc_min = n; }  // (tests shrink it)
    bool blocksOnDemand() const { return blocks_on_demand; }
    void ensureBlocks();
    void resetCaseCounts() { case_counts.fill(0); case_seconds.fill(0); }

    // ---------------- Downstream analytics (gir_query.cpp) ----------------
    // BFS levels (-1 = unreachable). bfs: direction-optimizing, bottom-up
    // steps scan only the vertices of SCCs reachable in the condensation DAG.
    QueryVec<int> bfs(int source);
    std::vector<char> reachableSccsPublic(int c) const { return reachableSccs(c); }
    // Paper Sec. 6 two-level model: wavefront over the condensation DAG
    // (independent SCCs as concurrent tasks), BC-tree block tasks inside the
    // source SCC, data parallelism inside large units.
    QueryVec<int> bfsTwoLevel(int source);
    QueryVec<int> bfsTaskParallel(int source);  // = bfs (kept for old callers)
    // SSSP with deterministic pseudo-random weights edgeWeight(u,v) in
    // [1, kMaxWeight]; UINT64_MAX = unreachable.
    static const uint32_t kMaxWeight = 100;
    static uint32_t edgeWeight(int u, int v);
    QueryVec<uint64_t> sssp(int source, uint64_t delta = 16);
    QueryVec<uint64_t> ssspTwoLevel(int source, uint64_t delta = 16);
    // Brandes betweenness accumulated over `sources` (unweighted), sources
    // processed as parallel tasks.
    QueryVec<double> bc(const std::vector<int>& sources);
    // Incremental single-source BFS (weighted = false) / SSSP (true), in the
    // style of KickStarter: trackIncremental computes the result for one
    // source once and from then on every edge change is logged; update
    // corrects the result for the logged changes (deleted tree edges
    // invalidate the subtree below them, which is re-seeded from its valid
    // in-neighbours; inserted edges relax; then parallel propagation).
    // incDist: distance (BFS: level) of v, UINT64_MAX = unreachable.
    void trackIncremental(int source, bool weighted);
    void updateIncremental(bool weighted);
    uint64_t incDist(bool weighted, int v) const {
        uint64_t k = inc[weighted].st[v];
        return k == ~0ULL ? UINT64_MAX : k >> 32;
    }
    // Granularity knobs (tests shrink them to exercise every code path):
    // frontiers below par_frontier run in the calling thread; SCCs/blocks at
    // least large_unit big use data parallelism, smaller ones are tasks, and
    // reachable regions at least large_unit big may use bottom-up BFS.
    void setQueryGranularity(size_t par_frontier, size_t large_unit) {
        q_par_frontier = par_frontier;
        q_large_unit = large_unit;
    }

private:
    std::string base_path_directed, base_path_undirected;
    graph* g_directed = nullptr;
    mutable graph_undirected* g_undirected = nullptr;
    int V = 0;

    // ---------------- Edge store: immutable base CSR + per-vertex deltas ----------------
    // The base CSR is sorted and deduplicated at load time, so presence tests
    // are binary searches. Insertions that are not base edges go into the
    // per-vertex delta buffers (paper Sec. 7, "Implementation"); deletions
    // of base edges are tombstoned.
    std::vector<std::vector<int>> extra_out_adj, extra_in_adj, extra_und_adj;
    long long extra_count = 0;  // number of directed delta entries
    // Presence index for delta edges (hub vertices accumulate long delta
    // lists, so a linear scan per presence test is too slow).
    // Each delta list is kept sorted: presence of a delta edge is a binary
    // search, and a batch merges into (or is subtracted from) each vertex's
    // list independently, so batches need no global index.
    // Tombstones for deleted base edges: one bit per base CSR entry (fw, bw,
    // symmetrized), so neighbour scans test a bit instead of hashing.
    std::vector<uint64_t> dead_fw, dead_bw, dead_sym;
    long long n_dead_dir = 0, n_dead_und = 0;
    static bool bitGet(const std::vector<uint64_t>& m, long long i) { return (m[(size_t)i >> 6] >> (i & 63)) & 1; }
    static void bitPut(std::vector<uint64_t>& m, long long i, bool on) {
        if (on) m[(size_t)i >> 6] |= 1ULL << (i & 63); else m[(size_t)i >> 6] &= ~(1ULL << (i & 63));
    }
    long long fwPos(int u, int v) const;   // index of u->v in fw_csr, or -1
    long long bwPos(int u, int v) const;   // index of u in bw list of v, or -1
    long long symPos(int u, int v) const;  // index of v in u's sym list, or -1
    bool deadD(int u, int v) const { if (!n_dead_dir) return false; long long p = fwPos(u, v); return p >= 0 && bitGet(dead_fw, p); }
    bool deadU(int u, int v) const { if (!n_dead_und) return false; long long p = symPos(u, v); return p >= 0 && bitGet(dead_sym, p); }
    // Symmetrized, sorted, deduplicated base adjacency (int64 offsets: the
    // symmetrized twitter graph has more than 2^31 entries).
    std::vector<long long> sym_beg;
    std::vector<int> sym_adj;

    static long long keyD(int u, int v) { return ((long long)u << 32) | (unsigned)v; }
    static long long keyU(int u, int v) { return u < v ? keyD(u, v) : keyD(v, u); }
    bool isDirectedEdgePresent(int u, int v) const;
    bool isUndirectedEdgePresent(int u, int v) const;
    void normalizeBaseCsr();
    void initAfterLoad();
    bool baseHasDirected(int u, int v) const;
    bool baseHasUndirected(int u, int v) const;
    void overlayAddDirected(int u, int v);
    void overlayRemoveDirected(int u, int v);
    void overlayAddUndirected(int u, int v);
    void overlayRemoveUndirected(int u, int v);
    // Apply an edge change to the edge store only. Returns false if it was
    // a no-op (already present / already absent). `und_changed` reports
    // whether the undirected edge {u,v} appeared / disappeared.
    bool storeInsert(int u, int v, bool& und_changed);
    bool storeDelete(int u, int v, bool& und_changed);
    // Insert distinct, currently absent directed edges at once (parallel,
    // grouped by vertex). Only for edges that change no layer.
    void storeInsertBulk(const std::vector<std::pair<int, int>>& es);
    // Delete distinct, currently present directed edges at once (parallel).
    void storeDeleteBulk(const std::vector<std::pair<int, int>>& es);

    template <class F> void forOut(int v, F&& f) const;
    template <class F> void forIn(int v, F&& f) const;
    template <class F> void forUnd(int v, F&& f) const;
    // Like forUnd, but stops as soon as f returns true; returns whether it did.
    template <class F> bool anyUnd(int v, F&& f) const;
    template <class F> bool anyOut(int v, F&& f) const;
    template <class F> bool anyIn(int v, F&& f) const;
    // Resumable neighbor cursors for the iterative DFS routines.
    struct Cursor { long long i, end; size_t j; };
    Cursor outBegin(int v) const;
    bool outNext(int v, Cursor& c, int& w) const;
    Cursor undBegin(int v) const;
    bool undNext(int v, Cursor& c, int& w) const;

    // ---------------- WCC: explicit labels, relabel smaller side ----------------
    std::vector<int> wcc_label;      // vertex -> label
    std::vector<int> wcc_size;       // label -> vertex count (0 = unused)
    void relabelComponent(int start, int from_label, int to_label);
    // Lockstep BFS from u and v over the undirected view. Returns true if
    // they are connected; otherwise `side` gets the side that was exhausted.
    bool undirectedConnected(int u, int v, std::vector<int>& side);

    // ---------------- SCC + condensation DAG ----------------
    std::vector<int> node_scc;
    std::vector<std::vector<int>> sccs;
    std::vector<char> scc_alive;
    std::vector<int> free_scc_ids;
    std::vector<std::unordered_map<int, int>> dag_out, dag_in;  // ref-counted
    std::vector<long long> scc_rank;
    std::map<long long, int> rank_to_scc;
    static const long long RANK_GAP = 1LL << 32;

    int newSccId();
    void retireScc(int id);
    void addDagEdge(int a, int b);
    void removeDagEdge(int a, int b);
    bool dagEdgeExists(int a, int b) const;
    void setRank(int scc, long long r);
    void renumberAllRanks();
    // Local order-maintenance respread: re-space a window of ranks around
    // key r (doubling the window until it is sparse enough) so that at least
    // `need` free ranks exist on each side of r. Falls back to renumbering
    // everything only if the window grows to the whole order.
    void respreadAround(long long r, long long need);
    // Position hints for O(1) removal from SCC / block vertex lists (a stale
    // hint falls back to a linear search).
    std::vector<int> scc_pos, blk_pos;
    void sccListRemove(int s, int x);
    void blkListRemove(int b, int x);
    void setBlkPos(int b);
    // Insert condensation edge a->b (different SCCs, same WCC). If it
    // violates the order, search forward from b (ranks <= rank a) and
    // backward from a (ranks >= rank b) in lockstep and stop at whichever
    // side is exhausted first; only that side is re-ranked, into the free
    // gap next to a (forward side) or b (backward side). Merges SCCs if a
    // cycle closes. Returns true if a merge happened.
    int insertDagEdgePK(int a, int b, bool cross_wcc = false);
    template <class InF, class InB>
    void hubMerge(int a, int b, int h, std::vector<int>& F, std::vector<int>& B, InF inF, InB inB);
    int chooseKeep(const std::vector<int>& M) const;
    std::vector<int> kahnOrder(const std::vector<int>& nodes);
    std::vector<int> pk_sf, pk_sb;  // search stamps over SCC ids
    int pk_ef = 0, pk_eb = 0;
    size_t hub_degree = 4096;       // condensation degree above which a node is a hub
    // Case 1 (a and b in different WCCs, so no cycle is possible): if the
    // edge violates the order, move only the smaller side's affected SCCs
    // (ancestors of a, or descendants of b) into a free rank gap next to the
    // other endpoint. Cost is bounded by the smaller WCC.
    void insertCrossWccEdge(int a, int b, bool move_a_side);
    // Merge the given SCC ids into `keep` (default: chooseKeep); returns it.
    int mergeSccs(const std::vector<int>& ids, int keep = -1);
    // Recompute SCCs inside one old SCC (after deletions). Returns true if
    // it split.
    bool resplitScc(int scc_id, bool incremental = false, std::pair<int, int> extra = {-1, -1});
    void detachPieces(int s, std::vector<std::vector<int>>& below, std::vector<std::vector<int>>& above,
                      std::vector<std::pair<int, int>> extras);
    void detachPieces(int s, std::vector<std::vector<int>>& below, std::vector<std::vector<int>>& above,
                      std::pair<int, int> extra) {
        detachPieces(s, below, above, extra.first >= 0 ? std::vector<std::pair<int, int>>{extra} : std::vector<std::pair<int, int>>{});
    }
    // Topological order of the SCCs of U (vertices marked node_scc == tmp,
    // taken out of SCC s) relative to the rest of s: pieces before the rest
    // go to `below`, after it to `above`.
    void orderPieces(int s, int tmp, const std::vector<int>& U, std::vector<std::vector<int>>& below,
                     std::vector<std::vector<int>>& above);
    // FMU: split large SCC s after a batch by trimming (see applyBatch).
    // `dels` are the batch's uncertified deletions inside s, `all_dels` every
    // deletion of the batch. Returns the trimmed vertices (empty: not applied).
    // `no_reach[i]`: dels[i] failed its reachability certificate (it may
    // split s); the others failed only the block certificate. `certified`:
    // the rest of s (s minus the result) was shown strongly connected
    // without any of its deletions, so they all keep the SCC.
    std::vector<int> bulkTrimSplit(int s, const std::vector<std::pair<int, int>>& dels,
                                   const std::vector<char>& no_reach,
                                   const std::vector<std::pair<int, int>>& all_dels,
                                   std::vector<std::pair<int, int>>& wcc_probe, bool& certified);
    std::pair<int, int> deletedEdgeAttachment(int u, int v);
    // Directed reachability u ~> v restricted to SCC `scc`, lockstep forward
    // from u / backward from v.
    // If not reachable and `side` is given, it receives the side the probe
    // exhausted: the forward set from u (*fwd = true; this is exactly the
    // new SCC of u) or the backward set to v (*fwd = false; exactly the new
    // SCC of v). See splitAfterDelete.
    bool reachableWithinScc(int u, int v, int scc, std::vector<int>* side = nullptr, bool* fwd = nullptr);
    // Single deletion u->v broke SCC s and the probe exhausted side X.
    // Fast path: if the rest R = s \ X is still strongly connected (checked
    // through X's few boundary vertices only), split into exactly {X, R}
    // without touching R's vertices, and remove X from R's blocks with path
    // certificates. Falls back to resplitScc otherwise.
    void splitAfterDelete(int s, int u, int v, std::vector<int>& X, bool X_is_forward);
    // Remove vertex set X (all non-APs) from the block-cut tree of their SCC,
    // which stays connected. Returns false (nothing changed) if a shortcut
    // doesn't apply; the caller then rebuilds.
    // `extra_att` = (x, y): the just-deleted undirected edge was an edge of
    // x's block, so y must count as an attachment vertex of x even though
    // the edge is no longer in the graph (-1 = none).
    bool removeVerticesFromBlocks(const std::vector<int>& X, const std::vector<std::pair<int, int>>& extra_att = {});

    // ---------------- Block-cut forest (per SCC) ----------------
    std::vector<int> v_block;          // vertex -> parent block (uf id)
    std::vector<int> v_bcnt;           // vertex -> #blocks containing it
    mutable std::vector<int> blk_uf;   // block id -> uf parent (path-compressed on lookup)
    std::vector<int> blk_parent_ap;    // root block id -> AP it hangs from, -1 = root
    std::vector<std::vector<int>> blk_verts;  // root block id -> vertices

    int findBlk(int b) const;
    int findBlkRO(int b) const;
    // Membership stamps for removed-vertex sets (two levels: split -> peel).
    std::vector<int> rm_stamp, rm2_stamp;
    int rm_epoch = 0, rm2_epoch = 0;
    int nextRmEpoch(std::vector<int>& st, int& ep);
    int commonBlockRO(int x, int y) const;
    int newBlock();
    bool inBlock(int x, int b) const;
    // Root id of the block containing both x and y, or -1.
    int commonBlock(int x, int y) const;
    // Installs `blocks` (vertex sets) as the new BiCC structure of a region.
    // base_cnt[x] (looked up via base_cnt_fn) is the number of blocks x keeps
    // outside the region; p is the AP the region's root hangs from (-1 = none).
    template <class BaseCnt>
    void installBlocks(std::vector<std::vector<int>>& blocks, BaseCnt base_cnt_fn, int p,
                       int reuse_id);
    // Full BiCC recompute of one SCC (used after a split).
    void rebuildSccBlocks(int scc_id);
    // Recompute one block after edge deletions inside it (SCC unchanged).
    void repairBlock(int b);
    // Single deletion of {u,v} inside block b (SCC unchanged). Tries two
    // sound shortcuts before falling back to repairBlock:
    //  (1) two internally vertex-disjoint u-v paths remain in b  =>  b is
    //      still biconnected (Menger), nothing to do;
    //  (2) u (or v) is left with a single neighbour w in b and two disjoint
    //      v-w paths avoid u  =>  b minus u is still biconnected and u hangs
    //      off w as a new bridge block.
    // skip_quick: the bounded two-path certificate is already known to fail.
    void repairBlockAfterDelete(int b, int u, int v, bool skip_quick = false);
    std::vector<int> blockPathOrdered(int b, int s, int t, long long budget,
                                      const std::unordered_set<int>* banned = nullptr);
    int mengerSplit(int b, int s, int t, std::vector<int>& chunk, int& cut, bool& from_s, long long budget,
                    const std::unordered_set<int>* banned = nullptr);
    // Carve the small side `chunk` (containing cut vertex `cut`) of a minimum
    // vertex cut out of block b after vertices `gone` leave it: recompute the
    // chunk's own blocks and hang them off `cut`. Returns false if b's parent
    // AP lies inside the chunk (caller then rebuilds).
    bool carveChunk(int b, const std::vector<int>& chunk, int cut, const std::unordered_set<int>& gone);
    void splitChunkOff(int b, const std::vector<int>& chunk, int cut);
    // Shared shortcut behind repairBlockAfterDelete / removeVerticesFromBlocks.
    // Block b loses the vertices `gone` (non-APs leaving the SCC) and/or one
    // deleted edge. Vertices left with a single neighbour are peeled off
    // (cascading) into bridge blocks; the remaining core is certified
    // 2-connected via attachment vertices (core vertices that lost an edge).
    // On success applies the change and returns true; otherwise changes
    // nothing and returns false.
    bool peelAndCertify(int b, const std::unordered_set<int>& gone, const std::vector<int>& seeds,
                        const std::vector<int>& extra_att);
    // Recompute block b minus `gone` from scratch (it stays connected).
    void recomputeBlockWithout(int b, const std::unordered_set<int>& gone);
    // Finds a path s..t inside block b (undirected view), avoiding `banned`
    // vertices and, if skip_direct, the direct edge {s,t}. Bidirectional
    // BFS; the search gives up after `budget` edge visits. Returns the
    // interior vertices in `interior`.
    bool blockPath(int b, int s, int t, const std::unordered_set<int>& banned, bool skip_direct,
                   long long budget, std::vector<int>& interior, std::vector<int>* ordered = nullptr);
    bool twoDisjointPaths(int b, int s, int t, const std::unordered_set<int>& banned, long long budget);
    // Thread-safe bounded searches (thread-local scratch, no shared stamps,
    // no union-find writes) for certificates checked in parallel. A false
    // answer may mean "budget exhausted"; callers fall back to exact searches.
    struct LocalSearch;
    enum class SearchMode { Directed, UndScc, UndBlock };
    long long undDeg(int x) const;  // undirected degree (base + delta)
    bool twoHop(int u, int v, int scc, int* mid = nullptr) const;  // directed u -> w -> v (w in scc if >= 0)
    // After deleting u->v inside SCC s: 0 = no u ~> v path found (s may
    // split), 1 = u ~> v found but the block certificate was not (or was not
    // asked for), 2 = u ~> v plus two disjoint undirected u-v paths in s.
    int certifyKeep(int u, int v, int s, bool need_block) const;
    // directed u ~> v (inside SCC scc if >= 0); on success *path (if given)
    // gets the interior vertices of one such path
    bool reachRO(int u, int v, int scc = -1, std::vector<int>* path = nullptr) const;
    // Hub anchors for bulkTrimSplit: the highest in- and out-degree vertices
    // of the base graph (picked once), so a probe x ~> anchor usually ends
    // after one or two hops.
    std::vector<int> hub_by_in, hub_by_out;
    // Anchor set for bulkTrimSplit, built once and kept: anc_list holds the
    // anchors; anc_in[w] = 1 + index of an anchor h with an edge w -> h (or
    // w == h), anc_out[w] = 1 + index of one with h -> w. The marks are only
    // hints: a probe that stops at w re-checks the edge and that h is still
    // in the SCC, so marks gone stale through later updates cost coverage,
    // never correctness.
    std::vector<int> anc_in, anc_out, anc_list;
    // Interior paths of passed anchor-pair probes (h ~> a0, a0 ~> h), kept
    // across batches; reused only after every vertex and edge is re-checked.
    std::unordered_map<unsigned long long, std::vector<int>> anc_pair_path;
    // Scratch for the parallel anchor search (V-sized): epoch stamp and BFS parent.
    std::vector<int> pa_stamp, pa_parent;
    int pa_epoch = 0;
    // WCC checks of a batch: per label, a ball around its representative
    // (ball_mark[v] = label's epoch) that every endpoint's search may stop at.
    std::vector<int> ball_mark;
    int ball_epoch = 0;
    // Level-synchronous BFS on all threads from x (dir 0: out-edges, 1:
    // in-edges, 2: undirected), inside SCC scc (any if < 0), until a vertex w
    // with hit(w). 1 = found (*path, if given: w, ..., x's successor), 0 =
    // exhausted. For the few searches that outgrow the sequential budgets.
    template <class Hit> int parBfs(int x, int dir, int scc, Hit hit, std::vector<int>* path = nullptr);
    // Incremental BFS/SSSP state: per vertex (dist << 32 | parent), all ones
    // = unreachable; `pending` collects edges changed since the last update
    // (net changes, logged where batches and single updates enter).
    struct IncState { int source = -1; QueryVec<uint64_t> st; std::vector<std::pair<int, int>> pending; };
    IncState inc[2];
    bool in_sync = false;  // sync() replays changes ingestBatch already logged
    void logChanges(const std::vector<std::pair<int, int>>& es) {
        if (in_sync) return;
        for (auto& s : inc) if (s.source >= 0) s.pending.insert(s.pending.end(), es.begin(), es.end());
    }
    std::vector<char> inc_flag;  // scratch for updateIncremental
    void pickHubCandidates();
    // Hubs a (in) and b (out) pay off as anchors of an SCC of size C only
    // when their neighbourhoods cover a good part of it (twitter: yes,
    // livejournal: no, there the local searches are cheaper).
    bool hubsCover(int a, int b, size_t C) const;
    // Bounded BFS inside scc from x along out-edges until a vertex with an
    // edge to h (fwd), or along in-edges until one with an edge from h.
    // 1 = reached (*path, if given: the vertices from x to the hit), 0 =
    // exhausted (no path), -1 = budget (edges scanned) hit.
    // Kept small: it only pays off where hubs are huge (twitter), and a miss
    // falls back to the general search.
    int hubProbe(int x, int h, bool fwd, int scc, long long budget, std::vector<int>* path = nullptr, bool big = false) const;
    bool reachUndRO(int u, int v) const;  // undirected, anywhere
    // Undirected search from x until a vertex of SCC `scc` (bounded,
    // read-only): x is then in that SCC's weak component.
    bool reachSccUndRO(int x, int scc) const;
    bool twoPathsRO(int s, int t, SearchMode mode, int filt, const std::unordered_set<int>* banned) const;
    std::vector<int> par_a, par_b;
    size_t fastpath_min_block = 2048;
    size_t q_par_frontier = 2048, q_large_unit = 4096;
    bool force_exact_cert = false, force_resplit = false;
    long long search_hub_deg = 2048;
    size_t hub_anchor_min = 1 << 16;
    bool force_anchor_set = false;
    long long seq_probe_budget = 1 << 14;
    bool blocks_on_demand = false;
    std::vector<int> blk_dirty_v;  // on demand: vertices whose SCC's blocks are stale
    void markBlocksDirty(int x) { if (blocks_on_demand) blk_dirty_v.push_back(x); }
    // Deferred maintenance: (edge key, present at the last sync?) for every
    // edge an ingested batch changed, in ingestion order.
    std::vector<std::pair<long long, char>> pending;
    // Last operation per edge wins; ins/del are the edges whose presence the
    // batch changes (first-appearance order). `before` (optional) gets
    // (key, presence before the batch) for each of them.
    void normalizeBatch(const std::vector<Update>& batch, std::vector<std::pair<int, int>>& ins,
                        std::vector<std::pair<int, int>>& del, long long& noops,
                        std::vector<std::pair<long long, char>>* before = nullptr);
    // Insert undirected edge {u,v} with u,v already in the same block-cut
    // tree but different blocks: collapse the tree path into one block.
    void mergeBlockPath(int u, int v);
    // Join two different block-cut trees with the new bridge {a,b}; b's tree
    // is re-rooted at b. Caller guarantees b's tree is the smaller one.
    void linkTrees(int a, int b);
    void rerootAt(int b);
    // Incremental BiCC for an SCC merge: the old SCCs' trees are joined by
    // processing every undirected edge between different old SCCs.
    void mergeBlockForests(const std::vector<int>& old_ids, int keep);

    // ---------------- Scratch space for traversals ----------------
    mutable std::vector<int> stamp_a, stamp_b;
    mutable int epoch_a = 0, epoch_b = 0;
    std::vector<int> scr_idx, scr_low;
    std::vector<char> scr_flag;
    int nextEpochA() const;
    int nextEpochB() const;

    // Iterative Tarjan SCC over `roots`, edges filtered by in_scope(w).
    template <class InScope, class Emit>
    void tarjanScc(const std::vector<int>& roots, InScope in_scope, Emit emit);
    // Iterative Hopcroft-Tarjan biconnected blocks (vertex sets), edges
    // filtered by keep_edge(v, w). Isolated vertices become singleton blocks.
    template <class KeepEdge, class Emit>
    void hopcroftTarjan(const std::vector<int>& roots, KeepEdge keep_edge, Emit emit);
    // Blocks of SCC s (its undirected view) in parallel: Tarjan-Vishkin over a
    // BFS spanning tree. Used by rebuildSccBlocks for SCCs of at least
    // par_bcc_min vertices.
    void parallelSccBlocks(int s);  // computes and installs s's block-cut forest
    std::vector<int> free_blk_ids;  // block ids retired by ensureBlocks, reused by newBlock
    size_t par_bcc_min = 1 << 16;

    // ---------------- Query helpers (gir_query.cpp) ----------------
    std::vector<char> reachableSccs(int src_scc) const;
    // Reachability cache for queries: valid while struct_version is
    // unchanged (bumped whenever an SCC or condensation edge appears or
    // disappears).
    unsigned long long struct_version = 1;
    unsigned long long block_version = 1;  // bumped by every block-cut forest mutation
    struct ReachCache {
        unsigned long long version = 0;
        int scc = -1;
        std::vector<char> sccs;     // reachable SCC mask
        std::vector<int> verts;     // reachable vertices
        long long edges = 0;        // their out-edges
        std::vector<int> pending;   // per SCC: #reachable predecessors (two-level driver)
    };
    struct KidsCache {
        unsigned long long sv = 0, bv = 0;
        int scc = -1;
        std::vector<std::pair<int, int>> kids;  // (parent AP, child block), sorted
    };
    mutable KidsCache kids_cache;
    std::vector<int> blk_stamp;  // per-vertex block membership stamp (scopedBfs on a block)
    int blk_epoch = 0;
    mutable ReachCache reach_cache;
    std::vector<int> round_mark;  // delta-stepping duplicate filter
    int round_epoch = 0;
    std::vector<int> settle_mark;  // delta-stepping: settled in the current bucket
    int settle_epoch = 0;
    const ReachCache& reachInfo(int src_scc) const;
    template <class Inner, class Relax> void dagWavefront(int source, Inner inner, Relax relax);
    void bfsInsideScc(int c, QueryVec<int>& dist);
    template <class InScope> void scopedBfs(const std::vector<int>& verts, InScope in_scope, QueryVec<int>& dist);
    std::vector<char> q_front;  // frontier bitmap for scopedBfs
    void bfsSourceSccByBlocks(int source, QueryVec<int>& dist);
    template <class InScope>
    void deltaStepping(QueryVec<uint64_t>& dist, std::vector<int> seeds, uint64_t delta, InScope in_scope, bool par);

    // ---------------- SPU dispatch ----------------
    void handleInsert(int u, int v, bool und_new);

    std::array<long long, NUM_CASES> case_counts{};
    std::array<double, NUM_CASES> case_seconds{};
public:
    // Why a single-deletion SCC split took the slow path (diagnostics):
    // [0] fast path used, [1] X too large, [2] boundary too large,
    // [3] rest not strongly connected, [4] block certificate failed (R rebuilt).
    std::array<long long, 5> split_stats{};
    std::array<double, 5> split_secs{};  // wall time per split path (same index as split_stats)
    // Order-maintenance path taken by violating condensation inserts:
    // [0] one-sided move, [1] hub merge, [2] full re-search, [3] cross-WCC fallback.
    std::array<long long, 4> order_stats{};
    // Sizes of the re-ranked sets: [0] total nodes moved, [1] max, [2] moves > 10k nodes.
    std::array<long long, 3> order_moved{};
    // Block repair after an intra-block deletion: [0] quick certificate,
    // [1] exact Menger certificate, [2] chain links split off, [3] fallback recompute.
    std::array<long long, 4> block_stats{};
    long long menger_nopath = 0, menger_nocut = 0, menger_budget = 0;
    long long carve_count = 0;
    long long anchor_set_calls = 0, hub_anchor_calls = 0;  // bulkTrimSplit rounds by anchor kind  // chunks carved out of a block when vertices leave it
private:
    void attributeTime(const std::array<long long, NUM_CASES>& before, double secs);
    // hint (from a batch's certification in the graph after the batch, which
    // is a subgraph of the current one): 0 none, 1 = u ~> v known (only the
    // block certificate failed), 2 = the bounded reach search failed,
    // 3 = u ~> v known (from the bulk split certificate; block not checked).
    void deleteEdgeImpl(int u, int v, int hint = 0);
};

#endif  // GIR_H
