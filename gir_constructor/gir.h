#ifndef GIR_H
#define GIR_H

#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>

class graph;               // 前置声明，用户已有 graph 类
class graph_undirected;    // 前置声明，用户已有 graph_undirected 类
class BlockCutTree;        // 前置声明，已有 BlockCutTree 类

// GIR maintains three structural layers (WCC / SCC+condensation-DAG / per-SCC
// BiCC+BC-tree) and keeps them consistent under single-edge insertions and
// deletions via the paper's Scope-Prune-Update (SPU) classification:
// see paper/maintenance.tex and paper/single_edge_update.tex for the five
// cases this class implements.
//
// Design notes (see docs/code-review.md and the design discussion that
// produced this rewrite):
//  - The base `graph`/`graph_undirected` CSR is immutable once loaded from
//    disk. Edge insertions/deletions are tracked in a small overlay
//    (extra_*_adj / removed_*) that every traversal helper consults
//    alongside the base CSR, so the base graph never needs to be resized.
//  - The vertex set is fixed at construction time (matches the paper's
//    stated scope: vertex insertions without edges don't affect
//    connectivity, and vertex deletions reduce to bulk edge deletions) --
//    every vertex is assumed to already have WCC/SCC/BiCC membership before
//    insertEdge/deleteEdge is ever called.
//  - WCC merge/split uses a real union-find (O(alpha(n)) merge); split is a
//    bounded lockstep BFS from both sides so cost is bounded by the smaller
//    resulting component.
//  - SCC layer keeps a topological "rank" per SCC with sparse spacing so
//    Case-2 cycle checks can be bounded to a rank interval; when a merge or
//    split needs more room than the interval currently has, the affected
//    interval (not the whole DAG) is locally re-topo-sorted and
//    renumbered. A full global renumber is a rare fallback, not the common
//    path -- see repairRankInterval/renumberAllRanks.
//  - BiCC/BC-tree is maintained *per SCC* (matching the paper's formal
//    model in task_parallel.tex: BiCC blocks are computed on the induced
//    undirected subgraph of each SCC, not on the whole graph). Repair is a
//    from-scratch recompute scoped to just the touched SCC's vertex set
//    (dynamic_repair::find_bicc_scoped), not a textbook fully-incremental
//    biconnectivity algorithm -- this is a deliberate trade-off: the higher
//    layers already bound the affected region, so a scoped recompute meets
//    the same O(|V_a|+|E_a|) bound with far less implementation risk.
class GIR {
public:
    GIR(const std::string& path_directed,
        const std::string& path_undirected);
    ~GIR();

    // 核心构建接口
    void buildWCC();
    void buildBCT();
    void buildSCC_DAG();

    // Getter (materialized on demand from the internal fast representation)
    std::map<int,int> getNodeToScc() const;
    std::map<int,int> getNodeToWcc() const;
    const std::vector<std::vector<int>>& getSCCs() const;
    std::vector<std::vector<int>> getDAG() const;
    std::vector<std::vector<int>> getWCCs() const;
    graph* getDirectedGraph() const;
    graph_undirected* getUndirectedGraph() const;
    // BC tree / BiCC blocks are per-SCC now (see class comment). Pass the
    // SCC id a vertex belongs to (getNodeToScc()[v]).
    BlockCutTree* getBCT(int scc_id) const;
    const std::vector<std::vector<int>>& getBCTComponents(int scc_id) const;

    // Single logical edge update: updates both the directed and the
    // undirected view of (u,v) atomically and classifies/repairs via SPU.
    // No-op if the edge is already present (insert) / already absent
    // (delete).
    void insertEdge(int u, int v);
    void deleteEdge(int u, int v);

    // Kept for older call sites; both route through insertEdge/deleteEdge
    // so the directed and undirected views can never desynchronize.
    void updateEdgeDirected(int u, int v, bool insert);
    void updateEdgeUndirected(int u, int v, bool insert);

    // Read-only access to the current (base CSR + overlay) graph view --
    // useful for tests/tools that need to independently recompute ground
    // truth to check GIR's incrementally-maintained state against.
    std::vector<int> getDirectedOutNeighbors(int v) const;
    std::vector<int> getDirectedInNeighbors(int v) const;
    std::vector<int> getUndirectedNeighborsPublic(int v) const;
    bool hasDirectedEdge(int u, int v) const;
    bool hasUndirectedEdge(int u, int v) const;
    int vertexCount() const;

    // GIR's two-level parallel BFS (paper/task_parallel.tex): Stage 1 prunes
    // to SCCs reachable from source via the condensation DAG; Stages 2/3 are
    // interleaved as a block-granularity wavefront (task parallelism across
    // blocks that become reachable in the same round, data parallelism via
    // OpenMP within each block's own local BFS); Stage 4 (global level =
    // entry level + local level) is folded directly into Stage 3 here.
    // Returns per-vertex BFS level, -1 for unreached vertices.
    std::vector<int> bfsTaskParallel(int source);

private:
    std::string base_path_directed;
    std::string base_path_undirected;

    // 图对象 (immutable base CSR, loaded once)
    graph* g_directed = nullptr;
    graph_undirected* g_undirected = nullptr;

    // ---------------- Dynamic edge overlay ----------------
    std::unordered_map<int, std::vector<int>> extra_out_adj, extra_in_adj;
    std::unordered_set<long long> removed_directed;
    std::unordered_map<int, std::vector<int>> extra_undirected_adj;
    std::unordered_set<long long> removed_undirected;
    // The base g_undirected CSR is not guaranteed to be symmetric in
    // practice (a known upstream bug: txt_to_bin_undirected.cpp never
    // symmetrizes edges, see docs/code-review.md). base_undirected_adj is
    // built once in buildBCT() by scanning the raw base CSR and adding both
    // directions of every edge, so every later query only has one place to
    // consult instead of re-deriving symmetry (expensively) on every call.
    std::vector<std::vector<int>> base_undirected_adj;

    static long long encodeDirected(int u, int v);
    static long long encodeUndirected(int u, int v);
    bool isDirectedEdgePresent(int u, int v) const;
    bool isUndirectedEdgePresent(int u, int v) const;
    std::vector<int> directedOutNeighbors(int v) const;
    std::vector<int> directedInNeighbors(int v) const;
    std::vector<int> undirectedNeighbors(int v) const;
    void addOverlayDirectedEdge(int u, int v);
    void removeOverlayDirectedEdge(int u, int v);
    void addOverlayUndirectedEdge(int u, int v);
    void removeOverlayUndirectedEdge(int u, int v);

    // ---------------- WCC layer: union-find ----------------
    std::vector<int> wcc_parent, wcc_size;
    int findWcc(int v);
    void unionWcc(int a, int b);

    // ---------------- SCC / condensation DAG layer ----------------
    std::map<int,int> nodeToScc;
    std::vector<std::vector<int>> sccs;
    std::vector<char> scc_alive;
    // scc id -> {successor/predecessor scc id -> number of underlying
    // vertex-level edges currently contributing to that condensation edge}.
    // Reference-counted rather than a plain set: two SCCs can be connected
    // by more than one real edge (common once enough random updates
    // accumulate), and deleting just one of them must not remove the
    // condensation-DAG edge while another still justifies it.
    std::vector<std::unordered_map<int, int>> dag_adj_set;
    std::vector<std::unordered_map<int, int>> dag_radj_set;
    std::vector<long long> scc_rank;
    std::map<long long,int> rank_to_scc;
    static const long long RANK_GAP = 1LL << 40;

    void addDagEdge(int from_scc, int to_scc);
    void removeDagEdge(int from_scc, int to_scc);
    void renumberAllRanks();
    void repairRankInterval(long long lo_rank, long long hi_rank);
    int newSccId();
    void shiftComponentBelow(int start_scc, long long must_be_below);
    void shiftComponentAbove(int start_scc, long long must_be_above);

    // ---------------- BiCC / BC-tree layer, per SCC ----------------
    std::vector<std::vector<std::vector<int>>> scc_bicc_blocks;
    std::vector<std::unordered_set<int>> scc_ap;
    std::vector<BlockCutTree*> scc_bct;

    void rebuildSccBicc(int scc_id);
    bool isBridgeWithinScc(int scc_id, int u, int v) const;
    bool inSameBiccBlock(int scc_id, int u, int v) const;

    // ---------------- SPU dispatch (five-case classifier) ----------------
    void handleInsert(int u, int v);
    void handleDelete(int u, int v);

    void mergeSccs(const std::vector<int>& scc_ids_to_merge);
    void splitScc(int scc_id);
    std::vector<int> sccVertices(int scc_id) const;

    // Lockstep bidirectional BFS over the undirected view; returns whether
    // u and v are still connected. If not, `smaller_side_out` receives the
    // exact vertex set of whichever side finished exploring first.
    bool probeUndirectedConnectivity(int u, int v, std::vector<int>& smaller_side_out);
    // Rebuilds union-find state for exactly `members` as one fresh
    // component (used after a confirmed WCC split).
    void rebuildWccComponent(const std::vector<int>& members);
    // Plain flood-fill over the condensation graph (dag_adj_set/dag_radj_set
    // treated as undirected) to check whether su can still reach sv.
    bool probeCondensationConnectivity(int su, int sv);

    void generate_and_save_dag(
        const std::map<int, int>& nodeToScc,
        const std::unordered_map<int, std::vector<int>>& adj_list,
        const std::string& output_filename);
};

#endif // GIR_H
