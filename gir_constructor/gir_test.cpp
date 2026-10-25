#include <queue>
#include <cmath>
#include "gir.h"
#include "graph.h"
#include "graph_undirected.h"
#include "dynamic_repair.h"
#include <iostream>
#include <random>
#include <vector>
#include <limits>
#include <chrono>
#include <type_traits>
#include <atomic>
#include <algorithm>
#include <omp.h>
#include <map>
#include <unordered_set>
#include <string>
#include <cstdlib>
#include <set>
#include <utility>
#include <queue>
#include <fstream>
//#include "fw_bw.h" // 引入 fw_bfs 函数


namespace ParallelBFS {

    template<typename T>
    static bool compare_and_swap(T &x, const T &old_val, const T &new_val) {
        #if defined _OPENMP && defined __GNUC__
            return __sync_bool_compare_and_swap(&x, old_val, new_val);
        #else
            bool result = false;
            #pragma omp critical
            {
                if (x == old_val) {
                    x = new_val;
                    result = true;
                }
            }
            return result;
        #endif
    }
    
    // --- 原子操作 ---
    template<typename T, typename U>
    static T fetch_and_add(T &x, U inc) {
        if constexpr (std::is_same_v<T, std::atomic<size_t>> ||
                      std::is_same_v<T, std::atomic<unsigned long>> ||
                      std::is_same_v<T, std::atomic<int>> ||
                      std::is_same_v<T, std::atomic<long long>>) {
            return x.fetch_add(inc, std::memory_order_relaxed);
        } else {
            return __sync_fetch_and_add(&x, inc);
        }
    }

    // --- Bitmap ---
    class Bitmap {
    public:
        explicit Bitmap(size_t size) : size_(size) {
            uint64_t num_words = (size + kBitsPerWord - 1) / kBitsPerWord;
            start_ = new uint64_t[num_words](); // 值初始化为0
            end_ = start_ + num_words;
        }
        ~Bitmap() { delete[] start_; }
        void reset() { std::fill(start_, end_, 0); }
        void set_bit(size_t pos) { start_[word_offset(pos)] |= (1l << bit_offset(pos)); }
        void set_bit_atomic(size_t pos) {
            uint64_t old_val, new_val;
            do {
                old_val = start_[word_offset(pos)];
                new_val = old_val | (1l << bit_offset(pos));
            } while (!compare_and_swap(start_[word_offset(pos)], old_val, new_val));
        }
        bool get_bit(size_t pos) const { return (start_[word_offset(pos)] >> bit_offset(pos)) & 1l; }
        void swap(Bitmap &other) { std::swap(start_, other.start_); std::swap(end_, other.end_); }
    private:
        uint64_t *start_;
        uint64_t *end_;
        size_t size_;
        static const uint64_t kBitsPerWord = 64;
        static uint64_t word_offset(size_t n) { return n / kBitsPerWord; }
        static uint64_t bit_offset(size_t n) { return n & (kBitsPerWord - 1); }
    };

    // --- SlidingQueue & QueueBuffer ---
    template <typename T> class QueueBuffer;

    template <typename T>
    class SlidingQueue {
        T *shared;
        std::atomic<size_t> shared_in;
        size_t shared_out_start;
        size_t shared_out_end;
        friend class QueueBuffer<T>;
    public:
        explicit SlidingQueue(size_t shared_size) {
            shared = new T[shared_size];
            reset();
        }
        ~SlidingQueue() { delete[] shared; }
        void push_back(T to_add) { shared[shared_in++] = to_add; }
        bool empty() const { return shared_out_start == shared_out_end; }
        void reset() { shared_out_start = 0; shared_out_end = 0; shared_in = 0; }
        void slide_window() { shared_out_start = shared_out_end; shared_out_end = shared_in; }
        T* begin() const { return shared + shared_out_start; }
        T* end() const { return shared + shared_out_end; }
        size_t size() const { return end() - begin(); }
    };

    template <typename T>
    class QueueBuffer {
        size_t in;
        T *local_queue;
        SlidingQueue<T> &sq;
        const size_t local_size;
    public:
        explicit QueueBuffer(SlidingQueue<T> &master, size_t given_size = 16384)
            : sq(master), local_size(given_size) {
            in = 0;
            local_queue = new T[local_size];
        }
        ~QueueBuffer() { delete[] local_queue; }
        void push_back(T to_add) {
            if (in == local_size) flush();
            local_queue[in++] = to_add;
        }
        void flush() {
            if (in == 0) return;
            size_t copy_start = fetch_and_add(sq.shared_in, in);
            std::copy(local_queue, local_queue + in, sq.shared + copy_start);
            in = 0;
        }
    };

    // ============================================================================
    // 并行BFS核心逻辑
    // ============================================================================

    // --- Top-Down Step ---
    template<typename GraphPtr>
    int64_t TDStep(GraphPtr g, SlidingQueue<int>& queue, std::vector<int>& distances,
                   const std::map<int, int>& nodeToComponent, int source_comp_id) {
        std::atomic<int64_t> scout_count(0);

        #pragma omp parallel
        {
            QueueBuffer<int> lqueue(queue);
            #pragma omp for nowait
            for (auto it = queue.begin(); it < queue.end(); ++it) {
                int u = *it;
                int u_dist = distances[u];
                if constexpr (std::is_same_v<GraphPtr, graph*>) { // 有向图
                    for (index_t i = g->fw_beg_pos[u]; i < g->fw_beg_pos[u + 1]; ++i) {
                        int v = g->fw_csr[i];
                        if (nodeToComponent.count(v) > 0 && nodeToComponent.at(v) == source_comp_id && distances[v] == -1) {
                            if (compare_and_swap(distances[v], -1, u_dist + 1)) {
                                lqueue.push_back(v);
                                scout_count++;
                            }
                        }
                    }
                } else if constexpr (std::is_same_v<GraphPtr, graph_undirected*>) { // 无向图
                    for (index_t i = g->beg_pos[u]; i < g->beg_pos[u + 1]; ++i) {
                        int v = g->csr[i];
                        if (nodeToComponent.count(v) > 0 && nodeToComponent.at(v) == source_comp_id && distances[v] == -1) {
                            if (compare_and_swap(distances[v], -1, u_dist + 1)) {
                                lqueue.push_back(v);
                                scout_count++;
                            }
                        }
                    }
                } else {
                    for (index_t i = g->fw_beg_pos[u]; i < g->fw_beg_pos[u + 1]; ++i) {
                        int v = g->fw_csr[i];
                        if (nodeToComponent.count(v) > 0 && nodeToComponent.at(v) == source_comp_id && distances[v] == -1) {
                            if (compare_and_swap(distances[v], -1, u_dist + 1)) {
                                lqueue.push_back(v);
                                scout_count++;
                            }
                        }
                    }
                }
            }
            lqueue.flush();
        }
        return scout_count;
    }

    // --- Bottom-Up Step ---
    template<typename GraphPtr>
    int64_t BUStep(GraphPtr g, Bitmap& front, Bitmap& next, std::vector<int>& distances, int level, const std::map<int, int>& nodeToComponent, int source_comp_id) {
        std::atomic<int64_t> awake_count(0);
        next.reset();

        #pragma omp parallel for schedule(dynamic, 1024)
        for (int u = 0; u < g->vert_count; ++u) {
            if (distances[u] == -1 && nodeToComponent.count(u) > 0 && nodeToComponent.at(u) == source_comp_id) {
                if constexpr (std::is_same_v<GraphPtr, graph*>) {
                    for (index_t i = g->bw_beg_pos[u]; i < g->bw_beg_pos[u + 1]; ++i) {
                        int v = g->bw_csr[i];
                        if (front.get_bit(v)) {
                            distances[u] = level + 1;
                            awake_count++;
                            next.set_bit(u);
                            break;
                        }
                    }
                } else {
                    for (index_t i = g->beg_pos[u]; i < g->beg_pos[u + 1]; ++i) {
                        int v = g->csr[i];
                        if (front.get_bit(v)) {
                            distances[u] = level + 1;
                            awake_count++;
                            next.set_bit(u);
                            break;
                        }
                    }
                }
            }
        }
        return awake_count;
    }

    // --- 队列与Bitmap转换 ---
    void QueueToBitmap(SlidingQueue<int>& queue, Bitmap& bm) {
        #pragma omp parallel for
        for (auto it = queue.begin(); it < queue.end(); ++it) {
            bm.set_bit(*it);
        }
    }

    void BitmapToQueue(Bitmap& bm, SlidingQueue<int>& queue, int vert_count) {
        #pragma omp parallel
        {
            QueueBuffer<int> lqueue(queue);
            #pragma omp for nowait
            for (int u = 0; u < vert_count; ++u) {
                if (bm.get_bit(u)) {
                    lqueue.push_back(u);
                }
            }
            lqueue.flush();
        }
        queue.slide_window();
    }

    // --- 主测试函数 ---
    void testParallelBFS(GIR& gir, int source, bool is_directed) {
        long long vert_count;
        long long edge_count;
        std::map<int, int> nodeToComponent;

        if (is_directed) {
            std::cout << "\n--- Running Parallel BFS Test on DIRECTED graph (pruned to SCC) ---" << std::endl;
            graph* g = gir.getDirectedGraph();
            if (!g) { std::cerr << "Test Error: Directed graph not loaded." << std::endl; return; }
            vert_count = g->vert_count;
            edge_count = g->edge_count;
            nodeToComponent = gir.getNodeToScc();
        } else {
            std::cout << "\n--- Running Parallel BFS Test on UNDIRECTED graph (pruned to WCC) ---" << std::endl;
            graph_undirected* g = gir.getUndirectedGraph();
            if (!g) { std::cerr << "Test Error: Undirected graph not loaded." << std::endl; return; }
            vert_count = g->vert_count;
            edge_count = g->edge_count;
            nodeToComponent = gir.getNodeToWcc();
        }

        if (source >= vert_count || nodeToComponent.empty() || nodeToComponent.find(source) == nodeToComponent.end() || nodeToComponent.at(source) == -1) {
            std::cerr << "Test Error: Source node " << source << " is invalid or not in any component." << std::endl;
            return;
        }
        int source_comp_id = nodeToComponent.at(source);
        std::cout << "Source: " << source << ", Component ID: " << source_comp_id << std::endl;

        auto start_time = std::chrono::high_resolution_clock::now();
        std::vector<int> distances(vert_count, -1);
        distances[source] = 0;

        SlidingQueue<int> queue(vert_count);
        queue.push_back(source);
        queue.slide_window();

        Bitmap front(vert_count);
        Bitmap next(vert_count);

        const int alpha = 15;
        const int beta = 18;

        int level = 0;
        long long edges_to_check = edge_count;
        long long scout_count = 0;
        if (is_directed) scout_count = gir.getDirectedGraph()->fw_beg_pos[source + 1] - gir.getDirectedGraph()->fw_beg_pos[source];
        else scout_count = gir.getUndirectedGraph()->beg_pos[source + 1] - gir.getUndirectedGraph()->beg_pos[source];

        while (!queue.empty()) {
            if (scout_count > (edges_to_check / alpha)) {
                QueueToBitmap(queue, front);
                queue.reset();
                long long awake_count;
                do {
                    std::cout << "  Level " << level << ": Bottom-Up Step" << std::endl;
                    if (is_directed) {
                        awake_count = BUStep(gir.getDirectedGraph(), front, next, distances, level, nodeToComponent, source_comp_id);
                    } else {
                        awake_count = BUStep(gir.getUndirectedGraph(), front, next, distances, level, nodeToComponent, source_comp_id);
                    }
                    front.swap(next);
                    level++;
                } while (awake_count > vert_count / beta);

                BitmapToQueue(front, queue, vert_count);
                scout_count = queue.size();

            } else {
                std::cout << "  Level " << level << ": Top-Down Step, frontier size = " << queue.size() << std::endl;
                edges_to_check -= scout_count;
                if (is_directed) {
                    scout_count = TDStep(gir.getDirectedGraph(), queue, distances, nodeToComponent, source_comp_id);
                } else {
                    scout_count = TDStep(gir.getUndirectedGraph(), queue, distances, nodeToComponent, source_comp_id);
                }
                queue.slide_window();
                level++;
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> bfs_time = end_time - start_time;

        long visited_count = 0;
        for(int d : distances) {
            if (d != -1) {
                visited_count++;
            }
        }

        std::cout << "BFS completed in: " << bfs_time.count() * 1000.0 << " ms" << std::endl;
        std::cout << "Nodes visited in component: " << visited_count << std::endl;
    }

//

} // namespace ParallelBFS

// ============================================================================
// EdgeUpdateTest -- differential validation of GIR's incremental maintenance
// (Scope-Prune-Update) against from-scratch recomputation, per
// paper/maintenance.tex's five-case classification. Unit-testing each case
// in isolation wouldn't catch cross-case interaction bugs, so this instead
// replays a random sequence of insertEdge/deleteEdge calls and, after each
// one, rebuilds WCC/SCC/BiCC from scratch on the *same* current graph view
// (base CSR + GIR's overlay, via the public getDirectedOutNeighbors/
// getDirectedInNeighbors/getUndirectedNeighborsPublic accessors) and diffs
// the two. See docs/code-review.md for why the previous update logic needed
// replacing instead of patching.
namespace EdgeUpdateTest {

    // Normalizes a partition (list of vertex groups) so two partitions of the
    // same underlying sets compare equal regardless of group/vertex order,
    // and regardless of empty "tombstoned" groups from retired SCC ids.
    std::vector<std::vector<int>> normalizePartition(std::vector<std::vector<int>> parts) {
        std::vector<std::vector<int>> result;
        for (auto& p : parts) {
            if (p.empty()) continue;
            std::sort(p.begin(), p.end());
            result.push_back(std::move(p));
        }
        std::sort(result.begin(), result.end());
        return result;
    }

    std::vector<std::vector<int>> groundTruthWcc(GIR& gir, int vert_count) {
        std::vector<int> comp(vert_count, -1);
        std::vector<std::vector<int>> groups;
        for (int s = 0; s < vert_count; ++s) {
            if (comp[s] != -1) continue;
            int gid = (int)groups.size();
            groups.push_back({});
            std::vector<int> stack = {s};
            comp[s] = gid;
            while (!stack.empty()) {
                int u = stack.back(); stack.pop_back();
                groups[gid].push_back(u);
                for (int w : gir.getDirectedOutNeighbors(u)) if (comp[w] == -1) { comp[w] = gid; stack.push_back(w); }
                for (int w : gir.getDirectedInNeighbors(u)) if (comp[w] == -1) { comp[w] = gid; stack.push_back(w); }
            }
        }
        return groups;
    }

    // Iterative Tarjan over the whole graph (ground truth SCC), mirroring
    // dynamic_repair::find_scc_scoped but over all vertices at once.
    std::vector<std::vector<int>> groundTruthScc(GIR& gir, int vert_count) {
        std::vector<int> index_of(vert_count, -1), low_of(vert_count, -1);
        std::vector<char> on_stack(vert_count, 0);
        std::vector<int> node_stack;
        std::vector<std::vector<int>> result;
        int next_index = 0;

        struct Frame { int v; std::vector<int> nbrs; size_t next_child; };

        for (int start = 0; start < vert_count; ++start) {
            if (index_of[start] != -1) continue;
            std::vector<Frame> call_stack;
            index_of[start] = low_of[start] = next_index++;
            node_stack.push_back(start);
            on_stack[start] = 1;
            call_stack.push_back({start, gir.getDirectedOutNeighbors(start), 0});

            while (!call_stack.empty()) {
                Frame& f = call_stack.back();
                if (f.next_child < f.nbrs.size()) {
                    int w = f.nbrs[f.next_child++];
                    if (index_of[w] == -1) {
                        index_of[w] = low_of[w] = next_index++;
                        node_stack.push_back(w);
                        on_stack[w] = 1;
                        call_stack.push_back({w, gir.getDirectedOutNeighbors(w), 0});
                    } else if (on_stack[w]) {
                        low_of[f.v] = std::min(low_of[f.v], index_of[w]);
                    }
                } else {
                    int v = f.v, v_low = low_of[v];
                    call_stack.pop_back();
                    if (!call_stack.empty()) {
                        int parent = call_stack.back().v;
                        low_of[parent] = std::min(low_of[parent], v_low);
                    }
                    if (v_low == index_of[v]) {
                        std::vector<int> comp;
                        while (true) {
                            int w = node_stack.back(); node_stack.pop_back();
                            on_stack[w] = 0;
                            comp.push_back(w);
                            if (w == v) break;
                        }
                        result.push_back(std::move(comp));
                    }
                }
            }
        }
        return result;
    }

    std::vector<std::vector<int>> groundTruthBiccForScc(GIR& gir, const std::vector<int>& verts) {
        // The undirected view is derived here from the *directed* edges
        // (out- plus in-neighbors, deduplicated), never from GIR's own
        // undirected overlay -- otherwise a bug in how GIR maintains that
        // overlay would be invisible to this check.
        std::unordered_set<int> vert_set(verts.begin(), verts.end());
        auto neighbor_fn = [&](int x) {
            std::set<int> nb;
            for (int w : gir.getDirectedOutNeighbors(x)) if (w != x && vert_set.count(w)) nb.insert(w);
            for (int w : gir.getDirectedInNeighbors(x)) if (w != x && vert_set.count(w)) nb.insert(w);
            return std::vector<int>(nb.begin(), nb.end());
        };
        auto result = dynamic_repair::find_bicc_scoped(verts, neighbor_fn);
        std::unordered_set<int> covered;
        for (auto& blk : result.blocks) for (int v : blk) covered.insert(v);
        for (int v : verts) if (!covered.count(v)) result.blocks.push_back({v});
        return result.blocks;
    }

    // Returns true (and prints nothing) if GIR's maintained state matches a
    // from-scratch rebuild; otherwise prints a diagnostic and returns false.
    bool checkConsistency(GIR& gir, int vert_count, const std::string& label) {
        gir.ensureBlocks();  // (blocks on demand: bring them up to date first)
        bool ok = true;

        // Undirected view must be exactly the symmetrization of the directed one.
        for (int a = 0; a < vert_count; ++a) {
            std::set<int> want, got;
            for (int w : gir.getDirectedOutNeighbors(a)) if (w != a) want.insert(w);
            for (int w : gir.getDirectedInNeighbors(a)) if (w != a) want.insert(w);
            for (int w : gir.getUndirectedNeighborsPublic(a)) got.insert(w);
            if (want != got) {
                std::cerr << "[MISMATCH][" << label << "] undirected view of vertex " << a
                          << " is not the symmetrized directed adjacency" << std::endl;
                ok = false;
                break;
            }
        }
        {
            std::string why;
            if (!gir.checkBlockForestInvariants(&why)) {
                std::cerr << "[MISMATCH][" << label << "] internal invariant: " << why << std::endl;
                ok = false;
            }
        }

        auto gt_wcc = normalizePartition(groundTruthWcc(gir, vert_count));
        auto gir_wcc = normalizePartition(gir.getWCCs());
        if (gt_wcc != gir_wcc) {
            std::cerr << "[MISMATCH][" << label << "] WCC partition differs from ground truth "
                      << "(gir=" << gir_wcc.size() << " groups, truth=" << gt_wcc.size() << " groups)" << std::endl;
            ok = false;
        }

        auto gt_scc = normalizePartition(groundTruthScc(gir, vert_count));
        std::vector<std::vector<int>> gir_scc_raw = gir.getSCCs();
        auto gir_scc = normalizePartition(gir_scc_raw);
        if (gt_scc != gir_scc) {
            std::cerr << "[MISMATCH][" << label << "] SCC partition differs from ground truth "
                      << "(gir=" << gir_scc.size() << " groups, truth=" << gt_scc.size() << " groups)" << std::endl;
            ok = false;
        }

        // Condensation-DAG edges: derive ground truth directly from real
        // vertex-level edges and compare against gir.getDAG() (keyed by the
        // *current* live scc ids, via nodeToScc). Catches the class of bug
        // where a DAG edge is dropped after deleting just one of several
        // vertex-level edges between the same pair of SCCs, even though
        // another still justifies it (see CHANGELOG.md).
        {
            std::map<int, int> nodeToSccGir0 = gir.getNodeToScc();
            std::set<std::pair<int, int>> gt_dag, gir_dag;
            for (int v = 0; v < vert_count; ++v) {
                int sv = nodeToSccGir0.at(v);
                for (int w : gir.getDirectedOutNeighbors(v)) {
                    int sw = nodeToSccGir0.at(w);
                    if (sv != sw) gt_dag.insert({sv, sw});
                }
            }
            auto dag = gir.getDAG();
            for (size_t s = 0; s < dag.size(); ++s)
                for (int t : dag[s]) gir_dag.insert({(int)s, t});
            if (gt_dag != gir_dag) {
                std::cerr << "[MISMATCH][" << label << "] condensation DAG edges differ from ground truth "
                          << "(gir=" << gir_dag.size() << " edges, truth=" << gt_dag.size() << " edges)"
                          << std::endl;
                ok = false;
            }
        }

        // BiCC is checked per (ground-truth) SCC: match each ground-truth SCC
        // to whichever live GIR scc id currently owns the same vertex set,
        // then compare block partitions.
        std::map<int, int> nodeToSccGir = gir.getNodeToScc();
        for (auto& gt_comp : gt_scc) {
            if (gt_comp.empty()) continue;
            int scc_id = nodeToSccGir.at(gt_comp[0]);
            auto gt_bicc = normalizePartition(groundTruthBiccForScc(gir, gt_comp));
            auto gir_bicc = normalizePartition(gir.getBCTComponents(scc_id));
            if (gt_bicc != gir_bicc) {
                std::cerr << "[MISMATCH][" << label << "] BiCC blocks differ from ground truth for SCC "
                          << scc_id << " (gir=" << gir_bicc.size() << " blocks, truth=" << gt_bicc.size()
                          << " blocks)" << std::endl;
                std::cerr << "  gt_comp verts:";
                for (int v : gt_comp) std::cerr << " " << v;
                std::cerr << "\n  gir blocks:";
                for (auto& b : gir_bicc) { std::cerr << " {"; for (int v : b) std::cerr << v << ","; std::cerr << "}"; }
                std::cerr << "\n  gt  blocks:";
                for (auto& b : gt_bicc) { std::cerr << " {"; for (int v : b) std::cerr << v << ","; std::cerr << "}"; }
                std::cerr << std::endl;
                ok = false;
            }
        }
        return ok;
    }

    // Runs `num_updates` random insert/delete operations (choosing whichever
    // is applicable for the sampled pair), checking full consistency after
    // every single update. Reports how many updates landed in each of the
    // five SPU cases (inferred from what actually changed) so a small test
    // run can confirm every case got exercised at least once.
    int runDifferentialTest(GIR& gir, int num_updates, unsigned seed) {
        int V = gir.vertexCount();
        if (V <= 1) { std::cout << "EdgeUpdateTest: graph too small, skipping." << std::endl; return 0; }
        std::srand(seed);

        std::cout << "\n=== EdgeUpdateTest: " << num_updates << " random updates on " << V << " vertices ===" << std::endl;

        int case_wcc_merge = 0, case_scc_merge_or_dag_edge = 0, case_bicc_change = 0,
            case_noop = 0, case_wcc_split = 0, case_scc_split = 0, mismatches = 0;

        if (!checkConsistency(gir, V, "initial")) {
            std::cerr << "EdgeUpdateTest: initial state already inconsistent with a from-scratch rebuild "
                         "-- aborting (see docs/code-review.md for known base-data caveats, e.g. "
                         "data/livejournal vs data/soc-LiveJournal)." << std::endl;
            return 1;
        }

        // Independent ledger of every edge that should currently exist, kept
        // entirely outside GIR's own overlay bookkeeping, so a mismatch can
        // be attributed to "the overlay itself lost track of an edge" vs
        // "the overlay is fine but a downstream layer misused it". Seeded
        // from the original base-graph edges (present since construction,
        // never routed through insertEdge/deleteEdge), then kept in sync
        // with every test-driven change below.
        std::set<std::pair<int,int>> ledger_directed, ledger_undirected;
        for (int a = 0; a < V; ++a) {
            for (int w : gir.getDirectedOutNeighbors(a)) ledger_directed.insert({a, w});
            for (int w : gir.getUndirectedNeighborsPublic(a)) ledger_undirected.insert({std::min(a, w), std::max(a, w)});
        }

        for (int iter = 0; iter < num_updates; ++iter) {
            int u = std::rand() % V;
            int v = std::rand() % V;
            // Half of the time, target an edge that exists (sampling random
            // pairs alone almost never hits one, so deletions -- and every
            // split case -- would go untested).
            if (std::rand() % 2) {
                auto out = gir.getDirectedOutNeighbors(u);
                if (!out.empty()) v = out[std::rand() % out.size()];
            }
            if (u == v) continue;

            bool present = gir.hasDirectedEdge(u, v);
            std::map<int, int> wcc_before = gir.getNodeToWcc();
            std::map<int, int> scc_before = gir.getNodeToScc();
            int wcc_groups_before = (int)gir.getWCCs().size();
            int scc_groups_before = (int)gir.getSCCs().size();

            std::string op;
            if (!present) {
                gir.insertEdge(u, v);
                op = "insert";
                ledger_directed.insert({u, v});
                ledger_undirected.insert({std::min(u,v), std::max(u,v)});
            } else {
                gir.deleteEdge(u, v);
                op = "delete";
                ledger_directed.erase({u, v});
                if (!ledger_directed.count({v, u}))
                    ledger_undirected.erase({std::min(u,v), std::max(u,v)});
            }

            int wcc_groups_after = (int)gir.getWCCs().size();
            int scc_groups_after = (int)gir.getSCCs().size();
            bool wcc_changed = wcc_before[u] != gir.getNodeToWcc()[u] || wcc_groups_after != wcc_groups_before;
            bool scc_changed = scc_before[u] != gir.getNodeToScc()[u] || scc_groups_after != scc_groups_before;

            if (op == "insert") {
                if (wcc_groups_after < wcc_groups_before) case_wcc_merge++;
                else if (scc_groups_after < scc_groups_before) case_scc_merge_or_dag_edge++;
                else if (scc_changed || wcc_changed) case_scc_merge_or_dag_edge++;
                else case_noop++;  // includes both true Case-4 no-ops and same-SCC BiCC-only merges
            } else {
                if (wcc_groups_after > wcc_groups_before) case_wcc_split++;
                else if (scc_groups_after > scc_groups_before) case_scc_split++;
                else case_bicc_change++;
            }

            bool ok = checkConsistency(gir, V, "iter " + std::to_string(iter) + " (" + op + " " +
                                             std::to_string(u) + "," + std::to_string(v) + ")");
            if (!ok) {
                mismatches++;
                // Cross-check the ledger against GIR's own view of edge
                // presence, to tell "overlay lost an edge" apart from
                // "overlay is fine but a downstream layer misused it".
                bool ledger_ok = true;
                for (int a = 0; a < V; ++a) {
                    for (int b = 0; b < V; ++b) {
                        if (a == b) continue;
                        bool led = ledger_directed.count({a, b}) > 0;
                        bool gir_has = gir.hasDirectedEdge(a, b);
                        if (led != gir_has) {
                            std::cerr << "  [LEDGER] directed(" << a << "," << b << "): ledger=" << led
                                      << " gir=" << gir_has << std::endl;
                            ledger_ok = false;
                        }
                    }
                }
                for (auto& p : ledger_undirected) {
                    if (!gir.hasUndirectedEdge(p.first, p.second)) {
                        std::cerr << "  [LEDGER] undirected(" << p.first << "," << p.second << "): ledger=1 gir=0" << std::endl;
                        ledger_ok = false;
                    }
                }
                for (int a = 0; a < V; ++a) {
                    for (int w : gir.getUndirectedNeighborsPublic(a)) {
                        if (a < w && !ledger_undirected.count({a, w})) {
                            std::cerr << "  [LEDGER] undirected(" << a << "," << w << "): ledger=0 gir=1" << std::endl;
                            ledger_ok = false;
                        }
                    }
                }
                if (ledger_ok) std::cerr << "  [LEDGER] edge presence matches ledger exactly -- bug is downstream of the overlay." << std::endl;
            }
        }

        std::cout << "EdgeUpdateTest done. Case tallies (heuristic, based on observed group-count changes): "
                  << "WCC-merge=" << case_wcc_merge
                  << " SCC-merge/DAG-edge=" << case_scc_merge_or_dag_edge
                  << " BiCC-only=" << case_bicc_change
                  << " no-op=" << case_noop
                  << " WCC-split=" << case_wcc_split
                  << " SCC-split=" << case_scc_split
                  << " | mismatches=" << mismatches << "/" << num_updates << std::endl;
        return mismatches;
    }

    // Deterministic, hand-verifiable check of Cases 5b/5c specifically
    // (random testing on a small dense graph rarely re-creates a genuine
    // bridge once enough edges accumulate, so it alone can't be trusted to
    // exercise splits -- see the case tallies from runDifferentialTest).
    // Meant to run on a *freshly built* GIR instance (paths.txt / the
    // example_graph's known structure), separate from any instance a random
    // test has already mutated.
    void runTargetedSplitTest(GIR& gir) {
        std::cout << "\n=== EdgeUpdateTest: targeted Case 5b/5c split check (example_graph) ===" << std::endl;
        int V = gir.vertexCount();
        bool all_ok = true;

        // --- Case 1 then Case 5c: bridge two WCCs, then split them back ---
        // Mirrors paper/maintenance.tex's own running example (there
        // 1-indexed vertices 4/13 across WCC1={1..12}/WCC2={13..18}; this
        // 0-indexed graph shifts that to vertices 3/12 across {0..11}/{12..17}).
        int wcc_before = (int)gir.getWCCs().size();
        gir.insertEdge(3, 12);
        int wcc_after_insert = (int)gir.getWCCs().size();
        bool case1_ok = (wcc_after_insert == wcc_before - 1) && checkConsistency(gir, V, "targeted case1");
        std::cout << "Case 1 (insert bridges two WCCs): " << wcc_before << " -> " << wcc_after_insert
                  << " WCCs " << (case1_ok ? "OK" : "FAILED") << std::endl;
        all_ok &= case1_ok;

        gir.deleteEdge(3, 12);
        int wcc_after_delete = (int)gir.getWCCs().size();
        bool case5c_ok = (wcc_after_delete == wcc_before) && checkConsistency(gir, V, "targeted case5c");
        std::cout << "Case 5c (delete that bridge splits the WCC back): " << wcc_after_insert << " -> "
                  << wcc_after_delete << " WCCs " << (case5c_ok ? "OK" : "FAILED") << std::endl;
        all_ok &= case5c_ok;

        // --- Case 5b: intra-SCC bridge deletion splits an SCC ---
        // The base graph already has both 2->9 and 9->2, merging {0,1,2}
        // and {9,10,11} (each their own 3-cycle) into one 6-vertex SCC via
        // that single bidirectional connection. Deleting just 9->2 removes
        // GIR's shared undirected {2,9} edge entirely (see gir.h's note on
        // deleteEdge treating directed+undirected as one logical edge), so
        // the only remaining path is 2->9 forward with no way back --
        // the SCC must split into {0,1,2} and {9,10,11}.
        int scc_before = (int)gir.getSCCs().size();
        std::map<int,int> nodeToScc_before = gir.getNodeToScc();
        bool same_scc_before = nodeToScc_before[0] == nodeToScc_before[9];
        gir.deleteEdge(9, 2);
        int scc_after = (int)gir.getSCCs().size();
        std::map<int,int> nodeToScc_after = gir.getNodeToScc();
        bool split_happened = nodeToScc_after[0] != nodeToScc_after[9] &&
                              nodeToScc_after[0] == nodeToScc_after[1] && nodeToScc_after[1] == nodeToScc_after[2] &&
                              nodeToScc_after[9] == nodeToScc_after[10] && nodeToScc_after[10] == nodeToScc_after[11];
        bool case5b_ok = same_scc_before && split_happened && (scc_after == scc_before + 1) &&
                         checkConsistency(gir, V, "targeted case5b");
        std::cout << "Case 5b (delete intra-SCC bridge 9->2 splits the SCC): " << scc_before << " -> " << scc_after
                  << " SCCs, {0,1,2}/{9,10,11} separated=" << split_happened << " "
                  << (case5b_ok ? "OK" : "FAILED") << std::endl;
        all_ok &= case5b_ok;

        std::cout << "Targeted split test: " << (all_ok ? "ALL OK" : "SOME FAILED") << std::endl;
    }

    // Filter-Merge-Update: random batches (with deliberate duplicates and
    // cancelling insert/delete pairs), full consistency check after every
    // batch. Returns the number of inconsistent batches.
    int runBatchDifferentialTest(GIR& gir, int num_batches, int batch_size, unsigned seed) {
        int V = gir.vertexCount();
        std::srand(seed);
        int bad = 0;
        for (int it = 0; it < num_batches; ++it) {
            std::vector<GIR::Update> batch;
            for (int k = 0; k < batch_size; ++k) {
                int u = std::rand() % V, v = std::rand() % V;
                int r = std::rand() % 4;
                if (r <= 1) {
                    auto out = gir.getDirectedOutNeighbors(u);
                    if (!out.empty()) { v = out[std::rand() % out.size()]; batch.push_back({u, v, false}); }
                } else if (r == 2) {
                    batch.push_back({u, v, true});
                } else if (!batch.empty()) {
                    GIR::Update prev = batch[std::rand() % batch.size()];
                    batch.push_back({prev.u, prev.v, !prev.insert});  // cancel / re-apply
                }
            }
            gir.applyBatch(batch);
            if (!checkConsistency(gir, V, "batch " + std::to_string(it))) ++bad;
        }
        std::cout << "Batch differential test: " << (num_batches - bad) << "/" << num_batches
                  << " batches consistent." << std::endl;
        return bad;
    }

    // Deferred maintenance: random batches go through ingestBatch; BFS (store
    // only) is checked while changes are pending, the layers after each sync.
    // Some rounds use applyBatch / single-edge updates instead, which must
    // sync on their own first.
    int runLazyDifferentialTest(GIR& gir, int num_batches, int batch_size, unsigned seed) {
        int V = gir.vertexCount();
        std::srand(seed);
        int bad = 0, checks = 0;
        auto plain_bfs = [&](int s) {
            std::vector<int> lv(V, -1);
            std::queue<int> q;
            lv[s] = 0; q.push(s);
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int w : gir.getDirectedOutNeighbors(u)) if (lv[w] < 0) { lv[w] = lv[u] + 1; q.push(w); }
            }
            return lv;
        };
        for (int it = 0; it < num_batches; ++it) {
            std::vector<GIR::Update> batch;
            for (int k = 0; k < batch_size; ++k) {
                int u = std::rand() % V, v = std::rand() % V;
                int r = std::rand() % 4;
                if (r <= 1) {
                    auto out = gir.getDirectedOutNeighbors(u);
                    if (!out.empty()) { v = out[std::rand() % out.size()]; batch.push_back({u, v, false}); }
                } else if (r == 2) {
                    batch.push_back({u, v, true});
                } else if (!batch.empty()) {
                    GIR::Update prev = batch[std::rand() % batch.size()];
                    batch.push_back({prev.u, prev.v, !prev.insert});
                }
            }
            int mode = std::rand() % 8;
            if (mode == 0) gir.applyBatch(batch);
            else if (mode == 1 && !batch.empty()) {
                gir.ingestBatch(std::vector<GIR::Update>(batch.begin(), batch.begin() + batch.size() / 2));
                auto& op = batch.back();
                op.insert ? gir.insertEdge(op.u, op.v) : gir.deleteEdge(op.u, op.v);
            } else gir.ingestBatch(batch);
            if (!gir.synced()) {
                int s = std::rand() % V;
                ++checks;
                if (gir.bfs(s) != plain_bfs(s)) {
                    std::cout << "[MISMATCH][lazy batch " << it << "] BFS on unsynced graph" << std::endl;
                    ++bad;
                }
            }
            // Reachability from the index (syncs) vs. a search on the store.
            for (int k = 0; k < 20; ++k) {
                int a = std::rand() % V, b = std::rand() % V;
                bool x = gir.reachableNoIndex(a, b), y = gir.reachable(a, b);
                ++checks;
                if (x != y) {
                    std::cout << "[MISMATCH][lazy batch " << it << "] reachable(" << a << "," << b << ") index " << y
                              << " store " << x << std::endl;
                    ++bad;
                }
            }
            if (it % 3 == 2 || it == num_batches - 1) {
                gir.sync();
                ++checks;
                if (!checkConsistency(gir, V, "lazy batch " + std::to_string(it))) ++bad;
            }
        }
        std::cout << "Lazy differential test: " << checks - bad << "/" << checks << " checks passed." << std::endl;
        return bad;
    }

} // namespace EdgeUpdateTest

// ============================================================================
// GirBfsTest -- correctness + timing harness for GIR::bfsTaskParallel
// (paper/task_parallel.tex's two-level parallel BFS: SCC pruning via the
// condensation DAG + BC-tree block-wavefront task parallelism + per-block
// local BFS). Validated the same way as the SPU rewrite: differentially,
// against a plain sequential BFS, rather than trusting the design by
// inspection alone.
namespace GirBfsTest {

    std::vector<int> sequentialBfs(GIR& gir, int source, int V) {
        std::vector<int> levels(V, -1);
        if (source < 0 || source >= V) return levels;
        levels[source] = 0;
        std::queue<int> q;
        q.push(source);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int w : gir.getDirectedOutNeighbors(u)) {
                if (levels[w] == -1) {
                    levels[w] = levels[u] + 1;
                    q.push(w);
                }
            }
        }
        return levels;
    }

    // Runs bfsTaskParallel from `num_sources` random sources and diffs every
    // vertex's level against a plain sequential BFS. Returns true iff all
    // sources matched exactly.
    bool runDifferentialBfsTest(GIR& gir, int num_sources, unsigned seed) {
        int V = gir.vertexCount();
        std::cout << "\n=== GirBfsTest: differential check of bfsTaskParallel vs sequential BFS "
                  << "(" << num_sources << " sources) ===" << std::endl;
        std::srand(seed);
        int mismatches = 0;
        for (int i = 0; i < num_sources; ++i) {
            int s = std::rand() % V;
            auto got = gir.bfsTaskParallel(s);
            auto want = sequentialBfs(gir, s, V);
            if (got != want) {
                mismatches++;
                int first_bad = -1;
                for (int v = 0; v < V; ++v) if (got[v] != want[v]) { first_bad = v; break; }
                std::cerr << "  [MISMATCH] source=" << s << " first differing vertex=" << first_bad
                          << " got=" << (first_bad >= 0 ? got[first_bad] : -999)
                          << " want=" << (first_bad >= 0 ? want[first_bad] : -999) << std::endl;
            }
        }
        std::cout << "GirBfsTest: " << (num_sources - mismatches) << "/" << num_sources
                  << " sources matched sequential BFS exactly." << std::endl;
        return mismatches == 0;
    }

    // Times bfsTaskParallel vs the existing data-parallel-only
    // ParallelBFS::testParallelBFS-style traversal (SCC/WCC-pruned top-down/
    // bottom-up BFS, no BC-tree block decomposition) across `num_sources`
    // sources, reporting total wall time for each so the benefit of task
    // parallelism (vs. data parallelism alone) is visible on real datasets.
    void benchmarkBfs(GIR& gir, int num_sources, unsigned seed, const std::string& dataset_label,
                       std::ostream& csv_out) {
        int V = gir.vertexCount();
        if (V == 0) return;
        std::srand(seed);
        std::vector<int> sources(num_sources);
        for (int i = 0; i < num_sources; ++i) sources[i] = std::rand() % V;

        auto t0 = std::chrono::high_resolution_clock::now();
        long long total_reached_task = 0;
        for (int s : sources) {
            auto levels = gir.bfsTaskParallel(s);
            for (int lv : levels) if (lv >= 0) total_reached_task++;
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double task_parallel_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        auto t2 = std::chrono::high_resolution_clock::now();
        long long total_reached_seq = 0;
        for (int s : sources) {
            auto levels = sequentialBfs(gir, s, V);
            for (int lv : levels) if (lv >= 0) total_reached_seq++;
        }
        auto t3 = std::chrono::high_resolution_clock::now();
        double sequential_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

        double speedup = sequential_ms > 0 ? sequential_ms / task_parallel_ms : 0.0;
        std::cout << "[BFS bench] " << dataset_label << ": " << num_sources << " sources, "
                  << "task-parallel=" << task_parallel_ms << " ms, sequential-baseline=" << sequential_ms
                  << " ms, speedup=" << speedup << "x "
                  << "(avg reached/source: task=" << (total_reached_task / num_sources)
                  << " seq=" << (total_reached_seq / num_sources) << ")" << std::endl;
        csv_out << dataset_label << "," << V << "," << num_sources << "," << task_parallel_ms << ","
                << sequential_ms << "," << speedup << "\n";
    }

} // namespace GirBfsTest
// ============================================================================
// QueryTest -- BFS / SSSP / BC engines against plain sequential references.
// ============================================================================
namespace QueryTest {

    std::vector<uint64_t> dijkstra(GIR& gir, int s) {
        int V = gir.vertexCount();
        std::vector<uint64_t> d(V, UINT64_MAX);
        using P = std::pair<uint64_t, int>;
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        d[s] = 0; pq.push({0, s});
        while (!pq.empty()) {
            auto [du, u] = pq.top(); pq.pop();
            if (du != d[u]) continue;
            for (int w : gir.getDirectedOutNeighbors(u)) {
                uint64_t nd = du + GIR::edgeWeight(u, w);
                if (nd < d[w]) { d[w] = nd; pq.push({nd, w}); }
            }
        }
        return d;
    }

    std::vector<double> brandes(GIR& gir, const std::vector<int>& sources) {
        int V = gir.vertexCount();
        std::vector<double> bc(V, 0.0);
        for (int s : sources) {
            std::vector<int> dist(V, -1), order;
            std::vector<double> sigma(V, 0), delta(V, 0);
            std::vector<std::vector<int>> pred(V);
            dist[s] = 0; sigma[s] = 1; order.push_back(s);
            for (size_t i = 0; i < order.size(); ++i) {
                int u = order[i];
                for (int w : gir.getDirectedOutNeighbors(u)) {
                    if (dist[w] < 0) { dist[w] = dist[u] + 1; order.push_back(w); }
                    if (dist[w] == dist[u] + 1) { sigma[w] += sigma[u]; pred[w].push_back(u); }
                }
            }
            for (size_t i = order.size(); i-- > 0;) {
                int w = order[i];
                for (int u : pred[w]) delta[u] += sigma[u] / sigma[w] * (1 + delta[w]);
                if (w != s) bc[w] += delta[w];
            }
        }
        return bc;
    }

    // Returns the number of mismatching checks.
    int run(GIR& gir, int num_sources, unsigned seed) {
        int V = gir.vertexCount(), bad = 0;
        std::srand(seed);
        std::vector<int> srcs;
        for (int i = 0; i < num_sources; ++i) srcs.push_back(std::rand() % V);
        for (int s : srcs) {
            auto want = GirBfsTest::sequentialBfs(gir, s, V);
            if (gir.bfs(s) != want) { ++bad; std::cerr << "[QUERY] bfs mismatch from " << s << std::endl; }
            if (gir.bfsTwoLevel(s) != want) { ++bad; std::cerr << "[QUERY] bfsTwoLevel mismatch from " << s << std::endl; }
            auto dw = dijkstra(gir, s);
            if (gir.sssp(s, 16) != dw) { ++bad; std::cerr << "[QUERY] sssp mismatch from " << s << std::endl; }
            if (gir.ssspTwoLevel(s, 16) != dw) { ++bad; std::cerr << "[QUERY] ssspTwoLevel mismatch from " << s << std::endl; }
        }
        auto got = gir.bc(srcs);
        auto want = brandes(gir, srcs);
        for (int v = 0; v < V; ++v)
            if (std::abs(got[v] - want[v]) > 1e-6 * (1 + std::abs(want[v]))) {
                ++bad; std::cerr << "[QUERY] bc mismatch at " << v << ": " << got[v] << " vs " << want[v] << std::endl; break;
            }
        return bad;
    }
    // Updates interleaved with queries from fixed sources: catches stale
    // query caches (the reachability cache must be invalidated by updates).
    int interleaved(GIR& gir, int rounds, unsigned seed) {
        int V = gir.vertexCount(), bad = 0;
        std::srand(seed);
        int s1 = std::rand() % V, s2 = std::rand() % V;
        for (int r = 0; r < rounds; ++r) {
            int u = std::rand() % V, v = std::rand() % V;
            if (std::rand() % 2) { auto o = gir.getDirectedOutNeighbors(u); if (!o.empty()) v = o[std::rand() % o.size()]; }
            if (u != v) { if (gir.hasDirectedEdge(u, v)) gir.deleteEdge(u, v); else gir.insertEdge(u, v); }
            for (int s : {s1, s2}) {
                auto want = GirBfsTest::sequentialBfs(gir, s, V);
                if (gir.bfs(s) != want) { ++bad; std::cerr << "[QUERY] stale bfs round " << r << std::endl; }
                if (gir.bfsTwoLevel(s) != want) { ++bad; std::cerr << "[QUERY] stale bfsTwoLevel round " << r << std::endl; }
                if (gir.ssspTwoLevel(s, 8) != dijkstra(gir, s)) { ++bad; std::cerr << "[QUERY] stale sssp round " << r << std::endl; }
            }
        }
        return bad;
    }
}  // namespace QueryTest

// Incremental BFS / SSSP against from-scratch results, under batches
// (applyBatch and ingestBatch) and single-edge updates.
namespace IncTest {
inline int run(GIR& gir, int rounds, unsigned seed) {
    const int V = gir.vertexCount();
    std::mt19937 rng(seed);
    int src = (int)(rng() % V), bad = 0;
    gir.trackIncremental(src, false);
    gir.trackIncremental(src, true);
    for (int r = 0; r < rounds && !bad; ++r) {
        std::vector<GIR::Update> batch;
        const int n = 1 + (int)(rng() % 40);
        for (int i = 0; i < n; ++i) {
            int u = (int)(rng() % V), v = (int)(rng() % V);
            if (u == v) continue;
            // Mostly deletions of existing edges (they exercise invalidation).
            if (rng() % 3) {
                auto out = gir.getDirectedOutNeighbors(u);
                if (!out.empty()) { batch.push_back({u, out[rng() % out.size()], false}); continue; }
            }
            batch.push_back({u, v, true});
        }
        const int how = (int)(rng() % 3);
        if (how == 0) gir.applyBatch(batch);
        else if (how == 1) gir.ingestBatch(batch);
        else for (auto& op : batch) op.insert ? gir.insertEdge(op.u, op.v) : gir.deleteEdge(op.u, op.v);
        gir.updateIncremental(false);
        gir.updateIncremental(true);
        auto lv = gir.bfs(src);
        auto d = gir.sssp(src);
        for (int v = 0; v < V; ++v) {
            uint64_t want_b = lv[v] < 0 ? UINT64_MAX : (uint64_t)lv[v];
            if (gir.incDist(false, v) != want_b || gir.incDist(true, v) != d[v]) {
                std::cerr << "[INC] round " << r << " (mode " << how << ") vertex " << v << ": bfs " << gir.incDist(false, v) << " vs " << want_b
                          << ", sssp " << gir.incDist(true, v) << " vs " << d[v] << std::endl;
                ++bad;
                break;
            }
        }
    }
    gir.sync();
    std::cout << "IncTest: " << (bad ? "MISMATCH" : "incremental BFS/SSSP match from-scratch") << std::endl;
    return bad;
}
}  // namespace IncTest
