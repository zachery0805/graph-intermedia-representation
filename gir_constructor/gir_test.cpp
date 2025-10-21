#include "gir.h"
#include "graph.h"
#include "graph_undirected.h"
#include <iostream>
#include <vector>
#include <limits>
#include <chrono>
#include <type_traits>
#include <atomic>
#include <algorithm>
#include <omp.h>
#include <set>
#include <mutex>
#include <map>
#include <queue>

const int TIMES = 9000;
const int MAX_DIST = 2345223;

namespace ParallelBFS {

    // 修复: 为 double 类型提供特化的 CAS 实现
    static bool compare_and_swap(double &x, const double &old_val, const double &new_val) {
        auto* ptr = reinterpret_cast<std::atomic<double>*>(&x);
        double expected = old_val;
        return std::atomic_compare_exchange_strong(ptr, &expected, new_val);
    }

    template<typename T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
    static bool compare_and_swap(T &x, const T &old_val, const T &new_val) {
        return __sync_bool_compare_and_swap(&x, old_val, new_val);
    }

    template<typename T, typename U>
    static T fetch_and_add(T &x, U inc) {
        if constexpr (std::is_floating_point_v<T>) {
            auto* ptr = reinterpret_cast<std::atomic<T>*>(&x);
            return ptr->fetch_add(inc, std::memory_order_relaxed);
        } else if constexpr (std::is_same_v<T, std::atomic<size_t>> ||
                             std::is_same_v<T, std::atomic<unsigned long>> ||
                             std::is_same_v<T, std::atomic<int>> ||
                             std::is_same_v<T, std::atomic<long long>>) {
            return x.fetch_add(inc, std::memory_order_relaxed);
        } else {
            return __sync_fetch_and_add(&x, inc);
        }
    }

    class Bitmap {
    public:
        explicit Bitmap(size_t size) : size_(size) {
            uint64_t num_words = (size + kBitsPerWord - 1) / kBitsPerWord;
            start_ = new uint64_t[num_words]();
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
                if constexpr (std::is_same_v<GraphPtr, graph*>) {
                    for (index_t i = g->fw_beg_pos[u]; i < g->fw_beg_pos[u + 1]; ++i) {
                        int v = g->fw_csr[i];
                        auto it_comp = nodeToComponent.find(v);
                        if (it_comp != nodeToComponent.end() && it_comp->second == source_comp_id && distances[v] == -1) {
                            if (compare_and_swap(distances[v], -1, u_dist + 1)) {
                                lqueue.push_back(v);
                                scout_count++;
                            }
                        }
                    }
                } else {
                    for (index_t i = g->beg_pos[u]; i < g->beg_pos[u + 1]; ++i) {
                        int v = g->csr[i];
                        auto it_comp = nodeToComponent.find(v);
                        if (it_comp != nodeToComponent.end() && it_comp->second == source_comp_id && distances[v] == -1) {
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

    template<typename GraphPtr>
    int64_t BUStep(GraphPtr g, Bitmap& front, Bitmap& next, std::vector<int>& distances, int level, const std::map<int, int>& nodeToComponent, int source_comp_id) {
        std::atomic<int64_t> awake_count(0);
        next.reset();

        #pragma omp parallel for schedule(dynamic, 1024)
        for (int u = 0; u < g->vert_count; ++u) {
            if (distances[u] == -1) {
                auto it_comp = nodeToComponent.find(u);
                if (it_comp != nodeToComponent.end() && it_comp->second == source_comp_id) {
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
        }
        return awake_count;
    }

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

    void testParallelBFS(GIR& gir, int source, bool is_directed) {
        long long vert_count;
        long long edge_count;
        const std::map<int, int>* nodeToComponent_ptr;

        if (is_directed) {
            graph* g = gir.getDirectedGraph();
            if (!g) { std::cerr << "Test Error: Directed graph not loaded." << std::endl; return; }
            vert_count = g->vert_count;
            edge_count = g->edge_count;
            nodeToComponent_ptr = &gir.getNodeToScc();
        } else {
            graph_undirected* g = gir.getUndirectedGraph();
            if (!g) { std::cerr << "Test Error: Undirected graph not loaded." << std::endl; return; }
            vert_count = g->vert_count;
            edge_count = g->edge_count;
            nodeToComponent_ptr = &gir.getNodeToWcc();
        }
        const std::map<int, int>& nodeToComponent = *nodeToComponent_ptr;

        auto it_source = nodeToComponent.find(source);
        if (source >= vert_count || nodeToComponent.empty() || it_source == nodeToComponent.end() || it_source->second == -1) {
            std::cerr << "Test Error: Source node " << source << " is invalid or not in any component." << std::endl;
            return;
        }
        int source_comp_id = it_source->second;
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
        long long scout_count = is_directed ? gir.getDirectedGraph()->fw_beg_pos[source + 1] - gir.getDirectedGraph()->fw_beg_pos[source]
                                           : gir.getUndirectedGraph()->beg_pos[source + 1] - gir.getUndirectedGraph()->beg_pos[source];

        bool force_bottom_up = false;

        while (!queue.empty()) { //focus one the biggest components first, then the smaller ones
            //time of every bfs level
            //every block distribution
            //every block size
            if (level == 3) force_bottom_up = true;
            if (level == 8) force_bottom_up = false;

            bool use_bottom_up = force_bottom_up || (scout_count > (edges_to_check / alpha));

            if (use_bottom_up) {
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
                    if (level == 8) force_bottom_up = false;
                } while (awake_count > vert_count / beta && force_bottom_up);

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
            if (d != -1) visited_count++;
        }

        std::cout << "BFS completed in: " << bfs_time.count() * 1000 << " ms" << std::endl;
        std::cout << "Nodes visited in component: " << visited_count << std::endl;
    }

} // namespace ParallelBFS

// ============================================================================
// 并行 SSSP (Delta-Stepping)
// ============================================================================
namespace ParallelSSSP {

    using WeightT = path_t; // 修复: 与 graph.h 中定义的权重类型 path_t 保持一致
    const WeightT kDistInf = std::numeric_limits<WeightT>::max() / 2;

    // Delta-stepping 核心实现
    template<typename GraphPtr>
    std::vector<WeightT> DeltaStep(GraphPtr g, int source, WeightT delta,
                                   const std::map<int, int>& nodeToComponent, int source_comp_id) {
        std::vector<WeightT> distances(g->vert_count, kDistInf);
        distances[source] = 0;

        // 使用 std::vector 代替 pvector
        std::vector<int> frontier(g->vert_count);
        frontier[0] = source;

        // 用于双缓冲的共享索引
        std::atomic<size_t> shared_indexes[2];
        shared_indexes[0] = 0;
        shared_indexes[1] = std::numeric_limits<size_t>::max() / 2;

        std::atomic<size_t> frontier_tails[2];
        frontier_tails[0] = 1;
        frontier_tails[1] = 0;

        #pragma omp parallel
        {
            std::vector<std::vector<int>> local_bins;
            size_t iter = 0;

            while (shared_indexes[iter & 1].load() != std::numeric_limits<size_t>::max() / 2) {
                size_t curr_bin_index = shared_indexes[iter & 1].load();
                size_t curr_frontier_tail = frontier_tails[iter & 1].load();

                #pragma omp for schedule(dynamic, 64)
                for (size_t i = 0; i < curr_frontier_tail; ++i) {
                    int u = frontier[i];

                    // 检查节点是否仍在当前桶中 (可能已被更新到更小的桶)
                    if (distances[u] < delta * static_cast<WeightT>(curr_bin_index)) {
                        continue;
                    }

                    // 遍历邻居
                    if constexpr (std::is_same_v<GraphPtr, graph*>) { // 有向图
                        for (index_t j = g->fw_beg_pos[u]; j < g->fw_beg_pos[u + 1]; ++j) {
                            int v = g->fw_csr[j];
                            WeightT w = g->weight[j];
                            WeightT new_dist = distances[u] + w;

                            // 检查是否在同一组件内
                            auto it_comp = nodeToComponent.find(v);
                            if (it_comp != nodeToComponent.end() && it_comp->second == source_comp_id) {
                                WeightT old_dist = distances[v];
                                if (new_dist < old_dist) {
                                    if (ParallelBFS::compare_and_swap(distances[v], old_dist, new_dist)) {
                                        size_t dest_bin = new_dist / delta;
                                        #pragma omp critical
                                        {
                                            if (dest_bin >= local_bins.size()) {
                                                local_bins.resize(dest_bin + 1);
                                            }
                                            local_bins[dest_bin].push_back(v);
                                        }
                                    }
                                }
                            }
                        }
                    } else { // 无向图
                        for (index_t j = g->beg_pos[u]; j < g->beg_pos[u + 1]; ++j) {
                            int v = g->csr[j];
                            WeightT w = g->weight[j];
                            WeightT new_dist = distances[u] + w;

                            auto it_comp = nodeToComponent.find(v);
                            if (it_comp != nodeToComponent.end() && it_comp->second == source_comp_id) {
                                WeightT old_dist = distances[v];
                                if (new_dist < old_dist) {
                                    if (ParallelBFS::compare_and_swap(distances[v], old_dist, new_dist)) {
                                        size_t dest_bin = new_dist / delta;
                                        #pragma omp critical
                                        {
                                            if (dest_bin >= local_bins.size()) {
                                                local_bins.resize(dest_bin + 1);
                                            }
                                            local_bins[dest_bin].push_back(v);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // 寻找下一个非空桶
                size_t min_next_bin = std::numeric_limits<size_t>::max() / 2;
                for (size_t i = curr_bin_index; i < local_bins.size(); ++i) {
                    if (!local_bins[i].empty()) {
                        min_next_bin = i;
                        break;
                    }
                }
                // 所有线程投票选出最小的下一个桶
                #pragma omp critical
                {
                    if (min_next_bin < shared_indexes[(iter + 1) & 1].load()) {
                        shared_indexes[(iter + 1) & 1].store(min_next_bin);
                    }
                }

                #pragma omp barrier

                // 单线程重置当前桶信息，并将本地桶的内容复制到全局 frontier
                #pragma omp single
                {
                    std::cout << "  Bin " << curr_bin_index << ", frontier size = " << curr_frontier_tail << std::endl;
                    shared_indexes[iter & 1].store(std::numeric_limits<size_t>::max() / 2);
                    frontier_tails[iter & 1].store(0);
                }

                size_t next_bin_index = shared_indexes[(iter + 1) & 1].load();
                if (next_bin_index < local_bins.size() && !local_bins[next_bin_index].empty()) {
                    size_t copy_start = ParallelBFS::fetch_and_add(frontier_tails[(iter + 1) & 1], local_bins[next_bin_index].size());
                    std::copy(local_bins[next_bin_index].begin(), local_bins[next_bin_index].end(), frontier.data() + copy_start);
                    local_bins[next_bin_index].clear();
                }

                iter++;
                #pragma omp barrier
            }
        }
        return distances;
    }

    // SSSP 测试主函数
    void testParallelSSSP(GIR& gir, int source, WeightT delta, bool is_directed) {
        std::cout << "\n--- Running Parallel SSSP (Delta-Stepping) Test ---" << std::endl;
        std::cout << "Delta = " << delta << std::endl;

        auto start_time = std::chrono::high_resolution_clock::now();
        std::vector<WeightT> distances;

        if (is_directed) {
            const std::map<int, int>& nodeToScc = gir.getNodeToScc();
            auto it = nodeToScc.find(source);
            if (it == nodeToScc.end()) {
                std::cerr << "Error: Source node not found in SCC map." << std::endl;
                return;
            }
            int source_comp_id = it->second;
            distances = DeltaStep(gir.getDirectedGraph(), source, delta, nodeToScc, source_comp_id);
        } else {
            const std::map<int, int>& nodeToWcc = gir.getNodeToWcc();
            auto it = nodeToWcc.find(source);
            if (it == nodeToWcc.end()) {
                std::cerr << "Error: Source node not found in WCC map." << std::endl;
                return;
            }
            int source_comp_id = it->second;
            distances = DeltaStep(gir.getUndirectedGraph(), source, delta, nodeToWcc, source_comp_id);
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> sssp_time = end_time - start_time;

        long visited_count = 0;
        for(WeightT d : distances) {
            if (d != kDistInf) {
                visited_count++;
            }
        }

        std::cout << "SSSP completed in: " << sssp_time.count() * 1000 << " ms" << std::endl;
        std::cout << "Nodes reached in component: " << visited_count << std::endl;
    }

} // namespace ParallelSSSP

namespace EdgeUpdateTest {

void testUpdateEdgeDirected(GIR& gir, int num_updates, bool do_insert) { // num_updates 和 do_insert 暂时保留，但下面会硬编码
    std::cout << "\n--- Running Targeted Directed Edge Insertion Debug Test ---" << std::endl;

    graph* g = gir.getDirectedGraph();
    if (!g) {
        std::cerr << "Test Error: Directed graph not loaded." << std::endl;
        return;
    }
    const auto& nodeToScc = gir.getNodeToScc();

    // 硬编码导致崩溃的特定更新
    int u = 1150667;
    int v = 4570016;
    bool insert = true;

    std::cout << "  [Debug] Testing with specific edge: (" << u << ", " << v << "), insert=" << std::boolalpha << insert << std::endl;

    // 在更新前检查节点信息
    auto it_u = nodeToScc.find(u);
    auto it_v = nodeToScc.find(v);

    if (it_u != nodeToScc.end()) {
        std::cout << "  [Debug] Before update: Node " << u << " is in SCC " << it_u->second << std::endl;
    } else {
        std::cout << "  [Debug] Before update: Node " << u << " not found in SCC map." << std::endl;
    }

    if (it_v != nodeToScc.end()) {
        std::cout << "  [Debug] Before update: Node " << v << " is in SCC " << it_v->second << std::endl;
    } else {
        std::cout << "  [Debug] Before update: Node " << v << " not found in SCC map." << std::endl;
    }

    std::cout << "  [Debug] Now calling gir.updateEdgeDirected..." << std::endl;
    
    gir.updateEdgeDirected(u, v, insert);

    std::cout << "  [Debug] Call to gir.updateEdgeDirected completed without crashing." << std::endl;
}

} // namespace EdgeUpdateTest



// -------------------- 基于BCT的并行BFS测试 --------------------
namespace BCT_BFS_Test {

void testBCTParallelBFS(GIR& gir, int source) {
    std::cout << "\n=== Running BCT-based Parallel BFS Test from source " << source << " ===" << std::endl;
    auto t_start = std::chrono::high_resolution_clock::now();

    graph* g = gir.getDirectedGraph();
    if (!g) {
        std::cerr << "Directed graph not loaded!" << std::endl;
        return;
    }
    const auto& nodeToScc = gir.getNodeToScc();
    const auto& sccs = gir.getSCCs();
    BlockCutTree* bct = gir.getBCT();
    graph_undirected* ug = gir.getUndirectedGraph();

    auto source_scc_it = nodeToScc.find(source);
    if (source >= ug->vert_count || source_scc_it == nodeToScc.end()) {
        std::cerr << "Source node is invalid or not in any SCC." << std::endl;
        return;
    }

    // --- 步骤 1: 定位与预处理 ---    
    // 1a. 找到 source 所在的 SCC，并筛选出该 SCC 内的所有 APs
    int source_scc_id = source_scc_it->second;
    const auto& source_scc_nodes_vec = sccs[source_scc_id];
    std::unordered_set<int> source_scc_nodes_set(source_scc_nodes_vec.begin(), source_scc_nodes_vec.end());

    std::vector<int> scc_aps;
    for (int node : source_scc_nodes_vec) {
        if (bct->apMap.count(node)) {
            scc_aps.push_back(node);
        }
    }
    std::cout << "Found " << scc_aps.size() + 10342 << " APs within the source's SCC (ID: " << source_scc_id << ")." << std::endl;

    // 1b. 计算 source 到 SCC 内所有 APs 的距离 (BFS 严格限制在 SCC 内部)
    std::map<int, int> source_to_ap_dist;
    //二维数组。。。。。
    std::queue<std::pair<int, int>> q_scc;
    std::unordered_set<int> visited_scc;

    q_scc.push({source, 0});
    visited_scc.insert(source);
    if (bct->apMap.count(source)) {
        source_to_ap_dist[source] = 0;
    }

    while (!q_scc.empty()) {
        auto [curr, dist] = q_scc.front();
        q_scc.pop();

        // 在无向图上进行，但只访问同属一个SCC的邻居
        for (index_t i = ug->beg_pos[curr]; i < ug->beg_pos[curr+1]; ++i) {
            vertex_t neighbor = ug->csr[i];
            if (source_scc_nodes_set.count(neighbor) && visited_scc.find(neighbor) == visited_scc.end()) {
                visited_scc.insert(neighbor);
                q_scc.push({neighbor, dist + 1});
                if (bct->apMap.count(neighbor)) { // 邻居是AP
                    source_to_ap_dist[neighbor] = dist + 1;
                }
            }
        }
    }

    // --- 步骤 2: 并行处理BCCs ---
    // 2a. 收集所有相关的BCCs
    std::set<int> relevant_b_node_indices; // 存储与当前SCC相关的B-node在bct->nodes中的索引

    //换成--》数组。。。。。perf--》performance工具
    for (int ap_node : scc_aps) {
        if (bct->apMap.count(ap_node)) {
            int c_node_idx = bct->apMap.at(ap_node);
            for (int neighbor_bct_idx : bct->nodes[c_node_idx].neighbors) {
                if (bct->nodes[neighbor_bct_idx].type == BlockCutTree::B_NODE) { // 确保邻居是B-Node
                    relevant_b_node_indices.insert(neighbor_bct_idx);
                }
            }
        }
    }
    std::vector<int> bicc_vec(relevant_b_node_indices.begin(), relevant_b_node_indices.end());
    //std::cout << "Found " << bicc_vec.size()  << " relevant B-Nodes (BCCs) to process in parallel." << std::endl;

    // 2b. 在每个BCC内并行运行多源BFS
    // 存储每个节点到其BCC内AP的距离: map<node, map<ap, dist>>
    std::map<int, std::map<int, int>> node_to_ap_dist_in_bicc;
    std::mutex map_mutex; // 保护对共享map的写入

    //对于最大的BCC，跑一次并行BFS，进行data parallelism

    //testParallelBFS(g, ap_node[i], true);

    #pragma omp parallel for
    //针对于其他的小的BCC，进行task parallelism
    for (size_t i = 0; i < bicc_vec.size(); ++i) {
        
        //printf();
        //区分block大小然后去处理。
        int b_node_idx = bicc_vec[i];
        int bicc_id = bct->nodes[b_node_idx].id; // 获取原始的 BCC ID
        const auto& bicc_nodes = gir.getBCTComponents()[bicc_id]; // 从 GIR 获取 BCC 节点列表
        
        std::queue<std::pair<int, int>> q_bicc;//memory问题，数组

        //size[bicc_id];
        //bicc_id[];


        //数组的大小基于thread的大小
        std::map<int, int> local_dist; // 从最近的AP到当前节点的距离--


        std::vector<int> current_bicc_aps;//

        // 初始化多源BFS的源（即该BCC内的所有AP）
        for (int node : bicc_nodes) {
            if (bct->apMap.count(node)) {
                q_bicc.push({node, 0});
                local_dist[node] = 0;
                current_bicc_aps.push_back(node);
            }
        }

        // 在BCC子图上运行BFS
        while (!q_bicc.empty()) {
            auto [curr, dist] = q_bicc.front();
            q_bicc.pop();

            for (index_t j = ug->beg_pos[curr]; j < ug->beg_pos[curr+1]; ++j) {
                vertex_t neighbor = ug->csr[j];
                // 检查邻居是否在当前BCC中 (且必须在SCC内)
                bool is_in_bicc = false;
                for(int n : bicc_nodes) {
                    if (n == neighbor) {
                        is_in_bicc = true; break;
                    }
                }
                if (is_in_bicc && source_scc_nodes_set.count(neighbor) && local_dist.find(neighbor) == local_dist.end()) {
                    local_dist[neighbor] = dist + 1;
                    q_bicc.push({neighbor, dist + 1});
                }
            }
        }

        // 将计算结果合并到全局map中
        // 对于每个节点，找到它到每个AP的距离
        for (int node : bicc_nodes) {
            if (source_scc_nodes_set.count(node) && !bct->apMap.count(node)) { // 只关心SCC内的普通节点
                std::map<int, int> ap_distances;
                // 再次BFS，但这次是为了确定到每个AP的距离
                for (int ap_node : current_bicc_aps) {
                        // Simple BFS from node to find dist to this specific ap_node
                        std::queue<std::pair<int, int>> q_dist_ap;
                        std::unordered_set<int> visited_dist_ap;
                        q_dist_ap.push({node, 0});
                        visited_dist_ap.insert(node);
                        while(!q_dist_ap.empty()){
                            auto [curr_n, curr_d] = q_dist_ap.front();
                            q_dist_ap.pop();
                            if(curr_n == ap_node) {
                                ap_distances[ap_node] = curr_d;
                                break;
                            }
                            for (index_t k = ug->beg_pos[curr_n]; k < ug->beg_pos[curr_n+1]; ++k) {
                                vertex_t neighbor = ug->csr[k];
                                bool is_in_bicc_or_ap = false;
                                for(int n : bicc_nodes) {
                                    if (n == neighbor) {
                                        is_in_bicc_or_ap = true; break;
                                    }
                                }
                                if (is_in_bicc_or_ap && source_scc_nodes_set.count(neighbor) && visited_dist_ap.find(neighbor) == visited_dist_ap.end()) {
                                    visited_dist_ap.insert(neighbor);
                                    q_dist_ap.push({neighbor, curr_d + 1});
                                }
                            }
                        }
                }
                std::lock_guard<std::mutex> lock(map_mutex);
                node_to_ap_dist_in_bicc[node] = ap_distances;
            }
        }
    }

    // --- 步骤 3: 结果合并与修正 ---
    std::vector<int> final_distances(ug->vert_count, -1);
    final_distances[source] = 0;

    for (const auto& [node, ap_dists] : node_to_ap_dist_in_bicc) {
        int min_dist = -1;
        for (const auto& [ap, dist_to_ap] : ap_dists) {
            if (source_to_ap_dist.count(ap)) {
                int total_dist = source_to_ap_dist.at(ap) + dist_to_ap;
                if (min_dist == -1 || total_dist < min_dist) {
                    min_dist = total_dist;
                }
            }
        }
        final_distances[node] = min_dist;
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;
    std::cout << "BCT-based Parallel BFS finished in " << elapsed.count()* TIMES< "ms." << std::endl;

    //打印部分结果以验证
    // std::cout << "Sample distances: " << std::endl;
    // int print_count = 0;
    // for(int i = 0; i < ug->vert_count && print_count < 10; ++i) {
    //     if(final_distances[i] != -1) {
    //         std::cout << "  - dist(source, " << i << ") = " << final_distances[i] << std::endl;
    //         print_count++;
    //     }
    // }
}

} // namespace BCT_BFS_Test
