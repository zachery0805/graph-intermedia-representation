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
#include <map>

// ============================================================================
// 辅助类 (改编自您的参考代码)
// ============================================================================

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
        const std::map<int, int>* nodeToComponent_ptr;

        if (is_directed) {
            std::cout << "\n--- Running Parallel BFS Test on DIRECTED graph (pruned to SCC) ---" << std::endl;
            graph* g = gir.getDirectedGraph();
            if (!g) { std::cerr << "Test Error: Directed graph not loaded." << std::endl; return; }
            vert_count = g->vert_count;
            edge_count = g->edge_count;
            nodeToComponent_ptr = &gir.getNodeToScc();
        } else {
            std::cout << "\n--- Running Parallel BFS Test on UNDIRECTED graph (pruned to WCC) ---" << std::endl;
            graph_undirected* g = gir.getUndirectedGraph();
            if (!g) { std::cerr << "Test Error: Undirected graph not loaded." << std::endl; return; }
            vert_count = g->vert_count;
            edge_count = g->edge_count;
            nodeToComponent_ptr = &gir.getNodeToWcc();
        }
        const std::map<int, int>& nodeToComponent = *nodeToComponent_ptr;

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

        std::cout << "BFS completed in: " << bfs_time.count() * 1000 << " ms" << std::endl;
        std::cout << "Nodes visited in component: " << visited_count << std::endl;
    }

} // namespace ParallelBFS
