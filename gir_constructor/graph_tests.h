#ifndef GRAPH_TESTS_H
#define GRAPH_TESTS_H
#include "gir.h"
#include "graph.h"
#include <iostream>
#include <vector>
#include <queue>
#include <chrono>
#include <limits>
#include <stack>   // 修复: 包含 std::stack 的定义
#include <cstdlib> // 包含 rand() 的定义

// -------------------- BFS Test --------------------
void testBFS(GIR& gir, int source) {
    graph* g = gir.getDirectedGraph();
    const auto& nodeToScc = gir.getNodeToScc();

    if (!g) {
        std::cerr << "BFS Test: Graph is not loaded." << std::endl;
        return;
    }
    if (source >= g->vert_count || nodeToScc.find(source) == nodeToScc.end() || nodeToScc.at(source) == -1) {
        std::cerr << "BFS Test: Source node " << source << " not in any SCC." << std::endl;
        return;
    }
    std::cout << "\n--- Running BFS Test from source " << source << " (pruned to its SCC) ---" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<int> distances(g->vert_count, -1);
    std::queue<int> q;
    long visited_count = 0;

    int source_scc_id = nodeToScc.at(source);
    distances[source] = 0;
    q.push(source);
    visited_count++;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (index_t i = g->fw_beg_pos[u]; i < g->fw_beg_pos[u + 1]; ++i) {
            vertex_t v = g->fw_csr[i];
            // 剪枝：只在同一个 SCC 内进行 BFS
            if (distances[v] == -1 && nodeToScc.count(v) > 0 && nodeToScc.at(v) == source_scc_id) {
                distances[v] = distances[u] + 1;
                q.push(v);
                visited_count++;
            }
        } // for
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> bfs_time = end - start;

    std::cout << "BFS completed in: " << bfs_time.count() << " seconds" << std::endl;
    std::cout << "Nodes visited: " << visited_count << "/" << g->vert_count << std::endl;
}

// -------------------- SSSP Test (Dijkstra) --------------------
// Note: Assumes uniform edge weight of 1 since weights are not loaded.
void testSSSP(GIR& gir, int source) {
    graph* g = gir.getDirectedGraph();
    const auto& nodeToScc = gir.getNodeToScc();

    if (!g) {
        std::cerr << "SSSP Test: Graph is not loaded." << std::endl;
        return;
    }
    if (source >= g->vert_count || nodeToScc.find(source) == nodeToScc.end() || nodeToScc.at(source) == -1) {
        std::cerr << "SSSP Test: Source node " << source << " not in any SCC." << std::endl;
        return;
    }
    std::cout << "\n--- Running SSSP (Dijkstra) Test from source " << source << " (pruned to its SCC) ---" << std::endl;
    std::cout << "(Assuming uniform edge weight of 1)" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    using WN = std::pair<float, int>; // Weight, Node
    std::priority_queue<WN, std::vector<WN>, std::greater<WN>> pq;
    std::vector<float> distances(g->vert_count, std::numeric_limits<float>::infinity());

    int source_scc_id = nodeToScc.at(source);
    distances[source] = 0;
    pq.push({0.0f, source});

    while (!pq.empty()) {
        float d = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (d > distances[u]) {
            continue;
        }

        for (index_t i = g->fw_beg_pos[u]; i < g->fw_beg_pos[u + 1]; ++i) {
            vertex_t v = g->fw_csr[i];
            // 剪枝：只在同一个 SCC 内进行 SSSP
            if (nodeToScc.count(v) > 0 && nodeToScc.at(v) == source_scc_id) {
                float weight = 1.0f; // 假设权重为 1
                if (distances[u] + weight < distances[v]) {
                    distances[v] = distances[u] + weight;
                    pq.push({distances[v], v});
                }
            }
        } // for
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sssp_time = end - start;
    std::cout << "SSSP completed in: " << sssp_time.count() << " seconds" << std::endl;
}

// -------------------- Betweenness Centrality (BC) Test --------------------
// Note: This is an approximation using a subset of nodes.
void testBC(GIR& gir, int num_samples) {
    graph* g = gir.getDirectedGraph();
    const auto& nodeToScc = gir.getNodeToScc();

    if (!g) {
        std::cerr << "BC Test: Graph is not loaded." << std::endl;
        return;
    }
    std::cout << "\n--- Running Approximate BC Test with " << num_samples << " samples (pruned to each source's SCC) ---" << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<double> bc(g->vert_count, 0.0);

    for (int i = 0; i < num_samples; ++i) {
        int source = rand() % g->vert_count;

        if (nodeToScc.find(source) == nodeToScc.end() || nodeToScc.at(source) == -1) {
            continue; // 跳过不在任何 SCC 中的采样点
        }
        int source_scc_id = nodeToScc.at(source);

        std::vector<int> dist(g->vert_count, -1);
        std::vector<double> sigma(g->vert_count, 0.0);
        std::vector<std::vector<int>> pred(g->vert_count);
        std::stack<int> s;
        std::queue<int> q;

        dist[source] = 0;
        sigma[source] = 1.0;
        q.push(source);

        while(!q.empty()){
            int v = q.front(); q.pop();
            s.push(v);
            for (index_t j = g->fw_beg_pos[v]; j < g->fw_beg_pos[v + 1]; ++j) {
                int w = g->fw_csr[j];
                // 剪枝：只在同一个 SCC 内进行
                if (nodeToScc.count(w) > 0 && nodeToScc.at(w) == source_scc_id) {
                    if (dist[w] < 0) {
                         dist[w] = dist[v] + 1;
                         q.push(w);
                    }
                    if (dist[w] == dist[v] + 1) {
                        sigma[w] += sigma[v];
                        pred[w].push_back(v);
                    }
                }
            } // for
        }

        std::vector<double> delta(g->vert_count, 0.0);
        while(!s.empty()){
            int w = s.top(); s.pop();
            if (sigma[w] > 0) { // 避免除以零
                for(int v : pred[w]){
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                }
            }
            if(w != source) bc[w] += delta[w];
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> bc_time = end - start;
    std::cout << "Approximate BC completed in: " << bc_time.count() << " seconds" << std::endl;
}

#endif // GRAPH_TESTS_H