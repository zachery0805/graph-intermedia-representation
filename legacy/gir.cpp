#include <iostream>
#include <vector>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>
#include <functional>
#include "graph.h" // 假设 graph.h 在当前目录下

using namespace std;

// --------------------- Tarjan Strongly Connected Components ---------------------
class TarjanSCC {
    int index = 0;
    unordered_map<int, int> indices, lowlink;
    stack<int> stk;
    unordered_set<int> onStack;
    vector<vector<int>> sccs;
    unordered_map<int, int> nodeToScc;
    unordered_map<int, vector<int>> graph;

public:
    TarjanSCC(const unordered_map<int, vector<int>>& g) : graph(g) {}

    void strongconnect(int v) {
        indices[v] = lowlink[v] = index++;
        stk.push(v);
        onStack.insert(v);

        if (graph.count(v)) { // 检查节点v是否存在于图中
            for (int w : graph[v]) {
                if (indices.find(w) == indices.end()) {
                    strongconnect(w);
                    lowlink[v] = min(lowlink[v], lowlink[w]);
                } else if (onStack.count(w)) {
                    lowlink[v] = min(lowlink[v], indices[w]);
                }
            }
        }

        if (lowlink[v] == indices[v]) {
            vector<int> scc;
            int w;
            do {
                w = stk.top(); stk.pop();
                onStack.erase(w);
                scc.push_back(w);
                nodeToScc[w] = sccs.size();
            } while (w != v);
            sccs.push_back(scc);
        }
    }

    void run() {
        // 遍历所有可能的节点，确保所有连通分量都被处理
        for (auto& p : graph) {
            int v = p.first;
            if (indices.find(v) == indices.end()) {
                strongconnect(v);
            }
        }
        // 对于那些没有出度的孤立点，它们可能不会被dfs遍历到，但是会存在于graph的key中
        // 确保所有节点都被尝试处理一次
        for (auto const& [node_id, _] : graph) {
            if (indices.find(node_id) == indices.end()) {
                strongconnect(node_id);
            }
        }
    }

    vector<vector<int>> getSCCs() const { return sccs; }
    unordered_map<int, int> getNodeToSccMap() const { return nodeToScc; }
};

// --------------------- WCC ---------------------
vector<vector<int>> computeWCC(const unordered_map<int, vector<int>>& graph) {
    unordered_map<int, vector<int>> undirected;
    // 构建无向图
    for (auto& [u, neighbors] : graph) {
        for (int v : neighbors) {
            undirected[u].push_back(v);
            undirected[v].push_back(u);
        }
    }
    // 确保所有在有向图中出现过的节点都存在于无向图中，即使它们没有出边或入边
    for (auto& [u, _] : graph) {
        if (undirected.find(u) == undirected.end()) {
            undirected[u] = {};
        }
    }


    unordered_set<int> visited;
    vector<vector<int>> wccs;

    for (auto& [node, _] : undirected) {
        if (visited.count(node)) continue;
        vector<int> comp;
        stack<int> stk;
        stk.push(node);
        visited.insert(node);
        while (!stk.empty()) {
            int curr = stk.top(); stk.pop();
            comp.push_back(curr);
            for (int nei : undirected[curr]) {
                if (!visited.count(nei)) {
                    visited.insert(nei);
                    stk.push(nei);
                }
            }
        }
        wccs.push_back(comp);
    }
    return wccs;
}

// --------------------- BCC ---------------------
class BCC {
public:
    unordered_map<int, vector<int>> g;
    unordered_map<int, int> disc, low;
    // vector<pair<int, int>> bridges; // 桥的计算在此代码中没有直接使用，但通常与BCC相关联
    vector<vector<int>> bccs;
    unordered_map<int, int> artPoints; // 记录割点，及其被发现的次数，理论上应该只被发现一次
    int time = 0;
    stack<pair<int, int>> edgeStack;

    BCC(const unordered_map<int, vector<int>>& subg) : g(subg) {}

    void dfs(int u, int parent) {
        disc[u] = low[u] = ++time;
        int children = 0;

        // 检查节点u是否存在于图中，防止访问空邻接列表
        if (g.count(u)) {
            for (int v : g[u]) {
                if (v == parent) continue; // 避免回溯到父节点

                if (!disc.count(v)) {
                    children++;
                    edgeStack.push({u, v});
                    dfs(v, u);
                    low[u] = min(low[u], low[v]);

                    // 检查割点和BCC的结束条件
                    if ((parent == -1 && children > 1) || (parent != -1 && low[v] >= disc[u])) {
                        artPoints[u]++; // 标记u为割点
                        vector<int> bcc;
                        // 从栈中弹出边直到遇到 (u,v)
                        pair<int, int> current_edge;
                        do {
                            current_edge = edgeStack.top(); edgeStack.pop();
                            bcc.push_back(current_edge.first);
                            bcc.push_back(current_edge.second);
                        } while (current_edge.first != u || current_edge.second != v); // 注意这里要完整匹配边

                        // 将u也加入到BCC中，因为它是割点
                        bcc.push_back(u); 
                        
                        sort(bcc.begin(), bcc.end());
                        bcc.erase(unique(bcc.begin(), bcc.end()), bcc.end());
                        bccs.push_back(bcc);
                    }
                } else if (disc[v] < disc[u]) { // 后向边
                    low[u] = min(low[u], disc[v]);
                    edgeStack.push({u, v});
                }
            }
        }
    }

    void run() {
        for (auto& [u, _] : g) {
            if (!disc.count(u)) {
                time = 0; // 重置时间戳，因为每个WCC内部的BCC是独立的
                // 清空栈，确保新WCC的BCC计算从干净的栈开始
                while(!edgeStack.empty()) edgeStack.pop(); 
                dfs(u, -1);

                // 处理DFS结束后栈中剩余的边，它们构成一个BCC
                if (!edgeStack.empty()) {
                    vector<int> bcc;
                    while (!edgeStack.empty()) {
                        bcc.push_back(edgeStack.top().first);
                        bcc.push_back(edgeStack.top().second);
                        edgeStack.pop();
                    }
                    sort(bcc.begin(), bcc.end());
                    bcc.erase(unique(bcc.begin(), bcc.end()), bcc.end());
                    bccs.push_back(bcc);
                }
            }
        }
        // 清理artPoints，因为同一个割点可能因为不同的子树而多次被标记
        // 实际上，只要被标记一次就是割点
        unordered_map<int, int> distinctArtPoints;
        for (auto& [ap, _] : artPoints) {
            distinctArtPoints[ap] = 1; 
        }
        artPoints = distinctArtPoints;
    }

    const vector<vector<int>>& getBCCs() const { return bccs; }
    const unordered_map<int, int>& getArticulationPoints() const { return artPoints; }
};

// 辅助函数：将CSR格式的图转换为邻接列表
unordered_map<int, vector<int>> convert_to_adj_list(const graph* g_csr) {
    unordered_map<int, vector<int>> adj_list;
    if (!g_csr) return adj_list;

    for (index_t u = 0; u < g_csr->vert_count; ++u) {
        // CSR 格式中，节点 u 的邻居从 fw_beg_pos[u] 到 fw_beg_pos[u+1]-1
        for (index_t i = g_csr->fw_beg_pos[u]; i < g_csr->fw_beg_pos[u+1]; ++i) {
            vertex_t v = g_csr->fw_csr[i];
            adj_list[u].push_back(v);
        }
        // 确保即使节点没有出边，它也存在于图中作为键
        if (g_csr->fw_beg_pos[u] == g_csr->fw_beg_pos[u+1] && adj_list.find(u) == adj_list.end()) {
             adj_list[u] = {}; // 添加空邻接列表
        }
    }
    return adj_list;
}

int main() {
    // 假设您的图文件存在，这里使用占位符
    
/*    
./scc_cpu /mnt/raid0_huge/yuede/data/friendster/fw_begin.bin /mnt/raid0_huge/yuede/data/friendster/fw_adjacent.bin /mnt/raid0_huge/yuede/data/friendster/bw_begin.bin /mnt/raid0_huge/yuede/data/friendster/bw_adjacent.bin 56 30 200 10 0.1 10
/home/cxl6029/projects/GIR/data/example_graph/fw_begin.bin /home/cxl6029/projects/GIR/data/example_graph/fw_adjacent.bin /home/cxl6029/projects/GIR/data/example_graph/bw_begin.bin /home/cxl6029/projects/GIR/data/example_graph/bw_adjacent.bin
*/
    const char* fw_beg_file = "/home/cxl6029/projects/GIR/data/example_graph/fw_begin.bin";
    const char* fw_csr_file = "/home/cxl6029/projects/GIR/data/example_graph/fw_adjacent.bin";
    const char* bw_beg_file = "/home/cxl6029/projects/GIR/data/example_graph/bw_begin.bin";
    const char* bw_csr_file = "/home/cxl6029/projects/GIR/data/example_graph/bw_adjacent.bin";

    double avg_time[15]; // graph_load 需要这个参数
    graph* my_graph_csr = new graph(); // 动态分配graph对象
    my_graph_csr = my_graph_csr->graph_load(fw_beg_file, fw_csr_file, bw_beg_file, bw_csr_file, avg_time);

    if (!my_graph_csr) {
        cerr << "Error: Failed to load graph." << endl;
        return 1;
    }

    // 将CSR格式的图转换为邻接列表格式
    unordered_map<int, vector<int>> graph_adj_list = convert_to_adj_list(my_graph_csr);

    // SCC
    TarjanSCC tarjan(graph_adj_list);
    tarjan.run();
    auto sccs = tarjan.getSCCs();
    auto nodeToScc = tarjan.getNodeToSccMap();

    cout << "SCCs:\n";
    for (size_t i = 0; i < sccs.size(); ++i) {
        cout << "SCC " << i << ": ";
        for (int node : sccs[i]) cout << node << " ";
        cout << "\n";
    }

    // SCC-DAG
    map<int, set<int>> sccGraph;
    for (auto& [u, neighbors] : graph_adj_list) {
        for (int v : neighbors) {
            // 确保节点在nodeToScc中存在，因为SCC算法可能不会处理孤立节点
            if (nodeToScc.count(u) && nodeToScc.count(v)) {
                int su = nodeToScc[u];
                int sv = nodeToScc[v];
                if (su != sv) sccGraph[su].insert(sv);
            }
        }
    }

    cout << "\nSCC-DAG:\n";
    for (auto& [u, neighbors] : sccGraph) {
        cout << "SuperNode " << u << " -> ";
        for (int v : neighbors) cout << v << " ";
        cout << "\n";
    }

    // WCC
    cout << "\nWCCs:\n";
    auto wccs = computeWCC(graph_adj_list);
    for (size_t i = 0; i < wccs.size(); ++i) {
        cout << "WCC " << i << ": ";
        for (int node : wccs[i]) cout << node << " ";
        cout << "\n";
    }

    // BCCs inside each SCC
    cout << "\nBCCs in each SCC:\n";
    for (size_t sccID = 0; sccID < sccs.size(); ++sccID) {
        unordered_map<int, vector<int>> subg;
        unordered_set<int> nodeSet(sccs[sccID].begin(), sccs[sccID].end());
        for (int u : sccs[sccID]) {
            // 检查graph_adj_list中是否存在u，即使没有出边也要考虑
            if (graph_adj_list.count(u)) {
                for (int v : graph_adj_list[u]) {
                    if (nodeSet.count(v)) {
                        // 在SCC内部构建无向子图
                        subg[u].push_back(v);
                        subg[v].push_back(u);
                    }
                }
            }
            // 确保SCC中的孤立节点也被加入到subg中，以便BCC算法能够处理
            if (subg.find(u) == subg.end()) {
                subg[u] = {}; 
            }
        }

        // 清理subg中的重复边（因为双向添加）
        for (auto& pair : subg) {
            sort(pair.second.begin(), pair.second.end());
            pair.second.erase(unique(pair.second.begin(), pair.second.end()), pair.second.end());
        }

        BCC bcc(subg);
        bcc.run();
        cout << "SCC " << sccID << ":\n";
        int bccID = 0;
        for (auto& comp : bcc.getBCCs()) {
            cout << "  BCC " << bccID++ << ": ";
            for (int node : comp) cout << node << " ";
            cout << "\n";
        }
        for (auto& [ap, _] : bcc.getArticulationPoints()) {
            cout << "  Articulation Point: " << ap << "\n";
        }
    }

    // 释放动态分配的graph对象
    delete my_graph_csr; 
    my_graph_csr = nullptr;

    return 0;
}

// 您可以移除或修改 construct_graph 函数，因为它现在不再需要原始的硬编码图
void construct_graph() {
    // 这个函数现在可以被移除，或者根据您的实际需求进行修改
    // 例如，如果您想在这里加载图并将其转换为邻接列表
}