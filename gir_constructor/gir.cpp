#include <stack>
#include "gir.h"
#include "scc_common.h"
#include "graph.h"
#include "graph_undirected.h"
#include "bicc_bfs.h"
#include "block_cut_tree.h" // 假设你有这个文件
#include <iostream>
#include <fstream>
#include <set>
#include <stack>
#include <set>
#include <map>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>



// -------------------- 构造/析构 --------------------
GIR::GIR(const std::string& path_directed,
         const std::string& path_undirected)
    : base_path_directed(path_directed),
      base_path_undirected(path_undirected) {}

GIR::~GIR() {
    if (g_directed) delete g_directed;
    if (g_undirected) delete g_undirected;
    if (bc_tree) delete bc_tree;
    if (ap) delete[] ap;
}

// -------------------- WCC --------------------
void GIR::buildWCC() {
    std::cout << "=== Build WCC ===" << std::endl;

    if (!g_directed) {
        std::string fw_bg_file = base_path_directed + "fw_begin.bin";
        std::string fw_aj_file = base_path_directed + "fw_adjacent.bin";
        std::string bw_bg_file = base_path_directed + "bw_begin.bin";
        std::string bw_aj_file = base_path_directed + "bw_adjacent.bin";
        g_directed = new graph(fw_bg_file.c_str(), fw_aj_file.c_str(),
                               bw_bg_file.c_str(), bw_aj_file.c_str());
    }

    // 转邻接表
    std::unordered_map<int,std::vector<int>> adj_list;
    for (index_t u=0; u<g_directed->vert_count; ++u) {
        for (index_t i=g_directed->fw_beg_pos[u]; i<g_directed->fw_beg_pos[u+1]; ++i)
            adj_list[u].push_back(g_directed->fw_csr[i]);
    }

    // 计算 WCC
    std::unordered_set<int> visited;
    wccs.clear();
    nodeToWcc.clear();

    for (auto& [node,_] : adj_list) {
        if (visited.count(node)) continue;
        std::vector<int> comp;
        std::stack<int> stk;
        stk.push(node);
        visited.insert(node);

        while (!stk.empty()) {
            int curr = stk.top(); stk.pop();
            comp.push_back(curr);
            for (int nei : adj_list[curr]) {
                if (!visited.count(nei)) {
                    visited.insert(nei);
                    stk.push(nei);
                }
            }
        }

        int wcc_id = wccs.size();
        for (int v : comp) nodeToWcc[v] = wcc_id;
        wccs.push_back(comp);
    }

    std::cout << "Found " << wccs.size() << " WCCs." << std::endl;
}

// -------------------- BCT --------------------
void GIR::buildBCT() {
    std::cout << "=== Build BCT ===" << std::endl;

    std::string undirected_bg_file = base_path_undirected + "undirected_begin.bin";
    std::string undirected_aj_file = base_path_undirected + "undirected_adjacent.bin";

    g_undirected = new graph_undirected(undirected_bg_file.c_str(),
                                        undirected_aj_file.c_str());
    vertex_t V = g_undirected->vert_count;
    ap = new bool[V];

    find_bicc_bfs(g_undirected, ap, bicc_components);

    bc_tree = new BlockCutTree();
    bc_tree->build(bicc_components, ap, g_undirected);
    bc_tree->printTree();
}

// -------------------- SCC + DAG --------------------
void GIR::buildSCC_DAG() {
    std::cout << "=== Build SCC + DAG ===" << std::endl;

    if (!g_directed) {
        std::string fw_bg_file = base_path_directed + "fw_begin.bin";
        std::string fw_aj_file = base_path_directed + "fw_adjacent.bin";
        std::string bw_bg_file = base_path_directed + "bw_begin.bin";
        std::string bw_aj_file = base_path_directed + "bw_adjacent.bin";
        g_directed = new graph(fw_bg_file.c_str(), fw_aj_file.c_str(),
                               bw_bg_file.c_str(), bw_aj_file.c_str());
    }

    double avg_time[15] = {0};
    int alpha=1, beta=1, gamma=1;
    double theta=0.5;
    index_t thread_count=1;

    scc_detection(g_directed, alpha, beta, gamma, theta, thread_count, avg_time, sccs, nodeToScc);

    // 转邻接表
    std::unordered_map<int,std::vector<int>> adj_list;
    for (index_t u=0; u<g_directed->vert_count; ++u)
        for (index_t i=g_directed->fw_beg_pos[u]; i<g_directed->fw_beg_pos[u+1]; ++i)
            adj_list[u].push_back(g_directed->fw_csr[i]);

    // DAG
    generate_and_save_dag(nodeToScc, adj_list, "dag.bin");
    buildDagAdjList(nodeToScc, adj_list);
}

// -------------------- 构建 DAG 邻接表 --------------------
void GIR::buildDagAdjList(const std::map<int,int>& nodeToScc,
                          const std::unordered_map<int,std::vector<int>>& adj_list) {
    scc_id_map.clear();
    int current_scc_id = 0;
    for (const auto& pair : nodeToScc)
        if (scc_id_map.find(pair.second) == scc_id_map.end())
            scc_id_map[pair.second] = current_scc_id++;

    int scc_count = current_scc_id;
    std::vector<std::set<int>> dag_adj_set(scc_count);

    for (auto& [u, neighbors] : adj_list) {
        if (!nodeToScc.count(u)) continue;
        int u_super = scc_id_map[nodeToScc.at(u)];
        for (int v : neighbors) {
            if (!nodeToScc.count(v)) continue;
            int v_super = scc_id_map[nodeToScc.at(v)];
            if (u_super != v_super) dag_adj_set[u_super].insert(v_super);
        }
    }

    dag_adj_list.assign(scc_count, {});
    for (int i=0; i<scc_count; ++i)
        dag_adj_list[i].assign(dag_adj_set[i].begin(), dag_adj_set[i].end());
}

void GIR::generate_and_save_dag(const std::map<int,int>& nodeToScc,
                                const std::unordered_map<int,std::vector<int>>& adj_list,
                                const std::string& output_filename) {
    // Map original SCC IDs to a contiguous range
    std::map<int,int> scc_id_map_local;
    int current_scc_id = 0;
    for (const auto& pair : nodeToScc) {
        if (scc_id_map_local.find(pair.second) == scc_id_map_local.end())
            scc_id_map_local[pair.second] = current_scc_id++;
    }
    int scc_count = current_scc_id;
    if (scc_count == 0) return;

    std::vector<std::set<int>> dag_adj_set(scc_count);
    for (const auto& pair : adj_list) {
        int u = pair.first;
        if (!nodeToScc.count(u)) continue;
        int u_super = scc_id_map_local[nodeToScc.at(u)];

        for (int v : pair.second) {
            if (!nodeToScc.count(v)) continue;
            int v_super = scc_id_map_local[nodeToScc.at(v)];
            if (u_super != v_super) dag_adj_set[u_super].insert(v_super);
        }
    }

    // 转 vector
    std::vector<std::vector<int>> dag_adj_list_local(scc_count);
    for (int i=0; i<scc_count; ++i)
        dag_adj_list_local[i].assign(dag_adj_set[i].begin(), dag_adj_set[i].end());

    // 保存到二进制文件
    std::ofstream outfile(output_filename, std::ios::out | std::ios::binary);
    if (!outfile) {
        std::cerr << "Error: cannot open file " << output_filename << std::endl;
        return;
    }

    outfile.write(reinterpret_cast<const char*>(&scc_count), sizeof(int));
    for (int i=0; i<scc_count; ++i) {
        int neighbor_count = dag_adj_list_local[i].size();
        outfile.write(reinterpret_cast<const char*>(&neighbor_count), sizeof(int));
        if (neighbor_count > 0)
            outfile.write(reinterpret_cast<const char*>(dag_adj_list_local[i].data()), neighbor_count*sizeof(int));
    }
    outfile.close();

    std::cout << "DAG saved to " << output_filename << std::endl;
}

// -------------------- Getter --------------------
const std::map<int,int>& GIR::getNodeToScc() const { return nodeToScc; }
const std::map<int,int>& GIR::getNodeToWcc() const { return nodeToWcc; }
const std::vector<std::vector<int>>& GIR::getSCCs() const { return sccs; }
const std::vector<std::vector<int>>& GIR::getDAG() const { return dag_adj_list; }
const std::vector<std::vector<int>>& GIR::getWCCs() const { return wccs; }
BlockCutTree* GIR::getBCT() const { return bc_tree; }
