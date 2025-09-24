#include <stack>
#include "gir.h"
#include "scc_common.h"
#include "graph.h"
#include "graph_undirected.h"
#include "wcc_cpu/wcc_common.h"
#include "bicc_bfs.h"
#include "block_cut_tree.h"
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <omp.h>

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

    vertex_t* wcc_id = new vertex_t[g_directed->vert_count];
    const index_t thread_count = 1;
    wcc_detection(g_directed, thread_count, wcc_id);

    std::map<int, std::vector<int>> tempWccs;
    wccs.clear();
    nodeToWcc.clear();

    for (vertex_t i = 0; i < g_directed->vert_count; ++i) {
        nodeToWcc[i] = wcc_id[i];
        tempWccs[wcc_id[i]].push_back(i);
    }
    for (const auto& pair : tempWccs) {
        wccs.push_back(pair.second);
    }
    delete[] wcc_id;

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
}

// -------------------- SCC + DAG（全并行版） --------------------
void GIR::buildSCC_DAG() {
    std::cout << "=== Build SCC + DAG (Parallel) ===" << std::endl;

    if (!g_directed) {
        std::string fw_bg_file = base_path_directed + "fw_begin.bin";
        std::string fw_aj_file = base_path_directed + "fw_adjacent.bin";
        std::string bw_bg_file = base_path_directed + "bw_begin.bin";
        std::string bw_aj_file = base_path_directed + "bw_adjacent.bin";
        g_directed = new graph(fw_bg_file.c_str(), fw_aj_file.c_str(),
                               bw_bg_file.c_str(), bw_aj_file.c_str());
    }

    double avg_time[15] = {0};
    int alpha = 1, beta = 1, gamma = 1;
    double theta = 0.5;
    index_t thread_count = omp_get_max_threads(); // 使用所有可用线程

    auto t0 = std::chrono::high_resolution_clock::now();
    scc_detection(g_directed, alpha, beta, gamma, theta, thread_count, avg_time, sccs, nodeToScc);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_time1 = t1 - t0;
    std::cout << "Time to find SCC: " << total_time1.count() << " seconds\n";

    // -------------------- 转邻接表（并行） --------------------
    std::unordered_map<int,std::vector<int>> adj_list;
    #pragma omp parallel
    {
        std::unordered_map<int,std::vector<int>> local_adj;
        #pragma omp for nowait
        for (index_t u = 0; u < g_directed->vert_count; ++u) {
            for (index_t i = g_directed->fw_beg_pos[u]; i < g_directed->fw_beg_pos[u+1]; ++i) {
                local_adj[u].push_back(g_directed->fw_csr[i]);
            }
        }
        // 合并线程本地数据
        #pragma omp critical
        for (auto& pair : local_adj) {
            adj_list[pair.first] = std::move(pair.second);
        }
    }

    // -------------------- 并行生成 DAG --------------------
    generate_and_save_dag_parallel(nodeToScc, adj_list, "dag.bin");

    // -------------------- 构建 DAG 邻接表（并行） --------------------
    buildDagAdjList_parallel(nodeToScc, adj_list);
}

// -------------------- 并行构建 DAG 邻接表 --------------------
void GIR::buildDagAdjList_parallel(const std::map<int,int>& nodeToScc,
                                   const std::unordered_map<int,std::vector<int>>& adj_list) {
    scc_id_map.clear();
    int current_scc_id = 0;
    for (const auto& pair : nodeToScc)
        if (scc_id_map.find(pair.second) == scc_id_map.end())
            scc_id_map[pair.second] = current_scc_id++;

    int scc_count = current_scc_id;
    std::vector<std::set<int>> dag_adj_set(scc_count);

    #pragma omp parallel for
    for (int idx = 0; idx < (int)adj_list.size(); ++idx) {
        auto it = std::next(adj_list.begin(), idx);
        int u = it->first;
        const auto& neighbors = it->second;

        if (!nodeToScc.count(u)) continue;
        int u_super = scc_id_map.at(nodeToScc.at(u));

        for (int v : neighbors) {
            if (!nodeToScc.count(v)) continue;
            int v_super = scc_id_map.at(nodeToScc.at(v));
            if (u_super != v_super) {
                #pragma omp critical
                dag_adj_set[u_super].insert(v_super);
            }
        }
    }

    dag_adj_list.assign(scc_count, {});
    #pragma omp parallel for
    for (int i = 0; i < scc_count; ++i)
        dag_adj_list[i].assign(dag_adj_set[i].begin(), dag_adj_set[i].end());
}

// -------------------- 并行生成并保存 DAG --------------------
void GIR::generate_and_save_dag_parallel(const std::map<int,int>& nodeToScc,
                                         const std::unordered_map<int,std::vector<int>>& adj_list,
                                         const std::string& output_filename) {
    std::map<int,int> scc_id_map_local;
    int current_scc_id = 0;
    for (const auto& pair : nodeToScc)
        if (scc_id_map_local.find(pair.second) == scc_id_map_local.end())
            scc_id_map_local[pair.second] = current_scc_id++;

    int scc_count = current_scc_id;
    if (scc_count == 0) return;

    std::vector<std::set<int>> dag_adj_set(scc_count);

    #pragma omp parallel for
    for (int idx = 0; idx < (int)adj_list.size(); ++idx) {
        auto it = std::next(adj_list.begin(), idx);
        int u = it->first;
        const auto& neighbors = it->second;

        if (!nodeToScc.count(u)) continue;
        int u_super = scc_id_map_local.at(nodeToScc.at(u));

        for (int v : neighbors) {
            if (!nodeToScc.count(v)) continue;
            int v_super = scc_id_map_local.at(nodeToScc.at(v));
            if (u_super != v_super) {
                #pragma omp critical
                dag_adj_set[u_super].insert(v_super);
            }
        }
    }

    std::vector<std::vector<int>> dag_adj_list_local(scc_count);
    #pragma omp parallel for
    for (int i = 0; i < scc_count; ++i)
        dag_adj_list_local[i].assign(dag_adj_set[i].begin(), dag_adj_set[i].end());

    std::ofstream outfile(output_filename, std::ios::out | std::ios::binary);
    if (!outfile) {
        std::cerr << "Error: cannot open file " << output_filename << std::endl;
        return;
    }

    outfile.write(reinterpret_cast<const char*>(&scc_count), sizeof(int));
    for (int i = 0; i < scc_count; ++i) {
        int neighbor_count = dag_adj_list_local[i].size();
        outfile.write(reinterpret_cast<const char*>(&neighbor_count), sizeof(int));
        if (neighbor_count > 0)
            outfile.write(reinterpret_cast<const char*>(dag_adj_list_local[i].data()), neighbor_count*sizeof(int));
    }
    outfile.close();
    std::cout << "DAG saved to " << output_filename << std::endl;
}
// -------------------- 动态更新接口 --------------------
void GIR::updateEdgeDirected(int u, int v, bool insert) {
    if (insert) {
        updateSCCInsert(u, v);
        updateWCCInsert(u, v);
    } else {
        updateSCCDelete(u, v);
        updateWCCDelete(u, v);
    }
}

void GIR::updateEdgeUndirected(int u, int v, bool insert) {
    if (insert) {
        updateBCCInsertInternal(u, v);
        updateWCCInsert(u, v);
    } else {
        updateBCCDeleteInternal(u, v);
        updateWCCDelete(u, v);
    }
}

// -------------------- WCC 删除更新 --------------------
void GIR::updateWCCDelete(int u, int v) {
    int wcc_id = nodeToWcc[u];
    if (wcc_id != nodeToWcc[v]) return;

    std::unordered_set<int> visited;
    std::stack<int> st;
    st.push(u);
    visited.insert(u);

    for (size_t idx=0; idx < wccs[wcc_id].size(); ++idx) {
        int curr = wccs[wcc_id][idx];
        if (!visited.count(curr)) visited.insert(curr);
    }

    // 拆分 WCC
    std::vector<int> new_component;
    std::vector<int> old_component;

    for (int node : wccs[wcc_id]) {
        if (visited.count(node)) old_component.push_back(node);
        else {
            new_component.push_back(node);
            nodeToWcc[node] = wccs.size();
        }
    }

    if (!new_component.empty()) wccs.push_back(new_component);
    wccs[wcc_id] = old_component;
}

// -------------------- SCC 删除更新 --------------------
void GIR::updateSCCDelete(int u, int v) {
    if (nodeToScc[u] != nodeToScc[v]) return;

    int scc_id = nodeToScc[u];
    std::vector<int> to_remove;

    for (int node : sccs[scc_id]) {
        // 仅根据 GIR 内部结构判定，不依赖 graph
        if (sccs[scc_id].size() <= 1) to_remove.push_back(node);
    }

    for (int node : to_remove) {
        sccs[scc_id].erase(std::remove(sccs[scc_id].begin(), sccs[scc_id].end(), node), sccs[scc_id].end());
        nodeToScc.erase(node);
    }
}

// -------------------- SCC 新增检测环 --------------------
bool GIR::isNewCycleFormed(int u, int v) {
    std::unordered_set<int> visited;
    std::stack<int> st;
    st.push(v);
    visited.insert(v);

    while (!st.empty()) {
        int curr = st.top(); st.pop();
        for (int nei : sccs[nodeToScc[curr]]) {
            if (nei == u) return true;
            if (!visited.count(nei)) {
                visited.insert(nei);
                st.push(nei);
            }
        }
    }
    return false;
}

// -------------------- BCC 删除内部更新 --------------------
void GIR::updateBCCDeleteInternal(int u, int v) {
    std::vector<int> queue = {u, v};
    std::unordered_set<int> visited;

    while (!queue.empty()) {
        int curr = queue.back(); queue.pop_back();
        if (visited.count(curr)) continue;
        visited.insert(curr);

        // 删除度 <= 1 的节点
        for (auto& bicc : bicc_components) {
            bicc.erase(std::remove(bicc.begin(), bicc.end(), curr), bicc.end());
        }

        for (auto& bicc : bicc_components) {
            for (int nei : bicc) {
                if (!visited.count(nei)) queue.push_back(nei);
            }
        }
    }
}

// -------------------- WCC 新增更新 --------------------
void GIR::updateWCCInsert(int u, int v) {
    int wcc_u = nodeToWcc.count(u) ? nodeToWcc[u] : -1;
    int wcc_v = nodeToWcc.count(v) ? nodeToWcc[v] : -1;

    if (wcc_u == -1 && wcc_v == -1) {
        // 两个节点都不在任何 WCC 中，新建一个
        std::vector<int> new_component = {u, v};
        wccs.push_back(new_component);
        int new_id = wccs.size() - 1;
        nodeToWcc[u] = nodeToWcc[v] = new_id;
    } else if (wcc_u != -1 && wcc_v == -1) {
        // v 加入 u 的 WCC
        wccs[wcc_u].push_back(v);
        nodeToWcc[v] = wcc_u;
    } else if (wcc_u == -1 && wcc_v != -1) {
        // u 加入 v 的 WCC
        wccs[wcc_v].push_back(u);
        nodeToWcc[u] = wcc_v;
    } else if (wcc_u != wcc_v) {
        // 合并两个 WCC
        auto& comp_u = wccs[wcc_u];
        auto& comp_v = wccs[wcc_v];

        comp_u.insert(comp_u.end(), comp_v.begin(), comp_v.end());
        for (int node : comp_v) nodeToWcc[node] = wcc_u;

        // 标记 comp_v 为已合并（可清空）
        comp_v.clear();
    }
}

// -------------------- SCC 新增更新 --------------------
void GIR::updateSCCInsert(int u, int v) {
    int scc_u = nodeToScc.count(u) ? nodeToScc[u] : -1;
    int scc_v = nodeToScc.count(v) ? nodeToScc[v] : -1;

    if (scc_u == -1 && scc_v == -1) {
        // 两个节点都不在任何 SCC 中，新建 SCC
        std::vector<int> new_scc = {u, v};
        sccs.push_back(new_scc);
        int new_id = sccs.size() - 1;
        nodeToScc[u] = nodeToScc[v] = new_id;
    } else if (scc_u != -1 && scc_v == -1) {
        sccs[scc_u].push_back(v);
        nodeToScc[v] = scc_u;
    } else if (scc_u == -1 && scc_v != -1) {
        sccs[scc_v].push_back(u);
        nodeToScc[u] = scc_v;
    } else if (scc_u != scc_v) {
        // 判断是否会形成环
        bool cycle = isNewCycleFormed(u, v);
        if (cycle) {
            // 合并两个 SCC
            auto& comp_u = sccs[scc_u];
            auto& comp_v = sccs[scc_v];

            comp_u.insert(comp_u.end(), comp_v.begin(), comp_v.end());
            for (int node : comp_v) nodeToScc[node] = scc_u;

            // 标记 comp_v 为已合并
            comp_v.clear();
        }
    }
}

// -------------------- BCC 内部新增更新 --------------------
void GIR::updateBCCInsertInternal(int u, int v) {
    std::vector<int> touched_nodes = {u, v};
    std::unordered_set<int> visited;

    for (int node : touched_nodes) {
        if (visited.count(node)) continue;
        visited.insert(node);

        bool added = false;
        for (auto& bicc : bicc_components) {
            if (std::find(bicc.begin(), bicc.end(), node) != bicc.end()) {
                added = true;
                break;
            }
        }

        if (!added) {
            bicc_components.push_back({node});
        }
    }
}

// -------------------- Getter --------------------
const std::map<int,int>& GIR::getNodeToScc() const { return nodeToScc; }
const std::map<int,int>& GIR::getNodeToWcc() const { return nodeToWcc; }
const std::vector<std::vector<int>>& GIR::getSCCs() const { return sccs; }
const std::vector<std::vector<int>>& GIR::getDAG() const { return dag_adj_list; }
const std::vector<std::vector<int>>& GIR::getWCCs() const { return wccs; }
BlockCutTree* GIR::getBCT() const { return bc_tree; }
