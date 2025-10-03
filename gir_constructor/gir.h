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

class GIR {
public:
    GIR(const std::string& path_directed,
        const std::string& path_undirected);
    ~GIR();

    // 核心构建接口
    void buildWCC();
    void buildBCT();
    void buildSCC_DAG();

    // Getter
    const std::map<int,int>& getNodeToScc() const;
    const std::map<int,int>& getNodeToWcc() const;
    const std::vector<std::vector<int>>& getSCCs() const;
    const std::vector<std::vector<int>>& getDAG() const;
    const std::vector<std::vector<int>>& getWCCs() const;
    graph* getDirectedGraph() const;
    graph_undirected* getUndirectedGraph() const;
    BlockCutTree* getBCT() const;
    void updateEdgeDirected(int u, int v, bool insert);
    void updateEdgeUndirected(int u, int v, bool insert);

private:
    std::string base_path_directed;
    std::string base_path_undirected;

    // 图对象
    graph* g_directed = nullptr;
    graph_undirected* g_undirected = nullptr;

    // BCT
    BlockCutTree* bc_tree = nullptr;
    bool* ap = nullptr;
    std::vector<std::vector<int>> bicc_components;

    // SCC & DAG
    std::map<int,int> nodeToScc;
    std::vector<std::vector<int>> sccs;
    std::vector<std::vector<int>> dag_adj_list;
    std::map<int,int> scc_id_map;

    // WCC
    std::map<int,int> nodeToWcc;
    std::vector<std::vector<int>> wccs;

private:
    void buildDagAdjList(const std::map<int,int>& nodeToScc,
                         const std::unordered_map<int,std::vector<int>>& adj_list);
    //void generate_and_save_dag_optimized(const std::unordered_map<int,int>& nodeToScc,
    //                                     const std::string& output_filename);
    void generate_and_save_dag(
    const std::map<int, int>& nodeToScc,
    const std::unordered_map<int, std::vector<int>>& adj_list,
    const std::string& output_filename);

        // WCC
    void updateWCCInsert(int u, int v);
    void updateWCCDelete(int u, int v);

    // SCC
    void updateSCCInsert(int u, int v);
    void updateSCCDelete(int u, int v);
    bool isNewCycleFormed(int u, int v);

    // BCC
    void updateBCCInsertInternal(int u, int v);
    void updateBCCDeleteInternal(int u, int v);
    void updateBCCInsertCross(int u, int v);
    void updateBCCDeleteCross(int u, int v);
};

#endif // GIR_H