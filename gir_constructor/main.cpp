#include "gir.h"
#include <iostream>
#include "block_cut_tree.h"

int main() {
    std::string path_undirected = "/home/cxl6029/projects/GIR/gir_constructor/";
    std::string path_directed   = "/home/cxl6029/projects/GIR/data/example_graph/";

    GIR gir(path_directed, path_undirected);

    // 构建 WCC
    gir.buildWCC();

    // 构建 BCT
    gir.buildBCT();
    // 构建 SCC + DAG
    gir.buildSCC_DAG();

    // --- 打印 WCC 映射 ---
    std::cout << "--- Node to WCC Mapping ---" << std::endl;
    for (auto& [node, wcc_id] : gir.getNodeToWcc())
        std::cout << "Node " << node << " -> WCC " << wcc_id << std::endl;
    std::cout << "---------------------------" << std::endl;

    // --- 打印 SCC-DAG ---
    std::cout << "--- SCC-DAG ---" << std::endl;
    const auto& dag = gir.getDAG();
    for (size_t i=0; i<dag.size(); ++i) {
        std::cout << "Supernode " << i << " -> ";
        for (int v : dag[i]) std::cout << v << " ";
        std::cout << std::endl;
    }
    std::cout << "----------------" << std::endl;

    // --- 打印 BCT ---
    std::cout << "--- Block-Cut Tree ---" << std::endl;
    if (auto bc_tree = gir.getBCT())
        bc_tree->printTree();
    else
        std::cout << "BCT not built." << std::endl;
    std::cout << "--------------------" << std::endl;

    return 0;
}
