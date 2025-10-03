#include "gir.h"
#include <iostream>
#include <chrono>
#include "block_cut_tree.h"
#include "graph_tests.h"
#include "gir_test.cpp"


int main() {
    std::string path_undirected = "/home/cxl6029/projects/GIR/data/web-Google_sym/";
    std::string path_directed   = "/home/cxl6029/projects/GIR/data/web-Google_sym/";
    GIR gir(path_directed, path_undirected);


    auto total_start = std::chrono::high_resolution_clock::now();

    // 构建 WCC
    auto t1 = std::chrono::high_resolution_clock::now();
    gir.buildWCC();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> wcc_time = t2 - t1;

    // 构建 BCT
    auto t3 = std::chrono::high_resolution_clock::now();
    gir.buildBCT();
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> bct_time = t4 - t3;

    // 构建 SCC + DAG
    auto t5 = std::chrono::high_resolution_clock::now();
    gir.buildSCC_DAG();
    auto t6 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> scc_dag_time = t6 - t5;

    auto total_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_time = total_end - total_start;

    // 输出时间统计
    std::cout << "Time for WCC:     " << wcc_time.count() << " seconds\n";
    std::cout << "Time for BCT:     " << bct_time.count() << " seconds\n";
    std::cout << "Time for SCC+DAG: " << scc_dag_time.count() << " seconds\n";
    std::cout << "Total time:       " << total_time.count() << " seconds\n";



    std::cout << "\n=== Running Parallel BFS Tests ===" << std::endl;
    ParallelBFS::testParallelBFS(gir, 0, true);  // 有向图测试
    ParallelBFS::testParallelBFS(gir, 0, false); // 无向图测试
    return 0;
}