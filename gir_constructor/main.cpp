#include "gir.h"
#include <iostream>
#include <chrono>
#include "block_cut_tree.h"
#include "graph_tests.h"
#include "gir_test.cpp"
#include <omp.h>
#include <cstdio>

void PrintOMPThreadInfo() {
    int max_threads = omp_get_max_threads();  // 理论最大可用线程
    int dynamic_threads = 0;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        {
            dynamic_threads = omp_get_num_threads();  // 实际运行线程数
            printf("[OpenMP] Available threads: %d\n", max_threads);
            printf("[OpenMP] Running with %d threads\n", dynamic_threads);
        }
    }
#else
    printf("[OpenMP] Not enabled in this build (compile with -fopenmp)\n");
#endif
}


int main() {
    std::string path_undirected = "/home/cxl6029/projects/GIR/data/livejournal/";
    std::string path_directed   = "/home/cxl6029/projects/GIR/data/livejournal/";
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
    PrintOMPThreadInfo();
    ParallelBFS::testParallelBFS(gir, 9, true);  // 有向图测试

    //BCT_BFS_Test::testBCTParallelBFS(gir, 9);
    // // 动态更新测试
    // int num_updates = 1000;
    // EdgeUpdateTest::testUpdateEdgeDirected(gir, num_updates, true); // 测试边插入
    // EdgeUpdateTest::testUpdateEdgeDirected(gir, num_updates, false); // 测试边删除



    // SSSP (Delta-Stepping) 测试
    /*
    std::cout << "\n=== Running Parallel SSSP (Delta-Stepping) Tests ===" << std::endl;
    long long delta = 1; // Delta-Stepping 的 delta 参数，对于权重为1的图，1是较优选择
    ParallelSSSP::testParallelSSSP(gir, 9, delta, true); // 有向图测试
    */
    // BCT-based 并行 BFS 测试





    return 0;
}