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

    // 构建 SCC + DAG (must run before BCT: BiCC is computed per-SCC, see
    // paper/task_parallel.tex and gir.h's class comment)
    auto t5 = std::chrono::high_resolution_clock::now();
    gir.buildSCC_DAG();
    auto t6 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> scc_dag_time = t6 - t5;

    // 构建 BCT (per SCC)
    auto t3 = std::chrono::high_resolution_clock::now();
    gir.buildBCT();
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> bct_time = t4 - t3;

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

    // 动态更新测试: differential check against from-scratch rebuild (see
    // docs/code-review.md -- the previous update functions were mostly
    // no-ops and this is the only thing that would have caught it).
    //
    // checkConsistency() rebuilds WCC/SCC/BiCC from scratch on the *entire*
    // graph after every single simulated update -- appropriate for a small
    // sanity-check graph (e.g. data/example_graph/, 18 vertices) but far too
    // slow to run unattended against a full production dataset like
    // livejournal/twitter. Point path_directed/path_undirected at
    // data/example_graph/ above to exercise this.
    if (gir.vertexCount() <= 10000) {
        EdgeUpdateTest::runDifferentialTest(gir, 3000, 42);

        // Targeted, hand-verifiable check of Cases 5b/5c on a fresh instance
        // (random testing alone rarely re-creates a genuine bridge once the
        // graph is well-connected, so it can't be trusted to exercise splits).
        GIR gir_split_test(path_directed, path_undirected);
        gir_split_test.buildWCC();
        gir_split_test.buildSCC_DAG();
        gir_split_test.buildBCT();
        EdgeUpdateTest::runTargetedSplitTest(gir_split_test);
    } else {
        std::cout << "\nSkipping EdgeUpdateTest (graph has " << gir.vertexCount()
                  << " vertices; the differential check rebuilds from scratch after every "
                     "update and is only intended for small sanity-check graphs, e.g. "
                     "data/example_graph/)." << std::endl;
    }



    // SSSP (Delta-Stepping) 测试
    /*
    std::cout << "\n=== Running Parallel SSSP (Delta-Stepping) Tests ===" << std::endl;
    long long delta = 1; // Delta-Stepping 的 delta 参数，对于权重为1的图，1是较优选择
    ParallelSSSP::testParallelSSSP(gir, 9, delta, true); // 有向图测试
    */
    // BCT-based 并行 BFS 测试





    return 0;
}