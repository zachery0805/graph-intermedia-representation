#include <unistd.h>
#include <stdlib.h>
#include "wtime.h"
#include "graph.h"
//#include "frontier_queue.h"
#include "wcc_common.h"
#include <omp.h>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <cstdio>
#include <cstring>



int main(int args, char **argv)
{
	if(args != 6)
    {
	    std::cout<<"Usage: ./bicc <beg_file> <csr_file> <thread_count> <run_times> <output_file>\n";
        exit(-1);
    }
	
	const char *beg_file = argv[1];
	const char *csr_file = argv[2];
	const vertex_t thread_count = atoi(argv[3]);
    const vertex_t run_times = atoi(argv[4]);
    const char *output_file = argv[5];
    
//    char cmd[50];
//    sprintf(cmd, "export OMP_NUM_THREADS=%d", thread_count);
//    printf("%s", cmd);
//    system(cmd);
    
//    double * avg_time = new double[15]; 
    ///step 1: load the graph
    graph_undirected *g = new graph_undirected(beg_file, csr_file);
    
    for(vertex_t i = 0; i < run_times; ++i)
    {
        printf("Test #%d\n", i);
        double time_begin = wtime();
        find_cc(g,
                thread_count,
                output_file);

    //    find_bicc_bfs(g, 
    //            thread_count,
    //            output_file); 
        double time_end = wtime();
        printf("\n");
    }

//    printf("args = %d\n", args);
//	if(args != 9)
//    {
//	    std::cout<<"Usage: ./wcc_cpu <fw_beg_file> <fw_csr_file> <bw_beg_file> <bw_csr_file> <thread_count> <run_times> <result_file_wcc> <result_file_time>\n";
//        exit(-1);
//    }
//	
//	const char *fw_beg_file = argv[1];
//	const char *fw_csr_file = argv[2];
//	const char *bw_beg_file = argv[3];
//	const char *bw_csr_file = argv[4];
//	const vertex_t thread_count=atoi(argv[5]);
//    const vertex_t run_times = atoi(argv[6]);
//    char *result_file_wcc = argv[7];
//    char *result_file_time = argv[8];
//
//	printf("Thread = %d, run_times = %d\n", thread_count, run_times);
//    
//    double * avg_time = new double[15]; 
//
//    for(vertex_t i = 0; i < 15; ++i)
//        avg_time[i] = 0.0;
//
//    ///step 1: load the graph
//    graph *g = new graph(fw_beg_file,
//                    fw_csr_file,
//                    bw_beg_file,
//                    bw_csr_file);
//    vertex_t i=0;
//    
//    ///step 2: detect scc
//    while(i++ < run_times)
//    {
//        printf("\nRuntime: %d\n", i);
//        wcc_detection(g,
//                thread_count,
//                avg_time,
//                result_file_wcc,
//                result_file_time);
////        sleep(2);
//
//    }

    ///step 3: print result
//    print_time_result(run_times,
//            avg_time);
//    delete[] avg_time;
    return 0;
}


void bfs_within_wcc_parallel(graph_undirected *g, vertex_t *wcc_id, vertex_t thread_count) {
    const vertex_t vert_count = g->vert_count;
    const index_t *beg_pos = g->beg_pos;
    const vertex_t *csr = g->csr;

    // Step 1: Group vertices by WCC ID
    std::map<vertex_t, std::vector<vertex_t>> wcc_groups;
    for (vertex_t v = 0; v < vert_count; ++v) {
        wcc_groups[wcc_id[v]].push_back(v);
    }

    printf("Number of WCCs: %lu\n", wcc_groups.size());

    // Step 2: BFS for each WCC in parallel
    #pragma omp parallel for num_threads(thread_count) schedule(dynamic)
    for (int i = 0; i < static_cast<int>(wcc_groups.size()); ++i) {
        auto it = std::next(wcc_groups.begin(), i);
        vertex_t wcc_root = it->first;
        const std::vector<vertex_t> &vertices = it->second;

        // Mark visited
        std::set<vertex_t> visited;
        std::queue<vertex_t> q;

        vertex_t start = vertices[0];
        q.push(start);
        visited.insert(start);

        while (!q.empty()) {
            vertex_t u = q.front();
            q.pop();

            index_t nbr_start = beg_pos[u];
            index_t nbr_end = beg_pos[u + 1];

            for (index_t i = nbr_start; i < nbr_end; ++i) {
                vertex_t v = csr[i];
                // Only visit if it's in the same WCC and not yet visited
                if (wcc_id[v] == wcc_root && visited.find(v) == visited.end()) {
                    visited.insert(v);
                    q.push(v);
                }
            }
        }
    }
}
