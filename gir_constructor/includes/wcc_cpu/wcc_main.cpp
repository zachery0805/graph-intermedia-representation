#include "wtime.h"
#include "graph.h"
//#include "frontier_queue.h"
#include "wcc_common.h"
#include <unistd.h>

int main(int args, char **argv)
{
    printf("args = %d\n", args);
	if(args != 9)
    {
	    std::cout<<"Usage: ./wcc_cpu <fw_beg_file> <fw_csr_file> <bw_beg_file> <bw_csr_file> <thread_count> <run_times> <result_file_wcc> <result_file_time>\n";
        exit(-1);
    }
	
	const char *fw_beg_file = argv[1];
	const char *fw_csr_file = argv[2];
	const char *bw_beg_file = argv[3];
	const char *bw_csr_file = argv[4];
	const vertex_t thread_count=atoi(argv[5]);
    const vertex_t run_times = atoi(argv[6]);
    char *result_file_wcc = argv[7];
    char *result_file_time = argv[8];

	printf("Thread = %d, run_times = %d\n", thread_count, run_times);
    
    double * avg_time = new double[15]; 

    for(vertex_t i = 0; i < 15; ++i)
        avg_time[i] = 0.0;

    ///step 1: load the graph
    graph *g = new graph(fw_beg_file,
                    fw_csr_file,
                    bw_beg_file,
                    bw_csr_file);
    vertex_t i=0;
    
    ///step 2: detect scc
    while(i++ < run_times)
    {
        printf("\nRuntime: %d\n", i);
        wcc_detection(g,
                thread_count,
                avg_time,
                result_file_wcc,
                result_file_time);
//        sleep(2);

    }

    ///step 3: print result
//    print_time_result(run_times,
//            avg_time);
//    delete[] avg_time;
    return 0;
}
