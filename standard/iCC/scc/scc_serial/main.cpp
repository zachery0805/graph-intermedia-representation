// A C++ program to find articulation points in an undirected graph 
#include <iostream> 
#include <list> 
#include <unistd.h>

#include "scc_dfs.h"
#include "wtime.h"
#include "graph.h"
#define NIL -1 
using namespace std; 


// Driver program to test above function 
int main(int args, char **argv)
{ 
	if(args != 5)
    {
	    std::cout<<"Usage: ./bicc <beg_file> <csr_file> <thread_count> <run_times>\n";
        exit(-1);
    }
	
	const char *beg_file = argv[1];
	const char *csr_file = argv[2];
	const index_t thread_count=atoi(argv[3]);
    const index_t run_times = atoi(argv[4]);
    
//    double * avg_time = new double[15]; 
    ///step 1: load the graph
    graph_undirected *g = new graph_undirected(beg_file, csr_file);
    
    double time_begin = wtime();
    find_scc(g); 
    double time_end = wtime();

    printf("Time(s),%.3lf\n", (time_end - time_begin));

    return 0; 
} 

