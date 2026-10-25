#ifndef __GRAPH_H__
#define __GRAPH_H__
#include "util.h"
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "wtime.h"
#include <vector>
#include <utility>
class graph
{
	public:
		index_t *fw_beg_pos;
		vertex_t *fw_csr;
		index_t *bw_beg_pos;
		vertex_t *bw_csr;
		path_t *weight;
		index_t vert_count;
		index_t edge_count;
		vertex_t *src_list;
		index_t src_count;

	public:
		graph() : fw_beg_pos(nullptr), fw_csr(nullptr), bw_beg_pos(nullptr), bw_csr(nullptr),
		          weight(nullptr), vert_count(0), edge_count(0), src_list(nullptr), src_count(0) {};
		~graph() { delete[] fw_beg_pos; delete[] fw_csr; delete[] bw_beg_pos; delete[] bw_csr; }
		// Builds forward and backward CSR from an in-memory edge list (used
		// by the benchmark's hold-out protocol). Caller owns the result.
		static graph* fromEdges(index_t V, const std::vector<std::pair<int, int>>& edges);
		graph(const graph&) = delete;
		graph& operator=(const graph&) = delete;
		graph(const char *fw_beg_file, 
				const char *fw_csr_file,
                const char *bw_beg_file,
                const char *bw_csr_file);
		void gen_src(){};
		void groupby(){};
		graph *graph_load(
        	const char *fw_beg_file, 
        	const char *fw_csr_file, 
        	const char *bw_beg_file, 
        	const char *bw_csr_file,
        	double *avg_time
        	)
		{
    		graph *g = new graph(fw_beg_file,
            		fw_csr_file,
            		bw_beg_file,
            		bw_csr_file);

    		for(index_t i=0; i<15; ++i)
        		avg_time[i] = 0.0;
    		return g;
		}

};
#endif
