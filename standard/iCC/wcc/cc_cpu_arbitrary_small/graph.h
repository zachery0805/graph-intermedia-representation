#ifndef __GRAPH_H__
#define __GRAPH_H__
#include "util.h"
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "wtime.h"
class graph_undirected
{
	public:
		index_t *beg_pos;
		vertex_t *csr;
		path_t *weight;
		index_t vert_count;
		index_t edge_count;
        vertex_t v_max_degree;
        vertex_t max_degree;
        vertex_t v_min_degree;
        vertex_t min_degree;

	public:
		graph_undirected(){};
		~graph_undirected(){};
		graph_undirected(const char *beg_file, 
				const char *csr_file);
};
#endif
