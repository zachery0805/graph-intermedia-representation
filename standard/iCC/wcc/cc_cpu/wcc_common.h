#ifndef __WCC_COMMON_H
#define __WCC_COMMON_H
#include <omp.h>
#include "graph.h"
#include "util.h"
#include "wtime.h"

void find_cc(
        graph_undirected *g,
        const vertex_t thread_count,
        const char *output_file
        );

//void wcc_detection(
//        graph *g,
//        const index_t thread_count,
//        double *avg_time,
//        char *result_file_wcc,
//        char *result_file_time
//        );

void get_wcc_result(
        index_t *wcc_id,
        index_t vert_count
        );

void print_time_result(
        index_t run_times,
        double *avg_time
        );

    
////step 2.1: pivot selection
//vertex_t pivot_selection(
//        graph *g,
//        vertex_t thread_count
//        )
//{
//    const index_t vert_count = g->vert_count;
//    const index_t edge_count = g->edge_count;
//    vertex_t *max_pivot_list = new vertex_t[thread_count];
//    vertex_t *max_degree_list = new vertex_t[thread_count];
//    vertex_t *fw_beg_pos = g->fw_beg_pos;
//    vertex_t *fw_csr = g->fw_csr;
//    vertex_t *bw_beg_pos = g->bw_beg_pos;
//    vertex_t *bw_csr = g->bw_csr;
//
//    vertex_t step = vert_count / thread_count;
//    if (vert_count != step * thread_count)
//    {
//        step += 1;
//    }
//
//    #pragma omp parallel \
//    num_threads(thread_count) 
//    {
//        const vertex_t tid = omp_get_thread_num();
//        vertex_t vert_beg = tid * step;
//        vertex_t vert_end = (tid == thread_count - 1 ? vert_count : vert_beg + step);
//    
//        vertex_t max_pivot_thread = 0;
//        vertex_t max_degree_thread = 0;
//
//        for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
//        {
//            index_t out_degree = fw_beg_pos[vert_id+1] - fw_beg_pos[vert_id];
//            index_t in_degree = bw_beg_pos[vert_id+1] - bw_beg_pos[vert_id];
//            index_t degree_mul = out_degree * in_degree;
//            
//            if(degree_mul > max_degree_thread)
//            {
//                max_degree_thread = degree_mul;
//                max_pivot_thread = vert_id;
//            }
//        }
//        max_pivot_list[tid] = max_pivot_thread;
//        max_degree_list[tid] = max_degree_thread;
//    }
//    
//    index_t max_pivot = 0;
//    index_t max_degree = 0;
//    
//    for(vertex_t i = 0; i < thread_count; ++i)
//    {
//        if(max_degree_list[i] > max_degree)
//        {
//            max_degree = max_degree_list[i];
//            max_pivot = max_pivot_list[i];
//        }
//    }
//    if(DEBUG)
//    {
//        printf("max_pivot, %d, max_degree, %d\n", max_pivot, max_degree);
//    }
//    return max_pivot;
//}

//inline index_t pivot_selection_from_fq(
//        index_t *scc_id,
//        index_t *fw_beg_pos,
//        index_t *bw_beg_pos,
//        index_t vert_beg,
//        index_t vert_end,
//        vertex_t *fw_csr,
//        vertex_t *bw_csr,
//        index_t *max_pivot_list,
//        index_t *max_degree_list,
//        index_t tid,
//        index_t thread_count,
//        index_t *small_queue
////        ,
////        index_t *mul_degree,
////        index_t *degree_prop
//        )
//{
//    index_t max_pivot_thread = 0;
//    index_t max_degree_thread = 0;
//    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
//    {
//        vertex_t vert_id = small_queue[fq_vert_id];
////        if(scc_id[vert_id] == 0)
////        {
//            index_t degree_mul = (fw_beg_pos[vert_id+1] - fw_beg_pos[vert_id]) * (bw_beg_pos[vert_id+1] - bw_beg_pos[vert_id]);
//            
////            mul_degree[vert_id] = degree_mul + fw_beg_pos[vert_id+1] - fw_beg_pos[vert_id];
////            degree_prop[vert_id] = mul_degree[vert_id];
//            if(degree_mul > max_degree_thread)
//            {
//                max_degree_thread = degree_mul;
//                max_pivot_thread = vert_id;
//            }
////        }
//    }
//    max_pivot_list[tid] = max_pivot_thread;
//    max_degree_list[tid] = max_degree_thread;
//
//    #pragma omp barrier
//    
//    index_t max_pivot = 0;
//    index_t max_degree = 0;
//    
//    for(index_t i=0; i<thread_count; ++i)
//    {
//        if(max_degree_list[i] > max_degree)
//        {
//            max_degree = max_degree_list[i];
//            max_pivot = max_pivot_list[i];
//        }
//    }
//    if(DEBUG)
//    {
//        if(tid == 0)
//            printf("max_pivot, %d, max_degree, %d\n", max_pivot, max_degree);
//    }
//    return max_pivot;
//}

#endif
