#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include <omp.h>
#include <algorithm>
#include "wtime.h"
#include "graph.h"
#include "scc_common.h"
#include "trim_1_gfq.h"
#include "trim_2_3.h"
#include "color_propagation.h"
#include "fw_bw.h"
#include "openmp_wcc.hpp"
#define INF -1
#define DEBUG 0

// 0 trim, 1 largest SCC, 2 small SCC, 3 total time
// 4 trim_size_1, 5 trim_size_2, 6 pivot_selection, 7 fw_bfs, 8 bw_bfs,
// 9 color propagation, 10 color identify, 11 color_init 

void scc_detection(
        const graph *g,
        const int alpha, 
        const int beta,
        const int gamma,
        const double theta,
        const index_t thread_count,
        double *avg_time,
        std::vector<std::vector<int>>& sccs,
        std::map<int, int>& nodeToScc)
{
    std::atomic<index_t> atomic_scc_label(2);
    const index_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    const double avg_degree = edge_count * 1.0 / vert_count;
    if(DEBUG)
        printf("vert_count = %d, edge_count = %d, avg_degree = %.3lf\n", vert_count, edge_count, avg_degree);

    index_t *fw_beg_pos = g->fw_beg_pos;
    vertex_t *fw_csr = g->fw_csr;
    index_t *bw_beg_pos = g->bw_beg_pos;
    vertex_t *bw_csr = g->bw_csr;

    if(VERBOSE)
    {
        for(int i=fw_beg_pos[vert_count]; i<fw_beg_pos[vert_count+1]; ++i)
            printf("%d\n", fw_csr[i]);
    }

    index_t *scc_id = new index_t[vert_count + 1];
    index_t *color = new index_t[vert_count + 1];
    index_t *max_pivot_list = new index_t[thread_count];
    index_t *max_degree_list = new index_t[thread_count];
	
    depth_t *fw_sa;
    depth_t *bw_sa;

    if(posix_memalign((void **)&fw_sa,getpagesize(),
        sizeof(depth_t)*(vert_count+1)))
        perror("posix_memalign");
	
    if(posix_memalign((void **)&bw_sa,getpagesize(),
        sizeof(depth_t)*(vert_count+1)))
        perror("posix_memalign");
    
    index_t *small_queue = new index_t[vert_count + 1];
    index_t *temp_queue = new index_t[vert_count + 1];
    index_t *inter_queue = new index_t[vert_count + 1];
    index_t *wcc_fq= new index_t[vert_count + 1];
    
    index_t *thread_bin = new index_t[thread_count];
    index_t *prefix_sum = new index_t[thread_count];
	
    index_t *front_comm=new index_t[thread_count];	
    index_t *work_comm=new index_t[thread_count];
    bool *color_change = new bool[thread_count];
    memset(color_change, 0, sizeof(bool) * thread_count);

    vertex_t wcc_fq_size = 0;

    if(DEBUG)
    {
        printf("Initialization\n");
    }
    for(index_t i=0; i<vert_count + 1; ++i)
    {
        color[i] = i;
        fw_sa[i] = -1;
        bw_sa[i] = -1;
        scc_id[i] = 0;
    }

    if(DEBUG)
    {
        printf("Parallel starts\n");
    }
    vertex_t vertex_fw = 0;
    vertex_t vertex_bw = 0;
    index_t size_3_1 = 0;
    index_t size_3_2 = 0;
    bool changed = false;
    double end_time;
    double start_time = wtime();
    #pragma omp parallel num_threads(thread_count) 
    {
        const index_t tid = omp_get_thread_num();
        index_t step = vert_count / thread_count;
        index_t vert_beg = tid * step;
        index_t vert_end = (tid == thread_count - 1 ? vert_count : vert_beg + step);
        double time_size_1_first;
        double time_size_1;
        double time_fw;
        double time_bw; 
        double time_size_2;
        double time_size_3;
        double time_gfq;
        double time_color_1;
        double time_color_2;
        double time_color;
        double pivot_time;
        double time_color_init;
        double time_wcc;
        double time_mice_fw_bw; 
        const vertex_t upper_bound = vert_count / thread_count * 5;
        vertex_t *thread_queue = new vertex_t[upper_bound];
        
        double time = wtime();
        index_t trim_times = 1;

        trim_1_first(
                atomic_scc_label,
                scc_id,
                fw_beg_pos,
                bw_beg_pos,
                vert_beg,
                vert_end);
        trim_times ++;
        #pragma omp barrier
        if(tid == 0)
        {
            time_size_1_first = wtime() - time;
        }
        time = wtime();

        generate_frontier_queue(vert_count,
                scc_id,
                thread_count,
                small_queue,
                thread_bin,
                prefix_sum,
                vert_beg,
                vert_end,
                tid);
        #pragma omp barrier

        if(tid == 0)
        {
            time_size_1 = wtime() - time;
        }
        index_t fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
        step = fq_size / thread_count;
        vert_beg = tid * step;
        vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
        
        if(DEBUG)
        {
            if(tid == 0)
            {
                printf("fq_size, %d\n", fq_size);
            }
        }
        
        time = wtime();
        vertex_t root = pivot_selection_from_fq(scc_id,
                        fw_beg_pos,
                        bw_beg_pos,
                        vert_beg,
                        vert_end,
                        fw_csr,
                        bw_csr,
                        max_pivot_list,
                        max_degree_list,
                        tid,
                        thread_count,
                        small_queue);
        #pragma omp barrier
        pivot_time = wtime() - time;

        time = wtime();
        fw_bfs_fq_queue(scc_id,
                fw_beg_pos,
                bw_beg_pos,
                vert_beg,
                vert_end,
                fw_csr,
                bw_csr,
                fw_sa,
                front_comm,
                work_comm,
                root,
                tid,
                thread_count,
                alpha,
                beta,
                gamma,
                small_queue,
                fq_size,
                avg_degree,
                vertex_fw,
                temp_queue,
                prefix_sum,
                upper_bound,
                thread_queue);
        #pragma omp barrier
        time_fw = wtime() - time;

        time = wtime();
        bw_gfq_from_fw(fw_sa, 
                thread_count,
                small_queue,
                thread_bin,
                prefix_sum,
                vert_beg,
                vert_end,
                tid,
                temp_queue);
        #pragma omp barrier
        time_bw = wtime() - time;

        index_t temp_fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
        index_t temp_step = temp_fq_size / thread_count;
        index_t temp_vert_beg = tid * temp_step;
        index_t temp_vert_end = (tid == thread_count - 1 ? temp_fq_size : temp_vert_beg + temp_step);
        
        if(DEBUG)
        {
            if(tid == 0)
            {
                printf("gfq bw, %.3lf, fq_size, %d\n", time_bw * 1000, temp_fq_size);
            }
        }
        #pragma omp barrier
        
        time = wtime();

        bw_bfs_fq_queue(scc_id,
                fw_beg_pos,
                bw_beg_pos,
                temp_vert_beg,
                temp_vert_end,
                fw_csr,
                bw_csr,
                fw_sa,
                bw_sa,
                front_comm,
                work_comm,
                root,
                tid,
                thread_count,
                alpha,
                beta,
                gamma,
                temp_queue,
                temp_fq_size,
                avg_degree,
                vertex_bw,
                inter_queue,
                prefix_sum,
                upper_bound,
                thread_queue);
        trim_times = 0;
        #pragma omp barrier
        time_bw += wtime() - time;
        delete[] thread_queue;
    
        vertex_t prev_fq_size = fq_size;

        while(prev_fq_size > 0 && trim_times < TRIM_TIMES)
        {
            time = wtime();
            trim_1_from_fq_gfq(
                    atomic_scc_label,
                    scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr,
                    thread_count,
                    small_queue,
                    thread_bin,
                    prefix_sum,
                    tid,
                    temp_queue);
            #pragma omp barrier
            
            if(tid == 0)
            {
                time_size_1 += wtime() - time;
            }

            fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
            
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("trim_1, 3rd, %.3lf, trimmed vertices, %d\n", time_size_1 * 1000, prev_fq_size - fq_size);
                }
            }
            step = fq_size / thread_count;
            vert_beg = tid * step;
            vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
            if(prev_fq_size - fq_size < prev_fq_size * theta)
                break;
            prev_fq_size = fq_size;
            trim_times += 1;
            #pragma omp barrier
        }
        if(DEBUG)
        {
            if(tid == 0)
                printf("trim-1 times before trim-2-3, %d\n", trim_times);
        }

        if(fq_size > 0)
        {
            fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
            step = fq_size / thread_count;
            vert_beg = tid * step;
            vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);

            if(DEBUG)
            {
                if(tid == 0)
                    printf("fq_size, %d\n", fq_size);
            }
            time = wtime();
            trim_2_from_fq(
                    atomic_scc_label,
                    scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr,
                    small_queue);
            #pragma omp barrier
            if(tid == 0)
                time_size_2 = wtime() - time;
            
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("time_size_2, %.3lf, fq_size, %d\n", time_size_2 * 1000, fq_size);
                }
            }
            trim_times = 0;
            #pragma omp barrier
            
            time = wtime();

            trim_3_1_from_fq(
                    atomic_scc_label,
                    scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr,
                    small_queue);
            #pragma omp barrier
            double time_size_3_1 = wtime() - time;

            if(DEBUG)
            {
                if(tid == 0)
                    printf("time_size_3_1, %.3lf\n", time_size_3_1 * 1000);
            }

            time = wtime();

            trim_3_2_from_fq(
                    atomic_scc_label,
                    scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    vert_beg,
                    vert_end,
                    fw_csr,
                    bw_csr,
                    small_queue);
            
            trim_times = 0;
            
            #pragma omp barrier
            double time_size_3_2 = wtime() - time;

            if(DEBUG)
            {
                if(tid == 0)
                    printf("time_size_3_2, %.3lf\n", time_size_3_2 * 1000);
            }
            time_size_3 = time_size_3_1 + time_size_3_2;
            
            time = wtime();
                
            prev_fq_size = fq_size;
            while(prev_fq_size > 0)
            {
                time = wtime();
                trim_1_from_fq_gfq(
                        atomic_scc_label,
                        scc_id,
                        fw_beg_pos,
                        bw_beg_pos,
                        vert_beg,
                        vert_end,
                        fw_csr,
                        bw_csr,
                        thread_count,
                        small_queue,
                        thread_bin,
                        prefix_sum,
                        tid,
                        temp_queue);
                #pragma omp barrier
                if(tid == 0)
                    time_size_1 += wtime() - time;

                fq_size = prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
                step = fq_size / thread_count;
                vert_beg = tid * step;
                vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
                if(DEBUG)
                {
                    if(tid == 0)
                    {
                        printf("trim_1, 4th, %.3lf, fq_size, %d, trimmed vertices, %d\n", time_size_1 * 1000, fq_size, prev_fq_size - fq_size);
                    }
                }
                if(prev_fq_size - fq_size < prev_fq_size * theta)
                    break;
                prev_fq_size = fq_size;
                #pragma omp barrier
                break;
            }
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("trim_1 times before Color, %d, fq_size, %d\n", trim_times, fq_size);
                }
            }
          
            time = wtime();
            coloring_wcc(fq_size,
                    scc_id,
                    thread_count,
                    small_queue,
                    vert_beg,
                    vert_end,
                    tid,
                    color,
                    color_change,
                    fw_beg_pos,
                    fw_csr,
                    bw_beg_pos,
                    bw_csr);
            #pragma omp barrier
            if(tid == 0)
            {
                time_wcc = wtime() - time;

            }

            time = wtime();
            init_fw_sa(vert_beg,
                    vert_end,
                    fw_sa,
                    small_queue,
                    fq_size,
                    wcc_fq,
                    tid,
                    color,
                    wcc_fq_size);
            #pragma omp barrier
            if(DEBUG)
            {
                if(tid == 0)
                {
                    printf("wcc_fq, %d\n", wcc_fq_size);
                }
            }
            mice_fw_bw(color,
                    scc_id,
                    fw_beg_pos,
                    bw_beg_pos,
                    fw_csr,
                    bw_csr,
                    fw_sa,
                    tid,
                    thread_count,
                    small_queue,
                    fq_size,
                    wcc_fq,
                    wcc_fq_size
                    ); 
            #pragma omp barrier
            if(tid == 0)
            {
                time_mice_fw_bw = wtime() - time;
            }
            
        }
        if(tid == 0)
        {
            avg_time[0] += time_size_1_first + time_size_1 + time_size_2 + time_size_3;
            avg_time[1] += time_fw + time_bw + pivot_time;
            avg_time[2] += time_wcc + time_mice_fw_bw;
            avg_time[4] += time_size_1_first + time_size_1;
            avg_time[5] += time_size_2;
            avg_time[6] += pivot_time;
            avg_time[7] += time_fw;
            avg_time[8] += time_bw;
            avg_time[9] += time_wcc;
            avg_time[10] += time_mice_fw_bw;
            avg_time[13] += time_size_3;
            avg_time[14] += time_gfq;
        }
        if(OUTPUT_TIME)
        {
            if(tid == 0)
            {
                printf("\ntime size_1_first, %.3lf\ntime size_1, %.3lf\ntime pivot, %.3lf\nlargest fw, %.3lf\nlargest bw, %.3lf\nlargest fw/bw, %.3lf\ntrim size_2, %.3lf\ntrim size_3, %.3lf\nwcc time, %.3lf\nmice fw-bw time, %.3lf\nmice scc time, %.3lf\ntotal time, %.3lf\n", time_size_1_first * 1000, time_size_1 * 1000, pivot_time * 1000, time_fw * 1000, time_bw * 1000, (pivot_time + time_fw + time_bw) * 1000, time_size_2 * 1000, time_size_3 * 1000, time_wcc * 1000, time_mice_fw_bw * 1000, (time_wcc + time_mice_fw_bw) * 1000, (time_size_1_first + time_size_1 + pivot_time + time_fw + time_bw + time_size_2 + time_size_3 + time_wcc + time_mice_fw_bw) * 1000);
            }
        }
        #pragma omp barrier
    }
    end_time = wtime() - start_time;
    avg_time[3] += end_time;
    if(DEBUG)
        printf("total time, %.3lf\n", end_time * 1000);

    get_scc_result(g,
            scc_id,
            vert_count);

  // Populate the output sccs and nodeToScc map
    std::map<int, std::vector<int>> tempSccs;
    sccs.clear();
    nodeToScc.clear();

    for (int i = 0; i < vert_count; ++i) {
        int current_scc_id = static_cast<int>(scc_id[i]);
        nodeToScc[i] = current_scc_id;
        // Assuming 0 is not a valid SCC ID for grouping, or represents trivial SCCs.
        if (current_scc_id != 0) {
            tempSccs[current_scc_id].push_back(i);
        }
    }

    if (DEBUG)
    {
        printf("--- Node to SCC ID Mapping ---\n");
        for (const auto& pair : nodeToScc)
        {
            printf("Node %d -> SCC ID %d\n", pair.first, pair.second);
        }
        printf("-----------------------------\n");
    }

    for (const auto& pair : tempSccs) {
        sccs.push_back(pair.second);
    }

    delete[] scc_id;
    delete[] color;
    delete[] max_pivot_list;
    delete[] max_degree_list;
    delete[] small_queue;
    delete[] temp_queue;
    delete[] thread_bin;
    delete[] prefix_sum;
    delete[] front_comm;
	delete[] work_comm;
    delete[] color_change;
//    delete[] wcc_color;
//    delete[] color_redirect;
//    delete[] is_redirect;

}

