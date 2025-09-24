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
#include <queue>
#include <map>
#include "wtime.h"
#include "graph.h"
////#include "frontier_queue.h"
#include "wcc_common.h"
////#include "trim_1.h"
//#include "trim_1_gfq.h"
//#include "trim_2_3.h"
//#include "fw_bw.h"
//#include "openmp_wcc.hpp"
#define INF 0x7fffffff

void bfs_sequential(graph *g, vertex_t root, vertex_t *wcc_id)
//void wcc_detect(graph *g, index_t source, index_t *wcc_color, index_t color_index, index_t *color_size)
{   
    std::queue<vertex_t> q; 
    q.push(root);
    wcc_id[root] = root;
    while(!q.empty())
    {
        vertex_t cur = q.front();
        q.pop();
        //fw
        for(index_t i = g->fw_beg_pos[cur]; i < g->fw_beg_pos[cur + 1]; ++i)
        {       
            vertex_t w = g->fw_csr[i];
            if(wcc_id[w] == INF)
            {
                q.push(w);
                wcc_id[w] = root;
            }
        } 
        //bw
        //Drawback: duplicate caused by edges a <--> b
        for(index_t i = g->bw_beg_pos[cur]; i < g->bw_beg_pos[cur + 1]; ++i)
        {       
            vertex_t w = g->bw_csr[i];
            if(wcc_id[w] == INF)
            {
                q.push(w);
                wcc_id[w] = root;
            }
        }   
    }
} 

void wcc_bfs_sequential(graph *g, vertex_t vert_count, vertex_t *wcc_id)
{
    for(vertex_t v = 0; v < vert_count; ++v)
    {
        if(wcc_id[v] == INF)
        {
//            printf("root = %d\n", v);
            bfs_sequential(g, v, wcc_id);
        }
    }
}

void get_queue(
        vertex_t *thread_queue,
        vertex_t *thread_bin,
        index_t *prefix_sum,
        index_t tid,
        vertex_t *temp_queue
        )
{
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
//    #pragma omp barrier
    for(index_t i = 0; i < tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }
//    #pragma omp barrier
    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = start_pos; vert_id < start_pos + thread_bin[tid]; ++vert_id)
    {
        temp_queue[vert_id] = thread_queue[vert_id - start_pos];
    }
//    #pragma omp barrier
}

////step 2.1: pivot selection
vertex_t pivot_selection(
        graph *g,
        vertex_t thread_count
        )
{
    const index_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    vertex_t *max_pivot_list = new vertex_t[thread_count];
    vertex_t *max_degree_list = new vertex_t[thread_count];
    vertex_t *fw_beg_pos = g->fw_beg_pos;
    vertex_t *fw_csr = g->fw_csr;
    vertex_t *bw_beg_pos = g->bw_beg_pos;
    vertex_t *bw_csr = g->bw_csr;

    vertex_t step = vert_count / thread_count;
    if (vert_count != step * thread_count)
    {
        step += 1;
    }

    #pragma omp parallel \
    num_threads(thread_count) 
    {
        const vertex_t tid = omp_get_thread_num();
        vertex_t vert_beg = tid * step;
        vertex_t vert_end = (tid == thread_count - 1 ? vert_count : vert_beg + step);
    
        vertex_t max_pivot_thread = 0;
        vertex_t max_degree_thread = 0;

        for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
        {
            vertex_t out_degree = fw_beg_pos[vert_id+1] - fw_beg_pos[vert_id];
//            vertex_t in_degree = bw_beg_pos[vert_id+1] - bw_beg_pos[vert_id];
 //           vertex_t degree_mul = out_degree * in_degree;
            
//            if(degree_mul > max_degree_thread)
            if(out_degree > max_degree_thread)
            {
                max_degree_thread = out_degree;
                max_pivot_thread = vert_id;
            }
        }
        max_pivot_list[tid] = max_pivot_thread;
        max_degree_list[tid] = max_degree_thread;
    }
    
    vertex_t max_pivot = 0;
    vertex_t max_degree = 0;
    
    for(vertex_t i = 0; i < thread_count; ++i)
    {
        if(max_degree_list[i] > max_degree)
        {
            max_degree = max_degree_list[i];
            max_pivot = max_pivot_list[i];
        }
    }
//    if(DEBUG)
//    {
//        printf("max_pivot, %d, max_degree, %d\n", max_pivot, max_degree);
//    }
    return max_pivot;
}
void bfs_concurrent(graph *g,
        vertex_t *wcc_id,
        depth_t *sa,
        vertex_t *frontier_queue,
        vertex_t *front_comm,
        vertex_t *work_comm,
        vertex_t root,
        vertex_t thread_count,
        vertex_t *temp_queue
        )
{
    const index_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    const double avg_degree = edge_count * 1.0 / vert_count;

    index_t *fw_beg_pos = g->fw_beg_pos;
    vertex_t *fw_csr = g->fw_csr;
    index_t *bw_beg_pos = g->bw_beg_pos;
    vertex_t *bw_csr = g->bw_csr;

    vertex_t *thread_bin = new vertex_t[thread_count];
    vertex_t *prefix_sum = new vertex_t[thread_count];
    bool *color_change = new bool[thread_count];
    memset(color_change, 0, sizeof(bool) * thread_count);
    
//    printf("ok");


    vertex_t step = vert_count / thread_count;
//    vertex_t fq_step = 0;
    vertex_t queue_size = 1;

    if (vert_count != step * thread_count)
    {
        step += 1;
    }
    const vertex_t upper_bound = vert_count / thread_count * 7;
    if(VERBOSE) 
        printf("root, %u, vert_count, %u, thread_count, %u, step, %d\n", root, vert_count, thread_count, step);
    
    double time_bfs = wtime(); 

    sa[root] = 0;
    wcc_id[root] = root;
    frontier_queue[0] = root;
   
    vertex_t my_beg = fw_beg_pos[root];
    vertex_t my_end = fw_beg_pos[root + 1];
    
    vertex_t index_temp = 0;
    for(; my_beg < my_end && index_temp < thread_count - 1; my_beg++, index_temp++)
    {
        vertex_t v = fw_csr[my_beg];
        sa[v] = 0;
        wcc_id[v] = root;
        frontier_queue[index_temp] = v;
    }


    
    #pragma omp parallel \
    num_threads(thread_count) 
    {
        const vertex_t tid = omp_get_thread_num();
        vertex_t vert_beg = tid * step;
        vertex_t vert_end = (tid == thread_count - 1 ? vert_count : vert_beg + step);
        vertex_t *thread_queue = new vertex_t[upper_bound];


        bool is_top_down = true;	
        bool is_top_down_queue = false;
        depth_t level = 0;
        vertex_t vertex_visited = 0;
        while(true)
        {
//            #pragma omp barrier
            double ltm = wtime();
            index_t front_count = 0;
//            index_t my_work_next = 0;
//            index_t my_work_curr = 0;
            vertex_t vertex_frontier = 0;

            if(is_top_down)
            {
//                for(vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++)
                vertex_t fq_step = queue_size / thread_count;
                 
                vertex_t fq_vert_beg = tid * fq_step;
                vertex_t fq_vert_end = (tid == thread_count - 1 ? queue_size : fq_vert_beg + fq_step);

//                printf("tid, %u, fq_step, %u, work_load, %u\n", tid, fq_step, fq_vert_end - fq_vert_beg);
//                #pragma omp barrier         
//                exit(EXIT_FAILURE);

                for(vertex_t index = fq_vert_beg; index < fq_vert_end; index++)
                {
                    vertex_t vert_id = frontier_queue[index];

//                    if(wcc_id[vert_id] == root && sa[vert_id] == level)
//                    {

                        //fw
                        index_t my_beg = fw_beg_pos[vert_id];
                        index_t my_end = fw_beg_pos[vert_id + 1];

                        for(; my_beg < my_end; my_beg++)
                        {
                            vertex_t nebr = fw_csr[my_beg];
                            //if(wcc_id[nebr] == INF && sa[nebr] == -1)
                            if(wcc_id[nebr] == INF)
                            {
                                sa[nebr] = level + 1;
                                wcc_id[nebr] = root;
//                                my_work_next += fw_beg_pos[nebr + 1] - fw_beg_pos[nebr];
                                thread_queue[front_count++] = nebr;
//                                front_count++;
                            }
                        }
                        //bw
                        my_beg = bw_beg_pos[vert_id];
                        my_end = bw_beg_pos[vert_id + 1];

                        for(; my_beg < my_end; my_beg++)
                        {
                            vertex_t nebr = bw_csr[my_beg];
                            //if(wcc_id[nebr] == INF && sa[nebr] == -1)
                            if(wcc_id[nebr] == INF)
                            {
                                sa[nebr] = level + 1;
                                wcc_id[nebr] = root;
//                                my_work_next += bw_beg_pos[nebr + 1] - bw_beg_pos[nebr];
                                thread_queue[front_count++] = nebr;

//                                front_count++;
                            }
                        }
//                    }
                }
//                work_comm[tid] = my_work_next;
//                front_comm[tid] = front_count; 
            }
            else
                if(!is_top_down_queue)
                {
                    for(vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++)
                    {
                        if(wcc_id[vert_id] == INF && sa[vert_id] == -1)
                        {
                            //fw
                            index_t my_beg = fw_beg_pos[vert_id];
                            index_t my_end = fw_beg_pos[vert_id + 1];
                            //my_work_curr += my_end - my_beg;

                            bool is_marked = false;

                            for(; my_beg < my_end; my_beg++)
                            {
                                vertex_t nebr = fw_csr[my_beg];
// Version 1: BFS
//                                if(sa[nebr] == level)
// Version 2: Rsync
                                if(sa[nebr] != -1)
                                {
                                    sa[vert_id] = level + 1;
                                    wcc_id[vert_id] = root;
                                    front_count++;
                                    is_marked = true;
                                    break;
                                }
                            }
                            if(is_marked)
                                continue;
                            //bw
                            my_beg = bw_beg_pos[vert_id];
                            my_end = bw_beg_pos[vert_id + 1];
                            //my_work_curr += my_end - my_beg;

                            for(; my_beg < my_end; my_beg++)
                            {
                                vertex_t nebr = bw_csr[my_beg];
// Version 1: BFS
//                                if(sa[nebr] == level)

// Version 2: Rsync
                                if(sa[nebr] != -1)
                                {
                                    sa[vert_id] = level + 1;
                                    wcc_id[vert_id] = root;
                                    front_count++;
                                    break;
                                }
                            }
                        }
                    }
//                    work_comm[tid] = my_work_curr;
                //    front_comm[tid] = front_count; 
                }
// Step 3: Async top-down, queue-based
                else
                {
                    vertex_t end_queue = upper_bound;
                    index_t head = 0;
                    index_t tail = 0;
                    //std::queue<index_t> q;
                    vertex_t step = queue_size / thread_count;
                    vertex_t queue_beg = tid * step;
                    vertex_t queue_end = (tid == thread_count - 1 ? queue_size: queue_beg + step);

                    //Option 1: put current level vertices into fq
                    for(vertex_t q_vert_id = queue_beg; q_vert_id < queue_end; q_vert_id++)
                    {
//                        thread_queue[tail] = temp_queue[q_vert_id];
                        thread_queue[tail] = frontier_queue[q_vert_id];
                        tail ++;
                    }
                    while(head != tail)
                    {
                        vertex_t temp_v = thread_queue[head++];
    //                    front_count ++;
                        if(head == end_queue)
                            head = 0;
                        // fw
                        index_t my_beg = fw_beg_pos[temp_v];
                        index_t my_end = fw_beg_pos[temp_v + 1];

                        for(; my_beg < my_end; ++my_beg)
                        {
                            vertex_t w = fw_csr[my_beg];
                            
                            if(wcc_id[w] == INF && sa[w] == -1)
                            {
                                thread_queue[tail++] = w;
                                if(tail == end_queue)
                                    tail = 0;
                                sa[w] = level + 1;
                                wcc_id[w] = root;
                            }
                        }
                        // bw
                        my_beg = bw_beg_pos[temp_v];
                        my_end = bw_beg_pos[temp_v + 1];

                        for(; my_beg < my_end; ++my_beg)
                        {
                            vertex_t w = bw_csr[my_beg];
                            
                            if(wcc_id[w] == INF && sa[w] == -1)
                            {
                                thread_queue[tail++] = w;
                                if(tail == end_queue)
                                    tail = 0;
                                sa[w] = level + 1;
                                wcc_id[w] = root;
                            }
                        }
                    }
                }
            front_comm[tid] = front_count;
//            work_comm[tid] = my_work_curr;
            #pragma omp barrier
            //printf("!!!!\n");
            front_count = 0;
            //my_work_next = 0;
    
            for(vertex_t i = 0; i < thread_count; ++i)
            {
                front_count += front_comm[i];
               // my_work_next += work_comm[i];
            }
                
            if(tid == 0)
                vertex_visited += front_count;

            if(VERBOSE)
            {
                double edge_frontier = (double)front_count * avg_degree;
                double edge_remaider = (double)(vert_count - vertex_visited) * avg_degree;
                if(tid==0 && level < 50) 
                    std::cout<<"Level-"<<(int)level<<" "
    //				<<"-frontier-time-visited:"
                    <<front_count<<" "
 //                   <<fq_size<<" "
//                    <<(double)(fq_size)/front_count<<" "
                    <<(wtime() - ltm) * 1000<<"ms "
                    <<vertex_visited<<" "
                    <<(int)edge_frontier<<" "
                    <<(int)edge_remaider<<" "
                    <<(int)edge_remaider/ALPHA<<" "
                    <<"\n";
            }
            
            if(front_count == 0) 
                break;
            
            if(is_top_down) 
            {
                double edge_frontier = (double)front_count * avg_degree;
                double edge_remainder = (double)(vert_count - vertex_visited) * avg_degree;
    //            printf("edge_remainder/alpha = %g, edge_froniter = %g\n", edge_remainder / alpha, edge_frontier);
                if(!is_top_down_queue && (edge_remainder / ALPHA) < edge_frontier)
                {
                    is_top_down = false;
                    if(VERBOSE)
                    {
                        if(tid==0)
                        {
    //                        double Nf = vertex_frontier;
    //                        double Nu = fq_size - vertex_visited;
    //                        double Mf = Nf * Nf / Nu + avg_degree * (Nu - Nf);
    //                        printf("mf=%.0lf, mu=%.0lf, alpha=%d, Mf=%.0lf\n", edge_frontier, edge_remainder, ALPHA, Mf);
                            std::cout<<"--->Switch to bottom up\n";
                        }
                    }
                }

            }
            else
                if((!is_top_down && !is_top_down_queue && (vert_count * 1.0 / BETA) > front_count) || (!is_top_down && !is_top_down_queue && level > GAMMA))
    //            if(level > 10)
                {
                    //if(!is_top_down_queue)
                    //
// Generate thread_queue for the final Async
                    double time_switch = wtime();
                    front_count = 0;
                    for(vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++)
                    {
//                        vertex_t vert_id = frontier_queue[fq_vert_id];
//                        vertex_t vert_id = [vert_id];
                        if(sa[vert_id] == level + 1)
                        {
                            thread_queue[front_count] = vert_id;
                            front_count ++;
                        }
                    }
                    front_comm[tid] = front_count;
                    
                    is_top_down = false;
                    is_top_down_queue = true;

                    if(VERBOSE)
                    {
                        if(tid==0)
                            std::cout<<"--->Switch to top down queue, "<<(wtime() - time_switch) * 1000<<"\n";
                    }
                }

//            #pragma omp barrier
            if(is_top_down || is_top_down_queue)
            //if(is_top_down_queue)
            {
                double queue_time = wtime();
                get_queue(thread_queue,
                        front_comm,
                        prefix_sum,
                        tid,
//                        temp_queue);
                        frontier_queue);
                
//                #pragma omp barrier
                if(tid == thread_count - 1)
                    queue_size = prefix_sum[thread_count - 1] + front_comm[thread_count - 1];
//                if(VERBOSE)
//                {
//                    if(tid == 0)
//                        printf("queue_size, %d, queue_time (ms), %.3lf\n", queue_size, (wtime() - queue_time) * 1000);
//                }
            }

            #pragma omp barrier
            level++;
//            if(tid == 0)
//            {
//                printf("root, %u, level, %d, front_count, %u, work_next, %u\n", root, level, front_count, my_work_next);
//            }
        }

        
    }
    if(VERBOSE)
        printf("BFS kernel time (ms), %.3lf\n", (wtime() - time_bfs) * 1000);
//    
//    vertex_t colored_num = 0;
//    for(vertex_t v = 0; v < vert_count; ++v)
//    {
//        if(wcc_id[v] != INF)
//            colored_num ++;
//    }
//    printf("root, %u, colored_num, %u\n", root, colored_num);
//    
//    exit(EXIT_FAILURE);
}
//void bfs_parallel(graph *g, vertex_t root, vertex_t *wcc_id)
void bfs_parallel(graph *g,
        vertex_t *wcc_id,
        depth_t *sa,
        vertex_t *frontier_queue,
        vertex_t *front_comm,
        vertex_t *work_comm,
        vertex_t root,
        vertex_t thread_count,
        vertex_t *temp_queue
        )
{
    const index_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    const double avg_degree = edge_count * 1.0 / vert_count;

    index_t *fw_beg_pos = g->fw_beg_pos;
    vertex_t *fw_csr = g->fw_csr;
    index_t *bw_beg_pos = g->bw_beg_pos;
    vertex_t *bw_csr = g->bw_csr;

    vertex_t *thread_bin = new vertex_t[thread_count];
    vertex_t *prefix_sum = new vertex_t[thread_count];
    bool *color_change = new bool[thread_count];
    memset(color_change, 0, sizeof(bool) * thread_count);
    
//    printf("ok");

    sa[root] = 0;
    wcc_id[root] = root;
    frontier_queue[0] = root;
//    fq_size = 1;

    vertex_t step = vert_count / thread_count;
//    vertex_t fq_step = 0;
    vertex_t queue_size = 1;

    if (vert_count != step * thread_count)
    {
        step += 1;
    }
    const vertex_t upper_bound = vert_count / thread_count * 7;
    if(VERBOSE) 
        printf("root, %u, vert_count, %u, thread_count, %u, step, %d\n", root, vert_count, thread_count, step);
    
    
    double time_bfs = wtime(); 
    #pragma omp parallel \
    num_threads(thread_count) 
    {
        const vertex_t tid = omp_get_thread_num();
        vertex_t vert_beg = tid * step;
        vertex_t vert_end = (tid == thread_count - 1 ? vert_count : vert_beg + step);
        vertex_t *thread_queue = new vertex_t[upper_bound];


        bool is_top_down = true;	
        bool is_top_down_queue = false;
        depth_t level = 0;
        vertex_t vertex_visited = 0;
        while(true)
        {
//            #pragma omp barrier
            double ltm = wtime();
            index_t front_count = 0;
//            index_t my_work_next = 0;
//            index_t my_work_curr = 0;
            vertex_t vertex_frontier = 0;

            if(is_top_down)
            {
//                for(vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++)
                vertex_t fq_step = queue_size / thread_count;

       //         if (queue_size != fq_step * thread_count)
       //         {
       //             fq_step += 1;
       //         }
                 
                vertex_t fq_vert_beg = tid * fq_step;
                vertex_t fq_vert_end = (tid == thread_count - 1 ? queue_size : fq_vert_beg + fq_step);

//                printf("tid, %u, fq_step, %u, work_load, %u\n", tid, fq_step, fq_vert_end - fq_vert_beg);
//                #pragma omp barrier         
//                exit(EXIT_FAILURE);

                for(vertex_t index = fq_vert_beg; index < fq_vert_end; index++)
                {
                    vertex_t vert_id = frontier_queue[index];

//                    if(wcc_id[vert_id] == root && sa[vert_id] == level)
//                    {

                        //fw
                        index_t my_beg = fw_beg_pos[vert_id];
                        index_t my_end = fw_beg_pos[vert_id + 1];

                        for(; my_beg < my_end; my_beg++)
                        {
                            vertex_t nebr = fw_csr[my_beg];
                            //if(wcc_id[nebr] == INF && sa[nebr] == -1)
                            if(wcc_id[nebr] == INF)
                            {
                                sa[nebr] = level + 1;
                                wcc_id[nebr] = root;
//                                my_work_next += fw_beg_pos[nebr + 1] - fw_beg_pos[nebr];
                                thread_queue[front_count++] = nebr;
//                                front_count++;
                            }
                        }
                        //bw
                        my_beg = bw_beg_pos[vert_id];
                        my_end = bw_beg_pos[vert_id + 1];

                        for(; my_beg < my_end; my_beg++)
                        {
                            vertex_t nebr = bw_csr[my_beg];
                            //if(wcc_id[nebr] == INF && sa[nebr] == -1)
                            if(wcc_id[nebr] == INF)
                            {
                                sa[nebr] = level + 1;
                                wcc_id[nebr] = root;
//                                my_work_next += bw_beg_pos[nebr + 1] - bw_beg_pos[nebr];
                                thread_queue[front_count++] = nebr;

//                                front_count++;
                            }
                        }
//                    }
                }
//                work_comm[tid] = my_work_next;
//                front_comm[tid] = front_count; 
            }
            else
                if(!is_top_down_queue)
                {
                    for(vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++)
                    {
                        if(wcc_id[vert_id] == INF && sa[vert_id] == -1)
                        {
                            //fw
                            index_t my_beg = fw_beg_pos[vert_id];
                            index_t my_end = fw_beg_pos[vert_id + 1];
                            //my_work_curr += my_end - my_beg;

                            bool is_marked = false;

                            for(; my_beg < my_end; my_beg++)
                            {
                                vertex_t nebr = fw_csr[my_beg];
// Version 1: BFS
//                                if(sa[nebr] == level)
// Version 2: Rsync
                                if(sa[nebr] != -1)
                                {
                                    sa[vert_id] = level + 1;
                                    wcc_id[vert_id] = root;
                                    front_count++;
                                    is_marked = true;
                                    break;
                                }
                            }
                            if(is_marked)
                                continue;
                            //bw
                            my_beg = bw_beg_pos[vert_id];
                            my_end = bw_beg_pos[vert_id + 1];
                            //my_work_curr += my_end - my_beg;

                            for(; my_beg < my_end; my_beg++)
                            {
                                vertex_t nebr = bw_csr[my_beg];
// Version 1: BFS
//                                if(sa[nebr] == level)

// Version 2: Rsync
                                if(sa[nebr] != -1)
                                {
                                    sa[vert_id] = level + 1;
                                    wcc_id[vert_id] = root;
                                    front_count++;
                                    break;
                                }
                            }
                        }
                    }
//                    work_comm[tid] = my_work_curr;
                //    front_comm[tid] = front_count; 
                }
// Step 3: Async top-down, queue-based
                else
                {
                    vertex_t end_queue = upper_bound;
                    index_t head = 0;
                    index_t tail = 0;
                    //std::queue<index_t> q;
                    vertex_t step = queue_size / thread_count;
                    vertex_t queue_beg = tid * step;
                    vertex_t queue_end = (tid == thread_count - 1 ? queue_size: queue_beg + step);

                    //Option 1: put current level vertices into fq
                    for(vertex_t q_vert_id = queue_beg; q_vert_id < queue_end; q_vert_id++)
                    {
//                        thread_queue[tail] = temp_queue[q_vert_id];
                        thread_queue[tail] = frontier_queue[q_vert_id];
                        tail ++;
                    }
                    while(head != tail)
                    {
                        vertex_t temp_v = thread_queue[head++];
    //                    front_count ++;
                        if(head == end_queue)
                            head = 0;
                        // fw
                        index_t my_beg = fw_beg_pos[temp_v];
                        index_t my_end = fw_beg_pos[temp_v + 1];

                        for(; my_beg < my_end; ++my_beg)
                        {
                            vertex_t w = fw_csr[my_beg];
                            
                            if(wcc_id[w] == INF && sa[w] == -1)
                            {
                                thread_queue[tail++] = w;
                                if(tail == end_queue)
                                    tail = 0;
                                sa[w] = level + 1;
                                wcc_id[w] = root;
                            }
                        }
                        // bw
                        my_beg = bw_beg_pos[temp_v];
                        my_end = bw_beg_pos[temp_v + 1];

                        for(; my_beg < my_end; ++my_beg)
                        {
                            vertex_t w = bw_csr[my_beg];
                            
                            if(wcc_id[w] == INF && sa[w] == -1)
                            {
                                thread_queue[tail++] = w;
                                if(tail == end_queue)
                                    tail = 0;
                                sa[w] = level + 1;
                                wcc_id[w] = root;
                            }
                        }
                    }
                }
            front_comm[tid] = front_count;
//            work_comm[tid] = my_work_curr;
            #pragma omp barrier
            //printf("!!!!\n");
            front_count = 0;
            //my_work_next = 0;
    
            for(vertex_t i = 0; i < thread_count; ++i)
            {
                front_count += front_comm[i];
               // my_work_next += work_comm[i];
            }
                
            if(tid == 0)
                vertex_visited += front_count;

            if(VERBOSE)
            {
                double edge_frontier = (double)front_count * avg_degree;
                double edge_remaider = (double)(vert_count - vertex_visited) * avg_degree;
                if(tid==0 && level < 50) 
                    std::cout<<"Level-"<<(int)level<<" "
    //				<<"-frontier-time-visited:"
                    <<front_count<<" "
 //                   <<fq_size<<" "
//                    <<(double)(fq_size)/front_count<<" "
                    <<(wtime() - ltm) * 1000<<"ms "
                    <<vertex_visited<<" "
                    <<(int)edge_frontier<<" "
                    <<(int)edge_remaider<<" "
                    <<(int)edge_remaider/ALPHA<<" "
                    <<"\n";
            }
            
            if(front_count == 0) 
                break;
            
            if(is_top_down) 
            {
                double edge_frontier = (double)front_count * avg_degree;
                double edge_remainder = (double)(vert_count - vertex_visited) * avg_degree;
    //            printf("edge_remainder/alpha = %g, edge_froniter = %g\n", edge_remainder / alpha, edge_frontier);
                if(!is_top_down_queue && (edge_remainder / ALPHA) < edge_frontier)
                {
                    is_top_down = false;
                    if(VERBOSE)
                    {
                        if(tid==0)
                        {
    //                        double Nf = vertex_frontier;
    //                        double Nu = fq_size - vertex_visited;
    //                        double Mf = Nf * Nf / Nu + avg_degree * (Nu - Nf);
    //                        printf("mf=%.0lf, mu=%.0lf, alpha=%d, Mf=%.0lf\n", edge_frontier, edge_remainder, ALPHA, Mf);
                            std::cout<<"--->Switch to bottom up\n";
                        }
                    }
                }

            }
            else
                if((!is_top_down && !is_top_down_queue && (vert_count * 1.0 / BETA) > front_count) || (!is_top_down && !is_top_down_queue && level > GAMMA))
    //            if(level > 10)
                {
                    //if(!is_top_down_queue)
                    //
// Generate thread_queue for the final Async
                    double time_switch = wtime();
                    front_count = 0;
                    for(vertex_t vert_id = vert_beg; vert_id < vert_end; vert_id++)
                    {
//                        vertex_t vert_id = frontier_queue[fq_vert_id];
//                        vertex_t vert_id = [vert_id];
                        if(sa[vert_id] == level + 1)
                        {
                            thread_queue[front_count] = vert_id;
                            front_count ++;
                        }
                    }
                    front_comm[tid] = front_count;
                    
                    is_top_down = false;
                    is_top_down_queue = true;

                    if(VERBOSE)
                    {
                        if(tid==0)
                            std::cout<<"--->Switch to top down queue, "<<(wtime() - time_switch) * 1000<<"\n";
                    }
                }

//            #pragma omp barrier
            if(is_top_down || is_top_down_queue)
            //if(is_top_down_queue)
            {
                double queue_time = wtime();
                get_queue(thread_queue,
                        front_comm,
                        prefix_sum,
                        tid,
//                        temp_queue);
                        frontier_queue);
                
//                #pragma omp barrier
                if(tid == thread_count - 1)
                    queue_size = prefix_sum[thread_count - 1] + front_comm[thread_count - 1];
//                if(VERBOSE)
//                {
//                    if(tid == 0)
//                        printf("queue_size, %d, queue_time (ms), %.3lf\n", queue_size, (wtime() - queue_time) * 1000);
//                }
            }

            #pragma omp barrier
            level++;
//            if(tid == 0)
//            {
//                printf("root, %u, level, %d, front_count, %u, work_next, %u\n", root, level, front_count, my_work_next);
//            }
        }

        
    }
    printf("BFS kernel time (ms), %.3lf\n", (wtime() - time_bfs) * 1000);
//    
//    vertex_t colored_num = 0;
//    for(vertex_t v = 0; v < vert_count; ++v)
//    {
//        if(wcc_id[v] != INF)
//            colored_num ++;
//    }
//    printf("root, %u, colored_num, %u\n", root, colored_num);
//    
//    exit(EXIT_FAILURE);
}

void wcc_bfs_concurrent(graph *g, vertex_t thread_count, vertex_t *wcc_id, depth_t *sa)
{
    // Step 1: Initilization for BFS
    
    const vertex_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    
    vertex_t *frontier_queue = new vertex_t[vert_count + 1];
    vertex_t *temp_queue = new vertex_t[vert_count + 1];
//    vertex_t *max_pivot_list = new vertex_t[thread_count];
    vertex_t *front_comm=new vertex_t[thread_count];	
	vertex_t *work_comm=new vertex_t[thread_count];
//    depth_t *sa = new depth_t[vert_count + 1];

//    #pragma omp parallel for
//    for(vertex_t i = 0; i < vert_count + 1; ++i)
//    {
//        sa[i] = -1;
////        bw_sa[i] = -1;
////        scc_id[i] = 0;
//    }

//    vertex_t pivot = pivot_selection(g, thread_count);
    vertex_t pivot = 1; 
    bfs_concurrent(g,
        wcc_id,
        sa,
        frontier_queue,
        front_comm,
        work_comm,
        pivot,
        thread_count,
        temp_queue);
    
    // Step 2: Compute WCC with parallel BFS
//    for(vertex_t v = 0; v < vert_count; ++v)
//    {
//        if(wcc_id[v] == INF)
//        {
//        }
//    }
    
}


void wcc_bfs_parallel(graph *g, vertex_t thread_count, vertex_t *wcc_id, depth_t *sa)
{
    // Step 1: Initilization for BFS
    
    const vertex_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    
    vertex_t *frontier_queue = new vertex_t[vert_count + 1];
    vertex_t *temp_queue = new vertex_t[vert_count + 1];
//    vertex_t *max_pivot_list = new vertex_t[thread_count];
    vertex_t *front_comm=new vertex_t[thread_count];	
	vertex_t *work_comm=new vertex_t[thread_count];
//    depth_t *sa = new depth_t[vert_count + 1];

//    #pragma omp parallel for
//    for(vertex_t i = 0; i < vert_count + 1; ++i)
//    {
//        sa[i] = -1;
////        bw_sa[i] = -1;
////        scc_id[i] = 0;
//    }
    
    // Step 2: Compute WCC with parallel BFS
    for(vertex_t v = 0; v < vert_count; ++v)
    {
        if(wcc_id[v] == INF)
        {
            bfs_parallel(g,
                wcc_id,
                sa,
                frontier_queue,
                front_comm,
                work_comm,
                v,
                thread_count,
                temp_queue);
        }
    }
    
}
vertex_t get_remain_frontier(graph *g, vertex_t *wcc_id, vertex_t *frontier_queue, vertex_t thread_count)
{
    const vertex_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    const vertex_t upper_bound = vert_count / thread_count;

    vertex_t *thread_bin = new vertex_t[thread_count];
    vertex_t *prefix_sum = new vertex_t[thread_count];

    vertex_t step = vert_count / thread_count;
    if (vert_count != step * thread_count)
    {
        step += 1;
    }

    for(vertex_t i = 0; i < thread_count; ++i)
    {
        thread_bin[i] = 0;
        prefix_sum[i] = 0;
    }
    
    double queue_time = wtime();
    #pragma omp parallel \
    num_threads(thread_count) 
    {
        vertex_t tid = omp_get_thread_num();
        vertex_t vert_beg = tid * step;
        vertex_t vert_end = (tid == thread_count - 1 ? vert_count : vert_beg + step);
        vertex_t *thread_queue = new vertex_t[upper_bound];
        
        //step 1: iterate the vertex list to find the remaining
        vertex_t index = 0;
        for(vertex_t v = vert_beg; v < vert_end; ++v)
        {
            if(wcc_id[v] == INF)
            {
                thread_queue[thread_bin[tid]++] = v;
//                thread_bin[tid]++;
                wcc_id[v] = v;
            }
        }
        //step 2: get the prefix_sum[*]
        #pragma omp barrier
        for(vertex_t i = 0; i < tid; ++i)
        {
            prefix_sum[tid] += thread_bin[i];
        }
        //step 3: write the vertices into fq
        vertex_t start_pos = prefix_sum[tid];
        for(vertex_t vert_id = start_pos; vert_id < start_pos + thread_bin[tid]; ++vert_id)
        {
            frontier_queue[vert_id] = thread_queue[vert_id - start_pos];
        }
    }
    if(VERBOSE)
        printf("Queue kernel time (ms), %.3lf\n", (wtime() - queue_time) * 1000);
    return prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];
}

void color_propagation_fq(graph *g, vertex_t *wcc_id, vertex_t thread_count, vertex_t fq_size, vertex_t *frontier_queue)
{
    depth_t depth = 0;
    const vertex_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    const double avg_degree = edge_count * 1.0 / vert_count;

    index_t *fw_beg_pos = g->fw_beg_pos;
    vertex_t *fw_csr = g->fw_csr;
    index_t *bw_beg_pos = g->bw_beg_pos;
    vertex_t *bw_csr = g->bw_csr;

    bool *color_change = new bool[thread_count];
    memset(color_change, 0, sizeof(bool) * thread_count);
    
    vertex_t step = fq_size / thread_count;
    if (fq_size != step * thread_count)
    {
        step += 1;
    }
    
    double time_color = wtime();
    #pragma omp parallel \
    num_threads(thread_count) 
    {
        depth_t level = 0;
        const vertex_t tid = omp_get_thread_num();
        vertex_t vert_beg = tid * step;
        vertex_t vert_end = (tid == thread_count - 1 ? fq_size : vert_beg + step);
//        printf("vert_beg, %u, vert_end, %u\n", vert_beg, vert_end);
        while(true)
        {
//            #pragma omp barrier
//            if(DEBUG)
//            {
//                depth += 1;
//            }
            
            bool color_changed = false;
            //option 1: bottom up
//            for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
//            {
//                vertex_t vert_id = frontier_queue[fq_vert_id];
////                if(wcc_id[vert_id] == 0)
////                {
//
//                    //1. using out edge
//                    vertex_t my_beg = fw_beg_pos[vert_id];
//                    vertex_t my_end = fw_beg_pos[vert_id + 1];
//
//                    for(; my_beg < my_end; ++my_beg)
//                    {
//                        vertex_t w = fw_csr[my_beg];
////                        if(vert_id == w)
////                            continue;
//    //                    if(scc_id[w] == 0 && color[w] > color[vert_id])
//                        if(vert_id != w && wcc_id[vert_id] > wcc_id[w])
//                        {
//                            wcc_id[vert_id] = wcc_id[w];
//                            if(!color_changed)
//                                color_changed = true;
//                        }
//                    }
//                    //2. using in edge
//                    my_beg = bw_beg_pos[vert_id];
//                    my_end = bw_beg_pos[vert_id + 1];
//
//                    for(; my_beg < my_end; ++my_beg)
//                    {
//                        vertex_t w = bw_csr[my_beg];
////                        if(vert_id == w)
////                            continue;
//    //                    if(scc_id[w] == 0 && color[w] > color[vert_id])
//                        if(vert_id != w && wcc_id[vert_id] > wcc_id[w])
//                        {
//                            wcc_id[vert_id] = wcc_id[w];
//                            if(!color_changed)
//                                color_changed = true;
//                        }
//                    }
////                }
//            }
            //option 2: top-down
            //Large color --> small color
            for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
            {
                vertex_t vert_id = frontier_queue[fq_vert_id];
//                if(wcc_id[vert_id] == 0)
//                {

                    //1. using out edge
                    vertex_t my_beg = fw_beg_pos[vert_id];
                    vertex_t my_end = fw_beg_pos[vert_id + 1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        vertex_t w = fw_csr[my_beg];
//                        if(vert_id == w)
//                            continue;
    //                    if(scc_id[w] == 0 && color[w] > color[vert_id])
                        if(vert_id != w && wcc_id[w] > wcc_id[vert_id])
                        {
                            wcc_id[w] = wcc_id[vert_id];
                            if(!color_changed)
                                color_changed = true;
                        }
                    }
                    //2. using in edge
                    my_beg = bw_beg_pos[vert_id];
                    my_end = bw_beg_pos[vert_id + 1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        vertex_t w = bw_csr[my_beg];
//                        if(vert_id == w)
//                            continue;
    //                    if(scc_id[w] == 0 && color[w] > color[vert_id])
                        if(vert_id != w && wcc_id[w] > wcc_id[vert_id])
                        {
                            wcc_id[w] = wcc_id[vert_id];
                            if(!color_changed)
                                color_changed = true;
                        }
                    }
//                }
            }
            color_change[tid] = color_changed;
//            printf("tid, %u, color_changed, %d\n", tid, color_changed);
            #pragma omp barrier
            
//            if(tid == 0)
//                printf("color propagation\n");
// Path compression                
// Version 1: reduce to depth 1 
            if(color_changed)
            {
                for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
                {
                    vertex_t vert_id = frontier_queue[fq_vert_id];
                    if(wcc_id[vert_id] != vert_id)
                    {
                        vertex_t root = wcc_id[vert_id];
                        depth_t c_depth = 0;
                        while(wcc_id[root] != root && c_depth < 50)// && degree_prop[root] <= degree_prop[color[root]])
                        {
                            root = wcc_id[root];
                            c_depth ++;
                        }
                        vertex_t v_id = vert_id;
                        c_depth = 0;
                        while(v_id != root && wcc_id[v_id] != root)// && degree_prop[v_id] < degree_prop[root])
                        {
                            vertex_t prev = v_id;
                            v_id = wcc_id[v_id];
                            wcc_id[prev] = root;
                            c_depth ++;
                        }
                    }
                }
            }
//            if(tid == 0)
//                printf("path compression\n");
//            else
//            {
//                color_change[tid] = false;
//            }

            bool final_color_change = false;
            for(vertex_t i = 0; i < thread_count; ++i)
            {
//                final_color_change |= color_change[i];
                if(color_change[i])
                {
                    final_color_change = true;
                    break;
                }
            }
//            if(tid == 0)
//                printf("level = %d\n", level);
            if(final_color_change == false)
            {
                break;
            }
            #pragma omp barrier
//            level++;
//            #pragma omp barrier
        }
    }
    if(VERBOSE)
        printf("Color kernel time (ms), %.3lf\n", (wtime() - time_color) * 1000);

}

void wcc_hybrid(graph *g, vertex_t thread_count, vertex_t *wcc_id, depth_t *sa)
{

    // Step 1: Initilization for BFS
    
    const vertex_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    
    vertex_t *frontier_queue = new vertex_t[vert_count + 1];
    vertex_t *temp_queue = new vertex_t[vert_count + 1];
    vertex_t *front_comm=new vertex_t[thread_count];	
	vertex_t *work_comm=new vertex_t[thread_count];

    double pivot_time = wtime();    
//    #pragma omp parallel for
//    for(vertex_t i = 0; i < vert_count + 1; ++i)
//    {
//        sa[i] = -1;
//    }
    
    vertex_t pivot = pivot_selection(g, thread_count);
    // Step 2: Compute WCC with parallel BFS
    // Assume vertex 1 as root
    //vertex_t pivot = 1;
    
    double start_time = wtime();
// Version 1: Parallel BFS, single root
//    bfs_parallel(g,
//        wcc_id,
//        sa,
//        frontier_queue,
//        front_comm,
//        work_comm,
//        pivot,
//        thread_count,
//        temp_queue);
// Version 2: Concurrent BFS, thread_count root
    bfs_concurrent(g,
        wcc_id,
        sa,
        frontier_queue,
        front_comm,
        work_comm,
        pivot,
        thread_count,
        temp_queue);
    double bfs_time = wtime();

    // Step 3: Get the frontier queue with the remaining vertices
    // wcc_id is initialized to own color

    // Merge with bfs_parallel, probably will get better performance
    vertex_t fq_size = get_remain_frontier(g, 
                            wcc_id,
                            frontier_queue,
                            thread_count);

//    std::sort(frontier_queue, frontier_queue + fq_size);
    double queue_time = wtime();
//    #pragma omp parallel for
//    for(vertex_t i = 0; i < fq_size; ++i)
//    {
//        vertex_t v = frontier_queue[i];
//        wcc_id[v] = v;
//    }
//    double init_time = wtime();
    
//    printf("# of remain vertex, %u\n", fq_size);
    // Step 4: Compute WCC with color propagation
    if(fq_size != 0)
    {
        color_propagation_fq(g,
                wcc_id,
                thread_count,
                fq_size,
                frontier_queue);
    }

    double color_time = wtime();

    printf("Pviot time (ms), %.3lf\nBFS time (ms), %.3lf\nQueue time (ms), %.3lf\nColor time (ms), %.3lf\n", (start_time - pivot_time) * 1000, (bfs_time - start_time) * 1000, (queue_time - bfs_time) * 1000, (color_time - queue_time) * 1000);
//    printf("BFS time (ms), %.3lf\nQueue time (ms), %.3lf\nInit time (ms), %.3lf\nColor time (ms), %.3lf\n", (bfs_time - start_time) * 1000, (queue_time - bfs_time) * 1000, (init_time - queue_time) * 1000, (color_time - init_time) * 1000);
    delete[] frontier_queue;
    delete[] temp_queue;
    delete[] front_comm;	
	delete[] work_comm;
}

void get_wcc_result(vertex_t *wcc_id, vertex_t vert_count, char *result_file_wcc)
{
// Step 1: get wcc_stat<key, value>, key is root_id, value is wcc_size
    std::map<vertex_t, vertex_t> wcc_stat;
    for(vertex_t v = 0; v < vert_count; ++v)
    {
        vertex_t id = wcc_id[v];
        if(wcc_stat.find(id) != wcc_stat.end())
        {
            wcc_stat[id] += 1;
        }
        else
        {
            wcc_stat[id] = 1;
        }
    }

// Step 2: get wcc_size<key, value>, key is wcc_size, value is # of WCCs with that size
    vertex_t wcc_num = 0;
    std::map<vertex_t, vertex_t, std::greater<vertex_t> > wcc_size;
    std::map<vertex_t, vertex_t>::iterator it;
    for(it = wcc_stat.begin(); it != wcc_stat.end(); it++)
    {
        vertex_t cur_size = it->second;
        if(wcc_size.find(cur_size) != wcc_size.end())
        {
            wcc_size[cur_size] += 1;
        }
        else
        {
            wcc_size[cur_size] = 1;
        }
        wcc_num++;
    }
    
    FILE *fp = fopen(result_file_wcc, "w");
    
    fprintf(fp, "wcc_number,%u\n", wcc_num);
    vertex_t i = 0;
    for(it = wcc_size.begin(); it != wcc_size.end(); it++)
    {
        if(i == 0)
        {
            printf("wcc_number, %u, largest wcc size, %u\n", wcc_num, it->first);
            i++;
        }
//        printf("%u, %u\n", it->first, it->second);
        fprintf(fp, "%u,%u\n", it->first, it->second);
    }
    fclose(fp);

}


void color_propagation(
        graph *g,
        vertex_t *wcc_id,
        vertex_t thread_count
        )
{
    depth_t depth = 0;
    const vertex_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    const double avg_degree = edge_count * 1.0 / vert_count;

    index_t *fw_beg_pos = g->fw_beg_pos;
    vertex_t *fw_csr = g->fw_csr;
    index_t *bw_beg_pos = g->bw_beg_pos;
    vertex_t *bw_csr = g->bw_csr;

    bool *color_change = new bool[thread_count];
    memset(color_change, 0, sizeof(bool) * thread_count);
    
    vertex_t step = vert_count / thread_count;
    if (vert_count != step * thread_count)
    {
        step += 1;
    }

    #pragma omp parallel \
    num_threads(thread_count) 
    {
        depth_t level = 0;
        const vertex_t tid = omp_get_thread_num();
        vertex_t vert_beg = tid * step;
        vertex_t vert_end = (tid == thread_count - 1 ? vert_count : vert_beg + step);
//        printf("vert_beg, %u, vert_end, %u\n", vert_beg, vert_end);
        while(true)
        {
            #pragma omp barrier
            if(DEBUG)
            {
                depth += 1;
            }
            
            bool color_changed = false;
            //option 1: bottom up
            for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
            {
//                vertex_t vert_id = small_queue[fq_vert_id];
//                if(wcc_id[vert_id] == 0)
//                {

                    //1. fw: using out edge
                    index_t my_beg = fw_beg_pos[vert_id];
                    index_t my_end = fw_beg_pos[vert_id + 1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        vertex_t w = fw_csr[my_beg];
//                        if(vert_id == w)
//                            continue;
    //                    if(scc_id[w] == 0 && color[w] > color[vert_id])
                        if(vert_id != w)
                        {
                            if(wcc_id[vert_id] > wcc_id[w])
                            {
                                wcc_id[vert_id] = wcc_id[w];
                                if(!color_changed)
                                    color_changed = true;
                            }
                        }
                    }
                    //2. bw: using in edge
                    my_beg = bw_beg_pos[vert_id];
                    my_end = bw_beg_pos[vert_id + 1];

                    for(; my_beg < my_end; ++my_beg)
                    {
                        vertex_t w = bw_csr[my_beg];
//                        if(vert_id == w)
//                            continue;
    //                    if(scc_id[w] == 0 && color[w] > color[vert_id])
                        if(vert_id != w)
                        {
                            if(wcc_id[vert_id] > wcc_id[w])
                            {
                                wcc_id[vert_id] = wcc_id[w];
                                if(!color_changed)
                                    color_changed = true;
                            }
                        }
                    }
//                }
            }
            #pragma omp barrier
            
// Path compression                
// Version 1: reduce to depth 1 
            if(color_changed)
            {
                for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
                {
//                    vertex_t vert_id = small_queue[fq_vert_id];
                    if(wcc_id[vert_id] != vert_id)
                    {
                        vertex_t root = wcc_id[vert_id];
                        depth_t c_depth = 0;
                        while(wcc_id[root] != root && c_depth < 100)// && degree_prop[root] <= degree_prop[color[root]])
                        {
                            root = wcc_id[root];
                            c_depth ++;
                        }
                        index_t v_id = vert_id;
                        while(v_id != root && wcc_id[v_id] != root)// && degree_prop[v_id] < degree_prop[root])
                        {
                            vertex_t prev = v_id;
                            v_id = wcc_id[v_id];
                            wcc_id[prev] = root;
                        }
                    }
                }
                color_change[tid] = true;
            }
            else
            {
                color_change[tid] = false;
            }
            #pragma omp barrier

            bool final_color_change = false;
            for(index_t i=0; i<thread_count; ++i)
            {
                if(color_change[i])
                {
                    final_color_change = true;
                    break;
                }
            }
            if(tid == 0)
                printf("level = %d\n", level);
            #pragma omp barrier
            if(final_color_change == false)
            {
                break;
            }
            level++;
        }
    }

}

void wcc_color_propagation(graph *g, vertex_t thread_count, vertex_t *wcc_id, depth_t *sa)
{
    const vertex_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    
//    vertex_t *frontier_queue = new vertex_t[vert_count + 1];
//    vertex_t *temp_queue = new vertex_t[vert_count + 1];
//    vertex_t *front_comm=new vertex_t[thread_count];	
//	vertex_t *work_comm=new vertex_t[thread_count];
//    depth_t *sa = new depth_t[vert_count + 1];

    // Step 1: Initilization for color propagation
    #pragma omp parallel for
    for(vertex_t i = 0; i < vert_count + 1; ++i)
    {
        wcc_id[i] = i;
    }
    
    // Step 2: Compute WCC with color propagation
    color_propagation(g,
        wcc_id,
        thread_count);

}

//0 trim, 1 largest SCC, 2 small SCC, 3 total time
//4 trim_size_1, 5 trim_size_2, 6 pivot_selection, 7 fw_bfs, 8 bw_bfs, 9 color propagation, 10 color identify, 11 color_init 

void wcc_detection(
        graph *g,
        const index_t thread_count,
        vertex_t *wcc_id
        )
{
    const vertex_t vert_count = g->vert_count;
    const index_t edge_count = g->edge_count;
    const double avg_degree = edge_count * 1.0 / vert_count;
    printf("vert_count = %d, edge_count = %u, avg_degree = %.3lf\n", vert_count, edge_count, avg_degree);

    depth_t *sa = new depth_t[vert_count + 1];

    for(vertex_t v = 0; v < vert_count; ++v)
    {
        wcc_id[v] = INF;
        sa[v] = -1;
    }
    double start_time = wtime();
// Version 1: Sequential BFS
//    wcc_bfs_sequential(g, vert_count, wcc_id);
// Version 2: Parallel BFS
//    wcc_bfs_parallel(g, thread_count, wcc_id, sa);
// Version 3: Parallel Color Propagation
//    wcc_color_propagation(g, thread_count, wcc_id, sa);
// Version 4: Concurrent BFS of first level
//    wcc_bfs_concurrent(g, thread_count, wcc_id, sa);
// Version 5: Hybrid of BFS + Color
    wcc_hybrid(g, thread_count, wcc_id, sa);
    double end_time = wtime();
    printf("wcc computation time, %.2lf (ms)\n", (end_time - start_time) * 1000);
    
    delete[] sa;
}
