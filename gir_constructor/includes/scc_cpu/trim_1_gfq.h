#ifndef TRIM_1_GFQ_H
#define TRIM_1_GFQ_H

#include "scc_common.h"
#include "util.h"
#include "wtime.h"
#include <omp.h>
#include <atomic> // Include the atomic header for thread-safe operations

// Macro for debug logging
#define DEBUG_LOG 0

inline void log_assign_scc(const char* caller, index_t scc_id, vertex_t vert_id, const char* reason)
{
#if DEBUG_LOG
    printf("[%s] Assigned scc_id %zu to vertex %zu (%s)\n", caller, scc_id, vert_id, reason);
#endif
}

// ------------------- trim_1_first -------------------
inline void trim_1_first(
        std::atomic<index_t>& current_scc_label,
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        const char* caller = "trim_1_first"
        )
{
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] != 0) continue;
        
        if(fw_beg_pos[vert_id+1] - fw_beg_pos[vert_id] == 0)
        {
            index_t id = current_scc_label.fetch_add(1, std::memory_order_relaxed);
            scc_id[vert_id] = id;
            log_assign_scc(caller, id, vert_id, "fw deg 0");
            continue;
        }
        else if(bw_beg_pos[vert_id+1] - bw_beg_pos[vert_id] == 0)
        {
            index_t id = current_scc_label.fetch_add(1, std::memory_order_relaxed);
            scc_id[vert_id] = id;
            log_assign_scc(caller, id, vert_id, "bw deg 0");
            continue;
        }
    }
}

// ------------------- trim_1_first_gfq -------------------
inline void trim_1_first_gfq(
        std::atomic<index_t>& current_scc_label,
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t tid,
        const char* caller = "trim_1_first_gfq"
        )
{
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] != 0) continue;

        if(fw_beg_pos[vert_id+1] - fw_beg_pos[vert_id] == 0)
        {
            index_t id = current_scc_label.fetch_add(1, std::memory_order_relaxed);
            scc_id[vert_id] = id;
            log_assign_scc(caller, id, vert_id, "fw deg 0");
            continue;
        }
        else if(bw_beg_pos[vert_id+1] - bw_beg_pos[vert_id] == 0)
        {
            index_t id = current_scc_label.fetch_add(1, std::memory_order_relaxed);
            scc_id[vert_id] = id;
            log_assign_scc(caller, id, vert_id, "bw deg 0");
            continue;
        }
        else
        {
            thread_bin_size++;
        }
    }

    thread_bin[tid] = thread_bin_size;

    // prefix sum
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
        prefix_sum[tid] += thread_bin[i];
    #pragma omp barrier

    // write to frontier queue
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
        if(scc_id[vert_id] == 0)
            frontier_queue[start_pos++] = vert_id;
}

// ------------------- trim_1_normal -------------------
inline void trim_1_normal(
        std::atomic<index_t>& current_scc_label,
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        const char* caller = "trim_1_normal"
        )
{
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] != 0) continue;

        index_t my_beg = fw_beg_pos[vert_id];
        index_t my_end = fw_beg_pos[vert_id+1];
        bool out_zero = true;
        for(; my_beg < my_end; ++my_beg)
        {
            index_t w = fw_csr[my_beg];
            if(scc_id[w] == 0 && w != vert_id)
            {
                out_zero = false;
                break;
            }
        }
        if(out_zero)
        {
            index_t id = current_scc_label.fetch_add(1, std::memory_order_relaxed);
            scc_id[vert_id] = id;
            log_assign_scc(caller, id, vert_id, "fw out-degree 0");
            continue;
        }

        my_beg = bw_beg_pos[vert_id];
        my_end = bw_beg_pos[vert_id+1];
        bool in_zero = true;
        for(; my_beg < my_end; ++my_beg)
        {
            index_t w = bw_csr[my_beg];
            if(scc_id[w] == 0 && w != vert_id)
            {
                in_zero = false;
                break;
            }
        }
        if(in_zero)
        {
            index_t id = current_scc_label.fetch_add(1, std::memory_order_relaxed);
            scc_id[vert_id] = id;
            log_assign_scc(caller, id, vert_id, "bw in-degree 0");
        }
    }
}

// ------------------- trim_1_from_fq -------------------
inline void trim_1_from_fq(
        std::atomic<index_t>& current_scc_label,
        index_t *scc_id,
        index_t *fw_beg_pos,
        index_t *bw_beg_pos,
        index_t vert_beg,
        index_t vert_end,
        vertex_t *fw_csr,
        vertex_t *bw_csr,
        index_t *small_queue,
        const char* caller = "trim_1_from_fq"
        )
{
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] != 0) continue;

        index_t my_beg = fw_beg_pos[vert_id];
        index_t my_end = fw_beg_pos[vert_id+1];
        bool out_zero = true;
        for(; my_beg < my_end; ++my_beg)
        {
            index_t w = fw_csr[my_beg];
            if(scc_id[w] == 0 && w != vert_id)
            {
                out_zero = false;
                break;
            }
        }
        if(out_zero)
        {
            index_t id = current_scc_label.fetch_add(1, std::memory_order_relaxed);
            scc_id[vert_id] = id;
            log_assign_scc(caller, id, vert_id, "fw out-degree 0");
            continue;
        }

        my_beg = bw_beg_pos[vert_id];
        my_end = bw_beg_pos[vert_id+1];
        bool in_zero = true;
        for(; my_beg < my_end; ++my_beg)
        {
            index_t w = bw_csr[my_beg];
            if(scc_id[w] == 0 && w != vert_id)
            {
                in_zero = false;
                break;
            }
        }
        if(in_zero)
        {
            index_t id = current_scc_label.fetch_add(1, std::memory_order_relaxed);
            scc_id[vert_id] = id;
            log_assign_scc(caller, id, vert_id, "bw in-degree 0");
        }
    }
}

// ------------------- trim_1_normal_only_size -------------------

// return the number of trimmed vertices
inline void trim_1_normal_only_size(
std::atomic<index_t>& current_scc_label,
index_t *scc_id,
index_t *fw_beg_pos,
index_t *bw_beg_pos,
index_t vert_beg,
index_t vert_end,
vertex_t *fw_csr,
vertex_t *bw_csr,
const index_t thread_count,
index_t *thread_bin,
index_t *prefix_sum,
index_t tid
)
{
index_t thread_bin_size = 0;
for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
{
if(scc_id[vert_id] != 0) continue; // ✅ Skip assigned vertices


index_t my_beg = fw_beg_pos[vert_id];
index_t my_end = fw_beg_pos[vert_id+1];
index_t out_degree = 0;
for(; my_beg < my_end; ++my_beg)
{
index_t w = fw_csr[my_beg];
if(scc_id[w] == 0 && w != vert_id)
{
out_degree = 1;
break;
}
}
if(out_degree == 0)
{
scc_id[vert_id] = current_scc_label.fetch_add(1, std::memory_order_relaxed);
thread_bin_size++;
continue;
}


index_t in_degree = 0;
my_beg = bw_beg_pos[vert_id];
my_end = bw_beg_pos[vert_id+1];
for(; my_beg < my_end; ++my_beg)
{
index_t w = bw_csr[my_beg];
if(scc_id[w] == 0 && w != vert_id)
{
in_degree = 1;
break;
}
}
if(in_degree == 0)
{
scc_id[vert_id] = current_scc_label.fetch_add(1, std::memory_order_relaxed);
thread_bin_size++;
}
}


thread_bin[tid] = thread_bin_size;


prefix_sum[tid] = 0;
#pragma omp barrier
for(index_t i=0; i<tid; ++i)
prefix_sum[tid] += thread_bin[i];
}

inline void trim_1_normal_gfq(
std::atomic<index_t>& current_scc_label,
index_t *scc_id,
index_t *fw_beg_pos,
index_t *bw_beg_pos,
index_t vert_beg,
index_t vert_end,
vertex_t *fw_csr,
vertex_t *bw_csr,
const index_t thread_count,
index_t *frontier_queue,
index_t *thread_bin,
index_t *prefix_sum,
index_t tid
)
{
index_t thread_bin_size = 0;
for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
{
if(scc_id[vert_id] != 0) continue;


index_t my_beg = fw_beg_pos[vert_id];
index_t my_end = fw_beg_pos[vert_id+1];
index_t out_degree = 0;
for(; my_beg < my_end; ++my_beg)
{
index_t w = fw_csr[my_beg];
if(scc_id[w] == 0 && w != vert_id)
{
out_degree = 1;
break;
}
}
if(out_degree == 0)
{
scc_id[vert_id] = current_scc_label.fetch_add(1, std::memory_order_relaxed);
continue;
}


index_t in_degree = 0;
my_beg = bw_beg_pos[vert_id];
my_end = bw_beg_pos[vert_id+1];
for(; my_beg < my_end; ++my_beg)
{
index_t w = bw_csr[my_beg];
if(scc_id[w] == 0 && w != vert_id)
{
in_degree = 1;
break;
}
}
if(in_degree == 0)
{
scc_id[vert_id] = current_scc_label.fetch_add(1, std::memory_order_relaxed);
continue;
}
thread_bin_size++;
}


thread_bin[tid] = thread_bin_size;


prefix_sum[tid] = 0;
#pragma omp barrier
for(index_t i=0; i<tid; ++i)
prefix_sum[tid] += thread_bin[i];


vertex_t start_pos = prefix_sum[tid];
for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
{
if(scc_id[vert_id] == 0)
frontier_queue[start_pos++] = vert_id;
}
}

inline void trim_1_from_fq_gfq(
std::atomic<index_t>& current_scc_label, // 🆕 Pass the atomic counter
index_t *scc_id,
index_t *fw_beg_pos,
index_t *bw_beg_pos,
index_t vert_beg,
index_t vert_end,
vertex_t *fw_csr,
vertex_t *bw_csr,
const index_t thread_count,
index_t *frontier_queue,
index_t *thread_bin,
index_t *prefix_sum,
index_t tid,
index_t *temp_queue
)
{
// step 1: get thread_bin_size
index_t thread_bin_size = 0;
for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
{
vertex_t vert_id = frontier_queue[fq_vert_id];
if(scc_id[vert_id] != 0) continue; // skip already assigned SCC


index_t my_beg = fw_beg_pos[vert_id];
index_t my_end = fw_beg_pos[vert_id+1];
index_t out_degree = 0;
for(; my_beg < my_end; ++my_beg)
{
index_t w = fw_csr[my_beg];
if(scc_id[w] == 0 && w != vert_id)
{
out_degree = 1;
break;
}
}
if(out_degree == 0)
{
scc_id[vert_id] = current_scc_label.fetch_add(1, std::memory_order_relaxed);
continue;
}


my_beg = bw_beg_pos[vert_id];
my_end = bw_beg_pos[vert_id+1];
index_t in_degree = 0;
for(; my_beg < my_end; ++my_beg)
{
index_t w = bw_csr[my_beg];
if(scc_id[w] == 0 && w != vert_id)
{
in_degree = 1;
break;
}
}
if(in_degree == 0)
{
scc_id[vert_id] = current_scc_label.fetch_add(1, std::memory_order_relaxed);
continue;
}
thread_bin_size++;
}
thread_bin[tid] = thread_bin_size;


// step 2: get prefix_sum for each thread
prefix_sum[tid] = 0;
#pragma omp barrier
for(index_t i = 0; i < tid; ++i)
{
prefix_sum[tid] += thread_bin[i];
}


// step 3: write the vertices into temp_queue
vertex_t start_pos = prefix_sum[tid];
for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
{
vertex_t vert_id = frontier_queue[fq_vert_id];
if(scc_id[vert_id] == 0)
{
temp_queue[start_pos++] = vert_id;
}
}


#pragma omp barrier
if(DEBUG && tid == 0)
{
printf("In normal trim, thread bin size: %d\n", prefix_sum[thread_count-1] + thread_bin[thread_count-1]);
}


// step 4: write back to frontier_queue
for(index_t i = prefix_sum[tid]; i < prefix_sum[tid] + thread_bin[tid]; ++i)
{
frontier_queue[i] = temp_queue[i];
}

}

inline static void get_queue(
        vertex_t *thread_queue,
        vertex_t *thread_bin,
        index_t *prefix_sum,
        index_t tid,
        vertex_t *temp_queue
        )
{
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }
    #pragma omp barrier
    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = start_pos; vert_id < start_pos + thread_bin[tid]; ++vert_id)
    {
        temp_queue[vert_id] = thread_queue[vert_id-start_pos];
    }
}

// Using prefix sum to generate frontier queue 
inline static void generate_frontier_queue(
        const index_t vert_count,
        index_t *scc_id,
        const index_t thread_count,
        index_t *frontier_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t vert_id = vert_beg; vert_id < vert_end; ++vert_id)
    {
        if(scc_id[vert_id] == 0)
        {
            frontier_queue[start_pos++] = vert_id;
        }
    }
//    #pragma omp barrier
//    return prefix_sum[thread_count - 1] + thread_bin[thread_count - 1];

}

inline static void gfq_from_queue(
        const index_t vert_count,
        index_t *scc_id,
        const index_t thread_count,
        index_t *small_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(scc_id[vert_id] == 0)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
    #pragma omp barrier
    //step 4: write back to small_queue
    for(index_t i=prefix_sum[tid]; i<prefix_sum[tid]+thread_bin[tid]; ++i)
    {
        small_queue[i] = temp_queue[i];
    }

}

inline static void bw_gfq_from_fw(
        index_t *fw_sa,
        const index_t thread_count,
        index_t *small_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(fw_sa[vert_id] != -1)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    #pragma omp barrier
    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(fw_sa[vert_id] != -1)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
//    #pragma omp barrier
//    //step 4: write back to small_queue
//    for(index_t i=prefix_sum[tid]; i<prefix_sum[tid]+thread_bin[tid]; ++i)
//    {
//        small_queue[i] = temp_queue[i];
//    }

}

inline static void gfq_fw_bw_from_queue(
        index_t *sa,
        const index_t thread_count,
        index_t *small_queue,
        index_t *thread_bin,
        index_t *prefix_sum,
        index_t vert_beg,
        index_t vert_end,
        index_t tid,
        index_t *temp_queue
        )
{
    //step 1: get thread_bin_size
    index_t thread_bin_size = 0;
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(sa[vert_id] == -1)
        {
            thread_bin_size ++;
        }
    }
    thread_bin[tid] = thread_bin_size;
    
    //step 2: get prefix_sum for each thread
    prefix_sum[tid] = 0;
    #pragma omp barrier
    for(index_t i=0; i<tid; ++i)
    {
        prefix_sum[tid] += thread_bin[i];
    }

    //step 3: write the vertices into fq
    vertex_t start_pos = prefix_sum[tid];
    for(vertex_t fq_vert_id = vert_beg; fq_vert_id < vert_end; ++fq_vert_id)
    {
        vertex_t vert_id = small_queue[fq_vert_id];
        if(sa[vert_id] == -1)
        {
            temp_queue[start_pos++] = vert_id;
        }
    }
}

#endif


