#ifndef __WCC_COMMON_H
#define __WCC_COMMON_H
#include <omp.h>
#include "graph.h"
#include "util.h"
#include "wtime.h"

#ifdef __cplusplus
extern "C" {
#endif

void wcc_detection(
        graph *g,
        const index_t thread_count,
        vertex_t *wcc_id
        );

#ifdef __cplusplus
}
#endif

void print_time_result(
        index_t run_times,
        double *avg_time
        );
    

#endif
