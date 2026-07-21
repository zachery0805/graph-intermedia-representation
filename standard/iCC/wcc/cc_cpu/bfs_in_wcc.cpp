#include <omp.h>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <cstdio>
#include <cstring>

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

        // Use first vertex as the root
        vertex_t start = vertices[0];
        q.push(start);
        visited.insert(start);

        printf("[Thread %d] BFS on WCC ID %u, size = %lu, start = %u\n",
               omp_get_thread_num(), wcc_root, vertices.size(), start);

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
