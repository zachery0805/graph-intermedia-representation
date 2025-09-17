// A C++ program to find articulation povertex_ts in an undirected graph 
#ifndef __BICC_BFS_H__
#define __BICC_BFS_H__

#include <iostream> 
#include <list> 
#include <unistd.h>
#include <queue>
#include <map>

#include "wtime.h"
#include "graph_undirected.h"
#define NIL -1 
using namespace std; 

vertex_t time_index = 0;

depth_t bfs_sequential(graph_undirected *g, 
        vertex_t root, 
        vertex_t *parent, 
        depth_t *level,
        vertex_t *node_type,
        vertex_t *lq_index,
        vertex_t *lq,
        vertex_t &lq_size)
{   
    std::queue<vertex_t> q; 
    q.push(root);
    parent[root] = root;

    vertex_t depth = 0;
    level[root] = 0;
    node_type[root] = 0;

//    vertex_t lq_size = 0;
    lq[lq_size++] = root;

    while(!q.empty())
    {
        //std::queue<vertex_t> q_temp;
        vertex_t queue_size = q.size();
        vertex_t index = 0; 

        ++depth;
        lq_index[depth] = lq_size;

        // Only inspect the vertices in current frontier
        while(index < queue_size)
        {
            vertex_t cur = q.front();
            vertex_t num_children = 0;
            ++index;
            q.pop();
            for(index_t i = g->beg_pos[cur]; i < g->beg_pos[cur + 1]; ++i)
            {       
                vertex_t w = g->csr[i];
                if(parent[w] == NIL)
                {
                    q.push(w);
                    parent[w] = cur;
                    level[w] = depth;
                    ++num_children;
                    lq[lq_size++] = w;
                }
            } 
            if(num_children == 0)
                node_type[cur] = 2;
            else
                node_type[cur] = 1;
        }
    }
    node_type[root] = 0;

    return depth;
} 

depth_t bfs_part(graph_undirected *g, 
        depth_t *level,
        vertex_t v,
        vertex_t p,
        bool *visited,
        vertex_t *track,
        vertex_t &index_track)
{
    std::queue<vertex_t> q; 
    q.push(v);

//    bool *visited = new bool[g->vert_count + 1];
//    for(vertex_t i = 0; i < g->vert_count; ++i)
//      k  visited[i] = false;
//    memset(visited, 0, sizeof(visited) * (g->vert_count + 1));

    visited[v] = true;
    visited[p] = true;
    
    track[index_track++] = v;
    track[index_track++] = p;

    while(!q.empty())
    {
        vertex_t cur = q.front();
        q.pop();
        for(index_t i = g->beg_pos[cur]; i < g->beg_pos[cur + 1]; ++i)
        {       
            vertex_t w = g->csr[i];
            if(!visited[w])
            {
                // Not AP
                if(level[w] < level[v])
                {
                    for(vertex_t j = 0; j < index_track; ++j)
                        visited[track[j]] = false;
                    
                    return level[w];
                }
                else
                {
                    q.push(w);
                    visited[w] = true;
                    track[index_track++] = w;
                }
            }
        } 
    }

//    for(vertex_t j = 0; j < index_track; ++j)
//        visited[track[j]] = false;

    return level[v];
} 


bool check_root(graph_undirected *g, 
        vertex_t root, 
        depth_t *level,
        vertex_t &num_total)
{
    
    vertex_t num_visited = 0;
    vertex_t root_degree = g->beg_pos[root + 1] - g->beg_pos[root];
    if(root_degree < 2)
    {
        num_visited = 1;
        return false;
    }

//    printf("%d, %d, %d, %d\n", root, root_degree, g->beg_pos[root], g->csr[g->beg_pos[root]]);

    vertex_t second_root = g->csr[g->beg_pos[root]];
//    printf("second_root, %d\n", second_root);

    std::queue<vertex_t> q; 
    q.push(second_root);
    num_visited = 1;
    num_total = 1;
// costly
    bool *visited = new bool[g->vert_count];
    memset(visited, 0, sizeof(bool) * g->vert_count);
//    for(vertex_t i = 0; i < g->vert_count; ++i)
//        visited[i] = false;

    visited[root] = true;
    visited[second_root] = true;
    num_total++;
//    vertex_t num_visited = 1;

    while(!q.empty())
    {
        vertex_t cur = q.front();
        q.pop();
        for(index_t i = g->beg_pos[cur]; i < g->beg_pos[cur + 1]; ++i)
        {       
            vertex_t w = g->csr[i];
            if(!visited[w])
            {
                q.push(w);
                visited[w] = true;
                num_total++;
                if(level[w] == 1)
                {
                    ++num_visited;
//                    printf("root = %d, num_visited = %d\n", root, num_visited);
//                    if(num_visited == root_degree)
//                    {
//                        return false;
//                    }
                }
            }
        } 
    }
    if(num_visited == root_degree)
    {
        return false;
    }
    return true;
}

void add_size(vertex_t cur_size, map<vertex_t, vertex_t, std::greater<vertex_t> > &size_map)
{
    if(size_map.find(cur_size) != size_map.end())
    {
        size_map[cur_size] += 1;
    }
    else
    {
        size_map[cur_size] = 1;
    }
}

vertex_t get_bicc(graph_undirected *g, 
        vertex_t root, 
        vertex_t *parent, 
        depth_t *level,
        vertex_t *node_type,
        bool *visited,
        bool *ap,
        vertex_t *track,
        bool *checked,
        vertex_t *low,
        vertex_t *par,
        vertex_t *lq_index,
        vertex_t *lq,
        vertex_t lq_size,
        vertex_t depth,
        map<vertex_t, vertex_t, std::greater<vertex_t> > &size_map,
        std::vector<std::vector<vertex_t>> &bicc_components)
{
    vertex_t num_bicc = 0;
    bool find_ap = false;
    // Step 1: check the non-root vertex filtering out leaf node 
    for(depth_t d = depth; d > 1; d--)
    {
        vertex_t l_begin = lq_index[d];
        vertex_t l_end = lq_index[d + 1];

        for(vertex_t i = l_begin; i < l_end; ++i)
        {
            vertex_t v = lq[i];
            //if(v == par[v])
    //        if(parent[v] == root)
    //            continue;
            if(!visited[v])
            {
                vertex_t p = parent[v];

//                if(v == p)
//                {
//                    continue;
//                }
                
                vertex_t index_track = 0;
                depth_t l = bfs_part(g, 
                                level, 
                                v, 
                                p,
                                visited,
                                track,
                                index_track);

                if(l > level[p])
                {
#ifdef DEBUG
                    printf("depth, %d, v, %d, p, %d\n", d, v, p);
#endif
                    add_size(index_track, size_map);
                    ap[p] = true;
                    find_ap = true;
                    // In case the traversed vertex is an AP
                    //for(vertex_t j = 0; j < index_track; ++j)
                    //{
                    //    vertex_t w = track[j];
                    //    if(ap[w])
                    //        visited[w] = false;
                    //}
    
                    visited[p] = false;
                    //visited[v] = true;
                    ++num_bicc;

                    // pusb back bicc componenet
                    std::vector<vertex_t> this_bicc(track, track + index_track);
                    bicc_components.push_back(this_bicc);

                    // printf("Track for BiCC starting at v=%d:\n", v);
                    // for (vertex_t i = 0; i < index_track; ++i) {
                    //     printf("%d ", track[i]);
                    // }
                    // printf("\n");

#ifdef DEBUG
                    for(vertex_t i = 0; i < index_track; ++i)
                    {
                        printf("%d ", track[i]);
                    }
                    printf("\n");
#endif
                }

            }

        }
    }
    

    // Step 2: check the root vertex
 //   vertex_t num_visited = 0;
 //   ap[root] = check_root(g,
 //                   root,
 //                   level,
 //                   num_visited);
 //   checked[root] = true;
//    printf("ap[root], %d, size, %d\n", ap[root], num_visited);

//    if(!ap[root])
    {
        /// Plus one is for root
//        add_size(num_visited, size_map);
//    else
        vertex_t my_beg = g->beg_pos[root];
        vertex_t my_end = g->beg_pos[root + 1];
        vertex_t num_bfs = 0;
        for(vertex_t i = my_beg; i < my_end; ++i)
        {
            vertex_t v = g->csr[i]; 
            if(!visited[v])
            {
                vertex_t index_track = 0;
                ++num_bfs;

                bfs_part(g, 
                        level, 
                        v, 
                        root,
                        visited,
                        track,
                        index_track);
                ++num_bicc;
                add_size(index_track, size_map);
                visited[root] = false;

                std::vector<vertex_t> this_bicc(track, track + index_track);
                bicc_components.push_back(this_bicc);

                // printf("2 Track for BiCC starting at v=%d:\n", v);
                // for (vertex_t i = 0; i < index_track; ++i) {
                //     printf("%d ", track[i]);
                // }
                // printf("\n");
#ifdef DEBUG
                    for(vertex_t i = 0; i < index_track; ++i)
                    {
                        printf("%d ", track[i]);
                    }
                    printf("\n");
#endif
            }
        }
        if(num_bfs > 1)
            ap[root] = true;
    }

//    printf("after, %d, %d\n", root, ap[root]);
//    printf("root, %d\n", root); 
    
    return num_bicc;
}


vertex_t bfs_bicc(graph_undirected *g, 
            vertex_t root, 
            bool *visited, 
            vertex_t *parent, 
            bool *ap,
            depth_t *level,
            vertex_t *node_type,
            vertex_t *track,
            bool *checked,
            vertex_t *low,
            vertex_t *par,
            vertex_t *lq_index,
            vertex_t *lq,
            std::map<vertex_t, vertex_t, std::greater<vertex_t> > &size_map,
            std::vector<std::vector<vertex_t>> &bicc_components)
{ 

    if(g->beg_pos[root + 1] == g->beg_pos[root])
    {
        add_size(1, size_map);
        return 1;
    }
#ifdef OUTPUT_TIME 
    double time_begin = wtime();
#endif
    vertex_t lq_size = 0;
// Step 1: Do BFS from the root vertex, record the parent, level, depth, and the vertices in each depth 
    depth_t depth = bfs_sequential(g, 
                            root, 
                            parent, 
                            level,
                            node_type,
                            lq_index,
                            lq,
                            lq_size);

#ifdef DEBUG
    printf("root, %d, lq_size, %d, depth, %d\n", root, lq_size, depth);
    
    for(vertex_t i = 0; i < depth + 1; ++i)
    {
        printf("%d ", lq_index[i]);
    }
    printf("\n");
    for(vertex_t i = 0; i < lq_size; ++i)
    {
        printf("%d ", lq[i]);
    }
    printf("\n");
#endif
    
#ifdef OUTPUT_TIME 
    double time_first = wtime();
    printf("BFS first, %.3lf\n", time_first - time_begin);
#endif

#ifdef DEBUG
// Debug, count the leaf and non-leaf node
    vertex_t non_leaf = 0;
    vertex_t leaf = 0;
    for(vertex_t i = 0; i < g->vert_count; ++i)
    {
        if(node_type[i] == 1)
            ++non_leaf;
        else
            if(node_type[i] == 2)
                ++leaf;
    }
    printf("non_leaf, %d, leaf, %d\n", non_leaf, leaf);
#endif
    
// Step 2: Check every vertex and see if it's AP
#ifdef OUTPUT_TIME 
    double time_check = wtime();
#endif
    vertex_t num_bicc = get_bicc(g,
                                root,
                                parent,
                                level,
                                node_type,
                                visited,
                                ap,
                                track,
                                checked,
                                low,
                                par,
                                lq_index,
                                lq,
                                lq_size,
                                depth,
                                size_map,
                                bicc_components);
#ifdef OUTPUT_TIME 
    double time_second = wtime();
    printf("BiCC check time, %.3lf\n", time_second - time_check);
#endif

    return num_bicc;
} 

void get_size_result(vertex_t num_bicc, map<vertex_t, vertex_t, std::greater<vertex_t> > &size_map, char *result_file)
{
    FILE *fp = fopen(result_file, "w");
    fprintf(fp, "bicc_number,%u\n", num_bicc);
    vertex_t i = 0;

    map<vertex_t, vertex_t>::iterator it;
    for(it = size_map.begin(); it != size_map.end(); it++)
    {
        if(i == 0)
        {
            printf("bicc_number, %u, largest bicc size, %u\n", num_bicc, it->first);
            i++;
        }
//        printf("%u, %u\n", it->first, it->second);
        fprintf(fp, "%u,%u\n", it->first, it->second);
    }
    fclose(fp);

}

void find_bicc_bfs(graph_undirected *g, bool *ap, std::vector<std::vector<vertex_t>>& bicc_components) 
{
    vertex_t V = g->vert_count;
    bool *visited = new bool[V]; 
    depth_t *level = new vertex_t[V]; 
    vertex_t *track = new vertex_t[V]; 
    vertex_t *parent = new vertex_t[V]; 
    vertex_t *low = new vertex_t[V];
    vertex_t *par = new vertex_t[V];
    bool *checked = new bool[V];
    vertex_t *node_type = new vertex_t[V];

// Level queue
    vertex_t *lq_index = new vertex_t[V];
    vertex_t *lq = new vertex_t[V];

// Statistical purpose
    map<vertex_t, vertex_t, std::greater<vertex_t> > size_map;
    // std::vector<std::vector<vertex_t>> bicc_components; 

    // bool *ap = new bool[V]; // To store articulation povertex_ts 

    // Initialize parent and visited, and ap(articulation povertex_t) arrays 
    for (vertex_t i = 0; i < V; i++) 
    { 
        parent[i] = NIL; 
        visited[i] = false; 
        ap[i] = false; 
        level[i] = NIL;
        node_type[i] = NIL;
        checked[i] = false;
        low[i] = i;
        par[i] = i;
        lq_index[i] = 0;
        lq[i] = NIL;
    } 

    // Call the recursive helper function to find articulation povertex_ts 
    // in DFS tree rooted with vertex 'i' 
    vertex_t num_bicc = 0;
    double time_begin = wtime();
    for (vertex_t i = 0; i < V; i++) 
    {
//        if (visited[i] == false) 
        if(node_type[i] == NIL) 
        {
            //bfs_ap(g, i, visited, disc, low, parent, ap); 
            vertex_t num = bfs_bicc (g, 
                                i, 
                                visited, 
                                parent, 
                                ap,
                                level,
                                node_type,
                                track,
                                checked,
                                low,
                                par,
                                lq_index,
                                lq,
                                size_map,
                                bicc_components);
            num_bicc += num;
        }
//        break;
    }
    double time_end = wtime();
    printf("Find BiCC time (s), %.3lf\n", (time_end - time_begin));
    
    printf("num_BiCC, %d\n", num_bicc);

// May change later
    const char *output_file = "result.log";
    get_size_result(num_bicc, size_map, const_cast<char*>(output_file));

    vertex_t num_ap = 0;
    // Now ap[] contains articulation povertex_ts, prvertex_t them 
    for (vertex_t i = 0; i < V; i++) 
    {
        if(ap[i] == true) 
        {
            num_ap ++;
            if(i < 100)
                cout << i << " "; 
            else
                if(i == 100)
                    cout << "\n" << num_ap;
        }
    }
    printf("\nnum_ap, %d\n", num_ap);

    // print bicc
    std::cout << "BICC results:" << std::endl;
     for (size_t i = 0; i < bicc_components.size(); ++i) {
         printf("BiCC %zu: ", i);
         for (vertex_t v : bicc_components[i]) {
             printf("%d ", v);
         }
         printf("\n");
   }
} 

#endif
