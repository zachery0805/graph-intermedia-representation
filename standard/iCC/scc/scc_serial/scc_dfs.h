// A C++ program to find articulation povertex_ts in an undirected graph 
#ifndef __SCC_DFS_H__
#define __SCC_DFS_H__

#include <iostream> 
#include <list> 
#include <unistd.h>
//#include <stack>
#include "wtime.h"
#include "graph.h"
#define NIL -1 
using namespace std; 

// A recursive function that find articulation povertex_ts using DFS traversal 
// u --> The vertex to be visited next 
// visited[] --> keeps tract of visited vertices 
// disc[] --> Stores discovery times of visited vertices 
// parent[] --> Stores parent vertices in DFS tree 
// ap[] --> Store articulation povertex_ts 

//vertex_t time_index = 0;
vertex_t num_bicc = 0;


//stack<vertex_t> st;

void dfs_scc(graph_undirected *g, vertex_t u, bool visited[], vertex_t disc[], vertex_t low[], vertex_t scc_id[], bool ap[], vertex_t &time_index, vertex_t *track, vertex_t &index_track) 
{ 

    // Mark the current node as visited 
    visited[u] = true; 

    // Initialize discovery time and low value 
    disc[u] = low[u] = ++time_index; 

    // Go through all vertices aadjacent to this 
     
    for(index_t i = g->beg_pos[u]; i < g->beg_pos[u + 1]; ++i)
    {
        vertex_t v = g->csr[i]; // v is current adjacent of u 
        // If v is not visited yet, then make it a child of u 
        // in DFS tree and recur for it 
        if (disc[v] == NIL)
        { 
//            st.push_back(v);
            dfs_scc(g, v, visited, disc, low, scc_id, ap, time_index, track, index_track); 
            // Check if the subtree rooted with v has a connection to 
            // one of the ancestors of u 
            low[u] = min(low[u], low[v]); 
        }
        else if (visited[v] == true)// && disc[v] < low[u]) 
        {
            low[u] = min(low[u], disc[v]); 
        }
    } 
    if(low[u] == disc[u])
    {
        scc_id[u] = u;
        while((index_track > 0) && (track[index_track - 1] != u))
        {
            vertex_t w = track[index_track - 1];
            scc_id[w] = u;
            index_track--;
        }
    }
} 

// The function to do DFS traversal. It uses recursive function APUtil() 
void find_scc(graph_undirected *g) 
{
    vertex_t V = g->vert_count;
    // Mark all the vertices as not visited 
    bool *visited = new bool[V]; 
    vertex_t *disc = new vertex_t[V]; 
    vertex_t *low = new vertex_t[V]; 
    vertex_t *scc_id = new vertex_t[V]; 
    bool *ap = new bool[V]; // To store articulation povertex_ts 
    
    vertex_t *track = new vertex_t[V];

//    list<vertex_t> st; // = new list<vertex_t>[g->edge_count];

    // Initialize parent and visited, and ap(articulation povertex_t) arrays 
    for (vertex_t i = 0; i < V; i++) 
    { 
        scc_id[i] = NIL; 
        low[i] = NIL;
        disc[i] = NIL;
        visited[i] = false; 
        ap[i] = false; 
    } 
    
//    vertex_t num_bicc = 0;
    // Call the recursive helper function to find articulation povertex_ts 
    // in DFS tree rooted with vertex 'i' 
    double dfs_time = wtime();
    for (vertex_t i = 0; i < V; i++) 
    {
        vertex_t index_track = 0;
        vertex_t time_index = 0;
        if (disc[i] == NIL) 
        {
//            printf("%d\n", i);
//            vertex_t prev_bicc = num_bicc;
            dfs_scc(g, i, visited, disc, low, scc_id, ap, time_index, track, index_track); 
//            vertex_t j = 0;
        //    break;
        }
//        if(index_track > 0)
//        {
////                j = 1;
////                st.pop_back();
//            num_bicc ++;
//
//        }
//        if(prev_bicc == num_bicc)
//            num_bicc ++;
    }
    printf("DFS Time, %.3lf\n", (wtime() - dfs_time) * 1000);
    
    vertex_t num_ap = 0;
} 

#endif
