#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

using namespace std;

typedef unsigned int vertex_t;
typedef unsigned int index_t;

inline off_t fsize(const char *filename) {
    struct stat st;
    if (stat(filename, &st) == 0)
        return st.st_size;
    return -1;
}

int main(int argc, char** argv){
    if(argc < 5){
        cout << "Usage: " << argv[0] << " <begin.bin> <adjacent.bin> <head.bin> <out_degree.bin>" << endl;
        return -1;
    }

    // === Step 1: load existing binary files ===
    int fd_begin = open(argv[1], O_RDONLY);
    int fd_adj   = open(argv[2], O_RDONLY);
    int fd_head  = open(argv[3], O_RDONLY);
    int fd_deg   = open(argv[4], O_RDONLY);

    if(fd_begin < 0 || fd_adj < 0 || fd_head < 0 || fd_deg < 0){
        perror("Error opening input files");
        return -1;
    }

    size_t begin_size = fsize(argv[1]);
    size_t adj_size   = fsize(argv[2]);
    size_t head_size  = fsize(argv[3]);
    size_t deg_size   = fsize(argv[4]);

    index_t* begin = (index_t*)mmap(NULL, begin_size, PROT_READ, MAP_PRIVATE, fd_begin, 0);
    vertex_t* adj  = (vertex_t*)mmap(NULL, adj_size, PROT_READ, MAP_PRIVATE, fd_adj, 0);
    vertex_t* head = (vertex_t*)mmap(NULL, head_size, PROT_READ, MAP_PRIVATE, fd_head, 0);
    index_t* degree= (index_t*)mmap(NULL, deg_size, PROT_READ, MAP_PRIVATE, fd_deg, 0);

    size_t vert_count = (begin_size / sizeof(index_t)) - 1;
    size_t edge_count = (adj_size / sizeof(vertex_t));

    cout << "Original vert_count: " << vert_count << endl;
    cout << "Original edge_count: " << edge_count << endl;

    // === Step 2: prepare new undirected graph storage ===
    size_t new_edge_count = edge_count * 2; // add reverse edges

    int fd4 = open("undirected_adjacent.bin", O_CREAT|O_RDWR, 0666);
    ftruncate(fd4, new_edge_count * sizeof(vertex_t));
    vertex_t* new_adj = (vertex_t*)mmap(NULL, new_edge_count * sizeof(vertex_t),
                                        PROT_READ|PROT_WRITE, MAP_SHARED, fd4, 0);

    int fd5 = open("undirected_head.bin", O_CREAT|O_RDWR, 0666);
    ftruncate(fd5, new_edge_count * sizeof(vertex_t));
    vertex_t* new_head = (vertex_t*)mmap(NULL, new_edge_count * sizeof(vertex_t),
                                         PROT_READ|PROT_WRITE, MAP_SHARED, fd5, 0);

    int fd2 = open("undirected_out_degree.bin", O_CREAT|O_RDWR, 0666);
    ftruncate(fd2, vert_count * sizeof(index_t));
    index_t* new_degree = (index_t*)mmap(NULL, vert_count * sizeof(index_t),
                                         PROT_READ|PROT_WRITE, MAP_SHARED, fd2, 0);

    int fd3 = open("undirected_begin.bin", O_CREAT|O_RDWR, 0666);
    ftruncate(fd3, (vert_count+1) * sizeof(index_t));
    index_t* new_begin = (index_t*)mmap(NULL, (vert_count+1) * sizeof(index_t),
                                        PROT_READ|PROT_WRITE, MAP_SHARED, fd3, 0);

    // === Step 3: compute new degrees (undirected, so add both ways) ===
    for(size_t i=0; i<vert_count; i++){
        new_degree[i] = 0;
    }

    for(size_t i=0; i<edge_count; i++){
        vertex_t u = head[i];
        vertex_t v = adj[i];
        new_degree[u]++;
        new_degree[v]++; // reverse edge
    }

    // === Step 4: compute begin array ===
    new_begin[0] = 0;
    for(size_t i=1; i<=vert_count; i++){
        new_begin[i] = new_begin[i-1] + new_degree[i-1];
        new_degree[i-1] = 0;
    }

    // === Step 5: fill new adjacency list ===
    for(size_t i=0; i<edge_count; i++){
        vertex_t u = head[i];
        vertex_t v = adj[i];

        new_adj[new_begin[u] + new_degree[u]] = v;
        new_head[new_begin[u] + new_degree[u]] = u;
        new_degree[u]++;

        new_adj[new_begin[v] + new_degree[v]] = u;
        new_head[new_begin[v] + new_degree[v]] = v;
        new_degree[v]++;
    }

    cout << "New undirected edge_count: " << new_edge_count << endl;

    // === Step 6: print first few results ===
    for(size_t i=0; i<min((size_t)5, vert_count); i++){
        cout << "Vertex " << i << ": ";
        for(index_t j=0; j<new_degree[i]; j++){
            cout << new_head[new_begin[i]+j] << "-" << new_adj[new_begin[i]+j] << " ";
        }
        cout << endl;
    }

    // === cleanup ===
    munmap(begin, begin_size);
    munmap(adj, adj_size);
    munmap(head, head_size);
    munmap(degree, deg_size);

    munmap(new_adj, new_edge_count * sizeof(vertex_t));
    munmap(new_head, new_edge_count * sizeof(vertex_t));
    munmap(new_begin, (vert_count+1) * sizeof(index_t));
    munmap(new_degree, vert_count * sizeof(index_t));

    close(fd_begin);
    close(fd_adj);
    close(fd_head);
    close(fd_deg);
    close(fd2);
    close(fd3);
    close(fd4);
    close(fd5);

    return 0;
}
