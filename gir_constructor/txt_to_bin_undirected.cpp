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

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_edge_file>" << endl;
        return 1;
    }

    int fd = open(argv[1], O_RDONLY);
    if (fd < 0) {
        perror("open input file failed");
        return 1;
    }

    size_t file_size = fsize(argv[1]);
    char* ss_head = (char*)mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
    if (ss_head == MAP_FAILED) {
        perror("mmap input file failed");
        return 1;
    }

    // Skip header lines
    size_t head_offset = 0;
    while (ss_head[head_offset] == '%') {
        while (ss_head[head_offset] != '\n') head_offset++;
        head_offset++;
    }
    char* ss = &ss_head[head_offset];
    file_size -= head_offset;

    size_t curr = 0, next = 0;
    size_t edge_count = 0;
    vertex_t v_max = 0, v_min = 0x7fffffff;
    vertex_t a;

    // Step 1: find max/min vertex ID and count edges
    while (next < file_size) {
        a = atoi(ss + curr);
        if (v_max < a) v_max = a;
        if (v_min > a) v_min = a;

        while (ss[next] != ' ' && ss[next] != '\n' && ss[next] != '\t') next++;
        while (ss[next] == ' ' || ss[next] == '\n' || ss[next] == '\t') next++;
        curr = next;

        a = atoi(ss + curr);
        if (v_max < a) v_max = a;
        if (v_min > a) v_min = a;

        while (ss[next] != ' ' && ss[next] != '\n' && ss[next] != '\t') next++;
        while (ss[next] == ' ' || ss[next] == '\n' || ss[next] == '\t') next++;
        curr = next;

        edge_count++;
    }

    vertex_t vert_count = v_max - v_min + 1;
    size_t total_edge_count = edge_count * 2;  // because it's undirected

    cout << "edge count (undirected): " << total_edge_count << endl;
    cout << "vertex count: " << vert_count << endl;

    // Create mmap output files
    int fd4 = open("undirected_adjacent.bin", O_CREAT | O_RDWR, 00666);
    int fd5 = open("undirected_head.bin", O_CREAT | O_RDWR, 00666);
    int fd2 = open("undirected_out_degree.bin", O_CREAT | O_RDWR, 00666);
    int fd3 = open("undirected_begin.bin", O_CREAT | O_RDWR, 00666);
    if (fd2 < 0 || fd3 < 0 || fd4 < 0 || fd5 < 0) {
        perror("open output file failed");
        return 1;
    }

    ftruncate(fd4, total_edge_count * sizeof(vertex_t));
    ftruncate(fd5, total_edge_count * sizeof(vertex_t));
    ftruncate(fd2, vert_count * sizeof(index_t));
    ftruncate(fd3, (vert_count + 1) * sizeof(index_t));

    vertex_t* adj = (vertex_t*)mmap(NULL, total_edge_count * sizeof(vertex_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd4, 0);
    vertex_t* head = (vertex_t*)mmap(NULL, total_edge_count * sizeof(vertex_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd5, 0);
    index_t* degree = (index_t*)mmap(NULL, vert_count * sizeof(index_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd2, 0);
    index_t* begin = (index_t*)mmap(NULL, (vert_count + 1) * sizeof(index_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd3, 0);

    // Step 3: count degrees
    memset(degree, 0, vert_count * sizeof(index_t));
    curr = next = 0;

    while (curr < file_size) {
        vertex_t u = atoi(ss + curr) - v_min;
        while (ss[next] != ' ' && ss[next] != '\n' && ss[next] != '\t') next++;
        while (ss[next] == ' ' || ss[next] == '\n' || ss[next] == '\t') next++;
        curr = next;

        vertex_t v = atoi(ss + curr) - v_min;
        while (ss[next] != ' ' && ss[next] != '\n' && ss[next] != '\t') next++;
        while (ss[next] == ' ' || ss[next] == '\n' || ss[next] == '\t') next++;
        curr = next;

        degree[u]++;
        degree[v]++;
    }

    // Step 4: compute begin[]
    begin[0] = 0;
    for (vertex_t i = 1; i <= vert_count; i++) {
        begin[i] = begin[i - 1] + degree[i - 1];
        degree[i - 1] = 0;  // reset for use as offset
    }

    // Step 5: fill adj/head
    curr = next = 0;
    while (curr < file_size) {
        vertex_t u = atoi(ss + curr) - v_min;
        while (ss[next] != ' ' && ss[next] != '\n' && ss[next] != '\t') next++;
        while (ss[next] == ' ' || ss[next] == '\n' || ss[next] == '\t') next++;
        curr = next;

        vertex_t v = atoi(ss + curr) - v_min;
        while (ss[next] != ' ' && ss[next] != '\n' && ss[next] != '\t') next++;
        while (ss[next] == ' ' || ss[next] == '\n' || ss[next] == '\t') next++;
        curr = next;

        adj[begin[u] + degree[u]] = v;
        head[begin[u] + degree[u]] = u;
        degree[u]++;

        adj[begin[v] + degree[v]] = u;
        head[begin[v] + degree[v]] = v;
        degree[v]++;
    }

    // Output debug info
    for (size_t i = 0; i < min((size_t)5, (size_t)vert_count); i++) {
        cout << "vertex " << i << ": ";
        for (index_t j = begin[i]; j < begin[i + 1]; j++) {
            cout << head[j] << "-" << adj[j] << " ";
        }
        cout << endl;
    }

    // Clean up
    munmap(ss_head, file_size + head_offset);
    munmap(adj, total_edge_count * sizeof(vertex_t));
    munmap(head, total_edge_count * sizeof(vertex_t));
    munmap(begin, (vert_count + 1) * sizeof(index_t));
    munmap(degree, vert_count * sizeof(index_t));

    close(fd);
    close(fd2);
    close(fd3);
    close(fd4);
    close(fd5);

    return 0;
}
