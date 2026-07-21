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

typedef int vertex_t;
typedef int index_t;

inline off_t fsize(const char *filename) {
    struct stat st; 
    if (stat(filename, &st) == 0)
        return st.st_size;
    return -1; 
}

int main(int argc, char** argv) {
    if(argc < 2) {
        cout << "Usage: " << argv[0] << " <fw_begin.bin>" << endl;
        return -1;
    }

    // 1. 读取顶点数量
    size_t file_size = fsize(argv[1]);
    if(file_size <= 0) {
        cerr << "Error reading file size." << endl;
        return -1;
    }

    int fd_begin = open(argv[1], O_RDONLY);
    if(fd_begin < 0) {
        perror("open fw_begin.bin");
        return -1;
    }

    index_t* begin = (index_t*) mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd_begin, 0);
    if(begin == MAP_FAILED) {
        perror("mmap fw_begin.bin");
        return -1;
    }

    size_t vert_count = file_size / sizeof(index_t) - 1; // begin数组长度为vert_count + 1
    cout << "Vertex count: " << vert_count << endl;

    // 2. 根据begin数组获取边数量
    size_t edge_count = begin[vert_count];
    cout << "Edge count: " << edge_count << endl;

    // 3. mmap 读取 fw_adjacent.bin 和 head.bin
    int fd_adj = open("fw_adjacent.bin", O_RDONLY);
    int fd_head = open("head.bin", O_RDONLY);
    if(fd_adj < 0 || fd_head < 0) {
        perror("open fw_adjacent.bin or head.bin");
        return -1;
    }

    vertex_t* adj = (vertex_t*) mmap(NULL, edge_count * sizeof(vertex_t), PROT_READ, MAP_PRIVATE, fd_adj, 0);
    vertex_t* head = (vertex_t*) mmap(NULL, edge_count * sizeof(vertex_t), PROT_READ, MAP_PRIVATE, fd_head, 0);
    if(adj == MAP_FAILED || head == MAP_FAILED) {
        perror("mmap adj or head");
        return -1;
    }

    // 4. mmap 读取 out_degree.bin
    int fd_degree = open("out_degree.bin", O_RDONLY);
    if(fd_degree < 0) {
        perror("open out_degree.bin");
        return -1;
    }
    index_t* degree = (index_t*) mmap(NULL, vert_count * sizeof(index_t), PROT_READ, MAP_PRIVATE, fd_degree, 0);
    if(degree == MAP_FAILED) {
        perror("mmap degree");
        return -1;
    }

    // 5. 输出前几个测试
    cout << "Printing first 8 vertices (forward edges):" << endl;
    for(size_t i = 0; i < vert_count && i < 8; i++) {
        cout << "Vertex " << i << ": begin=" << begin[i] << " degree=" << degree[i] << " neighbors: ";
        for(index_t j = 0; j < degree[i]; j++) {
            cout << adj[begin[i] + j] << "(" << head[begin[i] + j] << ") ";
        }
        cout << endl;
    }

    // 6. 可选：写入新的二进制文件 (bin_to_bin_copy_fw)
    int fd_adj_out = open("fw_adjacent_copy.bin", O_CREAT|O_RDWR, 0666);
    ftruncate(fd_adj_out, edge_count * sizeof(vertex_t));
    vertex_t* adj_out = (vertex_t*) mmap(NULL, edge_count * sizeof(vertex_t), PROT_READ|PROT_WRITE, MAP_SHARED, fd_adj_out, 0);
    memcpy(adj_out, adj, edge_count * sizeof(vertex_t));

    int fd_head_out = open("head_copy.bin", O_CREAT|O_RDWR, 0666);
    ftruncate(fd_head_out, edge_count * sizeof(vertex_t));
    vertex_t* head_out = (vertex_t*) mmap(NULL, edge_count * sizeof(vertex_t), PROT_READ|PROT_WRITE, MAP_SHARED, fd_head_out, 0);
    memcpy(head_out, head, edge_count * sizeof(vertex_t));

    int fd_begin_out = open("fw_begin_copy.bin", O_CREAT|O_RDWR, 0666);
    ftruncate(fd_begin_out, (vert_count+1) * sizeof(index_t));
    index_t* begin_out = (index_t*) mmap(NULL, (vert_count+1) * sizeof(index_t), PROT_READ|PROT_WRITE, MAP_SHARED, fd_begin_out, 0);
    memcpy(begin_out, begin, (vert_count+1) * sizeof(index_t));

    int fd_degree_out = open("out_degree_copy.bin", O_CREAT|O_RDWR, 0666);
    ftruncate(fd_degree_out, vert_count * sizeof(index_t));
    index_t* degree_out = (index_t*) mmap(NULL, vert_count * sizeof(index_t), PROT_READ|PROT_WRITE, MAP_SHARED, fd_degree_out, 0);
    memcpy(degree_out, degree, vert_count * sizeof(index_t));

    // 7. cleanup
    munmap(begin, file_size);
    munmap(adj, edge_count * sizeof(vertex_t));
    munmap(head, edge_count * sizeof(vertex_t));
    munmap(degree, vert_count * sizeof(index_t));
    munmap(adj_out, edge_count * sizeof(vertex_t));
    munmap(head_out, edge_count * sizeof(vertex_t));
    munmap(begin_out, (vert_count+1) * sizeof(index_t));
    munmap(degree_out, vert_count * sizeof(index_t));

    close(fd_begin); close(fd_adj); close(fd_head); close(fd_degree);
    close(fd_adj_out); close(fd_head_out); close(fd_begin_out); close(fd_degree_out);

    cout << "Forward-edge copy completed successfully." << endl;

    return 0;
}
