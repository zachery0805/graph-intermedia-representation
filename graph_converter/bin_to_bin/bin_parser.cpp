#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

using namespace std;

// 假设：与您现有程序一致，使用 unsigned int
typedef unsigned int vertex_t;
typedef unsigned int index_t;

// 请务必根据您的 grid_1000_100000_sym.bin 的实际内部结构调整以下常量！
const size_t VERTEX_COUNT = 1000 * 1000;  // V = 1,000,000 (来自文件名 grid_1000_100000)
const size_t EDGE_COUNT_GUESS = 100000 * 2; // E = 200,000 (如果是稀疏图，且是双向边，这只是一个猜测)

// !!! 警告：如果您不知道实际的 V 和 E，您无法继续。
// 在许多格式中，V 和 E 是文件最开始的两个 int/long 存储的。

int main(int argc, char** argv) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input_single_bin_file>" << endl;
        return -1;
    }

    const char* input_file = argv[1];

    // --- A. 假设读取 V 和 E ---
    // 在许多格式中，文件前几个字节存储 V 和 E
    ifstream fs(input_file, ios::binary);
    if (!fs) {
        perror("Error opening input file");
        return -1;
    }

    // 假设 V 和 E 是文件最开始的两个 index_t (4 字节无符号整数)
    index_t V_count, E_count;
    // 实际操作中，您应该知道 V 和 E 在哪里以及它们的类型。
    // 这里我们用硬编码的猜测代替，但在实际项目中这是错误的做法！
    // fs.read((char*)&V_count, sizeof(index_t));
    // fs.read((char*)&E_count, sizeof(index_t));
    
    // !!! 临时硬编码，如果您的文件结构不是这样，程序会崩溃或生成错误文件。
    V_count = VERTEX_COUNT;
    E_count = 20000000; // 假设一个更大的边数，比如 20M

    cout << "Parsing: V=" << V_count << ", E=" << E_count << endl;

    // --- B. 分配内存并读取数据 (按 CSR 结构猜测) ---
    // 假设文件内容顺序是：[BEGIN ARRAY] [ADJACENT ARRAY] [HEAD ARRAY] [DEGREE ARRAY]
    
    // 1. BEGIN ARRAY (V+1 个 index_t)
    size_t begin_bytes = (V_count + 1) * sizeof(index_t);
    vector<index_t> begin_vec(V_count + 1);
    fs.seekg(0, ios::beg);
    fs.read((char*)begin_vec.data(), begin_bytes);

    // 2. ADJACENT ARRAY (E 个 vertex_t)
    size_t adj_bytes = E_count * sizeof(vertex_t);
    vector<vertex_t> adj_vec(E_count);
    fs.seekg(begin_bytes, ios::beg);
    fs.read((char*)adj_vec.data(), adj_bytes);

    // 3. HEAD ARRAY (E 个 vertex_t)
    // 注意：许多 CSR 格式不包含 head 数组，这需要您在运行时生成。
    // 如果您的原始文件不包含 head 数组，请取消注释以下代码：
    /*
    size_t head_bytes = E_count * sizeof(vertex_t);
    vector<vertex_t> head_vec(E_count);
    fs.seekg(begin_bytes + adj_bytes, ios::beg);
    fs.read((char*)head_vec.data(), head_bytes);
    */
    
    // 4. DEGREE ARRAY (V 个 index_t)
    size_t degree_bytes = V_count * sizeof(index_t);
    vector<index_t> degree_vec(V_count);
    // 调整 seekg 位置，跳过 adj 和 head（如果存在）
    fs.seekg(begin_bytes + adj_bytes, ios::beg);
    fs.read((char*)degree_vec.data(), degree_bytes);
    
    fs.close();
    
    // --- C. 写入您的四个独立文件 ---
    ofstream out_begin("begin.bin", ios::binary);
    out_begin.write((char*)begin_vec.data(), begin_bytes);
    out_begin.close();

    ofstream out_adj("adjacent.bin", ios::binary);
    out_adj.write((char*)adj_vec.data(), adj_bytes);
    out_adj.close();
    
    // !!! 如果您的原始文件没有 head 数组，您需要从 begin 数组生成它 !!!
    // 此处假设您需要生成 head.bin
    vector<vertex_t> head_vec_generated(E_count);
    for(vertex_t u = 0; u < V_count; ++u) {
        for(index_t j = begin_vec[u]; j < begin_vec[u+1]; ++j) {
            head_vec_generated[j] = u;
        }
    }
    ofstream out_head("head.bin", ios::binary);
    out_head.write((char*)head_vec_generated.data(), E_count * sizeof(vertex_t));
    out_head.close();


    ofstream out_degree("out_degree.bin", ios::binary);
    out_degree.write((char*)degree_vec.data(), degree_bytes);
    out_degree.close();

    cout << "Successfully split input file into begin.bin, adjacent.bin, head.bin, out_degree.bin" << endl;
    return 0;
}