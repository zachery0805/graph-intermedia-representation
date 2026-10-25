// SNAP edge list (text, '#' comments, "u v" per line, optionally gzip'd by
// the caller via a pipe) -> GIR's binary CSR: fw_begin.bin / fw_adjacent.bin
// and bw_begin.bin / bw_adjacent.bin (4-byte ints; begin has V+1 entries).
// Vertex ids are remapped to 0..V-1 in increasing order of the original id;
// self loops and duplicate edges are dropped. Edges are taken as directed,
// in the direction listed.
// Usage: el_to_csr <edges.txt | -> <out dir>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <parallel/algorithm>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <omp.h>

static void writeArr(const std::string& path, const std::vector<int>& a) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) { std::perror(path.c_str()); std::exit(1); }
    std::fwrite(a.data(), sizeof(int), a.size(), f);
    std::fclose(f);
}

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: %s <edges.txt> <out dir>\n", argv[0]); return 2; }
    int fd = open(argv[1], O_RDONLY);
    if (fd < 0) { std::perror(argv[1]); return 1; }
    struct stat st; fstat(fd, &st);
    const size_t n = (size_t)st.st_size;
    const char* buf = (const char*)mmap(nullptr, n, PROT_READ, MAP_PRIVATE, fd, 0);
    // Parse in parallel chunks aligned to line starts.
    const int T = omp_get_max_threads();
    std::vector<size_t> cut(T + 1);
    cut[0] = 0; cut[T] = n;
    for (int t = 1; t < T; ++t) {
        size_t p = n * t / T;
        while (p < n && buf[p - 1] != '\n') ++p;
        cut[t] = p;
    }
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> part(T);
    #pragma omp parallel num_threads(T)
    {
        const int t = omp_get_thread_num();
        auto& out = part[t];
        const char* p = buf + cut[t];
        const char* e = buf + cut[t + 1];
        while (p < e) {
            if (*p == '#' || *p == '%') { while (p < e && *p != '\n') ++p; ++p; continue; }
            uint64_t x[2]; int k = 0;
            while (p < e && *p != '\n' && k < 2) {
                while (p < e && (*p == ' ' || *p == '\t' || *p == '\r')) ++p;
                if (p >= e || *p == '\n') break;
                uint64_t v = 0;
                while (p < e && *p >= '0' && *p <= '9') { v = v * 10 + (*p - '0'); ++p; }
                x[k++] = v;
            }
            while (p < e && *p != '\n') ++p;
            ++p;
            if (k == 2) out.push_back({x[0], x[1]});
        }
    }
    size_t m = 0;
    for (auto& v : part) m += v.size();
    std::vector<std::pair<uint64_t, uint64_t>> es;
    es.reserve(m);
    for (auto& v : part) { es.insert(es.end(), v.begin(), v.end()); std::vector<std::pair<uint64_t, uint64_t>>().swap(v); }
    std::fprintf(stderr, "read %zu edge lines\n", es.size());
    // Dense ids.
    std::vector<uint64_t> ids(2 * es.size());
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < es.size(); ++i) { ids[2 * i] = es[i].first; ids[2 * i + 1] = es[i].second; }
    __gnu_parallel::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    const size_t V = ids.size();
    std::vector<std::pair<int, int>> E(es.size());
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < es.size(); ++i)
        E[i] = {(int)(std::lower_bound(ids.begin(), ids.end(), es[i].first) - ids.begin()),
                (int)(std::lower_bound(ids.begin(), ids.end(), es[i].second) - ids.begin())};
    std::vector<std::pair<uint64_t, uint64_t>>().swap(es);
    if (std::getenv("EL_RANDOM_ORIENT")) {
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < E.size(); ++i) {
            uint64_t a = (uint64_t)std::min(E[i].first, E[i].second), b = (uint64_t)std::max(E[i].first, E[i].second);
            uint64_t x = a << 32 | b;
            x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
            E[i] = (x & 1) ? std::make_pair((int)a, (int)b) : std::make_pair((int)b, (int)a);
        }
        std::fprintf(stderr, "edges oriented by hash of their endpoints\n");
    }
    E.erase(std::remove_if(E.begin(), E.end(), [](const std::pair<int, int>& e) { return e.first == e.second; }), E.end());
    __gnu_parallel::sort(E.begin(), E.end());
    E.erase(std::unique(E.begin(), E.end()), E.end());
    if (E.size() >= (size_t)INT32_MAX) { std::fprintf(stderr, "too many edges for 32-bit offsets: %zu\n", E.size()); return 1; }
    std::fprintf(stderr, "V=%zu E=%zu (self loops and duplicates removed)\n", V, E.size());
    auto build = [&](bool rev, const char* bname, const char* aname) {
        std::vector<int> beg(V + 1, 0), adj(E.size());
        for (auto& e : E) beg[(rev ? e.second : e.first) + 1]++;
        for (size_t v = 0; v < V; ++v) beg[v + 1] += beg[v];
        std::vector<int> pos(beg.begin(), beg.end() - 1);
        for (auto& e : E) { int s = rev ? e.second : e.first, d = rev ? e.first : e.second; adj[pos[s]++] = d; }
        #pragma omp parallel for schedule(dynamic, 4096)
        for (size_t v = 0; v < V; ++v) std::sort(adj.begin() + beg[v], adj.begin() + beg[v + 1]);
        std::string dir = argv[2];
        writeArr(dir + "/" + bname, beg);
        writeArr(dir + "/" + aname, adj);
    };
    std::string mk = std::string("mkdir -p ") + argv[2];
    if (std::system(mk.c_str()) != 0) return 1;
    build(false, "fw_begin.bin", "fw_adjacent.bin");
    build(true, "bw_begin.bin", "bw_adjacent.bin");
    return 0;
}
