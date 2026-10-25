// Shared input reader / result writer for the baseline drivers.
//
// Inputs are produced by `bench_gir --export DIR` (gir_constructor/bench_main.cpp):
//   loaded.el    initial graph (50% of the edges), "u v" per line
//   full.el      the whole graph, for the query experiments
//   insert.ops / delete.ops / mixed.ops   "a u v" / "d u v" per line
//   sources.txt  BFS/SSSP sources; bc_sources.txt  BC sources
//   meta.txt     V, E, ops_per_phase, symmetric
// SSSP weights are computed with weight(u,v), identical to GIR::edgeWeight.
#pragma once
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace proto {

struct Op { char op; uint32_t u, v; };

inline double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Identical to GIR::edgeWeight (gir_constructor/gir_query.cpp), kMaxWeight = 100.
inline uint32_t weight(uint32_t u, uint32_t v) {
    uint64_t x = ((uint64_t)u << 32) | v;
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return 1 + (uint32_t)(x % 100);
}

struct Mapped {
    const char* p = nullptr; size_t n = 0; int fd = -1;
    explicit Mapped(const std::string& path) {
        fd = open(path.c_str(), O_RDONLY);
        if (fd < 0) { std::perror(path.c_str()); std::exit(1); }
        struct stat st; fstat(fd, &st); n = (size_t)st.st_size;
        p = n ? (const char*)mmap(nullptr, n, PROT_READ, MAP_PRIVATE, fd, 0) : nullptr;
    }
    ~Mapped() { if (p) munmap((void*)p, n); if (fd >= 0) close(fd); }
};

inline bool next_uint(const char*& p, const char* e, uint32_t& x) {
    while (p < e && (*p < '0' || *p > '9')) ++p;
    if (p >= e) return false;
    uint64_t v = 0;
    while (p < e && *p >= '0' && *p <= '9') v = v * 10 + (uint64_t)(*p++ - '0');
    x = (uint32_t)v;
    return true;
}

inline std::vector<std::pair<uint32_t, uint32_t>> read_el(const std::string& path) {
    Mapped m(path);
    std::vector<std::pair<uint32_t, uint32_t>> r;
    r.reserve(m.n / 16);
    const char *p = m.p, *e = m.p + m.n;
    uint32_t a, b;
    while (next_uint(p, e, a) && next_uint(p, e, b)) r.push_back({a, b});
    return r;
}

inline std::vector<Op> read_ops(const std::string& path) {
    Mapped m(path);
    std::vector<Op> r;
    const char *p = m.p, *e = m.p + m.n;
    while (p < e) {
        while (p < e && *p != 'a' && *p != 'd') ++p;
        if (p >= e) break;
        char c = *p++;
        uint32_t a, b;
        if (!next_uint(p, e, a) || !next_uint(p, e, b)) break;
        r.push_back({c, a, b});
    }
    return r;
}

inline std::vector<uint32_t> read_ints(const std::string& path) {
    Mapped m(path);
    std::vector<uint32_t> r;
    const char *p = m.p, *e = m.p + m.n;
    uint32_t a;
    while (next_uint(p, e, a)) r.push_back(a);
    return r;
}

inline long meta(const std::string& dir, const std::string& key) {
    std::ifstream f(dir + "/meta.txt");
    std::string k; long v;
    while (f >> k >> v) if (k == key) return v;
    std::fprintf(stderr, "meta key %s missing\n", key.c_str());
    std::exit(1);
}

// Appends one result row: dataset,system,metric,mode,value,unit,count,threads
inline void emit(const std::string& csv, const std::string& dataset, const std::string& system,
                 const std::string& metric, const std::string& mode, double value, const char* unit,
                 long long count, int threads) {
    FILE* f = std::fopen(csv.c_str(), "a");
    std::fprintf(f, "%s,%s,%s,%s,%.6g,%s,%lld,%d\n", dataset.c_str(), system.c_str(), metric.c_str(),
                 mode.c_str(), value, unit, count, threads);
    std::fclose(f);
    std::printf("  %-22s %-14s %.6g %s (n=%lld, threads=%d)\n", metric.c_str(), mode.c_str(), value, unit, count, threads);
}

// Splits a stream into consecutive batches of size B.
template <class F>
inline void batches(const std::vector<Op>& ops, size_t B, F f) {
    for (size_t i = 0; i < ops.size(); i += B) f(ops.data() + i, std::min(B, ops.size() - i));
}

}  // namespace proto
