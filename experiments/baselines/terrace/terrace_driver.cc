#include <cstdlib>
// Terrace (SIGMOD'21) on the GIR protocol: update phases + BFS/BC (unweighted
// build) or SSSP (weighted build, -DTERRACE_SSSP). See ../README.md.
#define ENABLE_LOCK 1
#ifdef TERRACE_SSSP
#define WEIGHTED 1
#else
#define WEIGHTED 0
#endif
#define VERIFY 0
#include <algorithm>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "graph.h"
#include "parallel.h"
#include "util.h"
#include "io_util.h"
#ifdef TERRACE_SSSP
#include "BellmanFord.h"
#else
#include "BFS.h"
#include "BC.h"
#endif
#include "../common/protocol.h"

using namespace graphstore;

// Bulk/batch insert: Terrace wants a batch sorted by (src, dst) without
// duplicates; it overwrites dests (and weights) in place, so pass copies.
static void insert_batch(Graph& G, std::vector<std::pair<uint32_t, uint32_t>>& es) {
    if (es.empty()) return;
    std::sort(es.begin(), es.end());
    es.erase(std::unique(es.begin(), es.end()), es.end());
    std::vector<uint32_t> s(es.size()), d(es.size());
    for (size_t i = 0; i < es.size(); ++i) { s[i] = es[i].first; d[i] = es[i].second; }
    std::vector<uint32_t> perm(es.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), std::mt19937_64(3));
#ifdef TERRACE_SSSP
    std::vector<uint32_t> w(es.size());
    for (size_t i = 0; i < es.size(); ++i) w[i] = proto::weight(s[i], d[i]);
    G.add_edge_batch(s.data(), d.data(), w.data(), (uint32_t)es.size(), perm);
#else
    G.add_edge_batch(s.data(), d.data(), (uint32_t)es.size(), perm);
#endif
}

static Graph* load(const std::string& path, uint32_t n) {
    auto es = proto::read_el(path);
    Graph* G = new Graph(n);
    // Load in chunks (add_edge_batch takes a 32-bit count).
    const size_t CH = 1u << 28;
    for (size_t i = 0; i < es.size(); i += CH) {
        std::vector<std::pair<uint32_t, uint32_t>> part(es.begin() + i, es.begin() + std::min(es.size(), i + CH));
        insert_batch(*G, part);
    }
    return G;
}

#ifndef TERRACE_SSSP
// Push-only BFS: Terrace's dense (pull) edgeMap reads out-edges as in-edges,
// which is only correct on symmetric graphs.
static int32_t* bfs_push(Graph& G, uint32_t src) {
    long n = G.get_num_vertices();
    int32_t* P = (int32_t*)malloc(n * sizeof(int32_t));
    parallel_for(long i = 0; i < n; i++) P[i] = -1;
    P[src] = src;
    VertexSubset fr(src, n);
    while (fr.not_empty()) { VertexSubset nx = edgeMap(G, fr, BFS_F(P), true, 1); fr.del(); fr = nx; }
    fr.del();
    return P;
}
#else
static int32_t* sssp_push(Graph& G, long start) {
    long n = G.get_num_vertices();
    int32_t* L = newA(int32_t, n);
    parallel_for(long i = 0; i < n; i++) L[i] = INT_MAX / 2;
    L[start] = 0;
    int* Vis = newA(int, n);
    parallel_for(long i = 0; i < n; i++) Vis[i] = 0;
    VertexSubset fr(start, n);
    while (fr.not_empty()) { VertexSubset o = edgeMap(G, fr, BF_F(L, Vis), true, 1); vertexMap(o, BF_Vertex_F(Vis)); fr.del(); fr = o; }
    fr.del();
    free(Vis);
    return L;
}
#endif

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: %s <protocol dir> <dataset label> <csv> [batch=10000]\n", argv[0]); return 2; }
    const std::string dir = argv[1], ds = argv[2], csv = argv[3];
    setvbuf(stdout, nullptr, _IOLBF, 0);
    const size_t B = argc > 4 ? std::stoul(argv[4]) : 10000;
    const uint32_t n = (uint32_t)proto::meta(dir, "V");
    const bool sym = proto::meta(dir, "symmetric") != 0;
    const int th = getWorkers();
    std::printf("Terrace: %s, workers %d, symmetric %d\n", ds.c_str(), th, (int)sym);

#ifndef TERRACE_SSSP
    // QUERIES_ONLY=1: skip the update phases (thread sweeps); NO_BC=1: skip BC.
    if (!std::getenv("QUERIES_ONLY")) {   // ---- updates on the 50% graph ----
        double t = proto::now();
        Graph* G = load(dir + "/loaded.el", n);
        std::printf("loaded 50%% graph in %.2f s (%u edges)\n", proto::now() - t, G->get_num_edges());
        auto run = [&](const std::string& f, const char* metric) {
            auto ops = proto::read_ops(dir + "/" + f);
            double t0 = proto::now();
            proto::batches(ops, B, [&](const proto::Op* o, size_t k) {
                std::vector<std::pair<uint32_t, uint32_t>> del, ins;
                for (size_t i = 0; i < k; ++i) (o[i].op == 'd' ? del : ins).push_back({o[i].u, o[i].v});
                // Terrace has no batch delete; its own benchmark deletes in a parallel loop.
                parallel_for(size_t i = 0; i < del.size(); i++) G->remove_edge(del[i].first, del[i].second);
                insert_batch(*G, ins);
            });
            double s = proto::now() - t0;
            proto::emit(csv, ds, "Terrace", metric, "batch_" + std::to_string(B), ops.size() / 1e6 / s, "MEPS", ops.size(), th);
        };
        run("insert.ops", "insert_MEPS");
        run("delete.ops", "delete_MEPS");
        run("mixed.ops", "mixed_MEPS");
        delete G;
    }
    {   // ---- queries on the full graph ----
        Graph* G = load(dir + "/full.el", n);
        auto src = proto::read_ints(dir + "/sources.txt");
        double tot = 0;
        long long reached0 = 0;
        for (size_t i = 0; i < src.size(); ++i) {
            double t0 = proto::now();
            int32_t* P = sym ? BFS_with_edge_map(*G, src[i]) : bfs_push(*G, src[i]);
            tot += proto::now() - t0;
            if (i == 0) for (uint32_t v = 0; v < n; ++v) reached0 += P[v] >= 0;
            free(P);
        }
        proto::emit(csv, ds, "Terrace", "bfs_ms", sym ? "stock" : "push-only", tot * 1000 / src.size(), "ms", src.size(), th);
        proto::emit(csv, ds, "Terrace", "bfs_reached_src0", "-", reached0, "vertices", 1, th);
        if (sym && !std::getenv("NO_BC")) {  // BC's backward pass follows out-edges: correct only on symmetric graphs
            auto bsrc = proto::read_ints(dir + "/bc_sources.txt");
            double t0 = proto::now();
            for (uint32_t s : bsrc) { auto* D = BC(*G, s); free(D); }
            proto::emit(csv, ds, "Terrace", "bc_ms_total", "stock", (proto::now() - t0) * 1000, "ms", bsrc.size(), th);
        }
        delete G;
    }
#else
    {
        Graph* G = load(dir + "/full.el", n);
        auto src = proto::read_ints(dir + "/sources.txt");
        double tot = 0;
        unsigned long long checksum = 0;
        for (size_t i = 0; i < src.size(); ++i) {
            double t0 = proto::now();
            int32_t* L = sym ? SSSP_BF(*G, src[i]) : sssp_push(*G, src[i]);
            tot += proto::now() - t0;
            if (i == 0) for (uint32_t v = 0; v < n; ++v) if (L[v] < INT_MAX / 2) checksum += (unsigned long long)L[v];
            free(L);
        }
        proto::emit(csv, ds, "Terrace", "sssp_ms", sym ? "bellman-ford" : "bellman-ford-push", tot * 1000 / src.size(), "ms", src.size(), th);
        proto::emit(csv, ds, "Terrace", "sssp_distsum_src0", "-", (double)checksum, "sum", 1, th);
        delete G;
    }
#endif
    return 0;
}
