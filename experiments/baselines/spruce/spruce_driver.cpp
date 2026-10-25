// Spruce (SIGMOD'24) on the GIR protocol. Spruce provides storage only; BFS,
// SSSP and BC are implemented here on top of its neighbour reads with the
// same algorithms GIR uses (level-synchronous BFS, Delta-stepping, Brandes
// with sources as parallel tasks). See ../README.md.
#include <omp.h>

#include <algorithm>
#include <atomic>
#include <climits>
#include <string>
#include <vector>

#include "spruce_transaction.h"
#include "../common/protocol.h"

std::atomic<long int> space_count_temp = 0;  // declared extern by Spruce; defined only in its demo

using G = SpruceTransVer;
using WeightedEdge = SpruceTransVer::WeightedEdge;
using WeightedOutEdgeSimple = SpruceTransVer::WeightedOutEdgeSimple;

static void load(G& g, const std::vector<std::pair<uint32_t, uint32_t>>& es) {
    #pragma omp parallel for schedule(dynamic, 4096)
    for (size_t i = 0; i < es.size(); ++i)
        G::InsertEdge(g, WeightedEdge{es[i].first, es[i].second, (double)proto::weight(es[i].first, es[i].second)});
}

static std::vector<int> bfs(G& g, uint32_t n, uint32_t s) {
    std::vector<int> lv(n, -1);
    lv[s] = 0;
    std::vector<uint32_t> fr{s};
    const int T = omp_get_max_threads();
    for (int L = 0; !fr.empty(); ++L) {
        std::vector<std::vector<uint32_t>> loc(T);
        #pragma omp parallel
        {
            std::vector<WeightedOutEdgeSimple> nb;
            auto& out = loc[omp_get_thread_num()];
            #pragma omp for schedule(dynamic, 256)
            for (size_t i = 0; i < fr.size(); ++i) {
                G::get_neighbours(g, fr[i], nb);
                for (auto& e : nb) {
                    int* p = &lv[e.des];
                    int expect = -1;
                    if (__atomic_load_n(p, __ATOMIC_RELAXED) == -1 &&
                        __atomic_compare_exchange_n(p, &expect, L + 1, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED))
                        out.push_back((uint32_t)e.des);
                }
            }
        }
        fr.clear();
        for (auto& l : loc) fr.insert(fr.end(), l.begin(), l.end());
    }
    return lv;
}

static std::vector<uint64_t> sssp(G& g, uint32_t n, uint32_t s, uint64_t delta) {
    std::vector<uint64_t> d(n, UINT64_MAX);
    std::vector<int> mark(n, 0);
    d[s] = 0;
    std::vector<std::vector<uint32_t>> bucket(1, {s});
    const int T = omp_get_max_threads();
    int round = 0;
    for (size_t k = 0; k < bucket.size(); ++k) {
        while (!bucket[k].empty()) {
            std::vector<uint32_t> fr;
            fr.swap(bucket[k]);
            ++round;
            std::vector<std::vector<uint32_t>> loc(T);
            #pragma omp parallel
            {
                std::vector<WeightedOutEdgeSimple> nb;
                auto& out = loc[omp_get_thread_num()];
                #pragma omp for schedule(dynamic, 256)
                for (size_t i = 0; i < fr.size(); ++i) {
                    uint32_t u = fr[i];
                    if (__atomic_exchange_n(&mark[u], round, __ATOMIC_RELAXED) == round) continue;
                    uint64_t du = __atomic_load_n(&d[u], __ATOMIC_RELAXED);
                    if (du / delta != k) continue;
                    G::get_neighbours(g, u, nb);
                    for (auto& e : nb) {
                        uint64_t nd = du + (uint64_t)e.weight;
                        uint64_t* p = &d[e.des];
                        uint64_t old = __atomic_load_n(p, __ATOMIC_RELAXED);
                        while (nd < old)
                            if (__atomic_compare_exchange_n(p, &old, nd, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) { out.push_back((uint32_t)e.des); break; }
                    }
                }
            }
            for (auto& l : loc)
                for (uint32_t w : l) {
                    size_t b = d[w] / delta;
                    if (bucket.size() <= b) bucket.resize(b + 1);
                    bucket[b].push_back(w);
                }
        }
    }
    return d;
}

static double bc_total(G& g, uint32_t n, const std::vector<uint32_t>& srcs) {
    std::vector<double> score(n, 0.0);
    double t0 = proto::now();
    #pragma omp parallel
    {
        std::vector<int> dist(n, -1);
        std::vector<double> sigma(n, 0), dep(n, 0);
        std::vector<uint32_t> order;
        std::vector<WeightedOutEdgeSimple> nb;
        #pragma omp for schedule(dynamic, 1)
        for (size_t si = 0; si < srcs.size(); ++si) {
            uint32_t s = srcs[si];
            order.clear();
            dist[s] = 0; sigma[s] = 1; order.push_back(s);
            for (size_t i = 0; i < order.size(); ++i) {
                uint32_t u = order[i];
                G::get_neighbours(g, u, nb);
                for (auto& e : nb) {
                    uint32_t w = (uint32_t)e.des;
                    if (dist[w] < 0) { dist[w] = dist[u] + 1; order.push_back(w); }
                    if (dist[w] == dist[u] + 1) sigma[w] += sigma[u];
                }
            }
            // Backward pass over out-edges (no reverse graph needed).
            for (size_t i = order.size(); i-- > 0;) {
                uint32_t v = order[i];
                G::get_neighbours(g, v, nb);
                for (auto& e : nb) {
                    uint32_t w = (uint32_t)e.des;
                    if (dist[w] == dist[v] + 1) dep[v] += sigma[v] / sigma[w] * (1 + dep[w]);
                }
                if (v != s) {
                    #pragma omp atomic
                    score[v] += dep[v];
                }
            }
            for (uint32_t x : order) { dist[x] = -1; sigma[x] = 0; dep[x] = 0; }
        }
    }
    return proto::now() - t0;
}

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: %s <protocol dir> <dataset label> <csv> [batch=10000]\n", argv[0]); return 2; }
    const std::string dir = argv[1], ds = argv[2], csv = argv[3];
    const size_t B = argc > 4 ? std::stoul(argv[4]) : 10000;
    const uint32_t n = (uint32_t)proto::meta(dir, "V");
    const int th = omp_get_max_threads();
    setvbuf(stdout, nullptr, _IOLBF, 0);
    std::printf("Spruce: %s, threads %d\n", ds.c_str(), th);
    // QUERIES_ONLY=1: skip the update phases; NO_SSSP=1 / NO_BC=1 skip those.
    if (!std::getenv("QUERIES_ONLY")) {
        G* g = new G();
        double t = proto::now();
        load(*g, proto::read_el(dir + "/loaded.el"));
        std::printf("loaded 50%% graph in %.2f s\n", proto::now() - t);
        auto run = [&](const std::string& f, const char* metric) {
            auto ops = proto::read_ops(dir + "/" + f);
            double t0 = proto::now();
            // Deletions of a batch first, then its insertions (each phase a
            // parallel loop). Running inserts and deletes of the mixed stream
            // fully concurrently livelocked Spruce (one thread spinning at
            // 100% for over an hour on livejournal), so the phases are split,
            // matching the in-batch order of Terrace, Aspen and DZiG.
            proto::batches(ops, B, [&](const proto::Op* o, size_t k) {
                #pragma omp parallel for schedule(dynamic, 64)
                for (size_t i = 0; i < k; ++i)
                    if (o[i].op == 'd') G::DeleteEdge(*g, o[i].u, o[i].v);
                #pragma omp parallel for schedule(dynamic, 64)
                for (size_t i = 0; i < k; ++i)
                    if (o[i].op == 'a') G::InsertEdge(*g, WeightedEdge{o[i].u, o[i].v, (double)proto::weight(o[i].u, o[i].v)});
            });
            double s = proto::now() - t0;
            proto::emit(csv, ds, "Spruce", metric, (std::string(metric) == "mixed_MEPS" ? "del-then-ins_" : "concurrent_") + std::to_string(B), ops.size() / 1e6 / s, "MEPS", ops.size(), th);
        };
        run("insert.ops", "insert_MEPS");
        run("delete.ops", "delete_MEPS");
        // Spruce bug: once DeleteEdge removes a vertex's last edge (Spruce
        // then deletes the vertex), a later InsertEdge at that vertex never
        // returns -- reproduced single-threaded on op 8,707 of the mixed
        // stream (vertex 518289). SPRUCE_SKIP_MIXED=1 records that instead.
        if (std::getenv("SPRUCE_SKIP_MIXED"))
            proto::emit(csv, ds, "Spruce", "mixed_MEPS", "HANGS:insert-at-emptied-vertex", -1, "MEPS", 0, th);
        else
            run("mixed.ops", "mixed_MEPS");
        // the store is not freed: Spruce has no destructor for its blocks
    }
    {
        G* g = new G();
        load(*g, proto::read_el(dir + "/full.el"));
        auto src = proto::read_ints(dir + "/sources.txt");
        double tb = 0, tsp = 0;
        long long reached0 = 0;
        unsigned long long sum0 = 0;
        for (size_t i = 0; i < src.size(); ++i) {
            double t0 = proto::now();
            auto lv = bfs(*g, n, src[i]);
            double t1 = proto::now();
            std::vector<uint64_t> d;
            if (!std::getenv("NO_SSSP")) d = sssp(*g, n, src[i], 16);
            double t2 = proto::now();
            tb += t1 - t0; tsp += t2 - t1;
            if (i == 0) {
                for (int x : lv) reached0 += x >= 0;
                for (auto x : d) if (x != UINT64_MAX) sum0 += x;
            }
        }
        proto::emit(csv, ds, "Spruce", "bfs_ms", "top-down", tb * 1000 / src.size(), "ms", src.size(), th);
        proto::emit(csv, ds, "Spruce", "bfs_reached_src0", "-", reached0, "vertices", 1, th);
        if (!std::getenv("NO_SSSP")) {
            proto::emit(csv, ds, "Spruce", "sssp_ms", "delta-stepping", tsp * 1000 / src.size(), "ms", src.size(), th);
            proto::emit(csv, ds, "Spruce", "sssp_distsum_src0", "-", (double)sum0, "sum", 1, th);
        }
        if (!std::getenv("NO_BC")) {
            auto bsrc = proto::read_ints(dir + "/bc_sources.txt");
            proto::emit(csv, ds, "Spruce", "bc_ms_total", "brandes", bc_total(*g, n, bsrc) * 1000, "ms", bsrc.size(), th);
        }
    }
    return 0;
}
