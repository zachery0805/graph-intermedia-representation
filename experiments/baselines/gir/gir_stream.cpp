// GIR in DZiG's streaming setting: load the 50% graph, then per batch of B
// operations apply the batch (FMU) and recompute BFS / SSSP from one source
// from scratch (GIR keeps no per-source state). Same streams and source as
// ../dzig/run_dzig.sh. Usage: gir_stream <protocol dir> <dataset> <csv> [batch=10000] [batches=10] [lazy]
// With "lazy", batches go through ingestBatch (edge store only; the layers
// are not updated) and rows are reported as system "GIR-lazy"; the layers'
// catch-up cost is reported once per stream as sync_*_ms.
// With "inc" / "lazy-inc", BFS and SSSP are corrected incrementally per batch
// (trackIncremental / updateIncremental, as DZiG/KickStarter do) instead of
// recomputed; rows go to "GIR-inc" / "GIR-lazy-inc".
#include <omp.h>
#include <cstdio>

#include <string>
#include <vector>

#include "gir.h"
#include "graph.h"
#include "../common/protocol.h"

int main(int argc, char** argv) {
    if (argc < 4) return 2;
    const std::string dir = argv[1], ds = argv[2], csv = argv[3];
    const size_t B = argc > 4 ? std::stoul(argv[4]) : 10000;
    const int NB = argc > 5 ? std::stoi(argv[5]) : 10;
    const std::string mode = argc > 6 ? argv[6] : "";
    const bool lazy = mode == "lazy" || mode == "lazy-inc";
    const bool incr = mode == "inc" || mode == "lazy-inc";
    const std::string sys = std::string("GIR") + (lazy ? "-lazy" : "") + (incr ? "-inc" : "");
    const int V = (int)proto::meta(dir, "V");
    const int th = omp_get_max_threads();
    const int src = (int)proto::read_ints(dir + "/sources.txt")[0];
    auto el = proto::read_el(dir + "/loaded.el");
    std::vector<std::pair<int, int>> edges(el.begin(), el.end());
    for (const char* s : {"insert", "delete", "mixed"}) {
        GIR g(graph::fromEdges(V, edges));
        g.buildWCC(); g.buildSCC_DAG(); g.buildBCT();
        double t0 = proto::now();
        if (incr) g.trackIncremental(src, false); else g.bfs(src);
        double t1 = proto::now();
        if (incr) g.trackIncremental(src, true); else g.sssp(src);
        double t2 = proto::now();
        if (std::string(s) == "insert") {
            proto::emit(csv, ds, sys, "bfs_initial_ms", "from-scratch", (t1 - t0) * 1000, "ms", 1, th);
            proto::emit(csv, ds, sys, "sssp_initial_ms", "from-scratch", (t2 - t1) * 1000, "ms", 1, th);
        }
        auto ops = proto::read_ops(dir + "/" + s + ".ops");
        double tu = 0, tb = 0, tsp = 0;
        int nb = 0;
        for (size_t i = 0; i < ops.size() && nb < NB; i += B, ++nb) {
            std::vector<GIR::Update> batch;
            for (size_t j = i; j < std::min(ops.size(), i + B); ++j) batch.push_back({(int)ops[j].u, (int)ops[j].v, ops[j].op == 'a'});
            double a = proto::now();
            if (lazy) g.ingestBatch(batch); else g.applyBatch(batch);
            double b = proto::now();
            if (incr) g.updateIncremental(false); else g.bfs(src);
            double c = proto::now();
            if (incr) g.updateIncremental(true); else g.sssp(src);
            double d = proto::now();
            tu += b - a; tb += c - b; tsp += d - c;
        }
        if (incr) {  // the corrected results must equal from-scratch ones
            auto lv = g.bfs(src);
            auto d = g.sssp(src);
            long long bad = 0;
            for (int v = 0; v < V; ++v) {
                uint64_t want = lv[v] < 0 ? UINT64_MAX : (uint64_t)lv[v];
                bool b = g.incDist(false, v) != want || g.incDist(true, v) != d[v];
                if (b && bad < 6) fprintf(stderr, "  v=%d bfs inc %llu true %llu | sssp inc %llu true %llu\n", v, (unsigned long long)g.incDist(false, v), (unsigned long long)want, (unsigned long long)g.incDist(true, v), (unsigned long long)d[v]);
                bad += b;
            }
            fprintf(stderr, "%s stream: incremental vs from-scratch mismatches: %lld\n", s, bad);
            proto::emit(csv, ds, sys, std::string("inc_mismatch_") + s, "final", (double)bad, "vertices", V, th);
        }
        const std::string st = s;
        const std::string mm = (lazy ? "ingest_batch_" : "fmu_batch_") + std::to_string(B);
        proto::emit(csv, ds, sys, st + "_mutation_MEPS", mm, nb * (double)B / 1e6 / tu, "MEPS", nb * (long long)B, th);
        proto::emit(csv, ds, sys, "bfs_" + st + "_recompute_ms", incr ? "incremental" : "from-scratch", tb * 1000 / nb, "ms", nb, th);
        proto::emit(csv, ds, sys, "bfs_" + st + "_endtoend_ms", "batch_" + std::to_string(B), (tu + tb) * 1000 / nb, "ms", nb, th);
        proto::emit(csv, ds, sys, "sssp_" + st + "_recompute_ms", incr ? "incremental" : "from-scratch", tsp * 1000 / nb, "ms", nb, th);
        proto::emit(csv, ds, sys, "sssp_" + st + "_endtoend_ms", "batch_" + std::to_string(B), (tu + tsp) * 1000 / nb, "ms", nb, th);
        if (lazy) {  // the layers' catch-up for the whole stream
            double a = proto::now(); g.sync(); double b = proto::now();
            proto::emit(csv, ds, sys, "sync_" + st + "_ms", "after_" + std::to_string(nb) + "_batches", (b - a) * 1000, "ms", nb, th);
        }
    }
    return 0;
}
