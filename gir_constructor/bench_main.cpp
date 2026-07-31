// First-version GIR-only benchmark harness (no external baselines -- see
// docs/ for why: Terrace/Aspen/Spruce comparisons and GraphBolt/DZiG are a
// separate, larger follow-up). Measures construction time, dynamic-update
// throughput (MEPS), and BFS runtime (GIR's two-level task-parallel engine
// vs. a plain sequential baseline) across whichever real datasets are
// available under data/, writing one CSV row per dataset to
// experiments/results/gir_benchmark.csv.
#include "gir.h"
#include "gir_test.cpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <functional>

struct DatasetSpec {
    std::string name;
    std::string path;
};

int main(int argc, char** argv) {
    bool twitter_only = argc > 1 && std::string(argv[1]) == "--twitter-only";

    std::vector<DatasetSpec> datasets = {
        {"example_graph", "/home/cxl6029/projects/GIR/data/example_graph/"},
        {"livejournal", "/home/cxl6029/projects/GIR/data/livejournal/"},
        {"web-Google", "/home/cxl6029/projects/GIR/data/web-Google_sym/"},
        {"RoadUSA", "/home/cxl6029/projects/GIR/data/RoadUSA_sym/"},
        {"grid_1000_100000", "/home/cxl6029/projects/GIR/data/grid_1000_100000/"},
    };
    if (argc > 1 && std::string(argv[1]) == "--with-twitter") {
        datasets.push_back({"twitter", "/home/cxl6029/projects/GIR/data/twitter/"});
    }
    if (twitter_only) {
        datasets = {{"twitter", "/home/cxl6029/projects/GIR/data/twitter/"}};
    }

    // --twitter-only appends one more row to the existing CSV (from the prior
    // 5-dataset run) instead of truncating it and re-running everything.
    std::string csv_path = "/home/cxl6029/projects/GIR/experiments/results/gir_benchmark.csv";
    std::ofstream csv(csv_path, twitter_only ? (std::ios::app) : std::ios::out);
    if (!twitter_only) {
        csv << "dataset,vertices,edges,load_s,wcc_build_s,sccdag_build_s,bct_build_s,"
               "num_scc,num_wcc,update_ops,insert_MEPS,delete_MEPS,mixed_MEPS,"
               "bfs_sources,bfs_taskparallel_ms_avg,bfs_sequential_ms_avg,bfs_speedup,bfs_diff_ok\n";
    }

    for (auto& ds : datasets) {
        std::cout << "\n########## Dataset: " << ds.name << " ##########" << std::endl;

        auto t0 = std::chrono::high_resolution_clock::now();
        GIR gir(ds.path, ds.path);
        auto t1 = std::chrono::high_resolution_clock::now();
        gir.buildWCC();
        auto t2 = std::chrono::high_resolution_clock::now();
        gir.buildSCC_DAG();
        auto t3 = std::chrono::high_resolution_clock::now();
        gir.buildBCT();
        auto t4 = std::chrono::high_resolution_clock::now();

        double load_s = std::chrono::duration<double>(t1 - t0).count();
        double wcc_s = std::chrono::duration<double>(t2 - t1).count();
        double scc_s = std::chrono::duration<double>(t3 - t2).count();
        double bct_s = std::chrono::duration<double>(t4 - t3).count();
        int V = gir.vertexCount();
        int num_scc = (int)gir.getSCCs().size();
        int num_wcc = (int)gir.getWCCs().size();

        std::cout << "V=" << V << " SCCs=" << num_scc << " WCCs=" << num_wcc << "\n"
                  << "load=" << load_s << "s wcc=" << wcc_s << "s sccdag=" << scc_s
                  << "s bct=" << bct_s << "s" << std::endl;

        // Update throughput. Capped per-phase operation count so this stays
        // tractable on very large graphs (the paper uses |E| ops per phase
        // and 10x|E| for mixed; we cap much lower here for a first pass and
        // note the methodology difference rather than hide it).
        //
        // Two real bugs were found and fixed while chasing multi-hour
        // stalls here (see CHANGELOG.md): Case 1 always shifted whichever
        // WCC contained the edge's source vertex instead of the smaller
        // one (catastrophic once a giant dominant WCC exists, as on
        // livejournal), and a DAG-edge-multiplicity bug. After both fixes,
        // individual Case-2 insertions between two numerically-distant but
        // topologically-unrelated SCCs (livejournal has 971K of them,
        // mostly near-singletons) can still cost 10-30+ seconds each: the
        // bounded rank-interval search has to explore every node in that
        // range to confirm "no cycle" even when the two endpoints have no
        // real relationship, and Kahn's-order rank assignment can place
        // unrelated components arbitrarily far apart in rank. This is a
        // genuine, not-yet-fixed scalability characteristic of the current
        // rank-interval design (see docs/code-review.md-style note in
        // CHANGELOG.md) rather than a bug being papered over -- op_count is
        // capped low enough for a first pass to actually finish. Scaled by
        // SCC count specifically (the actual driver of the pathology above,
        // confirmed by measurement: livejournal's 500-op insert phase alone
        // took ~82 minutes, ~10s/op average) rather than vertex count.
        long long op_count = (num_scc > 100000) ? 20 : std::min<long long>(V > 1 ? (long long)V : 1, 500);
        auto run_ops = [&](const std::vector<std::pair<int, int>>& ops, const char* label,
                           const std::function<void(int, int)>& fn) {
            auto t_start = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < ops.size(); ++i) {
                fn(ops[i].first, ops[i].second);
                if ((i + 1) % 50 == 0) {
                    double elapsed = std::chrono::duration<double>(
                        std::chrono::high_resolution_clock::now() - t_start).count();
                    std::cout << "  [" << label << "] " << (i + 1) << "/" << ops.size()
                              << " ops, " << elapsed << "s elapsed" << std::endl;
                }
            }
            return std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - t_start).count();
        };

        std::srand(1);
        std::vector<std::pair<int, int>> insert_pairs;
        insert_pairs.reserve(op_count);
        for (long long i = 0; i < op_count; ++i) {
            int u = std::rand() % V, v = std::rand() % V;
            if (u == v) { --i; continue; }
            insert_pairs.push_back({u, v});
        }

        double insert_s = run_ops(insert_pairs, "insert",
                                   [&](int u, int v) { gir.insertEdge(u, v); });
        double insert_meps = insert_s > 0 ? (op_count / 1e6) / insert_s : 0;

        double delete_s = run_ops(insert_pairs, "delete",
                                   [&](int u, int v) { gir.deleteEdge(u, v); });
        double delete_meps = delete_s > 0 ? (op_count / 1e6) / delete_s : 0;

        long long mixed_count = op_count;
        std::vector<std::pair<int, int>> mixed_ops;
        mixed_ops.reserve(mixed_count);
        for (long long i = 0; i < mixed_count; ++i) {
            int u = std::rand() % V, v = std::rand() % V;
            if (u == v) { --i; continue; }
            mixed_ops.push_back({u, v});
        }
        double mixed_s = run_ops(mixed_ops, "mixed", [&](int u, int v) {
            if (gir.hasDirectedEdge(u, v) || gir.hasUndirectedEdge(u, v))
                gir.deleteEdge(u, v);
            else
                gir.insertEdge(u, v);
        });
        double mixed_meps = mixed_s > 0 ? (mixed_count / 1e6) / mixed_s : 0;

        std::cout << "Insert MEPS=" << insert_meps << " Delete MEPS=" << delete_meps
                  << " Mixed MEPS=" << mixed_meps << " (op_count=" << op_count << ")" << std::endl;

        // BFS correctness (small sample) + timing.
        int diff_sources = std::min(V, 5);
        bool bfs_ok = GirBfsTest::runDifferentialBfsTest(gir, diff_sources, 7);

        int bfs_sources = std::min(V, 20);
        std::srand(7);
        std::vector<int> sources(bfs_sources);
        for (int i = 0; i < bfs_sources; ++i) sources[i] = std::rand() % V;

        auto tb0 = std::chrono::high_resolution_clock::now();
        for (int s : sources) gir.bfsTaskParallel(s);
        auto tb1 = std::chrono::high_resolution_clock::now();
        double task_ms_total = std::chrono::duration<double, std::milli>(tb1 - tb0).count();

        auto tb2 = std::chrono::high_resolution_clock::now();
        for (int s : sources) GirBfsTest::sequentialBfs(gir, s, V);
        auto tb3 = std::chrono::high_resolution_clock::now();
        double seq_ms_total = std::chrono::duration<double, std::milli>(tb3 - tb2).count();

        double task_ms_avg = task_ms_total / bfs_sources;
        double seq_ms_avg = seq_ms_total / bfs_sources;
        double bfs_speedup = task_ms_avg > 0 ? seq_ms_avg / task_ms_avg : 0;

        std::cout << "BFS: task-parallel avg=" << task_ms_avg << " ms, sequential avg=" << seq_ms_avg
                  << " ms, speedup=" << bfs_speedup << "x, diff_ok=" << bfs_ok << std::endl;

        csv << ds.name << "," << V << "," << gir.getDirectedGraph()->edge_count << "," << load_s << ","
            << wcc_s << "," << scc_s << "," << bct_s << "," << num_scc << "," << num_wcc << ","
            << op_count << "," << insert_meps << "," << delete_meps << "," << mixed_meps << ","
            << bfs_sources << "," << task_ms_avg << "," << seq_ms_avg << "," << bfs_speedup << ","
            << (bfs_ok ? 1 : 0) << "\n";
        csv.flush();
    }

    std::cout << "\nResults written to " << csv_path << std::endl;
    return 0;
}
