#include <cstdlib>
// Aspen (PLDI'19) on the GIR protocol, symmetric graphs only (Aspen requires
// them; its dense edge_map reads out-lists as in-lists). No SSSP: Aspen has
// no weighted graphs. See ../README.md.
#include "graph/api.h"
#include "algorithms/BFS.h"
#include "algorithms/BC.h"

#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

#include "../common/protocol.h"

using E = std::tuple<uintV, uintV>;

template <class Graph>  // BFS.h's BFS returns void; this copy returns parents
pbbs::sequence<uintV> bfs_parents(Graph& G, uintV src) {
    size_t n = G.num_vertices();
    auto vtxs = G.fetch_all_vertices();
    auto frontier = vertex_subset(n, src);
    timer a, b, c;
    auto parents = pbbs::sequence<uintV>(n, [&](size_t) { return UINT_V_MAX; });
    parents[src] = src;
    while (!frontier.is_empty()) {
        auto out = G.edge_map(frontier, BFS_F(parents.begin()), vtxs, a, b, c, stay_dense);
        frontier.del();
        frontier = out;
    }
    frontier.del();
    return parents;
}

static versioned_graph<treeplus_graph>* build_vg(size_t n, std::vector<std::pair<uint32_t, uint32_t>> el) {
    std::sort(el.begin(), el.end());
    el.erase(std::unique(el.begin(), el.end()), el.end());
    size_t m = el.size();
    uintE* off = pbbs::new_array_no_init<uintE>(n);
    uintV* ed = pbbs::new_array_no_init<uintV>(m);
    std::vector<size_t> deg(n, 0);
    for (auto& e : el) deg[e.first]++;
    size_t s = 0;
    for (size_t i = 0; i < n; ++i) { off[i] = s; s += deg[i]; }
    for (size_t i = 0; i < m; ++i) ed[i] = el[i].second;
    return new versioned_graph<treeplus_graph>(n, m, off, ed);
}

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: %s <protocol dir> <dataset label> <csv> [batch=10000]\n", argv[0]); return 2; }
    const std::string dir = argv[1], ds = argv[2], csv = argv[3];
    setvbuf(stdout, nullptr, _IOLBF, 0);
    const size_t B = argc > 4 ? std::stoul(argv[4]) : 10000;
    const size_t n = (size_t)proto::meta(dir, "V");
    if (!proto::meta(dir, "symmetric")) { std::fprintf(stderr, "Aspen needs a symmetric protocol directory\n"); return 2; }
    const int th = (int)num_workers();
    const size_t nn = 1UL << pbbs::log2_up(n + 1);
    std::printf("Aspen: %s, workers %d\n", ds.c_str(), th);
    // QUERIES_ONLY=1: skip the update phases (thread sweeps); NO_BC=1: skip BC.
    if (!std::getenv("QUERIES_ONLY")) {
        double t = proto::now();
        auto* VG = build_vg(n, proto::read_el(dir + "/loaded.el"));
        std::printf("loaded 50%% graph in %.2f s\n", proto::now() - t);
        auto run = [&](const std::string& f, const char* metric) {
            auto ops = proto::read_ops(dir + "/" + f);
            double t0 = proto::now();
            proto::batches(ops, B, [&](const proto::Op* o, size_t k) {
                std::vector<E> del, ins;
                for (size_t i = 0; i < k; ++i) (o[i].op == 'd' ? del : ins).push_back(E(o[i].u, o[i].v));
                if (!del.empty()) VG->delete_edges_batch(del.size(), del.data(), false, true, nn);
                if (!ins.empty()) VG->insert_edges_batch(ins.size(), ins.data(), false, true, nn);
            });
            double s = proto::now() - t0;
            proto::emit(csv, ds, "Aspen", metric, "batch_" + std::to_string(B), ops.size() / 1e6 / s, "MEPS", ops.size(), th);
        };
        run("insert.ops", "insert_MEPS");
        run("delete.ops", "delete_MEPS");
        run("mixed.ops", "mixed_MEPS");
    }
    {
        auto* VG = build_vg(n, proto::read_el(dir + "/full.el"));
        auto S = VG->acquire_version();
        auto src = proto::read_ints(dir + "/sources.txt");
        double tot = 0;
        long long reached0 = 0;
        for (size_t i = 0; i < src.size(); ++i) {
            double t0 = proto::now();
            auto p = bfs_parents(S.graph, src[i]);
            tot += proto::now() - t0;
            if (i == 0) for (size_t v = 0; v < p.size(); ++v) reached0 += p[v] != UINT_V_MAX;
        }
        proto::emit(csv, ds, "Aspen", "bfs_ms", "stock", tot * 1000 / src.size(), "ms", src.size(), th);
        proto::emit(csv, ds, "Aspen", "bfs_reached_src0", "-", reached0, "vertices", 1, th);
        if (!std::getenv("NO_BC")) {
            auto bsrc = proto::read_ints(dir + "/bc_sources.txt");
            double t0 = proto::now();
            for (uint32_t s : bsrc) { auto d = BC(S.graph, s); }
            proto::emit(csv, ds, "Aspen", "bc_ms_total", "stock", (proto::now() - t0) * 1000, "ms", bsrc.size(), th);
        }
        VG->release_version(std::move(S));
    }
    return 0;
}
