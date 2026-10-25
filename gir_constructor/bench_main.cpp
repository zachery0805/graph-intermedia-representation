// GIR benchmark following the paper's protocol (Sec. 7):
//
// Updates (DZiG-style hold-out): shuffle the edges, build GIR on 50% of them,
// then
//   insertion-only : stream K held-out edges back in,
//   deletion-only  : delete K edges sampled from the loaded half,
//   mixed          : K/2 fresh held-out insertions and K/2 fresh deletions
//                    from the loaded half, interleaved in one shuffled order.
// Each mode runs once with single-edge updates (SPU) and once on a fresh
// build with applyBatch (FMU, batch size B). Phases are capped by a time
// budget; MEPS is over the operations actually completed.
//
// Queries (on the full graph): BFS (direction-optimizing, two-level, and a
// sequential baseline), SSSP (Delta-stepping, two-level), BC (Brandes over k
// sources), each averaged over S random sources; a thread sweep for BFS.
//
// Usage: bench_gir [--datasets a,b] [--ops K] [--batch B] [--budget SEC]
//                  [--sources S] [--bc-sources k] [--threads 1,2,4,...]
//                  [--skip-updates] [--skip-queries] [--mode spu|fmu|lazy|none] [--workload] [--reach-queries N] [--blocks-on-demand] [--csv PATH]
//                  [--export DIR]  (write the protocol's edge lists, update streams
//                                   and query sources for the baseline systems)
#include "gir.h"
#include "gir_test.cpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <omp.h>
#include <cstdio>
#include <cstdlib>

#include "graph.h"

using Clock = std::chrono::steady_clock;
static double secs(Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); }

struct Phase {
    long long ops = 0;
    double s = 0;       // update time (for LAZY: ingest only)
    double sync_s = 0;  // LAZY: the sync() that brings the layers up to date at the end
    double blocks_s = 0;  // --blocks-on-demand: ensureBlocks() after the phase
    double meps() const { return s > 0 ? ops / 1e6 / s : 0; }
    double meps_total() const { return s + sync_s > 0 ? ops / 1e6 / (s + sync_s) : 0; }
};

// --verify: SCC and WCC partitions recomputed from scratch on the current
// edge store (iterative Tarjan, union-find) and compared with GIR's, plus
// rank(scc(u)) < rank(scc(v)) for every condensation edge u -> v.
static bool verifyStructure(GIR& g, const std::string& label) {
    auto t0 = Clock::now();
    const int n = g.vertexCount();
    std::vector<std::vector<int>> adj(n);
    #pragma omp parallel for schedule(dynamic, 4096)
    for (int v = 0; v < n; ++v) adj[v] = g.getDirectedOutNeighbors(v);
    bool ok = true;
    long long bad_rank = 0;
    #pragma omp parallel for schedule(dynamic, 4096) reduction(+ : bad_rank)
    for (int v = 0; v < n; ++v)
        for (int w : adj[v])
            if (g.sccOf(v) != g.sccOf(w) && g.sccRank(g.sccOf(v)) >= g.sccRank(g.sccOf(w))) ++bad_rank;
    if (bad_rank) { std::cout << "VERIFY " << label << ": " << bad_rank << " condensation edges against the order" << std::endl; ok = false; }
    // Tarjan, iterative.
    std::vector<int> idx(n, -1), low(n), comp(n, -1), st, cs;
    std::vector<size_t> it(n);
    int counter = 0, ncomp = 0;
    for (int r = 0; r < n; ++r) {
        if (idx[r] >= 0) continue;
        cs.push_back(r); idx[r] = low[r] = counter++; st.push_back(r); it[r] = 0;
        while (!cs.empty()) {
            int v = cs.back();
            if (it[v] < adj[v].size()) {
                int w = adj[v][it[v]++];
                if (idx[w] < 0) { idx[w] = low[w] = counter++; st.push_back(w); it[w] = 0; cs.push_back(w); }
                else if (comp[w] < 0) low[v] = std::min(low[v], idx[w]);
            } else {
                cs.pop_back();
                if (!cs.empty()) low[cs.back()] = std::min(low[cs.back()], low[v]);
                if (low[v] == idx[v]) {
                    int w;
                    do { w = st.back(); st.pop_back(); comp[w] = ncomp; } while (w != v);
                    ++ncomp;
                }
            }
        }
    }
    // Union-find for WCCs.
    std::vector<int> uf(n);
    for (int v = 0; v < n; ++v) uf[v] = v;
    auto find = [&](int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; };
    for (int v = 0; v < n; ++v) for (int w : adj[v]) { int a = find(v), b = find(w); if (a != b) uf[a] = b; }
    // Partition equality: truth id <-> gir id must be a bijection.
    auto same = [&](auto truth, auto gir, const char* what) {
        std::unordered_map<long long, long long> t2g, g2t;
        long long bad = 0;
        for (int v = 0; v < n; ++v) {
            long long a = truth(v), b = gir(v);
            auto p = t2g.emplace(a, b); auto q = g2t.emplace(b, a);
            if (p.first->second != b || q.first->second != a) ++bad;
        }
        if (bad) { std::cout << "VERIFY " << label << ": " << what << " partition differs at " << bad << " vertices (truth " << t2g.size() << " groups, gir " << g2t.size() << ")" << std::endl; ok = false; }
        return t2g.size();
    };
    size_t ns = same([&](int v) { return (long long)comp[v]; }, [&](int v) { return (long long)g.sccOf(v); }, "SCC");
    size_t nw = same([&](int v) { return (long long)find(v); }, [&](int v) { return (long long)g.wccOf(v); }, "WCC");
    std::cout << "VERIFY " << label << ": " << (ok ? "OK" : "FAILED") << " (" << ns << " SCCs, " << nw << " WCCs, " << secs(t0, Clock::now()) << " s)" << std::endl;
    return ok;
}

static Phase runSingle(GIR& g, const std::vector<GIR::Update>& ops, double budget) {
    Phase r;
    auto t0 = Clock::now();
    for (auto& op : ops) {
        if (op.insert) g.insertEdge(op.u, op.v); else g.deleteEdge(op.u, op.v);
        if ((++r.ops & 63) == 0 && secs(t0, Clock::now()) > budget) break;
    }
    r.s = secs(t0, Clock::now());
    return r;
}
static Phase runBatched(GIR& g, const std::vector<GIR::Update>& ops, size_t B, double budget) {
    Phase r;
    auto t0 = Clock::now();
    for (size_t i = 0; i < ops.size(); i += B) {
        std::vector<GIR::Update> b(ops.begin() + i, ops.begin() + std::min(ops.size(), i + B));
        auto tb = Clock::now();
        g.applyBatch(b);
        if (getenv("GIR_BATCH_TIMES")) std::cerr << "batch " << i / B << ": " << secs(tb, Clock::now()) * 1000 << " ms" << std::endl;
        r.ops += (long long)b.size();
        if (secs(t0, Clock::now()) > budget) break;
    }
    r.s = secs(t0, Clock::now());
    return r;
}

// Deferred maintenance: every batch goes to the edge store only; one sync()
// at the end of the phase updates the layers.
static Phase runLazy(GIR& g, const std::vector<GIR::Update>& ops, size_t B, double budget) {
    Phase r;
    auto t0 = Clock::now();
    for (size_t i = 0; i < ops.size(); i += B) {
        std::vector<GIR::Update> b(ops.begin() + i, ops.begin() + std::min(ops.size(), i + B));
        g.ingestBatch(b);
        r.ops += (long long)b.size();
        if (secs(t0, Clock::now()) > budget) break;
    }
    r.s = secs(t0, Clock::now());
    auto t1 = Clock::now();
    g.sync();
    r.sync_s = secs(t1, Clock::now());
    return r;
}

static GIR* build(int V, const std::vector<std::pair<int, int>>& edges, double& build_s) {
    auto t0 = Clock::now();
    GIR* g = new GIR(graph::fromEdges(V, edges));
    g->buildWCC();
    g->buildSCC_DAG();
    g->buildBCT();
    build_s = secs(t0, Clock::now());
    return g;
}

int main(int argc, char** argv) {
    std::vector<std::string> names = {"livejournal", "twitter"};
    long long K = 100000;
    size_t B = 10000;
    double budget = 120;
    int nsrc = 16, nbc = 64;
    std::vector<int> threads;
    bool do_updates = true, do_queries = true;
    int mode_lo = 0, mode_hi = 3;  // --mode spu|fmu|lazy: only that update mode
    bool blocks_od = false;
    int nreach = 10000;            // --reach-queries N: random u ~> v pairs (index vs. no index)
    long long hub_min = -1;        // --hub-anchor-min N: GIR::setHubAnchorMin (default: GIR's)
    bool verify = false;
    uint64_t sssp_delta = 16;      // --delta N: delta-stepping bucket width           // --verify: structure checked against a from-scratch recomputation after each phase
    bool workload = false;         // --workload: mixed stream batches each followed by reachability queries        // --blocks-on-demand: block layer rebuilt after each phase instead of per update
    std::string export_dir;  // --export DIR: write protocol inputs for baselines and exit
    std::string csv_path = "/home/cxl6029/projects/GIR/experiments/results/gir_benchmark.csv";
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(i + 1 < argc ? argv[++i] : ""); };
        auto list = [&](const std::string& s) { std::vector<std::string> r; std::stringstream ss(s); std::string t; while (std::getline(ss, t, ',')) r.push_back(t); return r; };
        if (a == "--datasets") names = list(next());
        else if (a == "--ops") K = std::stoll(next());
        else if (a == "--batch") B = std::stoul(next());
        else if (a == "--budget") budget = std::stod(next());
        else if (a == "--sources") nsrc = std::stoi(next());
        else if (a == "--bc-sources") nbc = std::stoi(next());
        else if (a == "--threads") for (auto& t : list(next())) threads.push_back(std::stoi(t));
        else if (a == "--skip-updates") do_updates = false;
        else if (a == "--skip-queries") do_queries = false;
        else if (a == "--blocks-on-demand") blocks_od = true;
        else if (a == "--reach-queries") nreach = std::stoi(next());
        else if (a == "--workload") workload = true;
        else if (a == "--verify") verify = true;
        else if (a == "--delta") sssp_delta = std::stoull(next());
        else if (a == "--hub-anchor-min") hub_min = std::stoll(next());
        else if (a == "--mode") {
            std::string m = next();
            mode_lo = m == "fmu" ? 1 : m == "lazy" ? 2 : 0;
            mode_hi = m == "none" ? mode_lo : mode_lo + 1;  // none: only --workload / queries
        }
        else if (a == "--export") export_dir = next();
        else if (a == "--csv") csv_path = next();
    }
    const std::string base = "/home/cxl6029/projects/GIR/data/";
    std::map<std::string, std::string> path = {{"example_graph", base + "example_graph/"},
                                               {"livejournal", base + "livejournal/"},
                                               {"twitter", base + "twitter/"},
                                               // the paper's other SNAP graphs (converted with el_to_csr)
                                               {"twitter_combined", base + "twitter_combined/"},
                                               {"web-Google", base + "web-Google/"},
                                               {"wiki-topcats", base + "wiki-topcats/"},
                                               {"cit-Patents", base + "cit-Patents/"},
                                               {"friendster", base + "friendster/"}};
    std::ofstream csv(csv_path);
    csv << "dataset,metric,mode,value,unit,count,threads\n";
    auto row = [&](const std::string& d, const std::string& m, const std::string& mode, double v, const char* unit,
                   long long n, int th) {
        csv << d << "," << m << "," << mode << "," << v << "," << unit << "," << n << "," << th << "\n";
        csv.flush();
    };
    const int maxth = omp_get_max_threads();

    for (auto& name : names) {
        // "<name>_sym": the symmetrized graph. Updates are undirected: each
        // held-out / deleted pair {u,v} becomes the two directed operations
        // u->v and v->u (counted as two edge operations, as Aspen/Terrace do).
        const bool sym = name.size() > 4 && name.compare(name.size() - 4, 4, "_sym") == 0;
        const std::string base_name = sym ? name.substr(0, name.size() - 4) : name;
        if (!path.count(base_name)) { std::cerr << "unknown dataset " << name << "\n"; continue; }
        std::cout << "\n########## " << name << " ##########" << std::endl;
        std::vector<std::pair<int, int>> all;
        int V;
        {
            const std::string p = path[base_name];
            graph g((p + "fw_begin.bin").c_str(), (p + "fw_adjacent.bin").c_str(), (p + "bw_begin.bin").c_str(),
                    (p + "bw_adjacent.bin").c_str());
            V = (int)g.vert_count;
            all.reserve(g.edge_count);
            for (int u = 0; u < V; ++u) {
                int last = -1;
                // drop self loops and duplicates (lists may be unsorted)
                std::vector<int> nb(g.fw_csr + g.fw_beg_pos[u], g.fw_csr + g.fw_beg_pos[u + 1]);
                std::sort(nb.begin(), nb.end());
                for (int w : nb) if (w != u && w != last) { all.push_back({u, w}); last = w; }
            }
        }
        if (sym) {  // unordered pairs u < v, deduplicated
            for (auto& e : all) if (e.first > e.second) std::swap(e.first, e.second);
            std::sort(all.begin(), all.end());
            all.erase(std::unique(all.begin(), all.end()), all.end());
        }
        auto both = [&](const std::vector<std::pair<int, int>>& es) {
            if (!sym) return es;
            std::vector<std::pair<int, int>> r;
            r.reserve(es.size() * 2);
            for (auto& e : es) { r.push_back(e); r.push_back({e.second, e.first}); }
            return r;
        };
        std::cout << "V=" << V << " E=" << all.size() << " (self-loops/duplicates removed)" << std::endl;
        std::mt19937_64 rng(2027);

        if (do_updates) {
            std::shuffle(all.begin(), all.end(), rng);
            const size_t half = all.size() / 2;
            std::vector<std::pair<int, int>> loaded_pairs(all.begin(), all.begin() + half);
            std::vector<std::pair<int, int>> loaded = both(loaded_pairs);
            // held-out stream: all[half ..]; deletions: sampled from loaded
            const size_t k = (size_t)std::min<long long>(K, (long long)half / 2);
            auto make = [&](auto first, auto last, bool ins) {
                std::vector<GIR::Update> v;
                for (auto it = first; it != last; ++it) {
                    v.push_back({it->first, it->second, ins});
                    if (sym) v.push_back({it->second, it->first, ins});
                }
                return v;
            };
            auto ins_ops = make(all.begin() + half, all.begin() + half + k, true);
            auto del_ops = make(loaded_pairs.begin(), loaded_pairs.begin() + k, false);
            std::vector<GIR::Update> mix = make(all.begin() + half + k, all.begin() + half + k + k / 2, true);
            auto mdel = make(loaded_pairs.begin() + k, loaded_pairs.begin() + k + k / 2, false);
            mix.insert(mix.end(), mdel.begin(), mdel.end());
            if (!sym) {
                std::shuffle(mix.begin(), mix.end(), rng);
            } else {  // shuffle whole pairs so u->v and v->u stay adjacent
                std::vector<size_t> idx(mix.size() / 2);
                for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
                std::shuffle(idx.begin(), idx.end(), rng);
                std::vector<GIR::Update> m2;
                m2.reserve(mix.size());
                for (size_t i : idx) { m2.push_back(mix[2 * i]); m2.push_back(mix[2 * i + 1]); }
                mix.swap(m2);
            }

            if (!export_dir.empty()) {
                // Identical inputs for the baseline systems (see experiments/baselines/).
                const std::string d = export_dir + "/" + name + "/";
                std::string mk = "mkdir -p " + d;
                if (std::system(mk.c_str()) != 0) { std::cerr << "cannot create " << d << "\n"; return 1; }
                auto write_el = [&](const std::string& f, const std::vector<std::pair<int, int>>& es, bool weighted) {
                    FILE* fp = std::fopen((d + f).c_str(), "w");
                    for (auto& e : es) {
                        if (weighted) std::fprintf(fp, "%d %d %u\n", e.first, e.second, GIR::edgeWeight(e.first, e.second));
                        else std::fprintf(fp, "%d %d\n", e.first, e.second);
                    }
                    std::fclose(fp);
                };
                auto write_ops = [&](const std::string& f, const std::vector<GIR::Update>& ops) {
                    FILE* fp = std::fopen((d + f).c_str(), "w");
                    for (auto& o : ops) std::fprintf(fp, "%c %d %d\n", o.insert ? 'a' : 'd', o.u, o.v);
                    std::fclose(fp);
                };
                write_el("loaded.el", loaded, false);
                if (!std::getenv("EXPORT_NO_WEIGHTED")) write_el("loaded_w.el", loaded, true);  // (no driver reads these; skipped for the billion-edge graphs)
                auto full = both(all);
                write_el("full.el", full, false);
                if (!std::getenv("EXPORT_NO_WEIGHTED")) write_el("full_w.el", full, true);  // (no driver reads these; skipped for the billion-edge graphs)
                write_ops("insert.ops", ins_ops);
                write_ops("delete.ops", del_ops);
                write_ops("mixed.ops", mix);
                // Query sources, chosen exactly as in the query phase below.
                std::vector<long long> outdeg(V, 0);
                for (auto& e : full) outdeg[e.first]++;
                std::mt19937 srng(7);
                std::vector<int> src;
                while ((int)src.size() < nsrc) { int s = (int)(srng() % V); if (outdeg[s] > 0) src.push_back(s); }
                std::mt19937 brng(11);
                FILE* fp = std::fopen((d + "sources.txt").c_str(), "w");
                for (int s : src) std::fprintf(fp, "%d\n", s);
                std::fclose(fp);
                fp = std::fopen((d + "bc_sources.txt").c_str(), "w");
                for (int i = 0; i < nbc; ++i) std::fprintf(fp, "%d\n", (int)(brng() % V));
                std::fclose(fp);
                fp = std::fopen((d + "meta.txt").c_str(), "w");
                std::fprintf(fp, "V %d\nE %zu\nloaded %zu\nops_per_phase %zu\nmax_weight %u\nsymmetric %d\n", V, full.size(), loaded.size(), ins_ops.size(), GIR::kMaxWeight, (int)sym);
                std::fclose(fp);
                std::cout << "exported protocol inputs to " << d << std::endl;
                continue;
            }

            // Update + query workload: each mixed batch is followed by Q random
            // u ~> v queries. Per batch we time the update and 10,000 queries;
            // the cost of a batch with Q queries is update + Q * mean query time.
            //   index:    applyBatch (blocks on demand) + reachable()
            //   deferred: ingestBatch + reachable() (the first query syncs)
            //   no index: ingestBatch only + bidirectional BFS on the store
            if (workload) {
                const int QS = 10000, NB = 10;
                const long long qlist[] = {0, 10, 100, 1000, 10000, 100000};
                for (int variant = 0; variant < 3; ++variant) {
                    const char* vn[] = {"index", "deferred", "no-index"};
                    double bs;
                    GIR* g = build(V, loaded, bs);
                    g->setBlocksOnDemand(true);
                    std::mt19937 qrng(17);
                    double upd = 0, qry = 0;  // totals over NB batches (query time per QS queries)
                    for (int b = 0; b < NB && (size_t)b * B < mix.size(); ++b) {
                        std::vector<GIR::Update> batch(mix.begin() + b * B, mix.begin() + std::min(mix.size(), (size_t)(b + 1) * B));
                        auto t0 = Clock::now();
                        if (variant == 0) g->applyBatch(batch); else g->ingestBatch(batch);
                        if (variant == 1) g->sync();  // charged to the update side
                        upd += secs(t0, Clock::now());
                        std::vector<std::pair<int, int>> qs(QS);
                        for (auto& q : qs) q = {(int)(qrng() % V), (int)(qrng() % V)};
                        auto t1 = Clock::now();
                        long long yes = 0;
                        for (auto& q : qs) yes += variant == 2 ? g->reachableNoIndex(q.first, q.second) : g->reachable(q.first, q.second);
                        qry += secs(t1, Clock::now());
                        (void)yes;
                    }
                    const double upd_ms = upd * 1000 / NB, q_us = qry * 1e6 / (NB * (double)QS);
                    std::cout << "workload " << vn[variant] << ": update " << upd_ms << " ms/batch, query " << q_us << " us; per batch with Q queries:";
                    for (long long q : qlist) {
                        double ms = upd_ms + q * q_us / 1000;
                        std::cout << " Q=" << q << ": " << ms << " ms";
                        row(name, "workload_ms_per_batch", std::string(vn[variant]) + "_Q" + std::to_string(q), ms, "ms", NB, maxth);
                    }
                    std::cout << std::endl;
                    row(name, "workload_update_ms", vn[variant], upd_ms, "ms", NB, maxth);
                    row(name, "workload_query_us", vn[variant], q_us, "us", NB * (long long)QS, maxth);
                    delete g;
                }
            }

            for (int batched = mode_lo; batched < mode_hi; ++batched) {
                double bs;
                GIR* g = build(V, loaded, bs);
                if (!batched) row(name, "build_50pct", "-", bs, "s", (long long)loaded.size(), maxth);
                const bool lazy = batched == 2;
                std::string mode = lazy ? "LAZY_B" + std::to_string(B) : batched ? "FMU_B" + std::to_string(B) : "SPU";
                if (blocks_od) mode += "_bod";
                g->setBlocksOnDemand(blocks_od);
                if (hub_min >= 0) g->setHubAnchorMin(hub_min);
                g->resetCaseCounts();
                auto run = [&](const std::vector<GIR::Update>& ops) {
                    Phase p = lazy ? runLazy(*g, ops, B, budget) : batched ? runBatched(*g, ops, B, budget) : runSingle(*g, ops, budget);
                    if (blocks_od) { auto t = Clock::now(); g->ensureBlocks(); p.blocks_s = secs(t, Clock::now()); }
                    if (verify) { if (lazy) g->sync(); verifyStructure(*g, mode); }
                    return p;
                };
                Phase pi = run(ins_ops), pd = run(del_ops), pm = run(mix);
                std::cout << mode << ": insert " << pi.ops << " ops " << pi.meps() << " MEPS | delete " << pd.ops
                          << " ops " << pd.meps() << " MEPS | mixed " << pm.ops << " ops " << pm.meps() << " MEPS"
                          << std::endl;
                std::cout << mode << " bulk split rounds: anchor set " << g->anchor_set_calls << ", hub " << g->hub_anchor_calls << std::endl;
                if (blocks_od) {
                    std::cout << mode << " block rebuild after phase: insert " << pi.blocks_s << " s | delete " << pd.blocks_s
                              << " s | mixed " << pm.blocks_s << " s" << std::endl;
                    row(name, "blocks_rebuild_s", mode + "_insert", pi.blocks_s, "s", pi.ops, maxth);
                    row(name, "blocks_rebuild_s", mode + "_delete", pd.blocks_s, "s", pd.ops, maxth);
                    row(name, "blocks_rebuild_s", mode + "_mixed", pm.blocks_s, "s", pm.ops, maxth);
                }
                if (lazy) {
                    std::cout << mode << " incl. sync: insert " << pi.meps_total() << " MEPS (sync " << pi.sync_s
                              << " s) | delete " << pd.meps_total() << " MEPS (sync " << pd.sync_s << " s) | mixed "
                              << pm.meps_total() << " MEPS (sync " << pm.sync_s << " s)" << std::endl;
                    mode = "LAZY_ingest_B" + std::to_string(B) + (blocks_od ? "_bod" : "");
                }
                row(name, "insert_MEPS", mode, pi.meps(), "MEPS", pi.ops, maxth);
                row(name, "delete_MEPS", mode, pd.meps(), "MEPS", pd.ops, maxth);
                row(name, "mixed_MEPS", mode, pm.meps(), "MEPS", pm.ops, maxth);
                if (lazy) {
                    std::string m2 = "LAZY_total_B" + std::to_string(B) + (blocks_od ? "_bod" : "");
                    row(name, "insert_MEPS", m2, pi.meps_total(), "MEPS", pi.ops, maxth);
                    row(name, "delete_MEPS", m2, pd.meps_total(), "MEPS", pd.ops, maxth);
                    row(name, "mixed_MEPS", m2, pm.meps_total(), "MEPS", pm.ops, maxth);
                    row(name, "sync_s", "insert", pi.sync_s, "s", pi.ops, maxth);
                    row(name, "sync_s", "delete", pd.sync_s, "s", pd.ops, maxth);
                    row(name, "sync_s", "mixed", pm.sync_s, "s", pm.ops, maxth);
                }
                const char* cn[] = {"C1", "C2-dag", "C2-merge", "C3", "C4", "D-inter", "D-wcc-split", "D-intra-keep", "D-scc-split"};
                auto& cc = g->caseCounts();
                auto& cs = g->caseSeconds();
                std::cout << "  cases:";
                for (int c = 0; c < GIR::NUM_CASES; ++c) {
                    std::cout << " " << cn[c] << "=" << cc[c];
                    if (!batched && cc[c]) std::cout << "(" << 1000.0 * cs[c] / cc[c] << "ms)";
                    row(name, std::string("case_") + cn[c], mode, (double)cc[c], "count", cc[c], maxth);
                }
                std::cout << "\n  order paths (one-sided, hub-merge, re-search, cross-wcc):";
                for (long long x : g->order_stats) std::cout << " " << x;
                std::cout << "\n  re-ranked nodes (total, max, moves>10k):";
                for (long long x : g->order_moved) std::cout << " " << x;
                std::cout << "\n  block repair (quick, menger, chain-splits, fallback):";
                for (long long x : g->block_stats) std::cout << " " << x;
                std::cout << " carved " << g->carve_count << " (menger gave up: nopath " << g->menger_nopath << ", nocut " << g->menger_nocut
                          << ", budget " << g->menger_budget << ")";
                std::cout << "\n  split paths (fast, X-big, boundary-big, R-not-sc, block-cert-fail):";
                for (long long x : g->split_stats) std::cout << " " << x;
                std::cout << "  seconds:";
                for (double x : g->split_secs) std::cout << " " << x;
                std::cout << std::endl;
                delete g;
            }
        }

        if (do_queries) {
            double bs;
            std::vector<std::pair<int, int>> full = both(all);
            std::vector<std::pair<int, int>>().swap(all);
            GIR* g = build(V, full, bs);
            std::vector<std::pair<int, int>>().swap(full);
            row(name, "build_full", "-", bs, "s", 0, maxth);
            std::cout << "built full graph in " << bs << " s" << std::endl;
            // Sources: random vertices with at least one out-edge.
            std::mt19937 srng(7);
            std::vector<int> src;
            while ((int)src.size() < nsrc) {
                int s = (int)(srng() % V);
                if (!g->getDirectedOutNeighbors(s).empty()) src.push_back(s);
            }
            auto timeit = [&](const char* metric, auto fn) {
                auto t0 = Clock::now();
                for (int s : src) fn(s);
                double ms = secs(t0, Clock::now()) * 1000 / src.size();
                std::cout << "  " << metric << ": " << ms << " ms/source" << std::endl;
                return ms;
            };
            // Structural queries: u ~> v for random pairs, answered from the
            // index and by bidirectional BFS on the edge store (no index).
            if (nreach > 0) {
                std::mt19937 qrng(13);
                std::vector<std::pair<int, int>> qs(nreach);
                for (auto& q : qs) q = {(int)(qrng() % V), (int)(qrng() % V)};
                std::vector<char> ai(nreach), an(nreach);
                auto t0 = Clock::now();
                for (int i = 0; i < nreach; ++i) ai[i] = g->reachable(qs[i].first, qs[i].second);
                double us_i = secs(t0, Clock::now()) * 1e6 / nreach;
                auto t1 = Clock::now();
                for (int i = 0; i < nreach; ++i) an[i] = g->reachableNoIndex(qs[i].first, qs[i].second);
                double us_n = secs(t1, Clock::now()) * 1e6 / nreach;
                long long agree = 0, yes = 0;
                for (int i = 0; i < nreach; ++i) { agree += ai[i] == an[i]; yes += ai[i]; }
                auto t2 = Clock::now();
                long long aps = 0;
                for (int i = 0; i < nreach; ++i) aps += g->articulationPoint(qs[i].first);
                double us_ap = secs(t2, Clock::now()) * 1e6 / nreach;
                std::cout << "  reachability (" << nreach << " random pairs, " << yes << " reachable): index " << us_i
                          << " us/query, no index (bidirectional BFS) " << us_n << " us/query, answers agree: "
                          << agree << "/" << nreach << "; articulation point " << us_ap << " us/query" << std::endl;
                row(name, "reach_us", "index", us_i, "us", nreach, maxth);
                row(name, "reach_us", "no-index-bibfs", us_n, "us", nreach, maxth);
                row(name, "reach_agree", "-", (double)agree, "count", nreach, maxth);
                row(name, "artpoint_us", "index", us_ap, "us", nreach, maxth);
            }
            // Correctness spot check on two sources.
            bool ok = true;
            for (int i = 0; i < 2 && i < (int)src.size(); ++i) {
                auto want = GirBfsTest::sequentialBfs(*g, src[i], V);
                ok &= g->bfs(src[i]) == want && g->bfsTwoLevel(src[i]) == want;
                ok &= g->sssp(src[i]) == g->ssspTwoLevel(src[i]);
            }
            std::cout << "  spot check (bfs/bfsTwoLevel vs sequential, sssp vs ssspTwoLevel): " << (ok ? "OK" : "MISMATCH") << std::endl;
            row(name, "query_spot_check_ok", "-", ok, "bool", 2, maxth);
            {   // the baselines report the same two checks for the first source
                auto lv = g->bfs(src[0]);
                long long reached = 0;
                for (int x : lv) reached += x >= 0;
                auto d = g->sssp(src[0], sssp_delta);
                unsigned long long sum = 0;
                for (auto x : d) if (x != UINT64_MAX) sum += x;
                row(name, "bfs_reached_src0", "-", (double)reached, "vertices", 1, maxth);
                row(name, "sssp_distsum_src0", "-", (double)sum, "sum", 1, maxth);
            }

            std::vector<int> ths = threads.empty() ? std::vector<int>{maxth} : threads;
            for (int th : ths) {
                omp_set_num_threads(th);
                std::cout << "threads=" << th << std::endl;
                row(name, "bfs_ms", "direction-optimizing", timeit("bfs (direction-optimizing)", [&](int s) { g->bfs(s); }), "ms", nsrc, th);
                row(name, "bfs_ms", "two-level", timeit("bfs (two-level, paper Sec. 6)", [&](int s) { g->bfsTwoLevel(s); }), "ms", nsrc, th);
                if (th == ths.front() || th == maxth) {
                    row(name, "bfs_ms", "sequential", timeit("bfs (sequential baseline)", [&](int s) { GirBfsTest::sequentialBfs(*g, s, V); }), "ms", nsrc, th);
                    row(name, "sssp_ms", "delta-stepping", timeit("sssp (delta-stepping)", [&](int s) { g->sssp(s, sssp_delta); }), "ms", nsrc, th);
                    row(name, "sssp_ms", "two-level", timeit("sssp (two-level)", [&](int s) { g->ssspTwoLevel(s); }), "ms", nsrc, th);
                    std::vector<int> bsrc;
                    std::mt19937 brng(11);
                    for (int i = 0; i < nbc; ++i) bsrc.push_back((int)(brng() % V));
                    auto t0 = Clock::now();
                    g->bc(bsrc);
                    double ms = secs(t0, Clock::now()) * 1000;
                    std::cout << "  bc (" << nbc << " sources): " << ms << " ms" << std::endl;
                    row(name, "bc_ms_total", "brandes", ms, "ms", nbc, th);
                }
            }
            omp_set_num_threads(maxth);
            delete g;
        }
    }
    std::cout << "\nResults written to " << csv_path << std::endl;
    return 0;
}
