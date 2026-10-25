// Correctness self-test for GIR's incremental maintenance. Generates small
// random directed graphs of several shapes, then for each one runs:
//   - single-edge differential test (insertEdge/deleteEdge, full rebuild
//     comparison after every update),
//   - batch differential test (applyBatch with duplicates and cancelling
//     pairs),
//   - BFS differential test against a plain sequential BFS,
// plus the hand-written Case 1/5b/5c check on data/example_graph.
//
// Usage: ./gir_selftest <scratch_dir> [example_graph_dir]
// Exit status is non-zero if anything mismatched.
#include "gir.h"
#include "gir_test.cpp"

#include <cstdio>
#include <fstream>
#include <random>
#include <set>
#include <sys/stat.h>

static void writeInts(const std::string& path, const std::vector<int>& v) {
    std::ofstream f(path, std::ios::binary);
    f.write(reinterpret_cast<const char*>(v.data()), (std::streamsize)(v.size() * sizeof(int)));
}

static void writeGraph(const std::string& dir, int n, const std::set<std::pair<int, int>>& edges) {
    mkdir(dir.c_str(), 0755);
    std::vector<std::vector<int>> out(n), in(n);
    for (auto& e : edges) { out[e.first].push_back(e.second); in[e.second].push_back(e.first); }
    auto dump = [&](const std::vector<std::vector<int>>& adj, const std::string& pre) {
        std::vector<int> beg{0}, csr;
        for (auto& a : adj) { csr.insert(csr.end(), a.begin(), a.end()); beg.push_back((int)csr.size()); }
        writeInts(dir + pre + "_begin.bin", beg);
        writeInts(dir + pre + "_adjacent.bin", csr);
    };
    dump(out, "fw");
    dump(in, "bw");
}

static std::set<std::pair<int, int>> randomEdges(std::mt19937& rng, int n, int m) {
    std::set<std::pair<int, int>> e;
    std::uniform_int_distribution<int> d(0, n - 1);
    while ((int)e.size() < m) { int u = d(rng), v = d(rng); if (u != v) e.insert({u, v}); }
    return e;
}

// Directed cycles of random length, chained by single edges and sprinkled
// with antiparallel pairs: lots of articulation points, bridges and
// SCCs that split on a single deletion.
static std::set<std::pair<int, int>> cycleChain(std::mt19937& rng, int n) {
    std::set<std::pair<int, int>> e;
    int start = 0;
    std::uniform_int_distribution<int> len(2, 6);
    while (start < n) {
        int L = std::min(len(rng), n - start);
        for (int i = 0; i < L; ++i) if (L > 1) e.insert({start + i, start + (i + 1) % L});
        if (start > 0) {
            int a = start - 1 - (int)(rng() % 3), b = start;
            if (a >= 0) { e.insert({a, b}); if (rng() % 3 == 0) e.insert({b, a}); }
        }
        start += L;
    }
    return e;
}

// Directed cycles glued at shared vertices (one big SCC with many blocks)
// plus random chords. Deleting a chord keeps the SCC but splits blocks, the
// case the Menger chain-splitting repair handles.
static std::set<std::pair<int, int>> cactusWithChords(std::mt19937& rng, int n, int chords) {
    std::set<std::pair<int, int>> e;
    std::vector<int> placed{0};
    int next = 1;
    std::uniform_int_distribution<int> len(3, 7);
    while (next < n) {
        int anchor = placed[rng() % placed.size()];
        int L = std::min(len(rng), n - next + 1);
        std::vector<int> cyc{anchor};
        for (int i = 1; i < L; ++i) { cyc.push_back(next); placed.push_back(next); ++next; }
        for (size_t i = 0; i < cyc.size(); ++i) e.insert({cyc[i], cyc[(i + 1) % cyc.size()]});
    }
    for (int i = 0; i < chords; ++i) {
        int u = (int)(rng() % n), v = (int)(rng() % n);
        if (u != v) e.insert({u, v});
    }
    return e;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: " << argv[0] << " <scratch_dir> [example_graph_dir]\n"; return 2; }
    std::string root = std::string(argv[1]) + "/";
    mkdir(root.c_str(), 0755);
    std::mt19937 rng(12345);

    struct Spec { std::string name; int n; std::set<std::pair<int, int>> edges; };
    std::vector<Spec> specs;
    specs.push_back({"sparse_60", 60, randomEdges(rng, 60, 80)});
    specs.push_back({"dense_40", 40, randomEdges(rng, 40, 160)});
    specs.push_back({"cycles_80", 80, cycleChain(rng, 80)});
    specs.push_back({"medium_200", 200, randomEdges(rng, 200, 320)});
    specs.push_back({"empty_30", 30, {}});
    // One large well-connected SCC: deletions mostly keep it intact, so the
    // block-repair shortcuts actually fire.
    specs.push_back({"dense_120", 120, randomEdges(rng, 120, 600)});
    specs.push_back({"cactus_150", 150, cactusWithChords(rng, 150, 40)});

    int failures = 0;
    for (int fast = 0; fast < 2; ++fast)
    for (size_t si = 0; si < specs.size(); ++si) {
        auto& s = specs[si];
        std::string dir = root + s.name + "/";
        writeGraph(dir, s.n, s.edges);
        std::cout << "\n##### " << s.name << " (" << s.n << " vertices, " << s.edges.size() << " edges)"
                  << (fast ? " [block fast path forced on]" : "") << std::endl;
        {
            GIR g(dir, dir); g.buildWCC(); g.buildSCC_DAG(); g.buildBCT();
            if (fast) { g.setFastPathMinBlock(0); g.setHubDegree(2); g.setForceExactBlockCert(true); g.setForceResplit(si % 2 == 1); g.setSearchHubDegree(3); }
            failures += EdgeUpdateTest::runDifferentialTest(g, 4000, 7) != 0;
            std::cout << "block repair (quick, menger, chain-splits, fallback):";
            for (long long x : g.block_stats) std::cout << " " << x;
            std::cout << " | carved " << g.carve_count << " | menger -1 reasons nopath/nocut/budget: " << g.menger_nopath << "/" << g.menger_nocut << "/" << g.menger_budget << std::endl;
            std::cout << "order paths (one-sided, hub-merge, re-search, cross-wcc):";
            for (long long x : g.order_stats) std::cout << " " << x;
            std::cout << std::endl;
            auto& cc = g.caseCounts();
            std::cout << "case counts:";
            for (long long x : cc) std::cout << " " << x;
            std::cout << std::endl;
            failures += !GirBfsTest::runDifferentialBfsTest(g, std::min(20, s.n), 3);
            // Queries on the mutated graph, with every granularity path forced:
            // default thresholds, then tiny ones (bottom-up BFS, data-parallel
            // blocks and SCCs, parallel frontiers all on).
            int qb = QueryTest::run(g, 8, 5);
            g.setQueryGranularity(1, 2);
            qb += QueryTest::run(g, 8, 6);
            qb += QueryTest::interleaved(g, 300, 9);
            qb += IncTest::run(g, 150, 31);
            std::cout << "QueryTest: " << (qb ? "MISMATCH" : "all BFS/SSSP/BC results match") << std::endl;
            failures += qb != 0;
        }
        {
            GIR g(dir, dir); g.buildWCC(); g.buildSCC_DAG(); g.buildBCT();
            if (fast) { g.setFastPathMinBlock(0); g.setHubDegree(2); g.setSearchHubDegree(3); g.setHubAnchorMin(2); }
            failures += EdgeUpdateTest::runBatchDifferentialTest(g, 300, 25, 11) != 0;
        }
        {
            GIR g(dir, dir); g.buildWCC(); g.buildSCC_DAG(); g.buildBCT();
            if (fast) { g.setFastPathMinBlock(0); g.setHubDegree(2); g.setSearchHubDegree(3); g.setHubAnchorMin(2); g.setForceAnchorSet(true); }
            failures += EdgeUpdateTest::runBatchDifferentialTest(g, 300, 40, 37) != 0;
            std::cout << "anchor-set rounds " << g.anchor_set_calls << ", hub rounds " << g.hub_anchor_calls << std::endl;
        }
        {
            GIR g(dir, dir); g.buildWCC(); g.buildSCC_DAG(); g.buildBCT();
            if (fast) { g.setFastPathMinBlock(0); g.setHubDegree(2); g.setSearchHubDegree(3); g.setHubAnchorMin(2); g.setForceAnchorSet(true); }
            failures += EdgeUpdateTest::runLazyDifferentialTest(g, 150, 25, 13) != 0;
        }
        // Blocks on demand: single-edge, batched and deferred updates, with
        // the blocks rebuilt before every check.
        {
            GIR g(dir, dir); g.buildWCC(); g.buildSCC_DAG(); g.buildBCT();
            if (fast) { g.setFastPathMinBlock(0); g.setHubDegree(2); g.setSearchHubDegree(3); g.setHubAnchorMin(2); g.setForceAnchorSet(si % 2 == 0); }
            g.setBlocksOnDemand(true);
            g.setParallelBlocksMin(2);  // rebuilds go through the parallel (Tarjan-Vishkin) path
            int bad = EdgeUpdateTest::runDifferentialTest(g, 1500, 17);
            bad += EdgeUpdateTest::runBatchDifferentialTest(g, 100, 25, 19);
            bad += EdgeUpdateTest::runLazyDifferentialTest(g, 100, 25, 23);
            std::cout << "Blocks-on-demand tests: " << (bad ? "MISMATCH" : "OK") << std::endl;
            failures += bad != 0;
        }
    }
    // A graph big enough for the batch split certificate's probe rounds to
    // run (on the graphs above the trimming gives up first): one giant SCC,
    // batches of mostly deletions, with each kind of anchor.
    {
        const int n = 5000;
        std::string dir = root + "giant_5000/";
        auto ge = randomEdges(rng, n, 6 * n);
        for (int k = 0; k < 300; ++k) {  // vertex 0 a hub (enough for the single-hub anchors)
            int w = 1 + (int)(rng() % (n - 1));
            ge.insert({0, w});
            ge.insert({w, 0});
        }
        writeGraph(dir, n, ge);
        std::cout << "\n##### giant_5000 (" << n << " vertices) [bulk split certificate]" << std::endl;
        for (int mode = 0; mode < 5; ++mode) {  // 0 local anchor, 1 hubs, 2 anchor set; 3, 4: 1, 2 with the parallel second pass
            GIR g(dir, dir); g.buildWCC(); g.buildSCC_DAG(); g.buildBCT();
            g.setFastPathMinBlock(0);
            g.setBlocksOnDemand(mode != 0);
            if (mode >= 1) g.setHubAnchorMin(2);
            if (mode == 2 || mode == 4) g.setForceAnchorSet(true);
            if (mode >= 3) g.setSeqProbeBudget(3);
            int bad = EdgeUpdateTest::runBatchDifferentialTest(g, 25, 400, 41 + mode);
            std::cout << "  mode " << mode << ": anchor-set rounds " << g.anchor_set_calls << ", hub rounds " << g.hub_anchor_calls
                      << (bad ? " MISMATCH" : " OK") << std::endl;
            failures += bad != 0;
        }
    }
    if (argc >= 3) {
        std::string ex = std::string(argv[2]) + "/";
        GIR g(ex, ex); g.buildWCC(); g.buildSCC_DAG(); g.buildBCT();
        EdgeUpdateTest::runTargetedSplitTest(g);
        GIR g2(ex, ex); g2.buildWCC(); g2.buildSCC_DAG(); g2.buildBCT();
        g2.deleteEdge(0, 1);  // non-bridge edge of the 3-cycle {0,1,2}: must split the SCC
        bool ok = g2.sccOf(0) != g2.sccOf(1) && EdgeUpdateTest::checkConsistency(g2, g2.vertexCount(), "delete 0->1");
        std::cout << "Regression (delete non-bridge 0->1 splits SCC): " << (ok ? "OK" : "FAILED") << std::endl;
        failures += !ok;
        GIR g3(ex, ex); g3.buildWCC(); g3.buildSCC_DAG(); g3.buildBCT();
        g3.insertEdge(1, 0); g3.deleteEdge(0, 1);
        ok = g3.hasUndirectedEdge(0, 1) && EdgeUpdateTest::checkConsistency(g3, g3.vertexCount(), "reverse kept");
        std::cout << "Regression (undirected {0,1} survives while 1->0 exists): " << (ok ? "OK" : "FAILED") << std::endl;
        failures += !ok;
        // Insertions only (delta edges, no tombstone bitmaps yet), then
        // bottom-up BFS steps: must not read the unallocated bitmaps.
        GIR g4(ex, ex); g4.buildWCC(); g4.buildSCC_DAG(); g4.buildBCT();
        for (int a = 0; a < g4.vertexCount(); ++a) g4.insertEdge(a, (a * 7 + 3) % g4.vertexCount());
        g4.setQueryGranularity(1, 2);
        ok = QueryTest::run(g4, 8, 4) == 0;
        std::cout << "Regression (queries after insertions only): " << (ok ? "OK" : "FAILED") << std::endl;
        failures += !ok;
        // Regression: both directions of 15<->299 deleted in one batch; the
        // SCC splits and the block {15,299} hangs off the departing
        // articulation point 15 only through the deleted edge.
        {
            // (Vertex ids as found by the stress test: they decide how the
            // block-cut forest is rooted, and the bug needs {15,299} to be a
            // child block of 15.)
            std::vector<std::pair<int, int>> E = {{15, 299}, {299, 15}, {299, 178}, {178, 299}, {15, 296}, {296, 15}};
            GIR g5(graph::fromEdges(300, E)); g5.buildWCC(); g5.buildSCC_DAG(); g5.buildBCT();
            g5.setFastPathMinBlock(0);
            g5.applyBatch({{15, 299, false}, {299, 15, false}});
            ok = EdgeUpdateTest::checkConsistency(g5, 300, "two-way bridge deletion");
            std::cout << "Regression (batch deletes both directions of a bridge at a departing AP): " << (ok ? "OK" : "FAILED") << std::endl;
            failures += !ok;
        }
    }
    std::cout << "\nSELFTEST " << (failures ? "FAILED" : "PASSED") << " (" << failures << " failing checks)" << std::endl;
    return failures ? 1 : 0;
}
