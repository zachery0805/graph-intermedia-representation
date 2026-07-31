#ifndef DYNAMIC_REPAIR_H
#define DYNAMIC_REPAIR_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <algorithm>

// Reusable "recompute-on-the-affected-region" primitives shared by GIR's
// initial construction (per-SCC BiCC, per-SCC/per-WCC condensation bridge
// info) and its incremental repair (Cases 2/3/5a/5b/5c). Both call sites
// hand in a small vertex subset plus a neighbor function scoped to whatever
// graph view they care about (the real graph + GIR's pending edge overlay);
// these algorithms don't know about GIR's data model at all.
//
// Iterative (not recursive) so repair on an unexpectedly large affected
// region can't blow the call stack.
namespace dynamic_repair {

// Tarjan's SCC algorithm, iterative. `out_neighbors(v)` should already be
// filtered to whatever universe the caller wants considered (e.g. only
// vertices within one SCC's induced subgraph, when re-checking a split).
inline std::vector<std::vector<int>> find_scc_scoped(
    const std::vector<int>& vertices,
    const std::function<std::vector<int>(int)>& out_neighbors) {
    std::unordered_map<int, int> index_of, low_of;
    std::unordered_set<int> on_stack;
    std::vector<int> node_stack;
    std::vector<std::vector<int>> result;
    int next_index = 0;

    struct Frame {
        int v;
        std::vector<int> nbrs;
        size_t next_child;
    };

    for (int start : vertices) {
        if (index_of.count(start)) continue;

        std::vector<Frame> call_stack;
        index_of[start] = low_of[start] = next_index++;
        node_stack.push_back(start);
        on_stack.insert(start);
        call_stack.push_back({start, out_neighbors(start), 0});

        while (!call_stack.empty()) {
            Frame& f = call_stack.back();
            if (f.next_child < f.nbrs.size()) {
                int w = f.nbrs[f.next_child++];
                if (!index_of.count(w)) {
                    index_of[w] = low_of[w] = next_index++;
                    node_stack.push_back(w);
                    on_stack.insert(w);
                    call_stack.push_back({w, out_neighbors(w), 0});
                } else if (on_stack.count(w)) {
                    low_of[f.v] = std::min(low_of[f.v], index_of[w]);
                }
            } else {
                int v = f.v;
                int v_low = low_of[v];
                call_stack.pop_back();
                if (!call_stack.empty()) {
                    int parent = call_stack.back().v;
                    low_of[parent] = std::min(low_of[parent], v_low);
                }
                if (v_low == index_of[v]) {
                    std::vector<int> component;
                    while (true) {
                        int w = node_stack.back();
                        node_stack.pop_back();
                        on_stack.erase(w);
                        component.push_back(w);
                        if (w == v) break;
                    }
                    result.push_back(std::move(component));
                }
            }
        }
    }
    return result;
}

struct BiccResult {
    std::vector<std::vector<int>> blocks;  // vertex sets, one per BiCC block
    std::unordered_set<int> articulation_points;
};

// Hopcroft-Tarjan biconnected components + articulation points, iterative.
// `neighbors(v)` must return a symmetric (undirected) adjacency scoped to
// whatever universe the caller wants considered. Vertices with no edges at
// all are not placed in any block here -- callers should add a singleton
// block for any vertex from `vertices` that isn't covered by the result
// (matching how a degree-0 vertex is still a trivial BiCC block of its own).
inline BiccResult find_bicc_scoped(
    const std::vector<int>& vertices,
    const std::function<std::vector<int>(int)>& neighbors) {
    BiccResult res;
    std::unordered_map<int, int> disc, low, parent;
    std::unordered_map<int, bool> parent_edge_used;
    int timer = 0;
    std::vector<std::pair<int, int>> edge_stack;

    struct Frame {
        int v;
        std::vector<int> nbrs;
        size_t next_child;
    };

    for (int root : vertices) {
        if (disc.count(root)) continue;

        std::vector<Frame> call_stack;
        disc[root] = low[root] = timer++;
        parent[root] = -1;
        parent_edge_used[root] = true;  // root has no parent edge to skip
        call_stack.push_back({root, neighbors(root), 0});
        int root_children = 0;

        while (!call_stack.empty()) {
            Frame& f = call_stack.back();
            int v = f.v;
            if (f.next_child < f.nbrs.size()) {
                int w = f.nbrs[f.next_child++];
                if (w == parent[v] && !parent_edge_used[v]) {
                    parent_edge_used[v] = true;
                    continue;
                }
                if (!disc.count(w)) {
                    disc[w] = low[w] = timer++;
                    parent[w] = v;
                    parent_edge_used[w] = false;
                    edge_stack.push_back({v, w});
                    if (v == root) root_children++;
                    call_stack.push_back({w, neighbors(w), 0});
                } else if (disc[w] < disc[v]) {
                    low[v] = std::min(low[v], disc[w]);
                    edge_stack.push_back({v, w});
                }
            } else {
                call_stack.pop_back();
                if (!call_stack.empty()) {
                    int p = call_stack.back().v;
                    low[p] = std::min(low[p], low[v]);
                    if (low[v] >= disc[p]) {
                        std::unordered_set<int> block_set;
                        while (!edge_stack.empty()) {
                            auto e = edge_stack.back();
                            edge_stack.pop_back();
                            block_set.insert(e.first);
                            block_set.insert(e.second);
                            if (e.first == p && e.second == v) break;
                        }
                        if (block_set.size() >= 2) {
                            res.blocks.emplace_back(block_set.begin(), block_set.end());
                        }
                        if (p != root) {
                            res.articulation_points.insert(p);
                        }
                    }
                }
            }
        }
        if (root_children > 1) {
            res.articulation_points.insert(root);
        }
    }
    return res;
}

}  // namespace dynamic_repair

#endif  // DYNAMIC_REPAIR_H
