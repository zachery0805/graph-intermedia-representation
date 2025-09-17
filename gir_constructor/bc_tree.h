#ifndef BLOCK_CUT_TREE_H
#define BLOCK_CUT_TREE_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <unordered_set>

#include <GraphRW.h>
#include <KSPGraph.h>
#include <PeeKAdaptiveWithEdgeSwap.h>
#include "graph.h"
#include "graph_undirected.h"
#include <limits>
#include <algorithm>
#include <omp.h>
#include <optional>

#ifdef KSP_PARALLEL_DEVIATIONS_L1
#define KSP_PD_L1 true
#else
#define KSP_PD_L1 false
#endif

// only support Num Threads > 1
#ifndef KSP_NUM_THREADS
#define KSP_NUM_THREADS 4
#endif

#if defined(KSP_PARALLEL_DEVIATIONS_L2) || defined(KSP_PARALLEL_DEVIATIONS_L1)
#define KSP_PD true
#else
#define KSP_PD false
#endif

#define NIL -1 \

#define K 3

using GraphType = BasicGraph<true, true>;
using IndexTuple = std::vector<int>;

// hash function
struct UnorderedPairHash {
    size_t operator()(const std::pair<int, int>& p) const {
        // order-independent hash: ensure (u, v) == (v, u)
        int a = std::min(p.first, p.second);
        int b = std::max(p.first, p.second);
        return std::hash<int>{}(a) ^ (std::hash<int>{}(b) << 1);
    }
};

// hash equal function
struct UnorderedPairEqual {
    bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
        return (a.first == b.first && a.second == b.second) ||
               (a.first == b.second && a.second == b.first);
    }
};

class BlockCutTree {
public:

    std::vector<std::vector<int>> tree_adj;

    // store edges that not included in sub blocks or pruned
    std::unordered_map<std::pair<index_t, index_t>, path_t, UnorderedPairHash, UnorderedPairEqual> bct_independent_edges;

    graph_undirected *g; // original graph

    bool* ap; // ap list

    std::unordered_map<int, int> node_to_ap;         // node -> AP
    std::unordered_map<int, int> ap_to_node;         // AP -> node ID in BCT

    std::vector<std::vector<int>> bicc_list;         // Each BiCC's vertex list
    std::unordered_map<int, int> node_to_block;      // node -> block

    std::unordered_map<index_t, std::shared_ptr<BasicGraph<true, true>>> sub_graph_mem; // block num -> sub graph (memory)
    std::unordered_map<std::pair<index_t, index_t>, std::vector<Path>, UnorderedPairHash, UnorderedPairEqual> sub_ksp_mem; // pair(ap1, ap2) -> sub ksp results (memory)
    std::unordered_map<std::pair<index_t, index_t>, bool, UnorderedPairHash, UnorderedPairEqual> sub_ksp_complete_mem; // track if the sub_ksp is completed in merge mode

    void validate_node_to_block(int total_nodes) {
        for (int i = 0; i < total_nodes; ++i) {
            if (node_to_block.find(i) == node_to_block.end()) {
                std::cerr << "Missing node_to_block mapping for node " << i << std::endl;
            }
        }
        std::cout << "Passed build test" << std::endl;
    }

    // build BC tree
    void build(const std::vector<std::vector<int>>& biccs, bool* ap, graph_undirected *g) {

        bicc_list = biccs;
        
        this->g = g;
        this->ap = ap;

        vertex_t V = this->g->vert_count;

        int m = biccs.size();
        
        int cur_node_id = m;

        // Map each articulation point to a new BCT node ID
        for (vertex_t v = 0; v < V; ++v) {
            if (ap[v]) {
                node_to_block[v] = cur_node_id; 
                ap_to_node[v] = cur_node_id;
                node_to_ap[cur_node_id] = v;
                cur_node_id++;
            }
        }       

        tree_adj.resize(cur_node_id);

        // Add BiCC <-> AP edges
        for (int i = 0; i < bicc_list.size(); ++i) {

            // block with size = 2
            if(biccs[i].size() == 2){
                
                // two node ids
                index_t node_id_1 = bicc_list[i][0];
                index_t node_id_2 = bicc_list[i][1];

                // add the weight to independent edges weight map
                path_t weight = get_edge_weight(node_id_1, node_id_2);
                bct_independent_edges[{node_id_1, node_id_2}] = weight;

                // std::cout << "Node weight stored to independent edges: " << node_id_1 << " " << node_id_2 << " " << bct_independent_edges[{node_id_2, node_id_1}] << std::endl;

                // first node is ap
                if(ap[node_id_1]){

                    // second node is also ap (bridge block): no need to this block to exist,
                    // directly connect those two aps
                    if(ap[node_id_2]){

                        int ap_id_1 = ap_to_node[node_id_1];
                        int ap_id_2 = ap_to_node[node_id_2];

                        tree_adj[ap_id_1].push_back(ap_id_2);
                        tree_adj[ap_id_2].push_back(ap_id_1);

                    }
                    // only first node is ap, prune the first node, connect block to ap
                    else{

                        int ap_id = ap_to_node[node_id_1];
                        bicc_list[i].erase(bicc_list[i].begin());
                        tree_adj[i].push_back(ap_id);
                        tree_adj[ap_id].push_back(i);

                        node_to_block[node_id_2] = i;
                        
                    }
                }
                // if first node is not ap, the second node must be ap, prune second ap
                else{

                    int ap_id = ap_to_node[node_id_2];
                    bicc_list[i].erase(bicc_list[i].begin() + 1);
                    tree_adj[i].push_back(ap_id);
                    tree_adj[ap_id].push_back(i);

                    node_to_block[node_id_1] = i;

                }

            }
            // general case (nodes >= 3)
            else{
                for (int v : biccs[i]) {
                    if (ap[v]) {
                        int ap_id = ap_to_node[v];
                        tree_adj[i].push_back(ap_id);
                        tree_adj[ap_id].push_back(i);
                    }
                    else{
                        node_to_block[v] = i;
                    }  
                }
            }          
        }
        print();

        // validate_node_to_block(g->vert_count);
    }  

    // get edge weight from two nodes in CSR format
    path_t get_edge_weight(int u, int v) {
        for (index_t i = g->beg_pos[u]; i < g->beg_pos[u + 1]; ++i) {
            if (g->csr[i] == v) {
                return g->weight[i];
            }
        }
        // edge doesn't exist
        return -1;
    }

    void print() const {
        for (int i = 0; i < tree_adj.size(); ++i) {
            std::cout << "BCT Node " << i << ": ";
            for (int nei : tree_adj[i]) {
                std::cout << nei << " ";
            }

            std::cout << "\n";

            // block node
            if(i < bicc_list.size()){
                std::cout << "Block node elements: ";

                for(int j = 0; j < bicc_list[i].size(); ++j){
                    std::cout << bicc_list[i][j] << " ";
                }

                std::cout << "\n\n";

            }
            // ap_node
            else{
                std::cout << "AP node in original graph: " << node_to_ap.at(i) << "\n\n";
            }
        }
    }

    // find bct path
    std::vector<int> find_path_bct(int source_bct, int target_bct) {
        std::vector<bool> visited(tree_adj.size(), false);
        std::unordered_map<int, int> parent;
        std::queue<int> q;

        visited[source_bct] = true;
        q.push(source_bct);

        // BFS
        while (!q.empty()) {
            int curr = q.front();
            q.pop();

            if (curr == target_bct) break;

            for (int neighbor : tree_adj[curr]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    parent[neighbor] = curr;
                    q.push(neighbor);
                }
            }
        }

        std::vector<int> path;
        if (!visited[target_bct]) {
            // std::cout << "No path found from " << source_bct << " to " << target_bct << ".\n";
            return path;  // empty path
        }

        // build path from target
        for (int v = target_bct; v != source_bct; v = parent[v]) {
            path.push_back(v);
        }
        path.push_back(source_bct);
        std::reverse(path.begin(), path.end());

        // print path
        // std::cout << "BCT Path from " << source_bct << " to " << target_bct << ": ";
        // for (int node : path) {
        //     std::cout << node << " ";
        // }
        // std::cout << "\n";

        return path;
    }

    // generate subgraph for each block in the path, store it for future use
    std::shared_ptr<BasicGraph<true, true>> gen_subgraph(std::vector<int> bicc, int bicc_index, int& src, int& tar){

        // change bicc to map
        std::unordered_map<int, int> bicc_index_map;
        for (int i = 0; i < bicc.size(); ++i) {
            bicc_index_map[bicc[i]] = i;
        }

        // new sub_graph object
        // graph_undirected subgraph;

        // CSR for subgraph
        std::vector<index_t> new_beg_pos;
        std::vector<vertex_t> new_csr;
        std::vector<path_t> new_weight;

        bool found_src = false;
        bool found_tar = false;

        // first element is 0 for beg
        new_beg_pos.push_back(0); 

        for (int n = 0; n < bicc.size(); n++) {

            int deg = 0;
            int node = bicc[n];

            // if source, store its position
            if(node == src && !found_src){
                src = n;
                found_src = true;
            }
            // if target, store its position
            else if(node == tar && !found_tar){
                tar = n;
                found_tar = false;
            }

            // is ap -> should delete edges that out of bicc
            if (ap[node]) {
                for (index_t i = g->beg_pos[node]; i < g->beg_pos[node + 1]; ++i) {
                    int neighbor = g->csr[i];
                    if (bicc_index_map.count(neighbor)) {
                        new_csr.push_back(bicc_index_map[neighbor]);
                        new_weight.push_back(g->weight[i]);
                        ++deg;
                    }
                }
            }
            // not ap -> store all edges 
            else {
                for (index_t i = g->beg_pos[node]; i < g->beg_pos[node + 1]; ++i) {
                    int neighbor = g->csr[i];
                    new_csr.push_back(bicc_index_map[neighbor]);
                    new_weight.push_back(g->weight[i]);
                    ++deg;
                }
            }

            new_beg_pos.push_back(new_beg_pos.back() + deg);
        }

        // subgraph.vert_count = new_beg_pos.size() - 1;
        // subgraph.edge_count = new_csr.size();

        // subgraph.beg_pos = new index_t[subgraph.vert_count + 1];
        // subgraph.csr = new vertex_t[subgraph.edge_count];
        // subgraph.weight = new path_t[subgraph.edge_count];

        // // copy contents
        // std::copy(new_beg_pos.begin(), new_beg_pos.end(), subgraph.beg_pos);
        // std::copy(new_csr.begin(), new_csr.end(), subgraph.csr);
        // std::copy(new_weight.begin(), new_weight.end(), subgraph.weight);

        // std::cout << "subgraph.vert_count = " << subgraph.vert_count << "\n";
        // std::cout << "subgraph.edge_count = " << subgraph.edge_count << "\n";
        // std::cout << "subgraph csr size = " << new_csr.size() << "\n";
        // std::cout << "subgraph weight size = " << new_weight.size() << "\n";

        // std::cout << "test gen_subgraph" << std::endl;
        // std::cout << "1st element of csr: " << subgraph.csr[0] << std::endl;

        // // beg_pos
        // std::cout << "subgraph.beg_pos:\n";
        // for (index_t i = 0; i <= subgraph.vert_count; ++i) {
        //     std::cout << i << ": " << subgraph.beg_pos[i] << "\n";
        // }

        // // csr & weight
        // std::cout << "subgraph edges (csr + weight):\n";
        // for (index_t i = 0; i < subgraph.vert_count; ++i) {
        //     for (index_t j = subgraph.beg_pos[i]; j < subgraph.beg_pos[i + 1]; ++j) {
        //         std::cout << "from " << i << " to " << subgraph.csr[j]
        //                 << " weight = " << subgraph.weight[j] << "\n";
        //     }
        // }

        const char* delta = "0.1";

        auto G = std::make_shared<BasicGraph<true, true>>(
            GraphRW::read_graph_bct<true, true, true>(
                new_beg_pos.data(), new_csr.data(), new_weight.data(),
                new_beg_pos.size() - 1, new_csr.size(), delta));

        return G;
    }

    // compute KSP for each sub graph
    std::vector<Path> compute_ksp(const int src, const int dest, const int k, const int bicc_index, BasicGraph<true, true>& G, path_t gap, bool& early_stopped){

        const char* delta = "0.1";

        // std::cout << "1st element of csr: " << g.csr[0] << std::endl;

        PeeKAdaptiveWithEdgeSwap<DeltaSteppingStatic, GraphType, KSP_NUM_THREADS, KSP_PD> peek(G, k);

        // double start_time = wtime();
        auto paths = peek.compute(src, dest, gap, early_stopped);
        // double end_time = wtime();

        // std::cout << "KSP src=" << src << ", dest=" << dest << ", k=" << k
        //         << ", time=" << (end_time - start_time) << "s\n";

        // std::cout << "Found " << paths.size() << " paths from " << src << " to " << dest << ":\n";

        for(auto& i : paths){
            // std::cout << "path length: " << i.length << std::endl;
            // std::cout << "path:";
            for(unsigned int j = 0; j < i.p.size(); j++){

                i.p[j] = bicc_list[bicc_index][i.p[j]];
                // std::cout << " " << i.p[j];
            }
            // std::cout << std::endl;               
        }

        return paths;
    }
    // joint sub paths using min-heap (base method)
    std::vector<std::vector<Path>> joint_ksp_minheap(std::vector<std::vector<Path>> &path_lists){

        using IndexTuple = std::vector<int>;

        struct HeapEntry {
            w_type total_len;
            IndexTuple indices;

            bool operator>(const HeapEntry& other) const {
                return total_len > other.total_len;
            }
        };

        const int num_blocks = path_lists.size();
        std::vector<std::vector<Path>> result;
        std::priority_queue<HeapEntry, std::vector<HeapEntry>, std::greater<HeapEntry>> min_heap;
        std::set<IndexTuple> visited;

        // initialization
        // shortest path (all 0)
        IndexTuple init_idx(num_blocks, 0);
        w_type init_len = 0;
        for (int i = 0; i < num_blocks; ++i)
            init_len += path_lists[i][0].length;

        min_heap.push({init_len, init_idx});
        visited.insert(init_idx);

        while (!min_heap.empty() && result.size() < K) {
            auto [cur_len, indices] = min_heap.top();
            min_heap.pop();

            std::vector<Path> current_combo;
            for (int i = 0; i < num_blocks; ++i) {
                current_combo.push_back(path_lists[i][indices[i]]);
            }
            result.push_back(std::move(current_combo));

            for (int i = 0; i < num_blocks; ++i) {
                IndexTuple next_idx = indices;
                next_idx[i]++;

                if (next_idx[i] < path_lists[i].size() && visited.find(next_idx) == visited.end()) {
                    w_type new_len = new_len = cur_len - path_lists[i][indices[i]].length + path_lists[i][next_idx[i]].length;
                    min_heap.push({new_len, next_idx});
                    visited.insert(next_idx);
                }
            }
        }

        return result;
    }

    // prune unnecessary paths
    void joint_ksp_minheap_prune(std::vector<std::vector<Path>> &path_lists){
        
        int combinations = 1;
        std::vector<path_t> gap_list;

        // initialization
        for(int i = 0; i < path_lists.size(); i++){
            // only one path found in a block, no increase for combinations, gap is 0
            if(path_lists[i].size() == 1){
                gap_list.push_back(0.0);
            }
            // general case
            else{
                int paths_num = path_lists[i].size();
                combinations *= paths_num;
                path_t gap = path_lists[i][paths_num - 1].length - path_lists[i][0].length;
                gap_list.push_back(gap);
            }     
        }

        // when there're no more combinations than required K, no need to prune
        if(combinations <= static_cast<double>(K)){
            return;
        }

        while(true){
            
            // find the max gap
            auto it = std::max_element(gap_list.begin(), gap_list.end());
            int max_index = std::distance(gap_list.begin(), it);

            // prune the useless paths
            int sub_size = path_lists[max_index].size();
            int temp_combinations = (combinations / sub_size) * (sub_size - 1);

            // can prune last path
            if(temp_combinations >= static_cast<double>(K) && sub_size > 1){
                // prune to only one path, gap becomes 0
                if(sub_size == 2){
                    gap_list[max_index] = 0.0;
                }
                // general case
                else{
                    gap_list[max_index] = path_lists[max_index][sub_size - 2].length - path_lists[max_index][0].length;
                }
                
                // prune path and update combinations
                path_lists[max_index].pop_back();
                combinations = temp_combinations;
            }
            // can't prune anymore, return
            else{
                break;
            }
        }

        return;
    }

    // given the sub ksp reults, combine them to the final result
    std::vector<Path> combine_ksp_result(std::vector<std::vector<Path>> &min_paths, std::vector<int> &skeleton_path, std::vector<int> &insertion_loc, path_t skeleton_weight){
        
        int num_blocks = insertion_loc.size();

        // std::cout << "num_blocks: " << std::endl;

        std::vector<Path> result_paths;

        for (int k = 0; k < min_paths.size(); ++k) {

            std::vector<NODE_ID> full_path = skeleton_path;

            // insert from the end to beginning
            for (int b = num_blocks - 1; b >= 0; --b) {

                const std::vector<NODE_ID> sub = min_paths[k][b].p;

                // std::cout << "insertion_loc[b]: " << insertion_loc[b] << std::endl;

                // first in bct path is block
                if(insertion_loc[b] == 0){
                    std::vector<NODE_ID> trimmed_sub(sub.begin(), sub.end() - 1);
                    int insert_pos = insertion_loc[b];
                    full_path.insert(full_path.begin(), trimmed_sub.begin(), trimmed_sub.end());
                }
                // last in bct path is block
                else if(insertion_loc[b] == skeleton_path.size()){
                    std::vector<NODE_ID> trimmed_sub(sub.begin() + 1, sub.end());
                    int insert_pos = insertion_loc[b];
                    full_path.insert(full_path.end(), trimmed_sub.begin(), trimmed_sub.end());
                }
                // general case
                else{
                    std::vector<NODE_ID> trimmed_sub(sub.begin() + 1, sub.end() - 1);
                    int insert_pos = insertion_loc[b];
                    full_path.insert(full_path.begin() + insert_pos, trimmed_sub.begin(), trimmed_sub.end());
                }
            }

            // construct complete paths
            Path combined;
            combined.parent_index = 0;
            combined.alpha = 0;
            combined.length = skeleton_weight;
            for (int b = 0; b < num_blocks; ++b) {
                combined.length += min_paths[k][b].length;
            }
            combined.p = std::move(full_path);
            result_paths.push_back(std::move(combined));
        }

        return result_paths;
    }

    // given two sub ksp reults, combine them to the merged result
    void combine_ksp_merge(std::vector<Path>& pre_result, std::vector<Path>& cur_result, bool last, path_t& gap){

        std::vector<std::vector<Path>> temp_path_list;

        temp_path_list.push_back(pre_result);
        temp_path_list.push_back(cur_result);

        std::vector<std::vector<Path>> temp_combo = joint_ksp_minheap(temp_path_list);

        std::vector<Path> result_paths;

        for (int k = 0; k < temp_combo.size(); ++k) {

            int insertion_loc = temp_combo[k][0].p.size();

            std::vector<NODE_ID> full_path = temp_combo[k][0].p;

            const std::vector<NODE_ID> sub = temp_combo[k][1].p;
            // std::cout << "insertion_loc[b]: " << insertion_loc[b] << std::endl;

            // no need to handle first block case here
            // last in bct path is block
            if(last){
                std::vector<NODE_ID> trimmed_sub(sub.begin() + 1, sub.end());
                full_path.insert(full_path.end(), trimmed_sub.begin(), trimmed_sub.end());
            }
            // general case
            else{
                std::vector<NODE_ID> trimmed_sub(sub.begin() + 1, sub.end() - 1);
                full_path.insert(full_path.begin() + insertion_loc, trimmed_sub.begin(), trimmed_sub.end());
            }

            // construct complete paths
            Path combined;
            combined.parent_index = 0;
            combined.alpha = 0;
            combined.length = temp_combo[k][0].length + temp_combo[k][1].length;
            combined.p = std::move(full_path);
            result_paths.push_back(std::move(combined));
        }

        gap = result_paths[result_paths.size() - 1].length - result_paths[0].length;

        // update pre_result to the merged result
        pre_result = result_paths;
    }

    // add skeleton path to sub_ksp
    void add_skeleton_merge(std::vector<int> skeleton_path, path_t skeleton_weight, std::vector<Path>& paths, path_t & gap, bool last){
        
        int size = paths.size();
    
        if(ap[skeleton_path[skeleton_path.size() - 1]]){
            // pop the duplicate articulation points
            skeleton_path.pop_back();
        }
        
        for(int i = 0; i < size; i++){

            if(!last){
                paths[i].p.pop_back();
            }
            
            paths[i].p.insert(paths[i].p.begin(), skeleton_path.begin(), skeleton_path.end());
            paths[i].length += skeleton_weight;
        }

        gap = paths[size - 1].length - paths[0].length;
    }
    

    void find_ksp(int source, int target, bool prune, bool merge){

        int source_bct, target_bct;

        // merge gap
        path_t inf_gap = std::numeric_limits<path_t>::infinity();
        path_t gap = std::numeric_limits<path_t>::infinity();

        bool early_stopped = false;

        std::vector<Path> merge_path;

        // store sub shortest paths
        std::vector<std::vector<Path>> path_lists;

        // source
        // if(ap_to_node.count(source)){
        //     source_bct = ap_to_node[source];
        // }
        // else{
        //     for(int i = 0; i < bicc_list.size(); i++){
        //         if(std::find(bicc_list[i].begin(), bicc_list[i].end(), source) != bicc_list[i].end()){
        //             source_bct = i;
        //         }
        //     }
        // }

        // target
        // if(ap_to_node.count(target)){
        //     target_bct = ap_to_node[target];
        // }
        // else{
        //     for(int i = 0; i < bicc_list.size(); i++){
        //         if(std::find(bicc_list[i].begin(), bicc_list[i].end(), target) != bicc_list[i].end()){
        //             target_bct = i;
        //         }
        //     }
        // }

        source_bct = node_to_block[source];
        target_bct = node_to_block[target];
        
        // bct path
        std::vector<int> path_bct = find_path_bct(source_bct, target_bct);

        // skeleton path on original graph, with insertion point
        std::vector<int> skeleton_path, insertion_loc;

        // skeleton weight for skeleton graph
        path_t skeleton_weight = 0;

        // traverse through the bct path we found
        for(int i = 0; i < path_bct.size(); i++){

            // if it is an ap
            if(path_bct[i] >= bicc_list.size()){

                if( !merge || (merge && merge_path.size() == 0)){

                    // add ap point to path
                    skeleton_path.push_back(node_to_ap[path_bct[i]]);

                    // if next val is avail and is ap, add weight to skeleton weight
                    if((i + 1) < path_bct.size()){
                        if(path_bct[i + 1] >= bicc_list.size()){

                            // std::cout << "Add weight to skeleton path weight: " << node_to_ap[path_bct[i]] << " " << node_to_ap[path_bct[i + 1]] << " " << bct_independent_edges[{node_to_ap[path_bct[i]], node_to_ap[path_bct[i + 1]]}] << std::endl;
                            skeleton_weight += bct_independent_edges[{node_to_ap[path_bct[i]], node_to_ap[path_bct[i + 1]]}];
                        }
                    }
                }
                else{

                    path_t temp_weight = 0.0;

                    // if next val is avail and is ap, add weight to merged path
                    if((i + 1) < path_bct.size()){
                        if(path_bct[i + 1] >= bicc_list.size()){
                           temp_weight = bct_independent_edges[{node_to_ap[path_bct[i]], node_to_ap[path_bct[i + 1]]}];
                        }
                    }

                    for(int j = 0; j < merge_path.size(); j++){
                        merge_path[j].p.push_back(node_to_ap[path_bct[i]]);
                        merge_path[j].length += temp_weight;
                    }

                }
            }
            // it is a block
            else{
                // if it's a block pruned with only one element, directly add to path
                if(bicc_list[path_bct[i]].size() == 1){

                    index_t node_id = bicc_list[path_bct[i]][0];

                    if(!merge || (merge && merge_path.size() == 0)){
                        
                        skeleton_path.push_back(node_id);

                        // add weight to skeleton weight
                        if((i - 1) >= 0){
                            skeleton_weight += bct_independent_edges[{node_id , node_to_ap[path_bct[i - 1]]}];
                        }
                        else if((i + 1) < path_bct.size()){
                            skeleton_weight += bct_independent_edges[{node_id, node_to_ap[path_bct[i + 1]]}];
                        }
                    }
                    else{

                        path_t temp_weight = 0.0;

                        if((i - 1) >= 0){
                            temp_weight = bct_independent_edges[{node_id , node_to_ap[path_bct[i - 1]]}];
                        }
                        else if((i + 1) < path_bct.size()){
                            temp_weight = bct_independent_edges[{node_id, node_to_ap[path_bct[i + 1]]}];
                        }

                        for(int j = 0; j < merge_path.size(); j++){
                            merge_path[j].p.push_back(bicc_list[path_bct[i]][0]);
                            merge_path[j].length += temp_weight;
                        }
                    }

                    
                }
                // blocks more than one element
                else{
                    // if it's the first block
                    if(i == 0){
                        // source and target in one block -> sub problem
                        if(path_bct.size() == 1){
                            
                            // insert insertion location
                            insertion_loc.push_back(skeleton_path.size());
                            // std::cout << "\nCalculate subpath from " << source << " to " << target << std::endl;

                            std::optional<std::shared_ptr<BasicGraph<true, true>>> sub_g;
                            index_t sub_source = source;
                            index_t sub_target = target;

                            // if there're sub_graph results, directly reuse it
                            // since we don't store sub ksp results in sub problem, we don't check avail sub_ksp
                            if(sub_graph_mem.find(path_bct[i]) != sub_graph_mem.end()){
                                
                                // std::cout << "Reuse subgraph" << std::endl;

                                // reuse sub_graph 
                                sub_g = sub_graph_mem[path_bct[i]];

                                bool found_src = false;
                                bool found_tar = false;

                                for(int j = 0; j < bicc_list[path_bct[i]].size(); ++j){
                                    if(bicc_list[path_bct[i]][j] == source){
                                        sub_source = j;
                                        found_src = true;
                                    }
                                    else if(bicc_list[path_bct[i]][j] == target){
                                        sub_target = j;
                                        found_tar = true;
                                    }

                                    if(found_src && found_tar){
                                        break;
                                    }
                                }
                            }
                            else{
                                // create sub graph and compute ksp, add to sub path lists
                                sub_g = gen_subgraph(bicc_list[path_bct[i]], path_bct[i], sub_source, sub_target);

                                // in this kind of sub problem, no need to store sub ksp, only store block's sub graph
                                sub_graph_mem[path_bct[i]] = *sub_g;
                            }

                            // compute ksp
                            std::vector<Path> paths = compute_ksp(sub_source, sub_target, K, path_bct[i], *sub_g.value(), inf_gap, early_stopped);

                            if(merge){
                                merge_path = paths;

                                // no need to update gap
                            }
                            else{
                                path_lists.push_back(paths); 
                            }

                        }
                        // general case: multiple elements in bct path
                        else{

                            // insert insertion location
                            insertion_loc.push_back(skeleton_path.size());
                            // std::cout << "\nCalculate subpath from " << source << " to " << node_to_ap[path_bct[i+1]] << std::endl;

                            std::optional<std::shared_ptr<BasicGraph<true, true>>> sub_g;
                            index_t sub_source = source;
                            index_t sub_target = node_to_ap[path_bct[i+1]];

                            // if there're sub_graph results, directly reuse it
                            // since we don't store sub ksp results in sub problem, we don't check avail sub_ksp
                            if(sub_graph_mem.find(path_bct[i]) != sub_graph_mem.end()){
                                
                                // std::cout << "Reuse subgraph" << std::endl;

                                // reuse sub_graph 
                                sub_g = sub_graph_mem[path_bct[i]];

                                bool found_src = false;
                                bool found_tar = false;

                                //std::cout << "Size: " << bicc_list[path_bct[i]].size() << std::endl;

                                for(int j = 0; j < bicc_list[path_bct[i]].size(); j++){
                                    if(bicc_list[path_bct[i]][j] == source){
                                        sub_source = j;
                                        found_src = true;
                                    }
                                    else if(bicc_list[path_bct[i]][j] == node_to_ap[path_bct[i+1]]){
                                        sub_target = j;
                                        found_tar = true;
                                    }

                                    if(found_src && found_tar){
                                        break;
                                    }
                                }
                            }
                            else{
                                // create sub graph and compute ksp, add to sub path lists
                                sub_g = gen_subgraph(bicc_list[path_bct[i]], path_bct[i], sub_source, sub_target);

                                // std::cout << "test1" << std::endl;
                                // std::cout << "1st element of csr: " << sub_g.csr[0] << std::endl;

                                // in this kind of sub problem, no need to store sub ksp, only store block's sub graph
                                sub_graph_mem[path_bct[i]] = *sub_g;
                            }

                            // std::cout << "test" << std::endl;
                            // std::cout << "1st element of csr: " << sub_g.csr[0] << std::endl;

                            // compute ksp
                            std::vector<Path> paths = compute_ksp(sub_source, sub_target, K, path_bct[i], *sub_g.value(), inf_gap, early_stopped);

                            // merge mode
                            if(merge){
                                merge_path = paths;
                                gap = merge_path[merge_path.size() - 1].length - merge_path[0].length;
                                for(Path& m_p : merge_path){
                                    m_p.p.pop_back();
                                }
                            }
                            // base & prune mode
                            else{
                                path_lists.push_back(paths);
                            }
                        }  
                    }
                    // last block
                    else if(i == path_bct.size() - 1){
                
                        // insert insertion location
                        insertion_loc.push_back(skeleton_path.size());
                        // std::cout << "\nCalculate subpath from " << node_to_ap[path_bct[i-1]] << " to " << target << std::endl;

                        std::optional<std::shared_ptr<BasicGraph<true, true>>> sub_g;
                        index_t sub_source = node_to_ap[path_bct[i-1]];
                        index_t sub_target = target;

                        // if there're sub_graph results, directly reuse it
                        // since we don't store sub ksp results in sub problem, we don't check avail sub_ksp
                        if(sub_graph_mem.find(path_bct[i]) != sub_graph_mem.end()){
                            
                            // std::cout << "Reuse subgraph" << std::endl;

                            // reuse sub_graph 
                            sub_g = sub_graph_mem[path_bct[i]];

                            bool found_src = false;
                            bool found_tar = false;

                            for(int j = 0; j < bicc_list[path_bct[i]].size(); ++j){
                                if(bicc_list[path_bct[i]][j] == node_to_ap[path_bct[i-1]]){
                                    sub_source = j;
                                    found_src = true;
                                }
                                else if(bicc_list[path_bct[i]][j] == target){
                                    sub_target = j;
                                    found_tar = true;
                                }

                                if(found_src && found_tar){
                                    break;
                                }
                            }
                        }
                        else{
                            // create sub graph and compute ksp, add to sub path lists
                            sub_g = gen_subgraph(bicc_list[path_bct[i]], path_bct[i], sub_source, sub_target);

                            // in this kind of sub problem, no need to store sub ksp, only store block's sub graph
                            sub_graph_mem[path_bct[i]] = *sub_g;
                        }

                        // compute ksp
                        std::vector<Path> paths = compute_ksp(sub_source, sub_target, K, path_bct[i], *sub_g.value(), gap, early_stopped);

                        // merge mode
                        if(merge){
                            // if not merged yet -> add skeleton path and skeleton weight 
                            if(merge_path.size() == 0){
                                add_skeleton_merge(skeleton_path, skeleton_weight, paths, gap, true);
                                merge_path = paths;
                            }
                            else{
                                combine_ksp_merge(merge_path, paths, true, gap);
                            }                          
                        }
                        // base & prune mode
                        else{
                            path_lists.push_back(paths);
                        }

                    }
                    // general case
                    else{
                        
                        // insert insertion location
                        insertion_loc.push_back(skeleton_path.size());
                        // std::cout << "\nCalculate subpath from " << node_to_ap[path_bct[i-1]] << " to " << node_to_ap[path_bct[i+1]] << std::endl;

                        std::optional<std::shared_ptr<BasicGraph<true, true>>> sub_g;
                        index_t sub_source = node_to_ap[path_bct[i-1]];
                        index_t sub_target = node_to_ap[path_bct[i+1]];

                        std::pair<index_t, index_t> key = {node_to_ap[path_bct[i-1]], node_to_ap[path_bct[i+1]]};

                        bool found_ksp = (sub_ksp_mem.find(key) != sub_ksp_mem.end());

                        bool ksp_valid_merge = false;

                        if(merge && found_ksp){

                            // if ksp already complete, it's valid
                            if(sub_ksp_complete_mem[key]){
                                ksp_valid_merge = true;
                            }
                            else{
                                int temp_size = sub_ksp_mem[key].size();
                                path_t temp_gap = sub_ksp_mem[key][temp_size - 1].length - sub_ksp_mem[key][0].length;

                                // key point: if current temp_gap is greater than merge gap, it's valid
                                // also current merge size needs to equal to K
                                if(temp_gap >= gap && merge_path.size() == K){
                                    ksp_valid_merge = true;
                                }
                            }
                            
                        }

                        // if there is sub_ksp results, reuse it
                        if( found_ksp && (!merge || (merge && ksp_valid_merge)) ){
                            
                            // std::cout << "Reuse subKSP" << std::endl;

                            if(sub_ksp_mem[key][0].p[0] == node_to_ap[path_bct[i-1]]){

                                std::vector<Path> paths = sub_ksp_mem[key];

                                if(merge){
                                    // if not merged yet -> add skeleton path and skeleton weight 
                                    if(merge_path.size() == 0){                      
                                        add_skeleton_merge(skeleton_path, skeleton_weight, paths, gap, false);
                                        merge_path = paths;                                    
                                    }
                                    else{
                                        combine_ksp_merge(merge_path, paths, false, gap);
                                    }  
                                }
                                else{
                                    path_lists.push_back(paths);
                                }
                                
                            }
                            // if reverse, have to reverse each path
                            else{
                                std::vector<Path> sub_paths = sub_ksp_mem[key];

                                for(auto& path : sub_paths){
                                    std::reverse(path.p.begin(), path.p.end());
                                }

                                if(merge){
                                    // if not merged yet -> add skeleton path and skeleton weight 
                                    if(merge_path.size() == 0){                      
                                        add_skeleton_merge(skeleton_path, skeleton_weight, sub_paths, gap, false);
                                        merge_path = sub_paths;                                    
                                    }
                                    else{
                                        combine_ksp_merge(merge_path, sub_paths, false, gap);
                                    }  
                                }
                                else{
                                    path_lists.push_back(sub_paths);
                                }

                            }
                        }
                        else{

                            if(sub_graph_mem.find(path_bct[i]) != sub_graph_mem.end()){
                            
                                // std::cout << "Reuse subgraph" << std::endl;

                                // reuse sub_graph 
                                sub_g = sub_graph_mem[path_bct[i]];

                                bool found_src = false;
                                bool found_tar = false;

                                for(int j = 0; j < bicc_list[path_bct[i]].size(); ++j){
                                    if(bicc_list[path_bct[i]][j] == node_to_ap[path_bct[i-1]]){
                                        sub_source = j;
                                        found_src = true;
                                    }
                                    else if(bicc_list[path_bct[i]][j] == node_to_ap[path_bct[i+1]]){
                                        sub_target = j;
                                        found_tar = true;
                                    }

                                    if(found_src && found_tar){
                                        break;
                                    }
                                }
                            }
                            else{
                                // create sub graph and compute ksp, add to sub path lists
                                sub_g = gen_subgraph(bicc_list[path_bct[i]], path_bct[i], sub_source, sub_target);

                                // in this kind of sub problem, no need to store sub ksp, only store block's sub graph
                                sub_graph_mem[path_bct[i]] = *sub_g;
                            }

                            // compute ksp

                            if(!merge || (merge && merge_path.size() != K)){

                                std::vector<Path> paths = compute_ksp(sub_source, sub_target, K, path_bct[i], *sub_g.value(), inf_gap, early_stopped);

                                // store ksp
                                sub_ksp_mem[key] = paths;

                                if(merge){

                                    early_stopped = false;
                                    sub_ksp_complete_mem[key] = true;

                                    if(merge_path.size() == 0){                      
                                        add_skeleton_merge(skeleton_path, skeleton_weight, paths, gap, false);
                                        merge_path = paths;                                    
                                    }
                                    else{
                                        combine_ksp_merge(merge_path, paths, false, gap);
                                    }

                                }
                                else{
                                    path_lists.push_back(paths);
                                }                            
                            }
                            else{
                                std::vector<Path> paths = compute_ksp(sub_source, sub_target, K, path_bct[i], *sub_g.value(), gap, early_stopped);

                                // store ksp
                                sub_ksp_mem[key] = paths;

                                if(merge){

                                    if(!early_stopped){
                                        sub_ksp_complete_mem[key] = true;
                                        early_stopped = false;
                                    }

                                    if(merge_path.size() == 0){                      
                                        add_skeleton_merge(skeleton_path, skeleton_weight, paths, gap, false);
                                        merge_path = paths;                                    
                                    }
                                    else{
                                        combine_ksp_merge(merge_path, paths, false, gap);
                                    }
                                }
                                else{
                                     path_lists.push_back(paths);
                                }
                            }                 
                        }
                    }
                }
            }

            // if(merge){
                
            //     std::cout << "\nBCT merged results Merge: " << std::endl;
            //     for(auto i : merge_path){
            //         i.print_path();
            //     }

            //     std::cout << "Current GAP: " << gap << std::endl;
            // }          
        }

        if(merge){
            
            // if there's no element in merge_path, skeleton path is the only path
            if(merge_path.size() == 0){
                Path temp_path;
                temp_path.p = skeleton_path;
                temp_path.length = skeleton_weight;
                merge_path.push_back(temp_path);
            }

            // std::cout << "\nBCT FINAL results Merge: " << std::endl;
            // for(auto i : merge_path){
            //     i.print_path();
            // }
        }
        // base + prune
        else{

             // path lists
            // std::cout << "\n\nPath lists before PRUNE:\n";
            // for (size_t i = 0; i < path_lists.size(); ++i) {
            //     std::cout << "== sub ksps in block " << i << " ==\n";
            //     for (size_t j = 0; j < path_lists[i].size(); ++j) {
            //         std::cout << "  -- Subpath from block " << j << ":\n";
            //         path_lists[i][j].print_path();
            //     }
            //     std::cout << std::endl;
            // }

            if(prune){
                // prune sub ksps first
                joint_ksp_minheap_prune(path_lists);
            }
            // path lists
            // std::cout << "\n\nPath lists AFTER PRUNE:\n";
            // for (index_t i = 0; i < path_lists.size(); ++i) {
            //     std::cout << "== sub ksps in block " << i << " ==\n";
            //     for (index_t j = 0; j < path_lists[i].size(); ++j) {
            //         std::cout << "  -- Subpath from block " << j << ":\n";
            //         path_lists[i][j].print_path();
            //     }
            //     std::cout << std::endl;
            // }

            // joint sub ksps
            std::vector<std::vector<Path>> min_paths = joint_ksp_minheap(path_lists);

            // combine results
            std::vector<Path> final_ksp = combine_ksp_result(min_paths, skeleton_path, insertion_loc, skeleton_weight);

            //print results
            // std::cout << "\nSkeleton Path: ";
            // for(int i = 0; i < skeleton_path.size(); i++){
            //     std::cout << skeleton_path[i] << " ";
            // }
            // std::cout << "\n";

            // std::cout << "Skeleton Path Weight: " << skeleton_weight << std::endl;

            // std::cout << "Insertion positions: ";
            // for(int i = 0; i < insertion_loc.size(); i++){
            //     std::cout << insertion_loc[i] << " ";
            // }
            // std::cout << "\n\n";

            // joint
            // std::cout << "Top " << min_paths.size() << " combined paths:\n";
            // for (size_t i = 0; i < min_paths.size(); ++i) {
            //     std::cout << "== Combined Path " << i << " ==\n";
            //     for (size_t j = 0; j < min_paths[i].size(); ++j) {
            //         std::cout << "  -- Subpath from block " << j << ":\n";
            //         min_paths[i][j].print_path();
            //     }
            //     std::cout << std::endl;
            // }

            // print final paths found
            // std::cout << "\nBCT final results: " << std::endl;
            // for(auto i : final_ksp){
            //     i.print_path();
            // }
        }       
    }






};

#endif  // BLOCK_CUT_TREE_H