#include <iostream>
#include <vector>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>
#include <functional>
#include <fstream>

// Assuming these header files exist based on the user's code

#include "scc_common.h"
#include "graph.h"
#include "graph_undirected.h"
#include "bicc_bfs.h"

using namespace std;


// --------------------- WCC ---------------------
vector<vector<int>> computeWCC(const unordered_map<int, vector<int>>& graph) {
    unordered_map<int, vector<int>> undirected;
    for (auto& [u, neighbors] : graph) {
        for (int v : neighbors) {
            undirected[u].push_back(v);
            undirected[v].push_back(u);
        }
    }
    for (auto& [u, _] : graph) {
        if (undirected.find(u) == undirected.end()) {
            undirected[u] = {};
        }
    }

    unordered_set<int> visited;
    vector<vector<int>> wccs;

    for (auto& [node, _] : undirected) {
        if (visited.count(node)) continue;
        vector<int> comp;
        stack<int> stk;
        stk.push(node);
        visited.insert(node);
        while (!stk.empty()) {
            int curr = stk.top(); stk.pop();
            comp.push_back(curr);
            for (int nei : undirected[curr]) {
                if (!visited.count(nei)) {
                    visited.insert(nei);
                    stk.push(nei);
                }
            }
        }
        wccs.push_back(comp);
    }
    return wccs;
}

// --------------------- BCC ---------------------
class BCC {
public:
    unordered_map<int, vector<int>> g;
    unordered_map<int, int> disc, low;
    vector<vector<int>> bccs;
    unordered_map<int, int> artPoints;
    int time = 0;
    stack<pair<int, int>> edgeStack;

    BCC(const unordered_map<int, vector<int>>& subg) : g(subg) {}

    void dfs(int u, int parent) {
        disc[u] = low[u] = ++time;
        int children = 0;

        if (g.count(u)) {
            for (int v : g[u]) {
                if (v == parent) continue;

                if (!disc.count(v)) {
                    children++;
                    edgeStack.push({u, v});
                    dfs(v, u);
                    low[u] = min(low[u], low[v]);

                    if ((parent == -1 && children > 1) || (parent != -1 && low[v] >= disc[u])) {
                        artPoints[u]++;
                        vector<int> bcc;
                        pair<int, int> current_edge;
                        do {
                            current_edge = edgeStack.top(); edgeStack.pop();
                            bcc.push_back(current_edge.first);
                            bcc.push_back(current_edge.second);
                        } while (current_edge.first != u || current_edge.second != v);

                        bcc.push_back(u);
                        
                        sort(bcc.begin(), bcc.end());
                        bcc.erase(unique(bcc.begin(), bcc.end()), bcc.end());
                        bccs.push_back(bcc);
                    }
                } else if (disc[v] < disc[u]) {
                    low[u] = min(low[u], disc[v]);
                    edgeStack.push({u, v});
                }
            }
        }
    }

    void run() {
        for (auto& [u, _] : g) {
            if (!disc.count(u)) {
                time = 0;
                while(!edgeStack.empty()) edgeStack.pop(); 
                dfs(u, -1);

                if (!edgeStack.empty()) {
                    vector<int> bcc;
                    while (!edgeStack.empty()) {
                        bcc.push_back(edgeStack.top().first);
                        bcc.push_back(edgeStack.top().second);
                        edgeStack.pop();
                    }
                    sort(bcc.begin(), bcc.end());
                    bcc.erase(unique(bcc.begin(), bcc.end()), bcc.end());
                    bccs.push_back(bcc);
                }
            }
        }
        unordered_map<int, int> distinctArtPoints;
        for (auto& [ap, _] : artPoints) {
            distinctArtPoints[ap] = 1; 
        }
        artPoints = distinctArtPoints;
    }

    const vector<vector<int>>& getBCCs() const { return bccs; }
    const unordered_map<int, int>& getArticulationPoints() const { return artPoints; }
};

unordered_map<int, vector<int>> convert_to_adj_list(const graph* g_csr) {
    unordered_map<int, vector<int>> adj_list;
    if (!g_csr) return adj_list;

    for (index_t u = 0; u < g_csr->vert_count; ++u) {
        for (index_t i = g_csr->fw_beg_pos[u]; i < g_csr->fw_beg_pos[u+1]; ++i) {
            vertex_t v = g_csr->fw_csr[i];
            adj_list[u].push_back(v);
        }
        if (g_csr->fw_beg_pos[u] == g_csr->fw_beg_pos[u+1] && adj_list.find(u) == adj_list.end()) {
             adj_list[u] = {};
        }
    }
    return adj_list;
}

// --------------------- Block-Cut Tree (bc_tree.h) ---------------------
class BlockCutTree {
public:
    enum NodeType { B_NODE, C_NODE };

    struct Node {
        NodeType type;
        int id; // For B-nodes, this is the index of the BCC. For C-nodes, this is the original vertex id.
        vector<int> neighbors;
        Node(NodeType t, int i) : type(t), id(i) {}
    };

    vector<Node> nodes;
    unordered_map<int, int> apMap; // Maps original vertex id to C-node id
    unordered_map<int, int> bccMap; // Maps BCC index to B-node id

    void build(const vector<vector<int>>& bicc_components, const bool* ap, const graph_undirected* ug) {
        // Step 1: Create C-nodes for each articulation point
        for (vertex_t i = 0; i < ug->vert_count; ++i) {
            if (ap[i]) {
                int c_node_id = nodes.size();
                nodes.emplace_back(C_NODE, i);
                apMap[i] = c_node_id;
            }
        }

        // Step 2: Create B-nodes for each BCC and link them to C-nodes
        for (size_t i = 0; i < bicc_components.size(); ++i) {
            const auto& bcc = bicc_components[i];
            int b_node_id = nodes.size();
            nodes.emplace_back(B_NODE, i);
            bccMap[i] = b_node_id;

            for (int node_in_bcc : bcc) {
                if (ap[node_in_bcc]) {
                    // This node is an articulation point, link B-node and C-node
                    int c_node_id = apMap[node_in_bcc];
                    nodes[b_node_id].neighbors.push_back(c_node_id);
                    nodes[c_node_id].neighbors.push_back(b_node_id);
                }
            }
        }
    }

    void printTree() {
        cout << "--- Block-Cut Tree Structure ---" << endl;
        for (size_t i = 0; i < nodes.size(); ++i) {
            const auto& node = nodes[i];
            if (node.type == B_NODE) {
                cout << "B-Node (BCC " << node.id << ") -> Neighbors: ";
            } else {
                cout << "C-Node (AP " << node.id << ") -> Neighbors: ";
            }
            for (int neighbor_id : node.neighbors) {
                const auto& neighbor = nodes[neighbor_id];
                if (neighbor.type == B_NODE) {
                    cout << "B-Node (BCC " << neighbor.id << ") ";
                } else {
                    cout << "C-Node (AP " << neighbor.id << ") ";
                }
            }
            cout << endl;
        }
        cout << "---------------------------------" << endl;
    }
};



// --------------------- DAG ---------------------
void generate_and_save_dag(
    const std::map<int, int>& nodeToScc,
    const std::unordered_map<int, std::vector<int>>& adj_list,
    const std::string& output_filename)
{
    // Map original SCC IDs to a contiguous range
    std::map<int, int> scc_id_map;
    int current_scc_id = 0;
    for (const auto& pair : nodeToScc) {
        if (scc_id_map.find(pair.second) == scc_id_map.end()) {
            scc_id_map[pair.second] = current_scc_id++;
        }
    }
    int scc_count = current_scc_id;

    if (scc_count == 0) {
        std::cerr << "No SCCs found to build DAG." << std::endl;
        return;
    }

    std::vector<std::set<int>> dag_adj_set(scc_count);

    // Build the DAG by iterating through original edges
    for (const auto& pair : adj_list) {
        int u = pair.first;
        if (nodeToScc.find(u) == nodeToScc.end()) continue;

        int u_scc_id = nodeToScc.at(u);
        int u_supernode_id = scc_id_map.at(u_scc_id);

        for (int v : pair.second) {
            if (nodeToScc.find(v) == nodeToScc.end()) continue;
            int v_scc_id = nodeToScc.at(v);

            // Add an edge if the source and destination are in different SCCs
            if (u_scc_id != v_scc_id) {
                int v_supernode_id = scc_id_map.at(v_scc_id);
                dag_adj_set[u_supernode_id].insert(v_supernode_id);
            }
        }
    }

    // Convert sets to vectors for easier use
    std::vector<std::vector<int>> dag_adj_list(scc_count);
    for (int i = 0; i < scc_count; ++i) {
        for (int neighbor : dag_adj_set[i]) {
            dag_adj_list[i].push_back(neighbor);
        }
    }

    // --- LOG: print the DAG ---
    std::cout << "--- DAG (Supernodes from SCCs) ---" << std::endl;
    for (int i = 0; i < scc_count; ++i) {
        std::cout << "Supernode " << i << " -> ";
        if (dag_adj_list[i].empty()) {
            std::cout << "(no outgoing edges)";
        } else {
            for (int v : dag_adj_list[i]) {
                std::cout << v << " ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << "---------------------------------" << std::endl;

    // Save the DAG to a binary file
    std::ofstream outfile(output_filename, std::ios::out | std::ios::binary);
    if (!outfile) {
        std::cerr << "Error: Could not open file " << output_filename << " for writing." << std::endl;
        return;
    }

    // Write the number of supernodes
    outfile.write(reinterpret_cast<const char*>(&scc_count), sizeof(int));

    // Write the DAG's adjacency list
    for (int i = 0; i < scc_count; ++i) {
        int neighbor_count = dag_adj_list[i].size();
        outfile.write(reinterpret_cast<const char*>(&neighbor_count), sizeof(int));
        if (neighbor_count > 0) {
            outfile.write(reinterpret_cast<const char*>(dag_adj_list[i].data()), neighbor_count * sizeof(int));
        }
    }

    outfile.close();
    std::cout << "DAG saved to " << output_filename << std::endl;
}

// --------------------- main ---------------------
int main() {
    // The main function is here to be a complete runnable example
    // Note: The file paths are specific to your environment and will need to be adjusted.
    const string string_path = "/home/cxl6029/projects/GIR/gir_constructor/";
    const string undirected_bg_file = string_path + "undirected_begin.bin";
    const string undirected_aj_file = string_path + "undirected_adjacent.bin";

    printf("start to calculate articulation points \n");
    graph_undirected *ug = new graph_undirected(undirected_bg_file.c_str(), undirected_aj_file.c_str());
    vertex_t V = ug->vert_count;
    bool *ap = new bool[V];
    std::vector<std::vector<int>> bicc_components;
    double time_begin = wtime();
    find_bicc_bfs(ug, ap, bicc_components);
    double time_end = wtime();
    printf("Time(s) to find AP in whole undirected graph, %.3lf\n", (time_end - time_begin));

    cout << "--- Building the Block-Cut Tree starting---" << endl;
    time_begin = wtime();
    BlockCutTree* bc_tree = new BlockCutTree();
    bc_tree->build(bicc_components, ap, ug);
    cout << "--- Block-Cut Tree built successfully ---." << endl;
    time_end = wtime();
    bc_tree->printTree(); // Add a print function to visualize the result
    printf("Time(s) to build bc tree, %.3lf\n", (time_end - time_begin));

    // --- New DAG generation part ---
    // You need a directed graph (graph type) to find SCCs.
    // Based on the graph constructor, you need four file paths.
    const string string_path2 = "/home/cxl6029/projects/GIR/data/example_graph/";
    const string directed_fw_bg_file = string_path2 + "fw_begin.bin"; // Assuming forward begin file
    const string directed_fw_aj_file = string_path2 + "fw_adjacent.bin"; // Assuming forward adjacent file
    const string directed_bw_bg_file = string_path2 + "bw_begin.bin"; // Assuming backward begin file
    const string directed_bw_aj_file = string_path2 + "bw_adjacent.bin"; // Assuming backward adjacent file

    // Instantiate the graph with four arguments
    graph *g = new graph(directed_fw_bg_file.c_str(), directed_fw_aj_file.c_str(), directed_bw_bg_file.c_str(), directed_bw_aj_file.c_str());

    std::map<int, int> nodeToScc;
    std::vector<std::vector<int>> sccs;
    double avg_time[15] = {0};

    // Call the scc_detection function
    int alpha = 1;
    int beta = 1;
    int gamma = 1;
    double theta = 0.5;
    index_t thread_count = 1;

    scc_detection(g, alpha, beta, gamma, theta, thread_count, avg_time, sccs, nodeToScc);

    // After scc_detection, nodeToScc is populated.
    // Now you can call the DAG generation function with the correct arguments.
    // First, get the adjacency list from the directed graph
    std::unordered_map<int, vector<int>> adj_list = convert_to_adj_list(g);
    
    // Then call the DAG generation function
    generate_and_save_dag(nodeToScc, adj_list, "dag.bin");

    // Clean up
    delete g;
    delete[] ap;
    delete ug;
    delete bc_tree;

    return 0;
}