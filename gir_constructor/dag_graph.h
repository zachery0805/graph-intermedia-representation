#ifndef DAG_GRAPH_H
#define DAG_GRAPH_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <stdexcept>

// Using a type alias for vertex IDs for clarity and easy modification.
// Your other files seem to use a similar convention.
using vertex_t = int;

/**
 * @struct SuperNode
 * @brief Represents a "super node" in the SCC-DAG, which corresponds to a
 *        Strongly Connected Component (SCC) from the original graph.
 */

struct SuperNode {
    int scc_id;
    std::vector<vertex_t> original_nodes;

    // Information about Bi-Connected Components (BCCs) within this SCC.
    std::vector<std::vector<vertex_t>> bccs;
    
    // Articulation points within this SCC.
    std::vector<vertex_t> articulation_points;
};

/**
 * @class dag_graph
 * @brief Represents a graph decomposed into a Directed Acyclic Graph (DAG) of
 * its Strongly Connected Components (SCCs).
 *
 * This class takes an original directed graph and computes its SCCs,
 * Weakly Connected Components (WCCs), and for each SCC, its
 * Bi-Connected Components (BCCs) and articulation points.
 *
 * The SCCs are treated as "super nodes" in a new DAG.
 */
class dag_graph {
public:
    /**
     * @brief Constructs the DAG graph representation by analyzing the original graph.
     * @param original_graph The original graph as an adjacency list.
     */
    explicit dag_graph(const std::unordered_map<vertex_t, std::vector<vertex_t>>& original_graph);

    // --- Accessors for computed graph properties ---

    /**
     * @brief Gets all super nodes (SCCs).
     * @return A constant reference to the vector of SuperNodes.
     */
    const std::vector<SuperNode>& get_super_nodes() const { return super_nodes_; }

    /**
     * @brief Gets the SCC-DAG.
     * @return A constant reference to the SCC-DAG adjacency list.
     *         Keys are SCC IDs, values are sets of successor SCC IDs.
     */
    const std::unordered_map<int, std::set<int>>& get_scc_dag() const { return scc_dag_; }

    /**
     * @brief Gets the SCC ID for a given original node.
     * @param node_id The ID of the node in the original graph.
     * @return The ID of the SCC containing the node.
     * @throws std::out_of_range if the node_id is not found.
     */
    int get_scc_id_for_node(vertex_t node_id) const { return node_to_scc_id_.at(node_id); }

    /**
     * @brief Gets the WCC ID for a given original node.
     * @param node_id The ID of the node in the original graph.
     * @return The ID of the WCC containing the node.
     * @throws std::out_of_range if the node_id is not found.
     */
    int get_wcc_id_for_node(vertex_t node_id) const { return node_to_wcc_id_.at(node_id); }

    /**
     * @brief Gets a specific SuperNode by its SCC ID.
     * @param scc_id The ID of the SCC.
     * @return A constant reference to the SuperNode.
     * @throws std::out_of_range if the scc_id is invalid.
     */
    const SuperNode& get_super_node(int scc_id) const { return super_nodes_.at(scc_id); }

    /**
     * @brief Gets all Weakly Connected Components.
     * @return A constant reference to the vector of WCCs.
     */
    const std::vector<std::vector<vertex_t>>& get_wccs() const { return wccs_; }
};