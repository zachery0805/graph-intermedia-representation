#ifndef BLOCK_CUT_TREE_H
#define BLOCK_CUT_TREE_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include "graph_undirected.h" // 你的无向图类
using namespace std;

// --------------------- Block-Cut Tree ---------------------
class BlockCutTree {
public:
    enum NodeType { B_NODE, C_NODE };

    struct Node {
        NodeType type;
        int id; // B-node: BCC索引, C-node: 原始顶点ID
        vector<int> neighbors;
        Node(NodeType t, int i) : type(t), id(i) {}
    };

    vector<Node> nodes;
    unordered_map<int, int> apMap;   // 原顶点ID -> C节点ID
    unordered_map<int, int> bccMap;  // BCC索引 -> B节点ID

    // 构建函数
    void build(const vector<vector<int>>& bicc_components, const bool* ap, const graph_undirected* ug) {
        // Step 1: 创建 C-node（关节点）
        for (vertex_t i = 0; i < ug->vert_count; ++i) {
            if (ap[i]) {
                int c_node_id = nodes.size();
                nodes.emplace_back(C_NODE, i);
                apMap[i] = c_node_id;
            }
        }

        // Step 2: 创建 B-node（每个 BCC）并与 C-node 相连
        for (size_t i = 0; i < bicc_components.size(); ++i) {
            const auto& bcc = bicc_components[i];
            int b_node_id = nodes.size();
            nodes.emplace_back(B_NODE, i);
            bccMap[i] = b_node_id;

            for (int node_in_bcc : bcc) {
                if (ap[node_in_bcc]) {
                    int c_node_id = apMap[node_in_bcc];
                    nodes[b_node_id].neighbors.push_back(c_node_id);
                    nodes[c_node_id].neighbors.push_back(b_node_id);
                }
            }
        }
    }

    // 打印函数
    void printTree() const {
        cout << "--- Block-Cut Tree ---" << endl;
        for (size_t i = 0; i < nodes.size(); ++i) {
            const auto& node = nodes[i];
            if (node.type == B_NODE)
                cout << "B-Node (BCC " << node.id << ") -> Neighbors: ";
            else
                cout << "C-Node (AP " << node.id << ") -> Neighbors: ";

            for (int neighbor_id : node.neighbors) {
                const auto& neighbor = nodes[neighbor_id];
                if (neighbor.type == B_NODE)
                    cout << "B-Node(" << neighbor.id << ") ";
                else
                    cout << "C-Node(" << neighbor.id << ") ";
            }
            cout << endl;
        }
        cout << "----------------------" << endl;
    }
};

#endif // BLOCK_CUT_TREE_H
