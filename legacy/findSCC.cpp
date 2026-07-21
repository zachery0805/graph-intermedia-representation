#include <iostream>
#include <vector>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>
using namespace std;

// --------------------- Tarjan Strongly Connected Components ---------------------
class TarjanSCC {
    int index = 0;
    unordered_map<int, int> indices, lowlink;
    stack<int> stk;
    unordered_set<int> onStack;
    vector<vector<int>> sccs;
    unordered_map<int, int> nodeToScc;
    unordered_map<int, vector<int>> graph;

public:
    TarjanSCC(const unordered_map<int, vector<int>>& g) : graph(g) {}

    void strongconnect(int v) {
        indices[v] = lowlink[v] = index++;
        stk.push(v);
        onStack.insert(v);

        for (int w : graph[v]) {
            if (indices.find(w) == indices.end()) {
                strongconnect(w);
                lowlink[v] = min(lowlink[v], lowlink[w]);
            } else if (onStack.count(w)) {
                lowlink[v] = min(lowlink[v], indices[w]);
            }
        }

        if (lowlink[v] == indices[v]) {
            vector<int> scc;
            int w;
            do {
                w = stk.top(); stk.pop();
                onStack.erase(w);
                scc.push_back(w);
                nodeToScc[w] = sccs.size();
            } while (w != v);
            sccs.push_back(scc);
        }
    }

    void run() {
        for (auto& p : graph) {
            int v = p.first;
            if (indices.find(v) == indices.end()) {
                strongconnect(v);
            }
        }
    }

    vector<vector<int>> getSCCs() const { return sccs; }
    unordered_map<int, int> getNodeToSccMap() const { return nodeToScc; }
};

// --------------------- AVL Tree for Forest Construction ---------------------
struct TreeNode {
    int key;
    vector<int> val; // original nodes in this SCC
    TreeNode* left = nullptr;
    TreeNode* right = nullptr;
    int height = 1;

    TreeNode(int k, vector<int> v) : key(k), val(v) {}
};

int height(TreeNode* n) {
    return n ? n->height : 0;
}

int getBalance(TreeNode* n) {
    return n ? height(n->left) - height(n->right) : 0;
}

TreeNode* rightRotate(TreeNode* y) {
    TreeNode* x = y->left;
    TreeNode* T2 = x->right;
    x->right = y;
    y->left = T2;
    y->height = max(height(y->left), height(y->right)) + 1;
    x->height = max(height(x->left), height(x->right)) + 1;
    return x;
}

TreeNode* leftRotate(TreeNode* x) {
    TreeNode* y = x->right;
    TreeNode* T2 = y->left;
    y->left = x;
    x->right = T2;
    x->height = max(height(x->left), height(x->right)) + 1;
    y->height = max(height(y->left), height(y->right)) + 1;
    return y;
}

TreeNode* insert(TreeNode* node, int key, vector<int> val) {
    if (!node) return new TreeNode(key, val);
    if (key < node->key)
        node->left = insert(node->left, key, val);
    else
        node->right = insert(node->right, key, val);

    node->height = 1 + max(height(node->left), height(node->right));
    int balance = getBalance(node);

    // Rebalance
    if (balance > 1 && key < node->left->key) return rightRotate(node);
    if (balance < -1 && key > node->right->key) return leftRotate(node);
    if (balance > 1 && key > node->left->key) {
        node->left = leftRotate(node->left);
        return rightRotate(node);
    }
    if (balance < -1 && key < node->right->key) {
        node->right = rightRotate(node->right);
        return leftRotate(node);
    }

    return node;
}

void printTree(TreeNode* node, int depth = 0) {
    if (!node) return;
    printTree(node->right, depth + 1);
    for (int i = 0; i < depth; ++i) cout << "    ";
    cout << "Node " << node->key << ": ";
    for (int n : node->val) cout << n << " ";
    cout << endl;
    printTree(node->left, depth + 1);
}

// --------------------- Main ---------------------
int main() {
    unordered_map<int, vector<int>> graph = {
        {1, {2, 4}},
        {2, {3, 4}},
        {3, {1, 4, 10}},
        {4, {5, 6, 8}},
        {5, {4}},
        {6, {7}},
        {7, {8}},
        {8, {9}},
        {9, {7}},
        {10, {3, 6, 11, 12}},
        {11, {12}},
        {12, {10}},
        {13, {14}},
        {14, {15, 17}},
        {15, {16}},
        {16, {13}},
        {17, {18}}
    };

    TarjanSCC tarjan(graph);
    tarjan.run();

    auto sccs = tarjan.getSCCs();
    auto nodeToScc = tarjan.getNodeToSccMap();

    cout << "SCCs:\n";
    for (size_t i = 0; i < sccs.size(); ++i) {
        cout << "SCC " << i << ": ";
        for (int node : sccs[i]) cout << node << " ";
        cout << "\n";
    }

    // 构建 SCC DAG
    map<int, set<int>> sccGraph;
    for (auto& [u, neighbors] : graph) {
        for (int v : neighbors) {
            int su = nodeToScc[u];
            int sv = nodeToScc[v];
            if (su != sv) sccGraph[su].insert(sv);
        }
    }

    cout << "\nSimplified DAG:\n";
    for (auto& [u, neighbors] : sccGraph) {
        cout << "SuperNode " << u << " -> ";
        for (int v : neighbors) cout << v << " ";
        cout << "\n";
    }

    // 构建森林
    unordered_set<int> visited;
    vector<TreeNode*> forest;

    function<TreeNode*(int, TreeNode*)> dfsBuildTree = [&](int node, TreeNode* root) -> TreeNode* {
        if (visited.count(node)) return root;
        visited.insert(node);
        root = insert(root, node, sccs[node]);
        for (int nei : sccGraph[node]) {
            root = dfsBuildTree(nei, root);
        }
        return root;
    };

    for (size_t node = 0; node < sccs.size(); ++node) {
        if (!visited.count(node)) {
            TreeNode* tree = dfsBuildTree(node, nullptr);
            forest.push_back(tree);
        }
    }

    cout << "\nBalanced Binary Forest:\n";
    for (size_t i = 0; i < forest.size(); ++i) {
        cout << "\nTree " << i << ":\n";
        printTree(forest[i]);
    }

    return 0;
}
