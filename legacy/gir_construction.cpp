#include <iostream>
#include <vector>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>
#include <functional>
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

// --------------------- WCC ---------------------
vector<vector<int>> computeWCC(const unordered_map<int, vector<int>>& graph) {
    unordered_map<int, vector<int>> undirected;
    for (auto& [u, neighbors] : graph) {
        for (int v : neighbors) {
            undirected[u].push_back(v);
            undirected[v].push_back(u);
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
    vector<pair<int, int>> bridges;
    vector<vector<int>> bccs;
    unordered_map<int, int> artPoints;
    int time = 0;
    stack<pair<int, int>> edgeStack;

    BCC(const unordered_map<int, vector<int>>& subg) : g(subg) {}

    void dfs(int u, int parent) {
        disc[u] = low[u] = ++time;
        int children = 0;

        for (int v : g[u]) {
            if (!disc.count(v)) {
                children++;
                edgeStack.push({u, v});
                dfs(v, u);
                low[u] = min(low[u], low[v]);

                if ((parent == -1 && children > 1) || (parent != -1 && low[v] >= disc[u])) {
                    artPoints[u]++;
                    vector<int> bcc;
                    while (!edgeStack.empty() && edgeStack.top() != make_pair(u, v)) {
                        bcc.push_back(edgeStack.top().first);
                        bcc.push_back(edgeStack.top().second);
                        edgeStack.pop();
                    }
                    bcc.push_back(u);
                    bcc.push_back(v);
                    edgeStack.pop();
                    sort(bcc.begin(), bcc.end());
                    bcc.erase(unique(bcc.begin(), bcc.end()), bcc.end());
                    bccs.push_back(bcc);
                }
            } else if (v != parent && disc[v] < disc[u]) {
                low[u] = min(low[u], disc[v]);
                edgeStack.push({u, v});
            }
        }
    }

    void run() {
        for (auto& [u, _] : g) {
            if (!disc.count(u)) {
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
    }

    const vector<vector<int>>& getBCCs() const { return bccs; }
    const unordered_map<int, int>& getArticulationPoints() const { return artPoints; }
};

// --------------------- Main ---------------------
int main() {
    unordered_map<int, vector<int>> graph = {
        {1, {2, 4}}, {2, {3, 4}}, {3, {1, 4, 10}},
        {4, {5, 6, 8}}, {5, {4}}, {6, {7}}, {7, {8}},
        {8, {9}}, {9, {7}}, {10, {3, 6, 11, 12}},
        {11, {12}}, {12, {10}}, {13, {14}}, {14, {15, 17}},
        {15, {16}}, {16, {13}}, {17, {18}}
    };

    // SCC
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

    // SCC-DAG
    map<int, set<int>> sccGraph;
    for (auto& [u, neighbors] : graph) {
        for (int v : neighbors) {
            int su = nodeToScc[u];
            int sv = nodeToScc[v];
            if (su != sv) sccGraph[su].insert(sv);
        }
    }

    cout << "\nSCC-DAG:\n";
    for (auto& [u, neighbors] : sccGraph) {
        cout << "SuperNode " << u << " -> ";
        for (int v : neighbors) cout << v << " ";
        cout << "\n";
    }

    // WCC
    cout << "\nWCCs:\n";
    auto wccs = computeWCC(graph);
    for (size_t i = 0; i < wccs.size(); ++i) {
        cout << "WCC " << i << ": ";
        for (int node : wccs[i]) cout << node << " ";
        cout << "\n";
    }

    // BCCs inside each SCC
    cout << "\nBCCs in each SCC:\n";
    for (size_t sccID = 0; sccID < sccs.size(); ++sccID) {
        unordered_map<int, vector<int>> subg;
        unordered_set<int> nodeSet(sccs[sccID].begin(), sccs[sccID].end());
        for (int u : sccs[sccID]) {
            for (int v : graph[u]) {
                if (nodeSet.count(v)) {
                    subg[u].push_back(v);
                    subg[v].push_back(u);
                }
            }
        }

        BCC bcc(subg);
        bcc.run();
        cout << "SCC " << sccID << ":\n";
        int bccID = 0;
        for (auto& comp : bcc.getBCCs()) {
            cout << "  BCC " << bccID++ << ": ";
            for (int node : comp) cout << node << " ";
            cout << "\n";
        }
        for (auto& [ap, _] : bcc.getArticulationPoints()) {
            cout << "  Articulation Point: " << ap << "\n";
        }
    }
    return 0;
}
