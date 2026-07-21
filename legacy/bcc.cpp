// A C++ program to find biconnected components in a given undirected graph
#include <iostream>
#include <list>
#include <stack> // Though we use list as a stack, <stack> header is good to include for conceptual clarity
#include <algorithm> // For std::min
#include <set> // To store unique nodes in a BCC

#define NIL -1
using namespace std;

int bcc_count = 0; // Renamed 'count' to 'bcc_count' to avoid potential conflicts

class Edge {
public:
    int u;
    int v;
    Edge(int u, int v);
};

Edge::Edge(int u, int v)
{
    this->u = u;
    this->v = v;
}

// A class that represents an undirected graph
class Graph {
    int V; // No. of vertices
    list<int>* adj; // A dynamic array of adjacency lists

    // A Recursive DFS based function used by BCC()
    void BCCUtil(int u, int disc[], int low[],
                 list<Edge>& edge_stack, int parent[]);

public:
    Graph(int V); // Constructor
    void addEdge(int v, int w); // function to add an edge to graph
    void BCC(); // prints biconnected components
};

Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
}

void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v); // For undirected graph
}

// A recursive function that finds and prints biconnected components
// u --> The vertex to be visited next
// disc[] --> Stores discovery times of visited vertices
// low[] -- >> earliest visited vertex (the vertex with minimum
// discovery time) that can be reached from subtree
// rooted with current vertex
// edge_stack -- >> To store visited edges
void Graph::BCCUtil(int u, int disc[], int low[], list<Edge>& edge_stack,
                     int parent[])
{
    // A static variable is used for simplicity, we can avoid use
    // of static variable by passing a pointer.
    static int time = 0;

    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
    int children = 0; // Count of children in DFS tree

    // Go through all vertices adjacent to this
    list<int>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i) {
        int v = *i; // v is current adjacent of 'u'

        // If v is not visited yet, then recur for it
        if (disc[v] == NIL) { // Use NIL for unvisited
            children++;
            parent[v] = u;
            // store the edge in stack
            edge_stack.push_back(Edge(u, v));
            BCCUtil(v, disc, low, edge_stack, parent);

            // Check if the subtree rooted with 'v' has a
            // connection to one of the ancestors of 'u'
            low[u] = min(low[u], low[v]);

            // If u is an articulation point,
            // pop all edges from stack till u -- v
            // A point 'u' is articulation point if following conditions are met
            // (1) u is root of DFS tree and has two or more children.
            // (2) u is not root of DFS tree and low[v] >= disc[u]
            if ((parent[u] == NIL && children > 1) || (parent[u] != NIL && low[v] >= disc[u])) {
                cout << "BCC " << ++bcc_count << ": ";
                set<int> current_bcc_nodes; // Use a set to store unique nodes

                while (true) {
                    Edge current_edge = edge_stack.back();
                    edge_stack.pop_back();
                    cout << current_edge.u << "--" << current_edge.v << " ";
                    current_bcc_nodes.insert(current_edge.u);
                    current_bcc_nodes.insert(current_edge.v);
                    if ((current_edge.u == u && current_edge.v == v) || (current_edge.u == v && current_edge.v == u)) {
                        break; // Found the (u,v) edge that triggered this BCC
                    }
                }
                cout << " (Nodes: ";
                for (int node : current_bcc_nodes) {
                    cout << node << " ";
                }
                cout << ")" << endl;
            }
        }
        // Update low value of 'u' only if 'v' is visited and not parent of 'u'.
        // (i.e. it's a back edge). Also, only consider edges that lead to an
        // ancestor (i.e., discovery_time[v] < discovery_time[u]).
        else if (v != parent[u]) {
            low[u] = min(low[u], disc[v]);
            // Only push back-edges that go up in the DFS tree, to avoid duplicates
            // and ensure correct BCC formation.
            // If v is visited and not parent of u, it must be an ancestor.
            if (disc[v] < disc[u]) { // This check is crucial for correctly identifying back-edges for BCCs
                edge_stack.push_back(Edge(u, v));
            }
        }
    }
}

// The function to do DFS traversal. It uses BCCUtil()
void Graph::BCC()
{
    int* disc = new int[V];
    int* low = new int[V];
    int* parent = new int[V];
    list<Edge> edge_stack; // Correctly declared as a single list

    // Initialize disc and low, and parent arrays
    for (int i = 0; i < V; i++) {
        disc[i] = NIL;
        low[i] = NIL;
        parent[i] = NIL;
    }

    // Call the recursive helper function for all vertices
    // so that all components are traversed
    for (int i = 1; i < V; i++) { // Starting from 1 for 1-based indexing
        if (disc[i] == NIL) {
            BCCUtil(i, disc, low, edge_stack, parent);

            // After a DFS tree is complete, if there are still edges in the stack,
            // it means they form a BCC. This typically happens for the very first
            // BCC found if the graph is fully connected and doesn't have an articulation point
            // that separates components during the main DFS run.
            if (!edge_stack.empty()) {
                cout << "BCC " << ++bcc_count << ": ";
                set<int> current_bcc_nodes;
                while (!edge_stack.empty()) {
                    Edge current_edge = edge_stack.back();
                    edge_stack.pop_back();
                    cout << current_edge.u << "--" << current_edge.v << " ";
                    current_bcc_nodes.insert(current_edge.u);
                    current_bcc_nodes.insert(current_edge.v);
                }
                cout << " (Nodes: ";
                for (int node : current_bcc_nodes) {
                    cout << node << " ";
                }
                cout << ")" << endl;
            }
        }
    }

    delete[] disc;
    delete[] low;
    delete[] parent;
}

// Driver program to test above function
int main()
{
    // The largest node number is 18, so we need V to be at least 19
    // to handle nodes 1 to 18 (node 0 will be unused)
    Graph g(19); // Changed V to 19 for 1-based indexing (nodes 1-18)

    // Add your graph edges
    g.addEdge(1, 2);
    g.addEdge(1, 3);
    g.addEdge(1, 4);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(3, 4);
    g.addEdge(3, 10);
    g.addEdge(4, 5);
    g.addEdge(4, 6);
    g.addEdge(4, 8);
    g.addEdge(6, 7);
    g.addEdge(6, 10);
    g.addEdge(7, 8);
    g.addEdge(7, 9);
    g.addEdge(8, 9);
    g.addEdge(10, 11);
    g.addEdge(10, 12);
    g.addEdge(11, 12);
    g.addEdge(13, 14);
    g.addEdge(13, 16);
    g.addEdge(14, 15);
    g.addEdge(14, 17);
    g.addEdge(15, 16);
    g.addEdge(17, 18);

    cout << "Biconnected Components in the graph:\n";
    g.BCC();
    cout << "\nTotal " << bcc_count << " biconnected components in graph.\n";

    return 0;
}