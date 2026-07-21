from collections import defaultdict, deque

def brandes_betweenness_centrality(graph):
    """
    Compute Betweenness Centrality for each node in the graph using Brandes' algorithm.

    Args:
        graph (dict): Adjacency list representation of the graph, e.g.:
                      {
                          'A': ['B', 'D'],
                          'B': ['A', 'C', 'E'],
                          'C': ['B', 'E'],
                          'D': ['A', 'E'],
                          'E': ['B', 'C', 'D'],
                          'F': ['I'],
                          'I': ['F'],
                          'G': ['H'],
                          'H': ['G'],
                          'J': ['K'],
                          'K': ['J']
                      }
                      Assumes an unweighted, undirected graph.

    Returns:
        dict: A dictionary where keys are nodes and values are their betweenness centrality.
    """
    nodes = list(graph.keys())
    num_nodes = len(nodes)
    bc = {node: 0.0 for node in nodes}  # Initialize BC for all nodes to 0

    print("--- Computing Betweenness Centrality using Brandes' Algorithm ---")
    print(f"Nodes in the graph: {nodes}\n")

    for s in nodes:  # Use each node as a source
        print(f"\n--- Starting Dijkstra from source node '{s}' ---")
        
        S = []  # Stack to record order of visits
        P = defaultdict(list)  # Predecessors on shortest paths
        dist = {node: -1 for node in nodes}  # Distance from source s
        sigma = {node: 0 for node in nodes}  # Number of shortest paths

        # Initialize BFS
        dist[s] = 0
        sigma[s] = 1
        queue = deque([s])

        while queue:
            v = queue.popleft()
            S.append(v)  # Push to stack for later back-propagation

            for w in graph[v]:  # Traverse neighbors
                if dist[w] == -1:  # If not yet visited
                    dist[w] = dist[v] + 1
                    queue.append(w)
                if dist[w] == dist[v] + 1:  # If w is a valid successor on a shortest path
                    sigma[w] += sigma[v]
                    P[w].append(v)

        print(f"  Shortest distances from source '{s}': {dist}")
        print(f"  Number of shortest paths from source '{s}' (sigma): {sigma}")
        print(f"  Predecessors for each node (P): {P}")
        print(f"  Visit order stack (S): {S}\n")

        delta = {node: 0.0 for node in nodes}  # Dependency values

        # Backward accumulation
        while S:
            w = S.pop()
            print(f"  --- Processing node '{w}' (popped from stack) ---")
            for v in P[w]:
                contribution = (sigma[v] / sigma[w]) * (1 + delta[w])
                delta[v] += contribution
                print(f"    Contribution from '{v}' to '{w}': sigma[{v}]={sigma[v]}, sigma[{w}]={sigma[w]}, delta[{w}]={delta[w]:.2f}")
                print(f"    Contribution = ({sigma[v]} / {sigma[w]}) * (1 + {delta[w]:.2f}) = {contribution:.2f}")
                print(f"    Updated delta[{v}] = {delta[v] - contribution:.2f} + {contribution:.2f} = {delta[v]:.2f}")
            
            if w != s:
                bc[w] += delta[w]
                print(f"    '{w}' is not the source '{s}', added delta[{w}] = {delta[w]:.2f} to BC[{w}]")
                print(f"    Updated BC[{w}] = {bc[w]:.2f}\n")
            else:
                print(f"    '{w}' is the source node '{s}', not added to BC.\n")

    # Normalization can be added here if needed (not performed in this version)

    print("\n--- Final Betweenness Centrality for All Nodes ---")
    for node, value in bc.items():
        print(f"Node '{node}': BC = {value:.2f}")

    return bc

# Example graph
graph_example = {
    'A': ['B', 'D'],
    'B': ['A', 'C', 'E'],
    'C': ['B', 'E'],
    'D': ['A', 'E'],
    'E': ['B', 'C', 'D'],
    'F': ['I'],
    'I': ['F'],
    'G': ['H'],
    'H': ['G'],
    'J': ['K'],
    'K': ['J']
}

# Run the algorithm
betweenness_scores = brandes_betweenness_centrality(graph_example)

# Print final result (already printed inside the function)
print("\n--- Final Betweenness Centrality Result (Dictionary Format) ---")
print(betweenness_scores)
