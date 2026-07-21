from collections import defaultdict, deque

def compute_wcc(graph):
    """
    Compute Weakly Connected Components (WCCs) of an undirected graph.
    
    Args:
        graph (dict): Adjacency list of the graph.

    Returns:
        list of list: List of WCCs, each is a list of nodes.
    """
    visited = set()
    wccs = []

    def bfs(start):
        queue = deque([start])
        component = []
        visited.add(start)
        while queue:
            node = queue.popleft()
            component.append(node)
            for neighbor in graph[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        return component

    for node in graph:
        if node not in visited:
            wcc = bfs(node)
            wccs.append(wcc)

    print(f"\n--- Detected {len(wccs)} Weakly Connected Components ---")
    for i, comp in enumerate(wccs):
        print(f"Component {i + 1}: {comp}")
    return wccs

def brandes_betweenness_centrality_optimized(graph, wccs):
    """
    Compute Betweenness Centrality using Brandes' algorithm,
    optimized by limiting to nodes within the same WCC.

    Args:
        graph (dict): Adjacency list of the full graph.
        wccs (list of list): List of WCCs, each is a list of nodes.

    Returns:
        dict: Betweenness centrality scores for all nodes.
    """
    bc = {node: 0.0 for node in graph}  # Initialize BC for all nodes to 0
    total_nodes = len(graph)

    print(f"\n--- Optimized Brandes' Algorithm with WCCs ---")
    print(f"Total number of nodes: {total_nodes}")
    print(f"Number of WCCs: {len(wccs)}")

    for idx, component in enumerate(wccs):
        print(f"\n>>> Processing WCC #{idx + 1} with {len(component)} nodes: {component}")

        # Create a subgraph view
        subgraph = {node: [nbr for nbr in graph[node] if nbr in component] for node in component}

        for s in component:  # Only loop inside this component
            print(f"\n--- Starting BFS from source node '{s}' in WCC #{idx + 1} ---")

            S = []
            P = defaultdict(list)
            dist = {node: -1 for node in component}
            sigma = {node: 0 for node in component}
            delta = {node: 0.0 for node in component}

            dist[s] = 0
            sigma[s] = 1
            queue = deque([s])

            while queue:
                v = queue.popleft()
                S.append(v)
                for w in subgraph[v]:
                    if dist[w] == -1:
                        dist[w] = dist[v] + 1
                        queue.append(w)
                    if dist[w] == dist[v] + 1:
                        sigma[w] += sigma[v]
                        P[w].append(v)

            print(f"  Visited nodes from '{s}': {S}")
            print(f"  Distance map: {dist}")
            print(f"  Path count (sigma): {sigma}")

            while S:
                w = S.pop()
                for v in P[w]:
                    delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w])
                if w != s:
                    bc[w] += delta[w]

        print(f"<<< Finished WCC #{idx + 1}, partial BC scores: {[f'{n}: {bc[n]:.2f}' for n in component]}")

    print("\n--- Final Betweenness Centrality for All Nodes ---")
    for node, value in bc.items():
        print(f"Node '{node}': BC = {value:.2f}")

    return bc

# ------------------------------
# Example Usage
# ------------------------------

if __name__ == "__main__":
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

    # Step 1: Compute Weakly Connected Components
    wccs = compute_wcc(graph_example)

    # Step 2: Run Optimized Brandes' Algorithm
    betweenness_scores = brandes_betweenness_centrality_optimized(graph_example, wccs)

    # Step 3: Final Output
    print("\n--- Final Betweenness Centrality Result (Dictionary Format) ---")
    print(betweenness_scores)
