# -*- coding: utf-8 -*-
from collections import defaultdict

class TarjanSCC:
    def __init__(self, graph):
        self.graph = graph
        self.index = 0
        self.dfn = {}         # discovery time of each node
        self.low = {}         # low-link value
        self.stack = []       # stack for Tarjan
        self.on_stack = set() # to check whether a node is in the current stack
        self.sccs = []        # list of SCCs

    def run(self):
        for node in self.graph:
            if node not in self.dfn:
                self._dfs(node)
        return self.sccs

    def _dfs(self, u):
        self.dfn[u] = self.low[u] = self.index
        self.index += 1
        self.stack.append(u)
        self.on_stack.add(u)

        for v in self.graph.get(u, []):
            if v not in self.dfn:
                self._dfs(v)
                self.low[u] = min(self.low[u], self.low[v])
            elif v in self.on_stack:
                self.low[u] = min(self.low[u], self.dfn[v])

        if self.dfn[u] == self.low[u]:
            scc = []
            while True:
                node = self.stack.pop()
                self.on_stack.remove(node)
                scc.append(node)
                if node == u:
                    break
            self.sccs.append(sorted(scc))

# 手动构建图（邻接表）
graph = {
    1: [2, 4],
    2: [3, 4],
    3: [1, 4, 10],
    4: [5, 6, 8],
    5: [4],
    6: [7],
    7: [8],
    8: [9],
    9: [7],
    10: [3, 6, 11],
    11: [12],
    12: [10],
    13: [14],
    14: [15, 17],
    15: [16],
    16: [13],
    17: [18],
    18: []
}

# 运行 Tarjan 算法
tarjan = TarjanSCC(graph)
scc_result = tarjan.run()

# 输出 SCCs
print("Strongly Connected Components (SCCs):")
for i, comp in enumerate(scc_result, 1):
    print(f"SCC {i}: {comp}")
