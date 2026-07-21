from collections import defaultdict

class Graph:
    def __init__(self, vertices):
        self.V = vertices  # 顶点数量
        self.graph = defaultdict(list)  # 邻接列表表示图
        self.time = 0  # 全局时间戳，用于记录发现时间

    # 添加边
    def add_edge(self, u, v):
        self.graph[u].append(v)
        self.graph[v].append(u)

    # 递归的 DFS 函数来寻找割点
    def find_articulation_points_util(self, u, visited, parent, discovery_time, low_link, articulation_points):
        visited[u] = True
        self.time += 1
        discovery_time[u] = self.time
        low_link[u] = self.time
        
        children = 0 # 记录当前节点的子节点数量（在DFS树中）

        # 遍历所有邻接节点
        for v in self.graph[u]:
            # 如果 v 没有被访问过
            if not visited[v]:
                children += 1
                parent[v] = u # 设置 v 的父节点为 u
                self.find_articulation_points_util(v, visited, parent, discovery_time, low_link, articulation_points)

                # 检查以 v 为根的子树是否能回溯到 u 或其祖先
                # 如果 low_link[v] >= discovery_time[u]，说明 u 是割点
                low_link[u] = min(low_link[u], low_link[v]) # 更新 u 的low_link
                if parent[u] == -1 and children > 1: # 情况 1: u 是DFS树的根节点且有多个子节点
                    articulation_points.add(u)
                elif parent[u] != -1 and low_link[v] >= discovery_time[u]: # 情况 2: u 不是根节点，且满足割点条件
                    articulation_points.add(u)
            elif v != parent[u]: # 如果 v 已经被访问过且不是 u 的父节点，说明是回边
                low_link[u] = min(low_link[u], discovery_time[v]) # 更新 u 的low_link

    # 主函数：找到图中的所有割点
    def find_articulation_points(self):
        visited = [False] * (self.V + 1) # 访问标记数组
        discovery_time = [-1] * (self.V + 1) # 发现时间
        low_link = [-1] * (self.V + 1) # 回溯最低时间
        parent = [-1] * (self.V + 1) # 父节点
        articulation_points = set() # 存储割点

        # 遍历所有节点，以防图不连通
        for i in range(1, self.V + 1): # 假设节点从1开始编号
            if not visited[i]:
                self.find_articulation_points_util(i, visited, parent, discovery_time, low_link, articulation_points)
        
        return articulation_points

# --- 示例用法 ---
if __name__ == "__main__":
    # 使用您提供的图数据来测试
    # 最大的节点号是18，所以顶点数量设置为18
    g = Graph(18) 
    g.add_edge(1, 2)
    g.add_edge(1, 3)
    g.add_edge(1, 4)
    g.add_edge(2, 3)
    g.add_edge(2, 4)
    g.add_edge(3, 4)
    g.add_edge(3, 10)
    g.add_edge(4, 5)
    g.add_edge(4, 6)
    g.add_edge(4, 8)
    g.add_edge(6, 7)
    g.add_edge(6, 10)
    g.add_edge(7, 8)
    g.add_edge(7, 9)
    g.add_edge(8, 9)
    g.add_edge(10, 11)
    g.add_edge(10, 12)
    g.add_edge(11, 12)
    g.add_edge(13, 14)
    g.add_edge(13, 16)
    g.add_edge(14, 15)
    g.add_edge(14, 17)
    g.add_edge(15, 16)
    g.add_edge(17, 18)

    print("图中的割点是:", g.find_articulation_points())