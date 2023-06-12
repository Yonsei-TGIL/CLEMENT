class Graph:
    def __init__(self):
        self.graph = {}
    def addVertex(self, v):
        if v not in self.graph.keys():
            self.graph[v] = []
    def getNode (self):
        return list(self.graph.keys())
    def getEdges(self):
        edges = []
        for v1 in self.graph.keys():
            for v2 in self.graph[v1]:
                edges.append ((v1,v2))
        return edges

class UnidirectedGraph(Graph):
    def addEdge(self, v1, v2):
        if v1 not in self.graph.keys():
            self.addVertex(v1)
        if v2 not in self.graph.keys():
            self.addVertex(v2)
        if v2 not in self.graph[v1]:          # v1이 보유하고 있는 자식 명단에 없다면 자식으로 추가시켜주자
            self.graph[v1].append(v2)
    def removeEdge(self, v1, v2):
        self.graph[v1].remove(v2)

    def findparent (self, v1):
        for someone in list(self.graph.keys()):
            if v1 in self.graph[someone] :       # 누군가의 자식 명단에 포함돼 있다면 someone을 부모라고 생각할 수 있다
                return someone
        return "None"

    def intervene (self, v1, v2_list):
        for previous_v1 in self.graph.keys():
            their_v2 = list (self.graph[previous_v1])       # 자식 명단
            if all(item in their_v2 for item in v2_list ) == True:       # 만약 previous_v1의 자식 명단이 내가 가진 자식 명단을 초과할 때
                print ("Grandparents : {},  Parents : {}, Child : {}".format(previous_v1, v1, v2_list))
                self.addEdge(previous_v1, v1)     #  할아버지 - 아버지 관계로 입적시키기
                for v2 in v2_list:
                    self.removeEdge (previous_v1, v2)    # 할아버지 -손자 관계는 직접적인 edge를 삭제하기
                break

    def _dfs (self, node, discovered, footage, PHYLOGENY_DIR):          # discovered : list
        discovered.append (node)
        footage.append (node)
        for w in self.graph[node]:     # 현재 자식을 다 돌자
            if w not in discovered:    # 현재까지 밝혀지지 않은 자식이라면
                self._dfs (w, discovered, footage, PHYLOGENY_DIR)

        if len(self.graph[node]) == 0:   # 말단 node일 때에만 출력
            print ("→".join([str(i) for i in footage]))
            with open (PHYLOGENY_DIR, "a", encoding = "utf8") as output_file:
                print ("→".join([str(i) for i in footage]), file = output_file)
        footage.remove (node)   # 갔다 온 node는 제거한다
    

    def dfs (self, root, PHYLOGENY_DIR):
        discovered = self._dfs(root, [], [], PHYLOGENY_DIR)         # 현재를 root로 삼아서 무한탐색을 돌자
        return discovered 

    def print (self):
        for v1 in self.graph.keys():
            for v2 in self.graph[v1]:
                print (v1, " - ", v2)

def main ():
    print ("Hi")