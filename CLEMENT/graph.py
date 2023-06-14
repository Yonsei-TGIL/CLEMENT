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
        if v2 not in self.graph[v1]:          # Enlist v2 the child list to v1  (self.graph[v1])
            self.graph[v1].append(v2)
    def removeEdge(self, v1, v2):
        self.graph[v1].remove(v2)

    def findparent (self, v1):
        for someone in list(self.graph.keys()):
            if v1 in self.graph[someone] :       # someone is parent if v1 is connected with someone (self.graph[someone])
                return someone
        return "None"

    def intervene (self, v1, v2_list):
        for previous_v1 in self.graph.keys():
            their_v2 = list (self.graph[previous_v1])       # List of child of previous_v1
            if all(item in their_v2 for item in v2_list ) == True:       # if the number of child of previous_v1 > the number of child (v2_list)
                print ("\t\tGrandparents : {},  Parents : {}, Child : {}".format(previous_v1, v1, v2_list))
                self.addEdge(previous_v1, v1)     #  add grandparent - parent relationship
                for v2 in v2_list:
                    self.removeEdge (previous_v1, v2)    # remove the edge that directly links grandparent - child
                break

    def _dfs (self, node, discovered, footage, PHYLOGENY_DIR):          # discovered : list
        discovered.append (node)
        footage.append (node)
        for w in self.graph[node]:     # Iterate the child list
            if w not in discovered:     # If undiscovered child yet,
                self._dfs (w, discovered, footage, PHYLOGENY_DIR)

        if len(self.graph[node]) == 0:   # Print only when  terminal node
            print (" → ".join([str(i) for i in footage]))
            with open (PHYLOGENY_DIR, "a", encoding = "utf8") as output_file:
                print (" → ".join([str(i) for i in footage]), file = output_file)
        footage.remove (node)   # Remove the node that had been
    

    def dfs (self, root, PHYLOGENY_DIR):
        discovered = self._dfs(root, [], [], PHYLOGENY_DIR)         # Present place is the root. Let's iterate.
        return discovered 

    def print (self):
        for v1 in self.graph.keys():
            for v2 in self.graph[v1]:
                print (v1, " - ", v2)

def main ():
    print ("Hi")