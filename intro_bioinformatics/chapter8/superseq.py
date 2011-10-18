# solve the smallest supersequence problem
# using a

import sys

# DAG
class Graph(object):

    def __init__(self, Vertex = {}, Edge = {}, DAG = False):
        self.Vertex = Vertex
        self.Edge = Edge
        self.DAG = DAG

    def hasVertex(self, vertex):
        return vertex in self.Vertex

    def hasEdge(self, v1, v2):
        ret = 0
        if v1 in Edge[v2]:
            ret = -1
        if v2 in Edge[v1]:
            ret = 1

        if self.DAG:
            return abs(ret)
        else:
            return ret

    def addVertex(self, vertex):
        if not hasVertex(vertex):
            self.Vertex[vertex] = 1
            self.Edge[vertex] = {}

    def addEdge(self, v1, v2):
        if not hasEdge(v1, v2):
            self.Edge[v1][v2] = 1
            if not DAG:
                self.Edge[v2][v1] = 1



# print all l-mers of a sequence
def print_lmer(seq, l):
    for i in range(0, len(seq)-l+1):
        print i*' '+str(seq[i:i+l])
    print




if __name__ == "__main__":

    # the shortest supersequence for all 'aab' like 3-mers
    s = 'abbbabaaab'
    t = 'abcabdabeab'

    # the assembly of the two strands gives an overall shortest supersequence
    # print_lmer(s + t[2:], 3)

    print_lmer('1011100010', 3)
    print_lmer('GTATGGGTG', 3)
    print_lmer('GTGGGTATG', 3)
