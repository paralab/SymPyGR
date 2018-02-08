##########################################################################
# module: dendro
# author: Pasindu Tennage
# email:  pasindu.13@cse.mrt.ac.lk
#
# python module to generate dependency graph
#
# (c) 2016 University of Utah, All rights reserved.
##########################################################################

from sympy import *
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering
from sklearn import metrics

def printGraph(G):
    nx.draw(G)
    plt.show()
    plt.savefig("dependency-graph.png")

def getReducedArguments(expression, result):
    args = expression.args
    func = expression.func
    if("grad" in str(func)):
        result.append(str(expression))
    else:
        for arg in args:
            if(len(arg.args)==0):
                if(str(arg) not in result):
                    result.append(str(arg))
            else:
                getReducedArguments(arg, result)
    return result

def makeDependencies(_v):
    G=nx.Graph()
    dependencies = {}
    counter = 1
    for i in _v[0]:
        result = []
        result = getReducedArguments(i[1], result)
        G.add_node(str(i[0]))
        for arg in result:
            G.add_node(str(arg))
            G.add_edge(str(i[0]),str(arg))
    for i in _v[1]:
        G.add_node("Equation" + str(counter))
        result = []
        result = getReducedArguments(i, result)
        for arg in result:
            G.add_node(str(arg))
            G.add_edge("Equation" + str(counter),str(arg))
        counter = counter+1
    return G

def getGraphLaplacian(G):
    return nx.normalized_laplacian_matrix(G)

def doSpectralClustering(graphLaplacian):
    # Cluster
    sc = SpectralClustering()
    sc.fit(graphLaplacian)
    return sc