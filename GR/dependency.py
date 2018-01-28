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
from graphviz import Digraph
import os

#Should install graphviz  and add the bin path to the environment path variable
os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'

variable_counter = 0
dot = Digraph(comment='Dependency Graph')
nodes = []

def generate_dependency_graph(expression):
    global dot
    global nodes
    global variable_counter
    args = expression.args
    func = expression.func
    variable_counter = variable_counter + 1
    node_name = "A" + str(variable_counter)
    dot.node( node_name, str(func))
    for arg in args:
        if(len(arg .args)==0):
            if(str(arg).startswith("DENDRO") and  (str(arg) in nodes)):
                dot.edge(node_name, str(arg))
            elif (str(arg).startswith("DENDRO") and (str(arg) not in nodes)):
                dot.node(str(arg), str(arg))
                nodes.append(str(arg))
                dot.edge(node_name, str(arg))
            else:
                variable_counter = variable_counter + 1
                dot.node(str(arg)+str(variable_counter), str(arg))
                dot.edge(node_name, str(arg))
        else:
            dot.edge(node_name, generate_dependency_graph(arg))
    return node_name

def drawgragh(_v):
    # Making the dependency graph

    for i in _v[0]:
        expr = Eq(i[0], i[1])
        generate_dependency_graph(expr)

    eqn_number = 1

    for i in _v[1]:
        expr = Eq(symbols("Equation" + str(eqn_number)), i)
        generate_dependency_graph(expr)
        eqn_number = eqn_number + 1

    global dot
    print(dot.source)
    dot.render('dependency.gv', view=True)



