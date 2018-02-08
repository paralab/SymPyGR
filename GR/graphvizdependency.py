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


def getAllArguments(expression, result):
    args = expression.args
    for arg in args:
        if(len(arg.args)==0):
            if(str(arg) not in result):
                result.append(str(arg))
        else:
            getAllArguments(arg, result)

    return result

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
                getAllArguments(arg, result)

    return result

def makeDependencies(_v):
    dependencies = {}
    counter = 1
    for i in _v[0]:
        result = []
        result = getAllArguments(i[1], result)
        dependencies[str(i[0])] = result
    for i in _v[1]:
        result = []
        result = getAllArguments(i, result)
        dependencies["Equation" + str(counter)+": "] = result
        counter= counter+1
    return dependencies

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]

def makeCompleteDependencies(_v):
    dependencies = {}
    counter = 1
    for i in _v[0]:
        result = []
        result = getReducedArguments(i[1], result)
        for j in result:
            if(str(j) in dependencies.keys()):
                externDepen = dependencies[str(j)]
                for k in externDepen:
                    if(str(k) not in result):
                        result.append(str(k))
                result = remove_values_from_list(result, str(j))
        dependencies[str(i[0])] = result
    for i in _v[1]:
        result = []
        result = getReducedArguments(i, result)
        for j in result:
            if (str(j) in dependencies.keys()):
                externDepen = dependencies[str(j)]
                for k in externDepen:
                    if (str(k) not in result):
                        result.append(str(k))
                result = remove_values_from_list(result, str(j))
        dependencies["Equation" + str(counter)+": "] = result
        counter= counter+1
    return dependencies

def printDependencies(dependencies):
    for i in dependencies.keys():
        print(str(i)+ str(dependencies[i]))
        print()


def getDependency(dependencies, var_name):
    return dependencies[var_name]