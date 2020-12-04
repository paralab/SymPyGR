# @author : Milinda Fernando (milinda@cs.utah.edu)
# School of Computing, University of Utah
# contain utility functions for optimized code generation for given sympy expressions.

import dendro
import sympy
import numpy as np
import networkx as nx
from sympy.printing.dot import dotprint
import matplotlib.pyplot as plt




"""
replace user defined functions with sympy symbols. 
mainly written to replace the derivative functions with corresponding symbol name. 

"""
def replace_userdef_funcs(expr):
    substitute_expr=dict()
    for item in expr.atoms(sympy.Function): ##in sympy.preorder_traversal(expr):
        if isinstance(item, sympy.core.function.AppliedUndef):
            sym_name=str(item.func)
            for a in item.args:
                sym_name = sym_name + '_' + str(a)
            #print(sym_name)
            #expr=expr.replace(item,sympy.Symbol(sym_name))
            substitute_expr.update({item:sympy.Symbol(sym_name)})
            #expr=expr.subs(item,sympy.Symbol(sym_name))

    #print(expr.free_symbols)
    #print(substitute_expr)
    #expr.subs(substitute_expr)
    #expr.replace()
    #print(expr)
    for k,v in substitute_expr.items():
        print("replacing function: %s with symbol %s: " %(k,v))
        expr=expr.replace(k,v)

    return expr

"""
advanced routine to free symbols in presence of user defined functions. 
@param: sympy expression. 
"""

def advanced_free_symbols(expr):
    sym_set=set()
    
    for item in expr.atoms(sympy.Function):
        sym_name=str(item.func)
        for a in item.args:
            sym_name = sym_name + '_' + str(a)

        print(item)
        sym_set.add(sympy.Symbol(sym_name))

    for item in expr.atoms(sympy.Symbol):
        print(item)
        sym_set.add(item)
    
    return sym_set
    

"""
Write sympy expression to a .dot(graphViz file format)
"""

def write_to_dot(outs, vnaems, suffix="[pp]", folder_ptah="."):
    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']
    idx=suffix
    
    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    for i, e in enumerate(outs):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                print("processing expr : %d var name %s[%s]" %(i,vnaems[i],str(j)))
                ev=replace_userdef_funcs(ev)
                d_str=str(dotprint(ev,labelfunc=sympy.srepr))
                gv_file = open(folder_ptah+"/"+vnaems[i]+"_"+str(j)+".dot",'w')
                gv_file.write(d_str)
                gv_file.close()
        elif type(e) == sympy.Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                print("processing expr : %d var name %s[%s]" %(i,vnaems[i],midx[j]))
                e[k]=replace_userdef_funcs(e[k])
                d_str=str(dotprint(e[k],labelfunc=sympy.srepr))
                gv_file = open(folder_ptah+"/"+vnaems[i]+"_"+str(midx[j])+".dot",'w')
                gv_file.write(d_str)
                gv_file.close()
                #exp_symbols = exp_symbols.union(replace_userdef_funcs(e[k]).free_symbols)
                #exp_symbols = exp_symbols.union(advanced_free_symbols(e[k]))
        else:
            num_e = num_e + 1
            print("processing expr : %d var name %s" %(i,vnaems[i]))
            e=replace_userdef_funcs(e)
            #d_str=str(dotprint(e,labelfunc=sympy.srepr))
            d_str=str(dotprint(e))
            gv_file = open(folder_ptah+"/"+vnaems[i]+".dot",'w')
            gv_file.write(d_str)
            gv_file.close()
            #exp_symbols = exp_symbols.union(replace_userdef_funcs(e).free_symbols)
            #exp_symbols = exp_symbols.union(advanced_free_symbols(e))


"""
Construct a networkX digraph from a dot file. 
"""
def construct_nx_digraph(file_name):
    nxgraph = nx.MultiDiGraph(nx.drawing.nx_pydot.read_dot(file_name))
    return nxgraph
    #print(nxgraph)
    #print(nx.linalg.graphmatrix.adj_matrix(nxgraph))
    #print(nxgraph.adj)
    #nx.draw(nxgraph,pos=nx.planar_layout(nxgraph))
    #nx.draw_networkx(nxgraph,pos=nx.planar_layout(nxgraph),font_size=8)
    #plt.show()
    #print(nx.linalg.graphmatrix.adj_matrix(nxgraph))

"""
construct the dependancy matrix for a given list of expressions
out     : symbolic expression list
vnames  : variable names for the symbolic expression
suffix  : suffix for append to vnames. 
"""
def construct_dep_matrix(outs, vnames, suffix="[pp]"):
    
    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']
    idx=suffix
    
    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp  = list()
    lname = list()
    exp_symbols = set()
    
    for i, e in enumerate(outs):
        print("processing expr : %d" %i)
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                #exp_symbols = exp_symbols.union(replace_userdef_funcs(ev).free_symbols)
                #print(exp_symbols)
                ev=replace_userdef_funcs(ev)
                #exp_symbols = exp_symbols.union(advanced_free_symbols(ev))
        elif type(e) == sympy.Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                e[k]=replace_userdef_funcs(e[k])
                #exp_symbols = exp_symbols.union(replace_userdef_funcs(e[k]).free_symbols)
                #exp_symbols = exp_symbols.union(advanced_free_symbols(e[k]))
        else:
            num_e = num_e + 1
            e=replace_userdef_funcs(e)
            #exp_symbols = exp_symbols.union(replace_userdef_funcs(e).free_symbols)
            #exp_symbols = exp_symbols.union(advanced_free_symbols(e))

    print(num_e)
    # for i,v in enumerate(lname):
    #     print("var : %s expr : %s" %(v,lexp[i]))
    

    print(exp_symbols)
    print(len(exp_symbols))


    
    



