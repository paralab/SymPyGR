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
        if isinstance(item.func, sympy.core.function.UndefinedFunction):
            sym_name=str(item.func)
            for a in item.args:
                sym_name = sym_name + '_' + str(a)
            #print(item)
            sym_set.add(sympy.Symbol(sym_name))

    for item in expr.atoms(sympy.Symbol):
        #print(item)
        sym_set.add(item)
    
    return sym_set
    

def graph_label_func(expr):
    return str(expr.func)


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
                d_str=str(dotprint(ev,labelfunc=sympy.srepr,repeat=False))
                gv_file = open(folder_ptah+"/"+vnaems[i]+"_"+str(j)+".dot",'w')
                gv_file.write(d_str)
                gv_file.close()
        elif type(e) == sympy.Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                print("processing expr : %d var name %s[%s]" %(i,vnaems[i],midx[j]))
                e[k]=replace_userdef_funcs(e[k])
                d_str=str(dotprint(e[k],labelfunc=sympy.srepr,repeat=False))
                gv_file = open(folder_ptah+"/"+vnaems[i]+"_"+str(midx[j])+".dot",'w')
                gv_file.write(d_str)
                gv_file.close()
                #exp_symbols = exp_symbols.union(replace_userdef_funcs(e[k]).free_symbols)
                #exp_symbols = exp_symbols.union(advanced_free_symbols(e[k]))
        else:
            num_e = num_e + 1
            print("processing expr : %d var name %s" %(i,vnaems[i]))
            e=replace_userdef_funcs(e)
            d_str=str(dotprint(e,labelfunc=sympy.srepr,repeat=False))
            gv_file = open(folder_ptah+"/"+vnaems[i]+".dot",'w')
            gv_file.write(d_str)
            gv_file.close()
            #exp_symbols = exp_symbols.union(replace_userdef_funcs(e).free_symbols)
            #exp_symbols = exp_symbols.union(advanced_free_symbols(e))


"""
Construct a networkX digraph from a dot file. 
"""
def construct_nx_digraph(file_name):
    G = nx.DiGraph(nx.drawing.nx_pydot.read_dot(file_name))
    return G
   

"""
draw networkX graph. 
"""

def draw_nx_graph(G,draw_labels=False):
    if( not draw_labels):
        nx.draw(G,pos=nx.circular_layout(G))
        plt.show()
    else:
        nx.draw_networkx(G,pos=nx.random_layout(G),font_size=8)
        plt.show()


"""
Breadth first traversal from the root expr.
"""
def bfs_traversal(G,g=None):
    nodes = iter(nx.nodes(g)) 
    root  = next(nodes) 
    #print(G.degree(root))

    #for n in nx.classes.function.all_neighbors(G,root):
    #    print(n)

    bfs_iter = dict(nx.bfs_successors(G,root))
    #print(bfs_iter)
    for n in bfs_iter:
        print("node %s has %d children" %(n,len(bfs_iter[n])))
        e=sympy.parse_expr(n)
        if(len(e.free_symbols)>1):
            print("node : %s \n \t\t has  dep %d  distinct symbols \n" %(n, len(e.free_symbols)))
        # for child in bfs_iter[n]:
        #     e=sympy.parse_expr(child)
        #     #print(e.free_symbols)
        #     if(len(e.free_symbols)>1):
        #         print("node : %s \n has a child node that violate the constraint dep %d \n child node :  %s\n" %(n, len(e.free_symbols), child))
            
        #for s in bfs_iter[n]:
        #    print("node %s has suc %s " %(n,s))



"""
extract all the expressions from quantities such as vectors and tensors
"""
def extract_expressions(outs,vnames,suffix="[pp]"):
    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']

    expr_dict=dict()
    num_e = 0

    for i, e in enumerate(outs):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                expr_name = vnames[i] + "_" + str(j) + suffix
                expr_dict[expr_name] = ev
        elif type(e) == sympy.Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                expr_name = vnames[i] + "_" +str(midx[j]) + suffix
                expr_dict[expr_name] = e[k]

        else:
            num_e = num_e + 1
            expr_name = vnames[i] + suffix
            expr_dict[expr_name] = e

    return expr_dict

