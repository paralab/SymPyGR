"""
@author: Milinda Fernando
@brief: graph utils for optimized code generation. 
"""

import matplotlib.pyplot as plt



"""
Computes the nodes and rank them based on their in node degree and out node degree threshold values. 
"""
def sorted_nodes_by_in_degree(expr_g,thresh_in,thresh_out=0):
    G=expr_g._G_
    n_list=list()
    for  n in G.nodes():
        if(G.in_degree(n)>=thresh_in  and G.out_degree(n) >thresh_out ):
            n_list.append((G.in_degree(n),n))
    
    n_list.sort(key=lambda  x : x[0],reverse=True)
    return n_list



            