'''
 BSSN core variables . 
'''
import sys as sys
import dendro
from sympy import *
import dendroutils
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import nxgraph
import gutils

###################################################################
# initialize
###################################################################

l1, l2, l3, l4, eta = symbols('LAMBDA0 LAMBDA1 LAMBDA2 LAMBDA3 ETA')
lf0, lf1 = symbols('LAMBDA_F0 LAMBDA_F1')

# Additional parameters for damping term
R0 = symbols('BSSN_ETA_R0')
ep1, ep2 = symbols('BSSN_ETA_POWER[0] BSSN_ETA_POWER[1]')

xi1, xi2, xi3 = symbols('BSSN_XI[0] BSSN_XI[1] BSSN_XI[2] ')

# declare variables
a   = dendro.scalar("alpha", "")
chi = dendro.scalar("chi", "")
K   = dendro.scalar("K", "")
Gt  = dendro.vec3("Gt", "")
b   = dendro.vec3("beta", "")
B   = dendro.vec3("B", "")
gt  = dendro.sym_3x3("gt", "")
At  = dendro.sym_3x3("At", "")
Gt_rhs  = dendro.vec3("Gt_rhs", "")

# Lie derivative weight
weight = -Rational(2,3)
weight_Gt = Rational(2,3)

# specify the functions for computing first and second derivatives

d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction
kod = dendro.set_kreiss_oliger_dissipation('kograd')

d2 = dendro.d2

dendro.set_metric(gt)
igt = dendro.get_inverse_metric()


eta_func = R0*sqrt(sum([igt[i,j]*d(i,chi)*d(j,chi) for i,j in dendro.e_ij]))/((1-chi**ep1)**ep2)

C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
C2_spatial = dendro.get_complete_christoffel(chi)
[R, Rt, Rphi, CalGt] = dendro.compute_ricci(Gt, chi)

a_rhs = l1*dendro.lie(b, a) - 2*a*K 

b_rhs = [(Rational(3,4) * (lf0 + lf1*a) * B[i] + l2 * dendro.vec_j_ad_j(b, b[i])) for i in dendro.e_i ] 
    
gt_rhs = dendro.lie(b, gt, weight) - 2*a*At 

chi_rhs = dendro.lie(b, chi, weight) + Rational(2,3) * (chi*a*K) 

AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])

At_rhs = dendro.lie(b, At, weight) + chi*dendro.trace_free( a*R - dendro.DiDj(a)) + a*(K*At - 2*AikAkj.reshape(3, 3)) 

K_rhs = dendro.lie(b, K) - dendro.laplacian(a,chi) + a*(K*K/3 + dendro.sqr(At)) 

At_UU = dendro.up_up(At)

Gt_rhs = Matrix([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i]) - \
        Matrix([sum(CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i]) + \
        Rational(2,3)*Matrix([ CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ]) + \
        Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i]) - \
        Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i]) + \
        Matrix([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
        Matrix([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i])
    

Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

B_rhs = [ (Gt_rhs[i] - eta * B[i] +
        l3 * dendro.vec_j_ad_j(b, B[i]) -
        l4 * dendro.vec_j_ad_j(b, Gt[i]))
        for i in dendro.e_i ]

###################################################################
# generate code
###################################################################

outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']

#outs = [a_rhs, b_rhs,gt_rhs]       #gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
#vnames = ['a_rhs', 'b_rhs','gt_rhs'] #, 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']


exp_graph=nxgraph.ExpressionGraph()
exp_graph.add_expressions(outs,vnames)
G=exp_graph.composed_graph(verbose=False)
#G_rev = G.reverse()
dd = nx.greedy_color(G)
#print(type(d))
colors=set()
for (k,v) in dd.items():
        colors.add(v)

print(colors)


#cse_nodes = gutils.sorted_nodes_by_in_degree(G,10)

# for (k,v) in cse_nodes:
#         print("node in degree %d out degree %d\n expression : %s\n\n" %(G.in_degree(v),G.out_degree(v),str(v)))

# for  n in G.nodes():
#         # if(G.in_degree(n)==0):
#         #         print(n)
#         if(G.in_degree(n)>10) and G.out_degree(n)>0:
#                print("node in degree %d out degree %d\n expression : %s\n\n" %(G.in_degree(n),G.out_degree(n),str(n)))
#                 #for v in G_rev.neighbors(n):
                #        print("used in %s " %v)
#nx.draw_networkx(G,pos=nx.random_layout(G),font_size=4)
#nx.draw(G,pos=nx.random_layout(G))
#plt.show()

#exp_graph.draw_graph("gt_rhs_00")
#exp_graph.plot_adjmatrix()

#G.add_expression(a_rhs,"a_rhs")

#dendroutils.construct_dep_matrix(outs,vnames,"[pp]")
#dendroutils.write_to_dot(outs,vnames,folder_ptah="../dot")
# a_G  = dendroutils.construct_nx_digraph("../dot/bkp/At_rhs_00.dot")
# a_G1 = exp_graph.get_graph("gt_rhs_00")

# print(nx.info(a_G))
# print(nx.info(a_G1))
# b_0_G  = dendroutils.construct_nx_digraph("../dot/b_rhs_0.dot")
# b_1_G  = dendroutils.construct_nx_digraph("../dot/b_rhs_1.dot")
# b_2_G  = dendroutils.construct_nx_digraph("../dot/b_rhs_2.dot")

# gt_00_G  = dendroutils.construct_nx_digraph("../dot/gt_rhs_00.dot")
# gt_01_G  = dendroutils.construct_nx_digraph("../dot/gt_rhs_01.dot")
# gt_02_G  = dendroutils.construct_nx_digraph("../dot/gt_rhs_02.dot")
# gt_11_G  = dendroutils.construct_nx_digraph("../dot/gt_rhs_11.dot")
# gt_12_G  = dendroutils.construct_nx_digraph("../dot/gt_rhs_12.dot")
# gt_22_G  = dendroutils.construct_nx_digraph("../dot/gt_rhs_22.dot")

# g_list=[a_G,b_0_G,b_1_G,b_2_G,gt_00_G,gt_01_G,gt_02_G,gt_11_G,gt_12_G,gt_22_G]
# G=nx.compose_all(g_list)
# dendroutils.bfs_traversal(G,gt_12_G)
# nx.write_multiline_adjlist(G.reverse(),'test_G_reverse.txt',delimiter="\n")
# print(nx.classes.function.info(G))
# A=nx.adj_matrix(G.reverse())
# plt.figure()
# plt.spy(A,markersize=3)
# plt.show()
#Gr = nx.classes.function.reverse_view(G)
#tc_Gr = nx.algorithms.dag.transitive_closure_dag(Gr)
#dendroutils.draw_nx_graph(G,draw_labels=False)

#print(type(set(nx.classes.function.nodes(G))))
#for n in nx.classes.function.nodes(G):
#        print(n)
#        print("neigh num: %s" % (nx.classes.function.neighbors(G,n)))


#for e in nx.classes.function.edges(G):
#        print(e)




