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
import re as regex

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
# a   = dendro.scalar("alpha", "")
# chi = dendro.scalar("chi", "")
# K   = dendro.scalar("K", "")
# Gt  = dendro.vec3("Gt", "")
# b   = dendro.vec3("beta", "")
# B   = dendro.vec3("B", "")
# gt  = dendro.sym_3x3("gt", "")
# At  = dendro.sym_3x3("At", "")
# Gt_rhs  = dendro.vec3("Gt_rhs", "")

a   = dendro.scalar("alpha"   , "")
chi = dendro.scalar("chi"     , "")
K   = dendro.scalar("K"       , "")
Gt  = dendro.vec3("Gt"        , "")
b   = dendro.vec3("beta"      , "")
B   = dendro.vec3("B"         , "")
gt  = dendro.sym_3x3("gt"     , "")
At  = dendro.sym_3x3("At"     , "")
Gt_rhs  = dendro.vec3("Gt_rhs", "")

# Lie derivative weight
weight = -Rational(2,3)
weight_Gt = Rational(2,3)

# specify the functions for computing first and second derivatives

d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('grad')  # first argument is direction
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

#outs = [a_rhs, b_rhs,gt_rhs,chi_rhs,K_rhs, Gt_rhs, B_rhs]                       #gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
#vnames = ['a_rhs', 'b_rhs','gt_rhs','chi_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs' ]     #, 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']


exp_graph=nxgraph.ExpressionGraph()
exp_graph.add_expressions(outs,vnames)
G=exp_graph.composed_graph(verbose=False)

#print(colors)
custom_functions = {'grad': 'grad', 'grad2': 'grad2', 'agrad': 'agrad', 'kograd': 'kograd'}
rename_dict={ 'alpha'   : 'alpha[pp]',
              'beta0'   : 'beta0[pp]',
              'beta1'   : 'beta1[pp]',
              'beta2'   : 'beta2[pp]',
              'gt0'     : 'gt0[pp]',
              'gt1'     : 'gt1[pp]',
              'gt2'     : 'gt2[pp]',
              'gt3'     : 'gt3[pp]',
              'gt4'     : 'gt4[pp]',
              'gt5'     : 'gt5[pp]',
              'At0'     : 'At0[pp]',
              'At1'     : 'At1[pp]',
              'At2'     : 'At2[pp]',
              'At3'     : 'At3[pp]',
              'At4'     : 'At4[pp]',
              'At5'     : 'At5[pp]',
              'B0'      : 'B0[pp]',
              'B1'      : 'B1[pp]',
              'B2'      : 'B2[pp]',
              'Gt0'     : 'Gt0[pp]',
              'Gt1'     : 'Gt1[pp]',
              'Gt2'     : 'Gt2[pp]',
              'chi'     : 'chi[pp]',
              'K'       : 'K[pp]'}



def change_deriv_names(str):
    c_str=str
    derivs=['agrad','grad','kograd']
    for deriv in derivs:
        key=deriv+'\(\d, \w+\[pp\]\)'
        slist=regex.findall(key,c_str)
        for s in slist:
            #print(s)
            w1=s.split('(')
            w2=w1[1].split(')')[0].split(',')
            #print(w1[0]+'_'+w2[0].strip()+'_'+w2[1].strip()+';')
            rep=w1[0]
            for v in w2:
                rep=rep+'_'+v.strip()
            #rep=rep+';'
            c_str=c_str.replace(s,rep)

    derivs2=['grad2']
    for deriv in derivs2:
        key=deriv+'\(\d, \d, \w+\[pp\]\)'
        slist=regex.findall(key,c_str)
        for s in slist:
            #print(s)
            w1=s.split('(')
            w2=w1[1].split(')')[0].split(',')
            #print(w1[0]+'_'+w2[0].strip()+'_'+w2[1].strip()+';')
            rep=w1[0]
            for v in w2:
                rep=rep+'_'+v.strip()
            #rep=rep+';'
            c_str=c_str.replace(s,rep)
    return c_str

gutils.generate_cpu_c_code(exp_graph,8,"bssn_eqs.cpp",rename_dict,change_deriv_names,custom_functions,None)


