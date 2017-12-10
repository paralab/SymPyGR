import dendro
from sympy import *

###################################################################
# initialize
###################################################################

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')

# declare variables
a   = dendro.scalar("alpha")
chi = dendro.scalar("chi")
K   = dendro.scalar("K")

Gt  = dendro.vec3("Gt")
b   = dendro.vec3("beta")
B   = dendro.vec3("B")

gt  = dendro.sym_3x3("gt")
At  = dendro.sym_3x3("At")

Gt_rhs  = dendro.vec3("Gt_rhs")

# Lie derivative weight
weight = -2/3
weight_Gt = 2/3

# specify the functions for computing first and second derivatives
d = dendro.set_first_derivative('grad')    # first argument is direction
d2 = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction

#f = Function('f')

# generate metric related quantities 
dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails 
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)
###################################################################
# evolution equations
###################################################################

a_rhs = l1*dendro.lie(b, a) - 2*a*K

#[ewh] In the had code, this is treated as an advective derivative.  
#      I think this should be:    
#         l2 * dendro.vec_j_del_j(b, b[i]) 
b_rhs = [ (3/4) * (lf0 + lf1*a) * B[i] +
        l2 * dendro.vec_j_ad_j(b, b[i]) 
         for i in dendro.e_i ]

gt_rhs = dendro.lie(b, gt, weight) - 2*a*At

chi_rhs = dendro.lie(b, chi, weight) + (2/3) * (chi*a*K)

AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])

#ewh2 At_rhs = dendro.lie(b, At, weight) + dendro.trace_free(chi*(dendro.DiDj(a) + a*R)) + a*(K*At - 2*AikAkj.reshape(3, 3))
At_rhs = dendro.lie(b, At, weight) + chi*dendro.trace_free( a*R - dendro.DiDj(a)) + a*(K*At - 2*AikAkj.reshape(3, 3))

#K_rhs = dendro.vec_k_del_k(b, K) - dendro.laplacian(a) + a*(1/3*K*K + dendro.sqr(At))
K_rhs = dendro.lie(b, K) - dendro.laplacian(a,chi) + a*(K*K/3 + dendro.sqr(At))

At_UU = dendro.up_up(At)

Gt_rhs = Matrix([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i]) - \
         Matrix([sum(CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i]) + \
         2/3*Matrix([ CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ]) + \
         Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i]) - \
         Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i]) + \
         Matrix([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
         Matrix([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + (4/3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i])

Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

B_rhs = [Gt_rhs[i] - eta * B[i] + 
         l3 * dendro.vec_j_ad_j(b, B[i]) - 
         l4 * dendro.vec_j_ad_j(b, Gt[i]) 
         for i in dendro.e_i]


#_I = gt*igt
#print(simplify(_I))

#_I = gt*dendro.inv_metric
#print(simplify(_I))

###################################################################
# generate code
###################################################################

outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']
#dendro.generate_debug(outs, vnames)
dendro.generate(outs, vnames)
