#file to generate cuda rhs code
import  dendro as dendro
from sympy import *


###################################################################
# initialize cuda variables
###################################################################

l1, l2, l3, l4, eta = symbols('lambda[0] lambda[1] lambda[2] lambda[3] eta')
lf0, lf1 = symbols('lambda_f[0] lambda_f[1]')

# declare variables
a   = dendro.scalar("dev_var_in","[alphaInt+pp]")
chi = dendro.scalar("dev_var_in","[chiInt+pp]")
K   = dendro.scalar("dev_var_in","[KInt+pp]")

Gt  = dendro.vec3_cuda("dev_var_in","Gt", "[Int+pp]")
b   = dendro.vec3_cuda("dev_var_in","beta", "[Int+pp]")
B   = dendro.vec3_cuda("dev_var_in","B", "[Int+pp]")

gt  = dendro.sym_3x3_cuda("dev_var_in","gt", "[Int+pp]")
At  = dendro.sym_3x3_cuda("dev_var_in","At", "[Int+pp]")

Gt_rhs  = dendro.vec3_cuda("dev_var_out","gt", "[Int+pp]")

# Lie derivative weight
weight = -Rational(2,3)
weight_Gt = Rational(2,3)

# specify the functions for computing first and second derivatives
d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction
kod = dendro.set_kreiss_oliger_dissipation('kograd')

d2 = dendro.d2

#f = Function('f')

#   metric related quantities
dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(Gt, chi)
###################################################################
# evolution equations
###################################################################

a_rhs = l1*dendro.lie(b, a) - 2*a*K + 0*dendro.kodiss(a)


b_rhs = [ S(3)/4 * (lf0 + lf1*a) * B[i] +
        l2 * dendro.vec_j_ad_j(b, b[i])
         for i in dendro.e_i ] + dendro.kodiss(b)

#fixing the issue of generating invalid b_rhs expressions         
del b_rhs[-3:]

gt_rhs = dendro.lie(b, gt, weight) - 2*a*At + 0*dendro.kodiss(gt)

chi_rhs = dendro.lie(b, chi, weight) + Rational(2,3) * (chi*a*K) + 0*dendro.kodiss(chi)

AikAkj = Matrix([sum([At[i, k] * sum([dendro.inv_metric[k, l]*At[l, j] for l in dendro.e_i]) for k in dendro.e_i]) for i, j in dendro.e_ij])

At_rhs = dendro.lie(b, At, weight) + chi*dendro.trace_free( a*R - dendro.DiDj(a)) + a*(K*At - 2*AikAkj.reshape(3, 3)) + 0*dendro.kodiss(At)

K_rhs = dendro.lie(b, K) - dendro.laplacian(a,chi) + a*(K*K/3 + dendro.sqr(At)) + 0*dendro.kodiss(K)

At_UU = dendro.up_up(At)

Gt_rhs = Matrix([sum(b[j]*ad(j,Gt[i]) for j in dendro.e_i) for i in dendro.e_i]) - \
         Matrix([sum(CalGt[j]*d(j,b[i]) for j in dendro.e_i) for i in dendro.e_i]) + \
         Rational(2,3)*Matrix([ CalGt[i] * sum(d(j,b[j]) for j in dendro.e_i)  for i in dendro.e_i ]) + \
         Matrix([sum([igt[j, k] * d2(j, k, b[i]) + igt[i, j] * d2(j, k, b[k])/3 for j, k in dendro.e_ij]) for i in dendro.e_i]) - \
         Matrix([sum([2*At_UU[i, j]*d(j, a) for j in dendro.e_i]) for i in dendro.e_i]) + \
         Matrix([sum([2*a*dendro.C2[i, j, k]*At_UU[j, k] for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
         Matrix([sum([a*(3/chi*At_UU[i,j]*d(j, chi) + Rational(4,3)*dendro.inv_metric[i, j]*d(j, K)) for j in dendro.e_i]) for i in dendro.e_i])

Gt_rhs = [item for sublist in Gt_rhs.tolist() for item in sublist]

B_rhs = [Gt_rhs[i] - eta * B[i] +
         l3 * dendro.vec_j_ad_j(b, B[i]) -
         l4 * dendro.vec_j_ad_j(b, Gt[i]) + 0*kod(i,B[i])
         for i in dendro.e_i]



###################################################################
# generate code
###################################################################
outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
arrayName ='dev_var_out'
arrIdxs = ['alpha', 'beta', 'gt', 'chi', 'At', 'K', 'Gt', 'B']
dendro.generate_cuda_code(outs,arrayName, arrIdxs, '[Int+pp]')

