import dendro
from sympy import *

###################################################################
# initialize
###################################################################
chi = dendro.scalar("chi", "[pp]")
K   = dendro.scalar("K", "[pp]")
Gt  = dendro.vec3("Gt", "[pp]")
gt  = dendro.sym_3x3("gt", "[pp]")
At  = dendro.sym_3x3("At", "[pp]")

#ham = dendro.scalar("Ham", "[pp]")
#mom = dendro.vec3("Mom", "[pp]")

d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
ad = dendro.set_advective_derivative('agrad')  # first argument is direction
kod = dendro.set_kreiss_oliger_dissipation('kograd')

d2 = dendro.d2

###################################################################
# Set metric
###################################################################
dendro.set_metric(gt)
igt = dendro.get_inverse_metric()

C1 = dendro.get_first_christoffel()
C2 = dendro.get_second_christoffel()
#what's this...tried to comment it out and python compilation fails
C2_spatial = dendro.get_complete_christoffel(chi)
R, Rt, Rphi, CalGt = dendro.compute_ricci(gt, chi)


###################################################################
# Constraint Equations
###################################################################

ham = sum(igt[j,k]*R[j,k] for j,k in dendro.e_ij) - dendro.sqr(At) + Rational(2/3)*K**2

mom = Matrix([sum([igt[j,k]*(  d(k,At[i,j]) - \
              sum(dendro.C2[m,k,i]*At[j,m] for m in dendro.e_i)) \
                  for j,k in dendro.e_ij]) for i in dendro.e_i]) - \
      Matrix([sum([Gt[j]*At[i,j] for j in dendro.e_i]) for i in dendro.e_i]) -\
      Rational(3,2)*Matrix([ \
            sum([igt[j,k]*At[k,i]*d(j,chi)/chi for j,k in dendro.e_ij])  \
            for i in dendro.e_i]) -\
      Rational(2,3)*Matrix([d(i,K) for i in dendro.e_i])
mom = [item for sublist in mom.tolist() for item in sublist]

###################################################################
# generate code
###################################################################

outs = [ham, mom] 
vnames = ['ham', 'mom']
dendro.generate(outs, vnames, '[pp]')


