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
# Calculate the Weyl scalar, Psi4 for graviational wave extraction
###################################################################

m_np_real = dendro.vec3("m_np_real")
m_np_img = dendro.vec3("m_np_img")
r_np = dendro.vec3("r_np")
psi4_real = dendro.scalar("psi4_real")
psi4_1_real = dendro.scalar("psi4_1_real")
psi4_2_real = dendro.scalar("psi4_2_real")
psi4_3_real = dendro.scalar("psi4_3_real")
psi4_4_real = dendro.scalar("psi4_4_real")
psi4_img = dendro.scalar("psi4_img")
psi4_1_img = dendro.scalar("psi4_1_img")
psi4_2_img = dendro.scalar("psi4_2_img")
psi4_3_img = dendro.scalar("psi4_3_img")
psi4_4_img = dendro.scalar("psi4_4_img")

#TODO : Need to define more variables

# HL : Formula for Weyl terms need to check it. I follow the 
#      python syntax but didn't check fully

#for i in range(1,3):  
#    for j in range(i,3):
#   MM[i,j] = m_np_real[i]*m_np_real[j] - m_np_img[i]*m_np_img[j] 
#   NN[i,j] = m_np_real[i]*m_np_img[j] + m_np_real[j]*m_np_img[i]  

#for i in range(1,3):
#    for j in range(i+1,3):
#   MM[j,i] = MM[i,j] 
#   NN[j,i] = NN[i,j] 

#for i in range(1,2):
#    for j in range(i+1,3):
#   MR[i,j] = m_np_real[i]*r_np[j] - m_np_real[j]*r_np[i] 
#   NR[i,j] = m_np_img[i]*r_np[j] - m_np_img[j]*r_np[i] 

#for i in range(1,3): 
#   MR[i,i] = 0 
#   NR[i,i] = 0  

#   A_vec[i] = sum([Atd[j,i]*r_np[j] for j in range(1,3)]) 


#for i in range(1,2): 
#    for j in range(i+1,3): 
#   MR[j,i] = - MR[i,j] 
#   NR[j,i] = - NR[i,j] 

#for a in range(1,3): 
#    for b in range(1,3): 
#   Uu[a,b] = sum([m_np_real[c] * (d_Atd[b,c,a] + 0.5 * sum([Ctd[i,c,a] * Atud[i,b] for i in range(1,3)])) for c in range(1,3)]) 
#   Vv[a,b] = sum([m_np_imag[c] * (d_Atd[b,c,a] + 0.5 * sum([Ctd[i,c,a] * Atud[i,b] for i in range(1,3)])) for c in range(1,3)]) 

#r_d_chi = sum([r_np[i] * d_chi[i] for i in range(1,3)]) 

#A_temp = inv_chi * inv_chi * ( sum([A_vec[i] * r_np[i] for i range(1,3)]) + 3 * trK * chi + 0.5 * r_d_chi ) 

#m_real_d_chi = sum([m_np_real[i] * d_chi[i] for i in range(1,3)])  
#m_imag_d_chi = sum([m_np_imag[i] * d_chi[i] for i in range(1,3)]) 

#m_real_A_vec = sum([m_np_real[i] *  A_vec[i] for i in range(1,3)]) 
#m_imag_A_vec = sum([m_np_imag[i] *  A_vec[i] for i in range(1,3)])  


#psi4_1_real = sum([(inv_chi * Rpd[i,i] + Rtd[i,i]) * MM[i,i] for i in range(1,3)]) 
#            + 2*sum([sum([(inv_chi * Rpd[i,j] + Rtd[i,j]) * MM[i,j] for j in range(i+1,3)]) for i range(1,2)]) 
#psi4_1_img = sum([( inv_chi * Rpd[i,i] + Rtd[i,i]) * NN[i,i] for i range(1,3)]) 
#           + 2*sum([sum([(inv_chi * Rpd[i,j] + Rtd[i,j]) * NN[i,j] for j in range(i+1,3)]) for i in range(1,2)]) 

#psi4_2_real = A_temp * (sum([Atd[i,i] * MM[i,i] for i in range(1,3)]) 
#psi4_2_img = A_temp * (sum([Atd[i,i] * NN[i,i] for i in range(1,3)]) 
#           + 2*sum([sum([Atd[i,j] * NN[i,j] for j in range(i+1,3)]) for i in range(1,2)]))  


#psi4_3_real = inv_chi * sum([sum([MR[a,b]* Uu[a,b] - NR[a,b]*Vv[a,b] for a in range(1,3)]) for b in range(1,3)])  
#psi4_3_img = inv_chi * sum([sum([NR[a,b]* Uu[a,b] + MR[a,b]*Vv[a,b] for a in range(1,3)]) for b in range(1,3)]) 

#psi4_4_real = inv_chi * inv_chi * (m_real_A_vec * (m_real_A_vec + 0.5 * m_real_d_chi)  
#                                   - m_imag_A_vec * (m_imag_A_vec + 0.5 * m_imag_d_chi))  
#psi4_4_img = inv_chi * inv_chi * (m_real_A_vec * (m_imag_A_vec - 0.5 * m_imag_d_chi ) 
#                                  + m_imag_A_vec * (m_real_A_vec - 0.5 * m_real_d_chi))  

#psi4_real = psi4_1_real + psi4_2_real - psi4_3_real - psi4_4_real

#psi4_img = - (psi4_1_img + psi4_2_imag - psi4_3_imag - psi4_4_imag)

# Output for this should be included psi4_real and psi4_img as double precision  

###################################################################
# Include the constraints
###################################################################

# For vacuum EFE, rho and J are zero but we need this eventually
rho_ADM = dendro.scalar("rho_ADM") # Energy density term
Jtd_ADM = dendro.vec3("Jtd_ADM") # Momentum density term

pi = 3.141592653589793
G = 6.676408*10**-11 #Graviational constant HL : Do we need this?
HamCon = dendro.scalar("HamCon") # For Hamiltonian constraint
MomCon = dendro.vec3("MonCon")   # For Momentum constraint

# HL : Formula for constraint equations need to check it. I follow the 
#      python syntax but didn't check fully

#HamCon = -16*pi*G * rho_ADM + sum([sum([(Rpd[i,j] + chi * Rtd[i,j]) * gtu[i,j] for i in range(1,3)]) for j in range(1,3)]) - sum([sum([Atd[i,j] * Atu[i,j] for i range (1,3)]) for j in range(1,3)]) + 2/3 * K * K

#MomCon = sum([sum([gtu[j,k] * (d_Atd[k,j,i] - 1/2 * sum([Ct[m,k,i]*Atd[m,j] for m in range(1,3)]))for j in range(1,3)]) for k in range(1,3)]) - sum([Gamt[j]*Atd[i,j] for j in range(1,3)]) - 3/2*sum([Atud[j,i] * d_chi[j] for j in range(1,3)])*inv_chi - 2/3*d_K[i] - 8*pi*G * Jtd_ADM[i] * inv_chi

#Rscalar = sum([sum([(Rpd[i,j]+chi*Rtd[i,j]*gtu[i,j] for i in range(1,3)]) for j in range(1,3)]) #Ricci scalar

# Some BSSN constraints that we should consider

#HL : Need different name?
gamt_con = dendro.vec3("gamt_con")
calgamt_con = dendro.vec3("calgamt_con")
trA = dendro.scalar("trA")
detgtm1 = dendro.scalar("detgtm1")

#trA = sum([Atud[i,i] for i in range(1,3)]) 
#detgtm1 = detgtd - 1.0 : 

#for i from 1 to 3 do 
#gamt_con[i] = sum([gtd[i,j] * Gamt[j] for j in range(1,3)]) - sum([sum([gtu[j,k] * d_gtd[j,k,i] for j in range(1,3)]) for k in range(1,3)]) 
#calgamt_con[i] = sum([gtd[i,j] * CalGamt[j] for j in range(1,3)]) - sum([sum([gtu[j,k] * d_gtd[j,k,i] for j in range(1,3)]) for k in range(1,3)]) 

# Output for this should be included HamCon, MomCon, Rscalar, trA, detgtm1, gamt_con, and calgamt_con

###################################################################
# generate code
###################################################################

outs = [a_rhs, b_rhs, gt_rhs, chi_rhs, At_rhs, K_rhs, Gt_rhs, B_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']
#dendro.generate_debug(outs, vnames)
dendro.generate(outs, vnames)
