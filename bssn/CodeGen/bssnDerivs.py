"""
Contains derivative computation for BSSN formulation of ET equations. 
"""

# first derivative
import cog

D = ["alpha", "beta0", "beta1", "beta2",
      "B0", "B1", "B2",
      "chi", "Gt0", "Gt1", "Gt2", "K",
      "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
      "At0", "At1", "At2", "At3", "At4", "At5" ]


# custom functions for code generation in cse.
custom_functions = {'grad': 'grad', 'grad2': 'grad2', 'agrad': 'agrad', 'kograd': 'kograd'}

# second derivs required for RHS
DD = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5", "chi",
       "alpha", "beta0", "beta1", "beta2" ]

# advective derivatives
AD = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
       "At0", "At1", "At2", "At3", "At4", "At5",
       "alpha", "beta0", "beta1", "beta2", "chi", "Gt0", "Gt1", "Gt2", "K",
       "B0", "B1", "B2"] 

KO=AD

# first derivs required for constraints--no gauge variables
CONSTRAINT_D = [ "chi", "Gt0", "Gt1", "Gt2", "K",
           "gt0", "gt1", "gt2", "gt3", "gt4", "gt5",
           "At0", "At1", "At2", "At3", "At4", "At5" ]

# second derivs required for constraints--no gauge variables
CONSTRAINT_DD = ["gt0", "gt1", "gt2", "gt3", "gt4", "gt5", "chi"]


PREFIX_D   = ["grad_0_", "grad_1_", "grad_2_"]
PREFIX_AD  = ["agrad_0_", "agrad_1_", "agrad_2_"]
PREFIX_KOD = ["kograd_0_", "kograd_1_", "kograd_2_"]
PREFIX_DD  = ["grad2_0_0_", "grad2_0_1_", "grad2_0_2_", "grad2_1_1_", "grad2_1_2_", "grad2_2_2_"]

# first derivative in i direction
FUNC_D_I=[]
for f in D:
    for p in PREFIX_D:
        FUNC_D_I.append(p+f)

# second derivative in ij direction
FUNC_D_IJ=[]
for f in DD:
    for p in PREFIX_DD:
        FUNC_D_IJ.append(p+f)

#advective derivative in i direction
FUNC_AD_I=[]
for f in AD:
    for p in PREFIX_AD:
        FUNC_AD_I.append(p+f)


#Kriess-Oliger derivative in i direction
FUNC_KOD_I=[]
for f in D:
    for p in PREFIX_KOD:
        FUNC_KOD_I.append(p+f)

FUNC_CONS=[]
for f in CONSTRAINT_D:
    for p in PREFIX_D:
        FUNC_CONS.append(p+f)
        
for f in CONSTRAINT_DD:
    for p in PREFIX_DD:
        FUNC_CONS.append(p+f)


def allocDerivMemory():
    
    for deriv in FUNC_D_I:
        cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")

    for deriv in FUNC_D_IJ:
        cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")

    for deriv in FUNC_AD_I:
        cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")
        
        
def computeRHSDerivs():
    
    for var in D:
        cog.outl("\t deriv_x(%s, %s, hx, sz, bflag);" %(PREFIX_D[0] + var ,var))
        cog.outl("\t deriv_y(%s, %s, hx, sz, bflag);" %(PREFIX_D[1] + var ,var))
        cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(PREFIX_D[2] + var ,var))

        if var in DD:
            cog.outl("\t deriv_xx(%s, %s, hx, sz, bflag);" %(PREFIX_DD[0] + var ,var))
            cog.outl("\t deriv_y(%s, %s, hx, sz, bflag);" %(PREFIX_DD[1] + var , PREFIX_D[0] + var ))
            cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(PREFIX_DD[2] + var , PREFIX_D[0] + var ))

            cog.outl("\t deriv_yy(%s, %s, hx, sz, bflag);" %(PREFIX_DD[3] + var ,var))
            cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(PREFIX_DD[4] + var , PREFIX_D[1] + var))

            cog.outl("\t deriv_zz(%s, %s, hx, sz, bflag);" %(PREFIX_DD[5] + var ,var))

        if var in AD:
            cog.outl("\t adv_deriv_x(%s, %s, hx, sz, bflag);" %(PREFIX_AD[0] + var ,var))
            cog.outl("\t adv_deriv_y(%s, %s, hx, sz, bflag);" %(PREFIX_AD[1] + var ,var))
            cog.outl("\t adv_deriv_z(%s, %s, hx, sz, bflag);" %(PREFIX_AD[2] + var ,var))



def deallocDerivMemory():
    
    for deriv in FUNC_D_I:
        cog.outl("\t free(%s);" %(deriv))

    for deriv in FUNC_D_IJ:
        cog.outl("\t free(%s);" %(deriv))

    for deriv in FUNC_AD_I:
        cog.outl("\t free(%s);" %(deriv))


