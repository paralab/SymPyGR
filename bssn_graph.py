
import sympy as sym
from bssn_eqns import dendroConfigs

idx = ""
iters = 1000

eqns, vnames = dendroConfigs.get_rhs_eqns_flat("evolution")


print(len(eqns), len(vnames))

print(vnames)