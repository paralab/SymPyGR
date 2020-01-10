from parameters import *
from decimal import Decimal
from sympy import *
from sympy.codegen.ast import String

parameters = Parameters()

parameters.add("__comment__(Evolution variable indices )", "U_CHI=0,U_PHI=1")
parameters.add("DENDRO_VERSION", 5.0)

parameters.setCategory("CONSTANTS")
parameters.add(RequiredParameter.ELE_ORDER, 4)

parameters.setCategory("IO")
parameters.add(RequiredParameter.RESTORE_SOLVER, 0)
parameters.add(RequiredParameter.IO_OUTPUT_FREQ, 10)
parameters.add(RequiredParameter.TIME_STEP_OUTPUT_FREQ,10)
parameters.add(RequiredParameter.REMESH_TEST_FREQ, 10)
parameters.add(RequiredParameter.CHECKPT_FREQ, 10000)
parameters.add(RequiredParameter.IO_OUTPUT_GAP, 1)

#for scratch output when running on the cluster
#parameters.add(RequiredParameter.VTU_FILE_PREFIX, "/scratch/kingspeak/serial/{uid}/nlsm_gr")
#parameters.add(RequiredParameter.CHKPT_FILE_PREFIX, "/scratch/kingspeak/serial/{uid}/nlsm_cp")
#parameters.add(RequiredParameter.PROFILE_FILE_PREFIX, "/scratch/kingspeak/serial/{uid}/nlsm_prof")
parameters.add(RequiredParameter.VTU_FILE_PREFIX, "nlsm_gr")
parameters.add(RequiredParameter.CHKPT_FILE_PREFIX, "nlsm_cp")
parameters.add(RequiredParameter.PROFILE_FILE_PREFIX, "nlsm_prof")
parameters.add(RequiredParameter.NUM_EVOL_VARS_VTU_OUTPUT, 2)
parameters.add(RequiredParameter.VTU_OUTPUT_EVOL_INDICES, [0, 1])

parameters.setCategory("LOAD BALANCING & MESH")
parameters.add(RequiredParameter.DENDRO_GRAIN_SZ, 100)
parameters.add(RequiredParameter.ASYNC_COMM_K, 2)
parameters.add(RequiredParameter.DENDRO_AMR_FAC, 1.0)
parameters.add(RequiredParameter.LOAD_IMB_TOL, .1)
parameters.add(RequiredParameter.DIM, 3)
parameters.add(RequiredParameter.MAXDEPTH, 10)
parameters.add(RequiredParameter.SPLIT_FIX, 2)

parameters.setCategory("WAVELET REFINEMENT")
parameters.add(RequiredParameter.NUM_REFINE_VARS, 2)
parameters.add(RequiredParameter.REFINE_VARIABLE_INDICES, [0,1])

parameters.setCategory("GRID")
parameters.add(RequiredParameter.GRID_MIN_X, -200.0)
parameters.add(RequiredParameter.GRID_MAX_X, 200.0)
parameters.add(RequiredParameter.GRID_MIN_Y, -200.0)
parameters.add(RequiredParameter.GRID_MAX_Y, 200.0)
parameters.add(RequiredParameter.GRID_MIN_Z, -200.0)
parameters.add(RequiredParameter.GRID_MAX_Z, 200.0)

parameters.setCategory("RK SOLVER")
parameters.add(RequiredParameter.RK45_TIME_BEGIN, 0)
parameters.add(RequiredParameter.RK45_TIME_END, 100)

parameters.add(RequiredParameter.CFL_FACTOR, .1)

compd_min = [parameters[RequiredParameter.GRID_MIN_X].value, parameters[RequiredParameter.GRID_MIN_Y].value, parameters[RequiredParameter.GRID_MIN_Z].value]
parameters.add(RequiredParameter.COMPD_MIN, compd_min)

compd_max = [parameters[RequiredParameter.GRID_MAX_X].value, parameters[RequiredParameter.GRID_MAX_Y].value, parameters[RequiredParameter.GRID_MAX_Z].value]
parameters.add(RequiredParameter.COMPD_MAX, compd_max)

octree_min = [0.0,0.0,0.0]
parameters.add(RequiredParameter.OCTREE_MIN, octree_min)

maxdepth = 1<<parameters[RequiredParameter.MAXDEPTH].value
octree_max = [maxdepth,maxdepth,maxdepth]
parameters.add(RequiredParameter.OCTREE_MAX, octree_max)

timestep = parameters[RequiredParameter.CFL_FACTOR].value * (parameters[RequiredParameter.COMPD_MAX].value[0] - parameters[RequiredParameter.COMPD_MIN].value[0]) / (1 << parameters[RequiredParameter.MAXDEPTH].value)
parameters.add(RequiredParameter.RK45_TIME_STEP_SIZE, timestep)
parameters.add(RequiredParameter.RK45_DESIRED_TOL, 1e-3)
parameters.add(RequiredParameter.KO_DISS_SIGMA, 1e-1)
parameters.add(RequiredParameter.RK_MIN_TOL, 1e-5)

parameters.setCategory("BLOCK Adaptivity (Not Recommended use AMR)")
parameters.add(RequiredParameter.ENABLE_BLOCK_ADAPTIVITY, 0)
parameters.add(RequiredParameter.BLK_MIN_X, -6.0)
parameters.add(RequiredParameter.BLK_MIN_Y, -6.0)
parameters.add(RequiredParameter.BLK_MIN_Z, -6.0)
parameters.add(RequiredParameter.BLK_MAX_X, 6.0)
parameters.add(RequiredParameter.BLK_MAX_Y, 6.0)
parameters.add(RequiredParameter.BLK_MAX_Z, 6.0)

namespace = "nlsm"

x = Symbol("x")
y = Symbol("y")
z = Symbol("z")

parameters.clearCategory()
amp = parameters.addInitialData("AMP", 1.3)
R = parameters.addInitialData("R", 8.0)
delta = parameters.addInitialData("DELTA", 3.0)
xc = parameters.addInitialData("XC", 0.0)
yc = parameters.addInitialData("YC", 0.0)
zc = parameters.addInitialData("ZC", 0.0)
epsx = parameters.addInitialData("EPSX", 1.0)
epsy = parameters.addInitialData("EPSY", 1.0)
epsz = parameters.addInitialData("EPSZ", 1.0)
nu = parameters.addInitialData("NU", 0.0)

waveSpeedX = parameters.addInitialData("WAVE_SPEED_X", 1.0)
waveSpeedY = parameters.addInitialData("WAVE_SPEED_Y", 0.0)
waveSpeedZ = parameters.addInitialData("WAVE_SPEED_Z", 0.0)

radiusExpr = sqrt( epsx*(x-xc)**2 + epsy*(y-yc)**2 + epsz*(z-zc)**2 )
radius = parameters.addInitialData("radius", radiusExpr)

chiInitialData = amp * exp(-(radius-R)**2/(delta**2))
chiRHS = Symbol("phi")
chi = parameters.addVar("chi", chiInitialData, chiRHS)

grad2_0_0_chi = parameters.addRhsPrecompute(PrecomputeFunc.deriv_xx, chi)
grad2_1_1_chi = parameters.addRhsPrecompute(PrecomputeFunc.deriv_yy, chi)
grad2_2_2_chi = parameters.addRhsPrecompute(PrecomputeFunc.deriv_zz, chi)

phiRHS = waveSpeedX*grad2_0_0_chi + waveSpeedY*grad2_1_1_chi + waveSpeedZ*grad2_2_2_chi;
parameters.addVar("phi", 0.0, phiRHS)