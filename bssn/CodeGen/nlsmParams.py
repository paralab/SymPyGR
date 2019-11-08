from parameters import *
from decimal import Decimal
from sympy import *

parameters = Parameters()

parameters.add("__comment__(Evolution variable indices )", "U_CHI=0,U_PHI=1")
parameters.add("DENDRO_VERSION", 5.0)

parameters.setCategory("CONSTANTS")
parameters.add(PresetParam.NUM_VARS, 2)
parameters.add(PresetParam.ELE_ORDER, 4)

parameters.setCategory("IO")
parameters.add(PresetParam.RESTORE_SOLVER, 0)
parameters.add(PresetParam.IO_OUTPUT_FREQ, 400)
parameters.add(PresetParam.REMESH_TEST_FREQ, 10)
parameters.add(PresetParam.CHECKPT_FREQ, 10000)
parameters.add(PresetParam.IO_OUTPUT_GAP, 1)
parameters.add(PresetParam.VTU_FILE_PREFIX, "nlsm_gr")
parameters.add(PresetParam.CHKPT_FILE_PREFIX, "nlsm_cp")
parameters.add(PresetParam.PROFILE_FILE_PREFIX, "nlsm_prof")
parameters.add(PresetParam.NUM_EVOL_VARS_VTU_OUTPUT, 2)
parameters.add(PresetParam.VTU_OUTPUT_EVOL_INDICES, [0, 1], "evolution variable ids")

parameters.setCategory("LOAD BALANCING & MESH")
parameters.add(PresetParam.DENDRO_GRAIN_SZ, 100)
parameters.add(PresetParam.ASYNC_COMM_K, 2)
parameters.add(PresetParam.DENDRO_AMR_FAC, 1.0)
parameters.add(PresetParam.LOAD_IMB_TOL, .1)
parameters.add(PresetParam.DIM, 3)
parameters.add(PresetParam.MAXDEPTH, 6)
parameters.add(PresetParam.SPLIT_FIX, 2)

parameters.setCategory("WAVELET REFINEMENT")
parameters.add(PresetParam.WAVELET_TOL, 1e-4)
parameters.add(PresetParam.NUM_REFINE_VARS, 2)
parameters.add(PresetParam.REFINE_VARIABLE_INDICES, [0,1])

parameters.setCategory("GRID")
parameters.add(PresetParam.GRID_MIN_X, -200.0)
parameters.add(PresetParam.GRID_MAX_X, 200.0)
parameters.add(PresetParam.GRID_MIN_Y, -200.0)
parameters.add(PresetParam.GRID_MAX_Y, 200.0)
parameters.add(PresetParam.GRID_MIN_Z, -200.0)
parameters.add(PresetParam.GRID_MAX_Z, 200.0)

parameters.setCategory("RK SOLVER")
parameters.add(PresetParam.RK45_TIME_BEGIN, 0)
parameters.add(PresetParam.RK45_TIME_END, 10)

parameters.add(PresetParam.CFL_FACTOR, .1)

compd_min = [parameters[PresetParam.GRID_MIN_X].value, parameters[PresetParam.GRID_MIN_Y].value, parameters[PresetParam.GRID_MIN_Z].value]
parameters.add(PresetParam.COMPD_MIN, compd_min)

compd_max = [parameters[PresetParam.GRID_MAX_X].value, parameters[PresetParam.GRID_MAX_Y].value, parameters[PresetParam.GRID_MAX_Z].value]
parameters.add(PresetParam.COMPD_MAX, compd_max)

timestep = parameters[PresetParam.CFL_FACTOR].value * (parameters[PresetParam.COMPD_MAX].value[0] - parameters[PresetParam.COMPD_MIN].value[0]) / (1 << parameters[PresetParam.MAXDEPTH].value)
parameters.add(PresetParam.RK45_TIME_STEP_SIZE, timestep)
parameters.add(PresetParam.RK45_DESIRED_TOL, 1e-3)
parameters.add(PresetParam.KO_DISS_SIGMA, 1e-1)

parameters.setCategory("BLOCK Adaptivity (Not Recommended use AMR)")
parameters.add(PresetParam.ENABLE_BLOCK_ADAPTIVITY, 0)
parameters.add(PresetParam.BLK_MIN_X, -6.0)
parameters.add(PresetParam.BLK_MIN_Y, -6.0)
parameters.add(PresetParam.BLK_MIN_Z, -6.0)
parameters.add(PresetParam.BLK_MAX_X, 6.0)
parameters.add(PresetParam.BLK_MAX_Y, 6.0)
parameters.add(PresetParam.BLK_MAX_Z, 6.0)

namespace = "nlsm"

x = symbols("x")
y = symbols("y")
z = symbols("z")

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

radiusExpr = sqrt( epsx*(x-xc)**2 + epsy*(y-yc)**2 + epsz*(z-zc)**2 )
radius = parameters.addInitialData("radius", radiusExpr)

chiExpr = amp * exp(-(radius-R)**2/(delta**2))
parameters.addVar("chi", chiExpr)
parameters.addVar("phi", 0.0)