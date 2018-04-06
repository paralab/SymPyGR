import os
import subprocess
import numpy as np

# GPU\CPU Performance test
# change test_param.h accordingly.
# test 1
# isGPU\isCPU 1

runtime = []
derivs = []
rhs = []
bdy = []

def print_report():
    runtimeNP = np.array(runtime)
    runtime_std = np.std(runtimeNP)
    runtime_avg = np.average(runtimeNP)
    runtime_max = np.max(runtimeNP)
    runtime_min = np.min(runtimeNP)

    derivsNP = np.array(derivs)
    derivs_std = np.std(derivsNP)
    derivs_avg = np.average(derivsNP)
    derivs_max = np.max(derivsNP)
    derivs_min = np.min(derivsNP)

    rhsNP = np.array(rhs)
    rhs_std = np.std(rhsNP)
    rhs_avg = np.average(rhsNP)
    rhs_max = np.max(rhsNP)
    rhs_min = np.min(rhsNP)

    bdyNP = np.array(bdy)
    bdy_std = np.std(bdyNP)
    bdy_avg = np.average(bdyNP)
    bdy_max = np.max(bdyNP)
    bdy_min = np.min(bdyNP)

    print("Runtime: \tavg=%f \tstd=%f \t max=%f \t min=%f"%(runtime_avg, runtime_std, runtime_max, runtime_min))
    print("DerivTime: \tavg=%f \tstd=%f \t max=%f \t min=%f"%(derivs_avg, derivs_std, derivs_max, derivs_min))
    print("RHSTime: \tavg=%f \tstd=%f \t max=%f \t min=%f"%(rhs_avg, rhs_std, rhs_max, rhs_min))
    print("BodyTime: \tavg=%f \tstd=%f \t max=%f \t min=%f"%(bdy_avg, bdy_std, bdy_max, bdy_min))
    print("\n")


os.chdir("../build") # Change to build folder

try:
    os.remove("performance.txt")
except:
    print("File was already deleted")

os.system("cmake ..")
os.system("make all") 


for i in range(100):
    output = str(subprocess.check_output("./computeBSSN 4 4 1", shell=True))
    output = output.split("\\n")

    run = float(output[-5][20:].strip())
    runtime.append(run)

    der = float(output[-4][20:].strip())
    derivs.append(der)

    rhstime = float(output[-3][20:].strip())
    rhs.append(rhstime)

    bdytime = float(output[-2][20:].strip())
    bdy.append(bdytime)

    if (i%20)==0:
        print("%dth checkpoint"%(i))
        print_report()

print("\nPrint final report")
print_report()