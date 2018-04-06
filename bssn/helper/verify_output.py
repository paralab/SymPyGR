import os
import subprocess
import numpy as np

# change malloc to calloc eg. "double *grad_1_alpha = (double *) calloc(n, sizeof(double));"
# Specify array in rhs_cuda line no. 127
# Specify array in rhs line no. 133

# change test_param.h accordingly.
# test 1
# isGPU 1
# isCPU 1

os.chdir("../build") # Change to build folder

os.system("cmake ..")
os.system("make all")

testcases = [
    "0 0 1", 
    "1 1 1", 
    "2 2 1", 
    "2 2 3", 
    "0 2 2",
    "0 2 1",
    ]

isCorrect = True
for test in testcases:
    subprocess.check_output("./computeBSSN " + test, shell=True)

    try:
        cpu = open("../build/output_cpu.txt", "r")
        gpu = open("../build/output_cuda.txt", "r")
    except:
        print("files not found")
    else:
        correct = True

        while(True):
            cpuline = cpu.readline()
            gpuline = gpu.readline()
            if cpuline!=gpuline:
                correct = False
                isCorrect = False
            if cpuline=="":
                if correct:
                    print("Outputs matched for %s"%(test))
                else:
                    print("Outputs did not matched for %s"%(test))
                break

        cpu.close()
        gpu.close()

if isCorrect:
    print("Output is correct")
else:
    print("Output is not correct")
