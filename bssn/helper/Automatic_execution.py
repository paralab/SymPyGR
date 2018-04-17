import os
import subprocess
import numpy as np

os.chdir("../build") # Change to build folder

os.system("clear")
os.system("cmake ..")
os.system("make all")

testcases = [
    "0 0 1"
    ]

for test in testcases:
    out = subprocess.check_output("./computeBSSN " + test, shell=True)
    print(out.decode("utf-8"))
    