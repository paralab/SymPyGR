import os
import subprocess
import numpy as np
import csv

def getRandomNumbersOnGausianDistribution(mean, std, count, lowerBound, upperBound):
    rand=np.random.normal(mean, std, count*2)
    rand=np.round(rand,0)
    randBound=[]

    for val in rand:
        if len(randBound)==count:
            break
        else:
            if val<=upperBound and val>=lowerBound:
                randBound.append(val)
    else:
        return -1

    return randBound

def get_two_csv_files(level, number_of_blocks):
    os.chdir("../build")
    os.system("cmake ..")
    os.system("make -j4 all")

    out1 = subprocess.check_output(flop_rate_command + " %d %d %d"%(level, level, number_of_blocks), shell=True)
    print("level=%d | blocks=%d flop_rate.csv generated..."%(level, number_of_blocks))

    out2 = subprocess.check_output(timing_command + " %d %d %d"%(level, level, number_of_blocks), shell=True)
    print("level=%d | blocks=%d timing.csv generated..."%(level, number_of_blocks))

def process_csvs():
    messurements = []
    units = []

    HtoD_time = 0
    HtoD_size = 0
    HtoD_bandwidth = 0

    DtoH_time = 0
    DtoH_size = 0
    DtoH_bandwidth = 0

    Overall_bandwidth = 0

    kernel_time = 0

    flop_dp_count = 0
    flop_dp_eff = 0

    flop_sp_count = 0
    flop_sp_eff = 0

    total_entries = 0

    flop_dp_rate = 0
    flop_sp_rate = 0

    flop_dp_avg_eff = 0
    flop_sp_avg_eff = 0

    with open("../build/timing.csv", 'r') as timing:
        timingCSV = csv.reader(timing, delimiter=',')
        for index, row in enumerate(timingCSV):
            # extract mesurement
            if index==3:
                messurements = (row[1], row[11], row[12], row[16])
            
            # extract units
            if index==4:
                units = (row[1], row[11], row[12], row[16])
            
            if index<5: continue

            # extract only important data
            duration = row[1].strip()
            size = row[11].strip()
            throughput = row[12].strip()
            name = row[16][:20].strip()

            # extract HtoD DtoH sizes and times
            if name=="[CUDA memcpy HtoD]":
                HtoD_time += float(duration)
                HtoD_size += float(size)
            
            if name=="[CUDA memcpy DtoH]":
                DtoH_time += float(duration)
                DtoH_size += float(size)
            
            if not(name in ["[CUDA memcpy HtoD]", "[CUDA memcpy DtoH]"]):
                kernel_time += float(duration)

    #bandwidth calc
    HtoD_bandwidth = (HtoD_size/1024)/(HtoD_time/1000)
    DtoH_bandwidth = (DtoH_size/1024)/(DtoH_time/1000)
    Overall_bandwidth = ((HtoD_size+DtoH_size)/1024)/((HtoD_time+DtoH_time)/1000)

    with open("../build/flop_rate.csv", 'r') as fr:
        flopCSV = csv.reader(fr, delimiter=',')
        for index, row in enumerate(flopCSV):
            if index<6: continue
            flop_dp = row[4].strip()
            flop_sp = row[5].strip()
            flop_dp_e = row[6].strip()
            flop_sp_e = row[7].strip()

            flop_dp_count += int(flop_dp)
            flop_sp_count += int(flop_sp)

            flop_dp_eff += float(flop_dp_e)
            flop_sp_eff += float(flop_sp_e)

            total_entries += 1
    
    #flop_rate calc
    flop_dp_rate = (flop_dp_count/1000000000)/(kernel_time/1000)
    flop_sp_rate = (flop_sp_count/1000000000)/(kernel_time/1000)

    flop_dp_avg_eff = flop_dp_eff/total_entries
    flop_sp_avg_eff = flop_sp_eff/total_entries

    print("DP-Rate:%.2fGigaFlops/s | HtoD:%.2fGB/s | DtoH:%.2fGB/s | Overall:%.2fGB/s"%(flop_dp_rate, HtoD_bandwidth, DtoH_bandwidth, Overall_bandwidth))

    return flop_dp_rate, HtoD_bandwidth, DtoH_bandwidth, Overall_bandwidth

#constants
cuda_path = "/usr/local/cuda-8.0/bin/nvprof"
flop_rate_command = "%s --unified-memory-profiling off --print-gpu-trace --metrics flop_count_dp,flop_count_sp,flop_dp_efficiency,flop_sp_efficiency --csv --log-file flop_rate.csv ./computeBSSN"%(cuda_path)
timing_command = "%s --unified-memory-profiling off --print-gpu-trace --csv --log-file timing.csv ./computeBSSN"%(cuda_path)

#distribution parameter
number_of_executions = 4
mean = 3
std = 0.5
lower_bound = 0
upper_bound = 4

number_of_blocks = 2

def main():
    flop_dp_rates = []
    HtoD_bandwidths = []
    DtoH_bandwidths = []
    Overall_bandwidths = []

    for ind, level in enumerate(getRandomNumbersOnGausianDistribution(mean, std, number_of_executions, lower_bound, upper_bound)):
    # for ind, level in enumerate([1, 1]):
        print("Execution no: %d"%(ind+1))
        get_two_csv_files(level, number_of_blocks)
        flop_dp_rate, HtoD_bandwidth, DtoH_bandwidth, Overall_bandwidth = process_csvs()
        print()

        flop_dp_rates.append(flop_dp_rate)
        HtoD_bandwidths.append(HtoD_bandwidth)
        DtoH_bandwidths.append(DtoH_bandwidth)
        Overall_bandwidths.append(Overall_bandwidth)
    
    print("DP-Rate:%.2fGigaFlops/s | HtoD:%.2fGB/s | DtoH:%.2fGB/s | Overall:%.2fGB/s"%(np.mean(flop_dp_rates), np.mean(HtoD_bandwidth), np.mean(DtoH_bandwidth), np.mean(Overall_bandwidth)))


main()













# os.chdir("../build")
# os.system("cmake ..")
# os.system("make -j4 all")



# isCorrect = True
# for test in testcases:
#     out = subprocess.check_output("./computeBSSN " + test, shell=True)
#     # print(out.decode("utf-8"))

#     try:
#         cpu = open("../build/output_cpu.txt", "r")
#         gpu = open("../build/output_cuda.txt", "r")
#     except:
#         print("files not found")
#     else:
#         correct = True
#         count_line = 0
#         error_count = 0
#         roundingTo = 9
#         while(True):
#             count_line += 1
#             cpuline = cpu.readline().strip()
#             gpuline = gpu.readline().strip()
#             if cpuline!=gpuline:
#                 if round(float(cpuline), roundingTo)==round(float(gpuline), roundingTo): continue # special requirement
#                 error_count += 1
#                 print("line-%d \t| CPU:%s \t| GPU:%s \t| DIFF:%.13f"%(count_line, cpuline, gpuline, round(float(cpuline), roundingTo)-round(float(gpuline), roundingTo)))
#                 correct = False
#                 isCorrect = False
#                 # break # if you want to print only the first error, uncomment this line
                
#             if cpuline=="":
#                 if correct:
#                     print("Outputs matched for %s"%(test))
#                 else:
#                     print("Outputs did not matched for %s"%(test))
#                 break

#         cpu.close()
#         gpu.close()

# if isCorrect:
#     print("Output is correct")
# else:
#     print("Output is not correct")
#     print("Error count", error_count)
