try:
    cpu = open("../build/CPU_CHECK.txt", "r")
    gpu = open("../build/GPU_CHECK.txt", "r")
except:
    print("files not found")
else:
    correct = True
    count_line = 0
    while(True):
        count_line += 1
        cpul = cpu.readline().strip().split("--")
        gpul = gpu.readline().strip().split("--")

        if len(cpul)<2:
            if correct: 
                print("Values matched")
            else:
                print("Values are not matching")
            break


        cpuline = cpul[0].strip()
        gpuline = gpul[0].strip()
        dendroVar = cpul[1].strip()

        if cpuline!=gpuline:
            if round(float(cpuline), 20)==round(float(gpuline), 20): continue # special requirement
            print("line-%d \t| CPU-%s \t| GPU-%s | %s | dif-%.20f"%(count_line, cpuline, gpuline, dendroVar, abs(float(cpuline)-float(gpuline))))
            correct = False
            # break # if you want to print only the first error, uncomment this line

    cpu.close()
    gpu.close()