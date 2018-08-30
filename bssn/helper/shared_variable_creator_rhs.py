import re

fileName = "../include/cuda_bssneqs.h"
f1 = open(fileName, "r")
lines = f1.readlines()
variablesCreated = []
sharedVariables = []
correctedLines = []
valueAssignment = []
count = 0
blockSize = 250
for word in lines:
    if "dev_var_in" in word:
        #get the term
        terms = re.findall(r'dev_var_in\[.*?\]', word)
        newWord = word
        for var in terms:
            if var not in variablesCreated:
                variableName = "var_" + str(count)
                dev_var_in_term = "__shared__ double "+ variableName + "[" + str(blockSize) + "]"
                sharedVariables.append(dev_var_in_term)
                variablesCreated.append(var)
                valueAssignment.append(variableName+'[threadIdx.x] = '+var)
                count += 1
            newWord = newWord.replace(var, sharedVariables[variablesCreated.index(var)].split()[2].split('[')[0]+'[threadIdx.x]')
        correctedLines.append(newWord)

    else:
        correctedLines.append(word)


newFile = open(fileName, 'w+')
for word in sharedVariables:
    newFile.write(word+';\n')

for word in valueAssignment:
    newFile.write(word+';\n')

for word in correctedLines:
    newFile.write(word)
