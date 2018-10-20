##########################################################################
# author: Milinda Fernando
# email:  milinda@cs.utah.edu
# date: 10/1/2018
#
# Manages the GPU shared memory
# Note that currently memory manager use sequential allocation not similar to 
# heap like allocations. (We can add those if it is needed)
#
##########################################################################

import sys as sys
import math as math

class MemoryManager:

    # shared memory size in bytes in the device
    __memMaxSz=48*1024

    # shared usable memory in bytes
    __memUsable=48*1024

    # allocated memory
    __memAlloc=0

    # memory manager
    __memMap=dict()

    # basic data type for the memory allocation
    __basicType={'char':1,'bool':1,'int':4,'float':4,'double':8}

    # default prefix value for base declaration
    __decPrefix="__shared__ "

    # variable type of the mem manager
    __varType=""

    # maximum number of elements/ base size. 
    __maxElem=0

    # current location for new variable allocations
    __currentLoc=0

    __baseName=""

    def __init__(self, maxMemSz,memUsable,cout,baseName="__sm_base",varType="double"):
        self.__memMaxSz = maxMemSz
        self.__memUsable = memUsable
        self.__memAlloc=0
        self.__memMap=dict()
        self.__baseName=baseName
        
        if varType in self.__basicType:
            maxElems=int(math.floor((self.__memUsable)/(self.__basicType[varType])))
            self.__maxElem=maxElems
            #print(self.__maxElem)
            self.__varType=varType
            print(self.__decPrefix+" "+varType+" "+baseName+"["+str(maxElems)+"];",file=cout)
            
        else:
            print("MemManager Error: Invalid data type\n")
            sys.exit(0)



        self.__currentLoc=0

    def malloc(self,varName,numElements,cout,prefix=""):
        if varName in self.__memMap:
            print("MemManager Error: Duplicate declaration of "+varName+" \n")
            sys.exit(0)
            #self.__memMap[varName]=[self.__currentLoc,numElements]
            #print(prefix+varName+" = "+self.__baseName+" + "+str(self.__currentLoc)+";",file=cout)
            #self.__currentLoc=self.__currentLoc+numElements
        else:
            if self.__currentLoc + numElements > self.__maxElem :
                print("MemManager Error: Max element count reached while trying to malloc "+varName)
                sys.exit(0)
            self.__memMap.update({varName:[self.__currentLoc,numElements]})
            print(prefix+self.__varType+" * "+varName+" = "+self.__baseName+" + "+str(self.__currentLoc)+";",file=cout)
            self.__currentLoc=self.__currentLoc+numElements
    
    def deallocAll(self):
        self.__currentLoc=0
        #self.__memMap.clear()

    def clearScopeVariables(self):
        self.__memMap=dict()

    def getMemUsable(self):
        return self.__memUsable


'''
    def mfree(self,varName):
        if varName  in self.__memMap:
            print("MemManager Error: currently not available \n")
            sys.exit(0)
        else:
            print("MemManager Error: variable "+varName+" is not found in mem. map \n")
            sys.exit(0)
'''

        


