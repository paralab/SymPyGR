##########################################################################
# module: dendro
# author: Pasindu Tennage
# email:  pasindu.13@cse.mrt.ac.lk
#
# python module to generate dependency graph
#
# (c) 2018 University of Utah, All rights reserved.
##########################################################################
from sympy import *
import csv
import numpy as np
from sklearn.cluster import *

def getArguments(expression, result):
    #Returns the terms in the expression with duplicates removed
    args = expression.args
    func = expression.func
    if("grad" in str(func)):
        result.append(str(expression))
    elif("Symbol" in str(func)):
        result.append(str(expression))
    else:
        for arg in args:
            if(len(arg.args)==0):
                if((str(arg) not in result) and not(sympify(arg).is_real)):
                    result.append(str(arg))
            else:
                getArguments(arg, result)
    return result

def getAllOriginalVariables(dependencies):
    # Returns a list of all unique ground variables that are present
    originalVariables = []
    for i in dependencies.keys():
        values = dependencies[i]
        for j in values:
            originalVariables.append(j)
    originalVariables = list(set(originalVariables))
    return originalVariables

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]

def makeCompleteDependencies(_v):
    # Returns a dictionary of dependencies
    dependencies = {}
    counter = 1
    for i in _v[0]:
        result = []
        result = getArguments(i[1], result)
        for j in result:
            if(str(j) in dependencies.keys()):
                externDepen = dependencies[str(j)]
                for k in externDepen:
                    if(str(k) not in result):
                        result.append(str(k))
                result = remove_values_from_list(result, str(j))
        dependencies[str(i[0])] = list(set(result))
    for i in _v[1]:
        result = []
        result = getArguments(i, result)
        for j in result:
            if (str(j) in dependencies.keys()):
                externDepen = dependencies[str(j)]
                for k in externDepen:
                    if (str(k) not in result):
                        result.append(str(k))
                result = remove_values_from_list(result, str(j))
        dependencies["Equation" + str(counter)+": "] = list(set(result))
        counter= counter+1
    return dependencies

def getFeatureVectors(dependencies, originalVariables):
    #Return an an array that is used for clustering
    pointNames = dependencies.keys()
    featureVectors = []
    for i in pointNames:
        dependencyList = dependencies[i]
        featureVector = []
        for j in originalVariables:
            if(j in dependencyList):
                featureVector.append(1)
            else:
                featureVector.append(0)
        featureVectors.append(featureVector)
    return featureVectors

def writeFeatureVectorstoCSV(originalVariables, featureVectors):
    data = []
    data.append(originalVariables)
    for featureVector in featureVectors:
        data.append(featureVector)
    myFile = open('featureVectors.csv', 'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerows(data)


def kMeansClusterPoints(dependencies, featureVectors):
    #Using K means clustering, cluster the points
    X = np.array(featureVectors)
    kmeans = KMeans().fit(X)
    return kmeans

def miniBatchKMeansClusterPoints(dependencies, featureVectors):
    #Using K means clustering, cluster the points
    X = np.array(featureVectors)
    minibkmeans = MiniBatchKMeans().fit(X)
    return minibkmeans

def affinityPropagationClusterPoints(dependencies, featureVectors):
    #Using K means clustering, cluster the points
    X = np.array(featureVectors)
    affinityPropagation = AffinityPropagation().fit(X)
    return affinityPropagation

def meanShiftClusterPoints(dependencies, featureVectors):
    #Using K means clustering, cluster the points
    X = np.array(featureVectors)
    meanShift = MeanShift().fit(X)
    return meanShift

def spectralClusteringClusterPoints(dependencies, featureVectors):
    #Using K means clustering, cluster the points
    X = np.array(featureVectors)
    spectralClustering = SpectralClustering().fit(X)
    return spectralClustering

def agglomerativeClusteringClusterPoints(dependencies, featureVectors):
    #Using K means clustering, cluster the points
    X = np.array(featureVectors)
    agglomerative = AgglomerativeClustering().fit(X)
    return agglomerative

def dBSCANClusteringClusterPoints(dependencies, featureVectors):
    #Using K means clustering, cluster the points
    X = np.array(featureVectors)
    dBSCAN = DBSCAN().fit(X)
    return dBSCAN

def birchClusteringClusterPoints(dependencies, featureVectors):
    #Using K means clustering, cluster the points
    X = np.array(featureVectors)
    birch = Birch().fit(X)
    return birch

def writeClusteredResultstoCSV(name, method, dependencies):
    #Wrte cluster labels and dependencies to a csv file
    labels = method.labels_.tolist()
    data = []
    data.append(["Name", "Num_Dependecies", "Label"])
    k = 0
    for i in dependencies.keys():
        row = [i, len(dependencies[i]),labels[k]]
        k=k+1
        data.append(row)
    myFile = open(name+".csv", 'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerows(data)
    return


