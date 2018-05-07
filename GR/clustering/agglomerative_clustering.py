from sympy import *
import csv
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

class Cluster:
  index = 0
  left_child = None
  right_child = None
  ancestor_index = None
  item_list = None
  level = None

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

def makeCompleteDependencies(_v, lname):
    # Returns a dictionary of dependencies
    dependencies = {}
    counter = 0
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
        dependencies[lname[counter]] = list(set(result))
        counter= counter+1
    return dependencies

def getFeatureVectors(dependencies, originalVariables):
    #Return an an array that is used for clustering
    featureVectors = []
    for i in dependencies.keys():
        featureVector = []
        dependency_list = dependencies[i]
        for j in originalVariables:
            if(j in dependency_list):
                featureVector.append(1)
            else:
                featureVector.append(0)
        featureVectors.append(featureVector)
    return featureVectors

def writeFeatureVectorstoCSV(row_names,column_names,featureVectors):
    data = []
    data.append(["Variable Name"]+list(column_names))
    i = 0
    for featureVector in featureVectors:
        data.append([row_names[i]]+featureVector)
        i = i+1
    myFile = open('featureVectors.csv', 'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerows(data)

def cluster(feature_vectors):
    X = np.asarray(feature_vectors)
    Z = linkage(X, method='average', metric='cosine')
    # calculate full dendrogram
    plt.figure()
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    )
    #plt.show()
    return Z

def printClusterResults(z):
    print(type(z))
    np.set_printoptions(threshold=np.nan)
    print(z)

def createVariableClusterGraphList(z, dendro_and_equations, threshold):
    # z has the format of [[idx1, idx2, dist, sample_count], [...]]
    # All indices idx >= len(X) actually refer to the cluster formed in Z[idx - len(X)]
    print()
    z = z.tolist()
    np.set_printoptions(threshold=np.nan)
    newz= []
    substitutions = {}
    substitutions_counter = 0
    #dealing with the bottom layer variables
    for i in dendro_and_equations:
        s = Cluster()
        s.left_child=None
        s.right_child=None
        s.ancestor_index = None
        s.item_list = [i]
        s.index = dendro_and_equations.index(i)
        newz.append(s)
        s.level = 0

    index = len(dendro_and_equations)
    for i in z:
        s = Cluster()
        left_child = int(i[0])
        right_child = int(i[1])
        UpdateChildAncestor(left_child, newz, index)
        UpdateChildAncestor(right_child, newz, index)
        s.left_child = left_child
        s.right_child = right_child
        s.level = max(getClusterLevel(newz, left_child),getClusterLevel(newz, right_child))+1
        s.ancestor_index = None
        s.index = index
        left_dependency_list = getClusterItemList(newz, left_child)
        right_dependency_list = getClusterItemList(newz, right_child)

        i_item_list = list(set(left_dependency_list+right_dependency_list))

        if(len(i_item_list)<threshold):
            s.item_list = i_item_list
        elif(len(i_item_list)==threshold):
            substitutions["Reduction_"+str(substitutions_counter)] =i_item_list
            s.item_list = ["Reduction_"+str(substitutions_counter)]
            #print("Since item at index "+str(index)+"  has "+str(threshold)+
                  #" number of dependencies, it's item list is renamed as  "+"Reduction_"+str(substitutions_counter))
            substitutions_counter = substitutions_counter+1
        elif (len(i_item_list) > threshold):
            if(len(left_dependency_list)>1):
                substitutions["Reduction_" + str(substitutions_counter)] = left_dependency_list
                left_dependency_list = ["Reduction_" + str(substitutions_counter)]
                updateChildDependencyList(newz, left_child,"Reduction_" + str(substitutions_counter))
                #print("Since item at index " + str(index) + "  has greater than " + str(threshold) +
                      #" number of dependencies, it's left child's (which is at index "+
                      #str(left_child)+ ") item list is renamed as  " +
                      #"Reduction_" + str(substitutions_counter))
                substitutions_counter = substitutions_counter + 1

            if (len(right_dependency_list) > 1):
                substitutions["Reduction_" + str(substitutions_counter)] = right_dependency_list
                right_dependency_list = ["Reduction_" + str(substitutions_counter)]
                updateChildDependencyList(newz, right_child, "Reduction_" + str(substitutions_counter))
                #print("Since item at index " + str(index) + "  has greater than " + str(threshold) +
                      #" number of dependencies, it's right child's (which is at index " + str(right_child) +
                      #") item list is renamed as  " + "Reduction_" + str(substitutions_counter))
                substitutions_counter = substitutions_counter + 1

                s.item_list = left_dependency_list+right_dependency_list

        newz.append(s)
        index=index+1

    return (newz,substitutions)

def getMaxLevel(newz):
    max_level = 0
    for i in newz:
        if (i.level>max_level):
            max_level = i.level
    return max_level

def getClustersByLevel(newz, level):
    clusters = []
    for i in newz:
        if (i.level==level):
            clusters.append(i)
    return clusters

def isReducedCLuster(cluster):
    if(len(cluster.item_list)==1):
        return True
    else:
        return False


def writeClusterListToFile(newz):
    data = []
    data.append(["Index","Item List","Left_Child_Index","Right_Child_Index","Parent_Index","Level"])
    k = 0
    for i in newz:
        data.append([str(i.index),str(i.item_list),str(i.left_child),str(i.right_child),str(i.ancestor_index),str(i.level)])
        k = k + 1
    myFile = open('cluster_results.csv', 'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerows(data)

def getClusterItemList(newz, index):
    for i in newz:
        if(i.index==index):
            return i.item_list
    return None

def getClusterItemListUsingCluster(cluster):
    return cluster.item_list



def getClusterLevel(newz, index):
    for i in newz:
        if(i.index==index):
            return i.level
    return None

def updateChildDependencyList(newz, index, new_name):
    for i in newz:
        if(i.index==index):
            i.item_list = [new_name]
            break

def UpdateChildAncestor(index, newz, parent_index):
    for i in newz:
        if(i.index==index):
            i.ancestor_index = parent_index
            break

def createClusteringGraph(newz):
    G = nx.DiGraph()
    for cluster in newz:
        G.add_node(str(cluster.index)+": "+str(cluster.item_list))
        G.add_node(str(cluster.ancestor_index)+": "+str(getClusterItemList(newz, cluster.ancestor_index)))
        G.add_edge(str(cluster.index)+": "+str(cluster.item_list), str(cluster.ancestor_index)+": "+str(getClusterItemList(newz, cluster.ancestor_index)))

    pos = nx.spring_layout(G, scale=10)  # double distance between all nodes
    nx.draw(G, pos=pos, with_labels=True)
    plt.show()
