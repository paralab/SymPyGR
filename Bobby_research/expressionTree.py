#This will be the graph class
import networkx as nx
import pydot
import graphviz
import disjoint_set as ds
from collections import deque

class expressionTree:
    def __init__(self):
        self.idDisplayMap = {}
        self.leavesNodeMap = {}
        self.graph = nx.DiGraph()
        self.sources = set()
        self.mult_sources = set()

    def jaccard_similarity(self, nodeId, otherExpressionTree, otherId):
        ''' similiarity based on the leaf dependecies using the jaccard similarity'''

        count = self.jaccard_similiarity_count(nodeId, otherExpressionTree, otherId)
        return count / (self.getNumLeafDependents(nodeId) + otherExpressionTree.getNumLeafDependents(otherId) - count)

    def jaccard_similiarity_count(self, nodeId, otherExpressionTree, otherId):
        ''' returns the number of leaf dependecies that are the same between two trees'''

        leaf_values = self.getLeafDependents(nodeId)
        other_leaf_values = otherExpressionTree.getLeafDependents(otherId)

        count = 0
        for leaf_value in leaf_values:
            #only check for similarity of id is in both trees
            if leaf_value in other_leaf_values:
                    count = count + 1
        return count        

    def createCodeOutput(self, nodeId, globals=set()):
        if(not self.hasNode(nodeId)):
            raise ValueError("the nodeId " + str(nodeId) + "is not in the graph")

        node = self.getNode(nodeId)

        if node['numChildren'] == 0:
            if not self.isOperator(node['value']):
                if node['value'].startswith("_") or node['value'].startswith("DENDRO"):
                    output = node['value']
                    if nodeId in globals:
                        output = output + '[pp]'
                else:
                    output =  str(node['value'])
            else:
                output = node['nodeID']
        else:
            if node['numChildren'] != 2:
                raise ValueError("Node has a a number of children that is not 2")
            else:
                if node['value'] == 'sqrt':
                    output = self.sqrtCode(nodeId, globals)
                elif node['value'] == 'pow':
                    output = self.powCode(nodeId, globals)
                elif node['value'] == '+' or node['value'] == '-':
                    output = self.addSubCode(nodeId, globals)
                elif node['value'] == '*' or node['value'] == '/':
                    output = self.multDivCode(nodeId, globals)
                else:
                    raise ValueError('unRecognized value: ' + node['value'])
        return output

    def addSubCode(self, nodeId, globals):

        node = self.getNode(nodeId)
        if node['leftChild']['value'] == 'Negate':
            raise ValueError('cant have negate in addition or subtraction')
        if node['rightChild']['value'] == 'Negate':
            raise ValueError('cant have negate in addition or subtraction')

        return self.createCodeOutput(node['leftChild']['nodeID'], globals) + node['value'] + self.createCodeOutput(node['rightChild']['nodeID'], globals)
    
    def multDivCode(self, nodeId, globals):

        node = self.getNode(nodeId)
        if node['leftChild']['value'] == 'Negate' and node['rightChild']['value'] == 'Negate':   
            raise ValueError('cannot have double negate')
        elif node['leftChild']['value'] == 'Negate':
            return '-(' + self.createCodeOutput(node['rightChild']['nodeID'], globals) + ')'
        elif node['rightChild']['value'] == 'Negate':
            return '-(' + self.createCodeOutput(node['leftChild']['nodeID'], globals) + ')'
        else:
            left = self.createCodeOutput(node['leftChild']['nodeID'], globals)
            if node['leftChild']['value'] == '+' or node['leftChild']['value'] == '-':
                left = '(' + left + ')'
            right = self.createCodeOutput(node['rightChild']['nodeID'], globals)
            if node['rightChild']['value'] == '+' or node['rightChild']['value'] == '-':
                right = '(' + right + ')'
            if node['value'] == '*':
                return left + node['value'] + right
            else:
                return '(' + left + node['value'] + right + ')'
    
    def powCode(self, nodeId, globals):

        node = self.getNode(nodeId)
        leftCode = self.createCodeOutput(node['leftChild']['nodeID'], globals)
        rightCode = self.createCodeOutput(node['rightChild']['nodeID'], globals)

        return 'pow('+leftCode + ',' + rightCode + ')'

    def sqrtCode(self, nodeId, globals):

        node = self.getNode(nodeId)
        return 'sqrt(' + self.createCodeOutput(node['leftChild']['nodeID'], globals) + ')'

    def nonUniformStaging(self, disjoint_set = None):
        "stages a tree comprised of both multiplication and addition nodes."

        multiplyTree = self.splitAddMuliplyTrees()
        additionTree = self

        if disjoint_set != None:
            for id in multiplyTree.idDisplayMap.keys():
                disjoint_set.add(id)
        else:
            disjoint_set = ds.disjoint_set(multiplyTree.idDisplayMap.keys())

        multiplyTree.uniformStaging(disjoint_set)
        additionTree.uniformStaging(disjoint_set)
        additionTree.remergeMultiplicationTree(multiplyTree)

    def splitAddMuliplyTrees(self):
        "splits of the lower multiplication tree. Lower part of tree, multiplication symbols and upper tree, addition symbols. The source multiplies become leaf nodes in addition tree"

        additionTree = self
        multiplyTree = expressionTree()

        while len(additionTree.mult_sources) != 0:
            reductionNodes = list(additionTree.mult_sources)
            for reductionNodeId in reductionNodes:
                multiplyTree.addExpresionTree(additionTree.createSubTree(reductionNodeId))

        return multiplyTree

    def remergeMultiplicationTree(self, mult_tree):
        "used by the nonUniform Staging method to reinsert the multiplication tree into the addition tree"

        if(not isinstance(mult_tree, expressionTree)):
            raise ValueError("the passed object is not an expression tree")

        leafIndexes = sorted(mult_tree.leavesNodeMap.keys())
        for leafIndex in leafIndexes:
            nodeIds = mult_tree.leavesNodeMap[leafIndex]
            for nodeId in nodeIds:
                if(leafIndex == 1):
                    if(not self.hasNode(nodeId)):
                        self.addLeafNode(nodeId, mult_tree.getNodeValue(nodeId), mult_tree.idDisplayMap[nodeId])
                    else:
                        if(not self.are_identical(self.getNode(nodeId), mult_tree.getNode(nodeId))):
                            raise ValueError("The new tree has a leaf node id/value conflict at id: " + str(nodeId)) 
                else:
                    #handle adding nonLeaf nodes
                    if self.hasNode(nodeId):
                        mult_node = mult_tree.getNode(nodeId)
                        self.addChild(nodeId, mult_node['leftChild']['nodeID'])
                        self.addChild(nodeId, mult_node['rightChild']['nodeID'])
                        self.getNode(nodeId)['value'] = mult_node['value']
                    else:
                        node = mult_tree.getNode(nodeId)
                        self.addNonLeafNode(nodeId, node['value'], node['leftChild']['nodeID'], node['rightChild']['nodeID'])   
   
    def registerAdaptTrees(self, registerSize):
        "reduces a tree into smaller trees. Choosing nodes with the largest number of parent to extract first. This does change the tree"
        subtrees = []

        
        for i in range(registerSize):
            reductionNodeId = ''
            reduction_node_indegree = 0

            for source_id in self.sources:
                best_id, best_in_degree = self.findMaxInDegree(source_id,cache={}, include_leaves=False)
                if best_in_degree > reduction_node_indegree:
                    reductionNodeId = best_id
                    reduction_node_indegree = best_in_degree
                    
            subtrees.append(self.createSubTree(reductionNodeId))
        return subtrees

    def cacheAdaptTrees(self, cacheSize):
        "reduces trees into smaller trees that are smaller than or equal to the cache Size. This algorithm return a list of trees in the order they should be evaluated. Does not change current tree"

        subtrees = []
        parentTree = expressionTree()
        parentTree.addExpresionTree(self)
            
        while parentTree.largestLeafDependecy() > cacheSize:
            reductionNodeId = parentTree.findReductionNode(cacheSize)
            subtrees.append(parentTree.createSubTree(reductionNodeId))
        subtrees.append(parentTree)

        return subtrees

    def findReductionNode(self, dependencyTarget):
        "finds the best node in the tree that must be reduced to hit the dependency target. Returns the node Id of the necessary target."

        closestLength = -1
        closestId = -1

        best_id = None
        best_distance = dependencyTarget
        best_in_degree = -1

        #cache = {}
        #check over all source
        for sourceId in self.sources:
            target_id, target_distance, target_in_degree = self.findClosestLeafDependecy(dependencyTarget, sourceId)

            if target_distance < best_distance:
                best_id = target_id
                best_distance = target_distance
                best_in_degree = target_in_degree
            elif target_distance == best_distance and target_in_degree > best_in_degree:
                best_id = target_id
                best_distance = target_distance
                best_in_degree = target_in_degree

        return best_id

    def findMaxInDegree(self, nodeId, cache ={}, include_leaves = True):
        "finds the node with the largets in degree that also fulfills the dependecy target"

        if(not self.hasNode(nodeId)):
            raise ValueError("the nodeId is not in the graph")

        if(nodeId in cache):
            bestId = cache[nodeId]
            if self.getNode(bestId)['numChildren'] == 0:
                return bestId, -1
            return cache[nodeId], len(self.getParents(cache[nodeId]))

        

        node = self.getNode(nodeId)

        if not include_leaves and node['numChildren'] == 0:
            cache[nodeId] = nodeId
            return nodeId, -1


        currentBestInDegree = len(self.getParents(nodeId))
        bestNodeId = nodeId        

        if node['leftChild'] != None:
            leftId, leftDegree = self.findMaxInDegree(node['leftChild']["nodeID"], cache, include_leaves)
            if leftDegree > currentBestInDegree:
                    currentBestInDegree = leftDegree
                    bestNodeId = leftId

        if node['rightChild'] != None:
            rightId, rightDegree = self.findMaxInDegree(node['rightChild']["nodeID"], cache, include_leaves)
            if rightDegree > currentBestInDegree:
                    currentBestInDegree = rightDegree
                    bestNodeId = rightId

        cache[nodeId] = bestNodeId

        return bestNodeId, currentBestInDegree
            
    def findClosestLeafDependecyGreedy(self, target):

        
        nodeIds = set()

        bestDegree = target

        while len(nodeIds) == 0:
            while bestDegree not in self.leavesNodeMap:
                bestDegree  = bestDegree -1
            nodeIds = self.leavesNodeMap[bestDegree]
            bestDegree = bestDegree -1

        bestId = -1
        bestInDegree =-1

        for nodeId in nodeIds:
            if bestId == -1:
                bestId = nodeId
                bestInDegree = len(self.getParents(nodeId))
            else:
                currentDegree = len(self.getParents(nodeId))
                if currentDegree>bestInDegree:
                    bestId = nodeId
                    bestInDegree = currentDegree
        if bestId == -1:
            raise ValueError('something went wrong when finding the best node')
        numLeafDependents = self.getNumLeafDependents(bestId)
        return bestId, target - numLeafDependents

    def findClosestLeafDependecy(self, depedencyTarget, nodeId):
        "finds the node with the closest Leaf Dependency. The search starts at nodeId. returns the nodeId of the node and how far off from the dependency target"

        if(not self.hasNode(nodeId)):
            raise ValueError("the nodeId is not in the graph")

        stack = []        
        stack.append(nodeId)

        
        best_id = None
        best_num_dependents = -1
        best_in_degree = None

        while len(stack) != 0:
            challenger_id = stack.pop()
            num_dependents = self.getNumLeafDependents(challenger_id)
            in_degree = len(self.getParents(challenger_id))

            if num_dependents > depedencyTarget:
                node = self.getNode(challenger_id)
                stack.append(node['leftChild']['nodeID'])
                stack.append(node['rightChild']['nodeID'])
            elif num_dependents > best_num_dependents:
                best_id = challenger_id
                best_num_dependents = num_dependents
                best_in_degree = in_degree
            elif num_dependents == best_num_dependents and in_degree > best_in_degree:
                best_id = challenger_id
                best_num_dependents = num_dependents
                best_in_degree = in_degree

        return best_id, depedencyTarget - best_num_dependents, best_in_degree

    def largestLeafDependecy(self):
        "finds the node with the largest leaf dependecies. Only considers the number of leaves and not the number of times each leaf is needed. Returns the largest number of leaf dependecies."
        
        target = -1
        for sourceId in self.sources:
            numLeafDependents = self.getNumLeafDependents(sourceId)
            if(numLeafDependents>target):
                target = numLeafDependents
        return target

    def createSubTree(self, nodeId):
        "copies a a part of the expression tree. Starts at nodeId and copies everything down. Returns a new tree with the copied chunk."

        if(not self.hasNode(nodeId)):
            raise ValueError("The node does not exist in the expression tree")

        node = self.getNode(nodeId)
        subtree = expressionTree()
        subtree.createSubTreeHelper(nodeId, self)

        #used at the end to remove children if possible
        leftChildId =  node['leftChild']["nodeID"]
        rightChildId = node['rightChild']["nodeID"]

        if(node['value'] == "*" and nodeId in self.mult_sources):
            self.remove_mult_source(nodeId)

        #remove connection to both children. Remove left child twice since the right child will be moved to the left child after first removal
        self.removeChild(nodeId, node['leftChild']["nodeID"])
        self.removeChild(nodeId, node['rightChild']["nodeID"])

        #update node to become a leaf node
        node['value'] = nodeId
        self.initializeAttributes(node, nodeId, node['value'])      
        node['numLeafDependents'] = 1
        node['leafDependecies'][node['value']] = 1
        for parentId in self.getParents(nodeId):   
            self.updateParentDependeciesHelperAdd(parentId,nodeId)

        #remove children if necesary
        if(len(self.getParents(leftChildId))==0):
            self.removeNode(leftChildId)
        elif(self.getNodeValue(leftChildId) == "*" and len(self.getMultiplyParents(leftChildId)) == 0):
            self.add_mult_source(leftChildId)
        if(len(self.getParents(rightChildId))==0):
            self.removeNode(rightChildId)
        elif(self.getNodeValue(rightChildId) == "*" and len(self.getMultiplyParents(rightChildId)) == 0):
            self.add_mult_source(rightChildId)

        return subtree        

    def updateParentDependeciesHelperRemove(self, nodeId, childId):
        "updates dependecies when a child is removed from the graph. The algorithm updates the given node and recurses upwards to update all parent."

        childNode = self.getNode(childId)
        node = self.getNode(nodeId)

        #handle edge case with a leaf node
        if(node["numLeafDependents"] in self.leavesNodeMap):
            self.leavesNodeMap[node["numLeafDependents"]].remove(nodeId)
            if(len(self.leavesNodeMap[node['numLeafDependents']]) == 0 ):
                del self.leavesNodeMap[node['numLeafDependents']]
                
        self.removeChildDependecies(node, childNode) 
        self.addLeafMapping(node["numLeafDependents"],nodeId)

        '''
        for parentId in self.getParents(nodeId):     
            self.updateParentDependeciesHelperRemove(parentId, childId)
        '''

        for ancestorId in self.get_allAncestors(nodeId):
            node = self.getNode(ancestorId)

            #handle edge case with a leaf node
            if(node["numLeafDependents"] in self.leavesNodeMap):
                self.leavesNodeMap[node["numLeafDependents"]].remove(ancestorId)
                if(len(self.leavesNodeMap[node['numLeafDependents']]) == 0 ):
                    del self.leavesNodeMap[node['numLeafDependents']]
                    
            self.removeChildDependecies(node, childNode) 
            self.addLeafMapping(node["numLeafDependents"],ancestorId)  

    def updateParentDependeciesHelperAdd(self, nodeId, childId):
        "updates dependecies when a child is added to the graph. The algorithm updates the given node and recurses upwards to update all parent."

        childNode = self.getNode(childId)
        node = self.getNode(nodeId)

        #handle edge case with a leaf node
        if(node["numLeafDependents"] in self.leavesNodeMap):
            self.leavesNodeMap[node["numLeafDependents"]].remove(nodeId)
            if(len(self.leavesNodeMap[node['numLeafDependents']]) == 0 ):
                del self.leavesNodeMap[node['numLeafDependents']]
                
        self.addLeafDependecies(node, childNode) 
        self.addLeafMapping(node["numLeafDependents"],nodeId)

        '''
        for parentId in self.getParents(nodeId):     
            self.updateParentDependeciesHelperAdd(parentId, childId) 
        '''

        #do the same work while recursing up
        for ancestorId in self.get_allAncestors(nodeId):
            node = self.getNode(ancestorId)

            #handle edge case with a leaf node
            if(node["numLeafDependents"] in self.leavesNodeMap):
                self.leavesNodeMap[node["numLeafDependents"]].remove(ancestorId)
                if(len(self.leavesNodeMap[node['numLeafDependents']]) == 0 ):
                    del self.leavesNodeMap[node['numLeafDependents']]
                    
            self.addLeafDependecies(node, childNode) 
            self.addLeafMapping(node["numLeafDependents"],ancestorId)      
                
    def createSubTreeHelper(self, nodeId, oldTree):
        "This is a recursive helper method for the subtree"

        if(not oldTree.hasNode(nodeId)):
            raise ValueError("The node does not exist in the old expression tree")

        if(not self.hasNode(nodeId)):
            node = oldTree.getNode(nodeId)
            nodeDisplay = oldTree.idDisplayMap[nodeId]

            if(node['numChildren'] == 0):
                self.addLeafNode(nodeId, node['value'], nodeDisplay)
            else:
                #recursively add left child
                leftChildNode = oldTree.getNode(node['leftChild']["nodeID"])
                self.createSubTreeHelper(leftChildNode["nodeID"], oldTree)

                #recursively add right child
                rightChildNode = oldTree.getNode(node['rightChild']["nodeID"])
                self.createSubTreeHelper(rightChildNode["nodeID"], oldTree)

                self.addNonLeafNode(nodeId, node['value'], leftChildNode["nodeID"], rightChildNode["nodeID"], nodeDisplay,)          

    def evaluateNode(self, nodeId, symbolMapping = None, cache = None):
        "evaluates the value of the node in the expression tree. Symbol Mapping allows for nondecimal leaf nodes to be assigned a value"

        if(not self.hasNode(nodeId)):
            raise ValueError("The node does not exist in the expression tree")
        if(cache != None and not isinstance(cache, dict)):
            raise ValueError("The cache must be a dictionary")
        if(symbolMapping != None and not isinstance(symbolMapping, dict)):
            raise ValueError("The symbol mapping must be a dictionary")
        
        if(cache != None and nodeId in cache):
            return cache[nodeId]

        if(cache == None):
            cache = {}

        node = self.getNode(nodeId)
        value = None
        if(node['numChildren'] == 0):
            if(symbolMapping != None and node['value'] in symbolMapping):
                value = symbolMapping[node['value']]
            else:
                try:
                    value = float(node['value'])
                except:
                    raise TypeError("The leaf node value could not be parsed into a float and does not appear in the symbol map")
        else:
            if(node['value'] == '+'):
                value = self.evaluateNode(node['leftChild']['nodeID'], symbolMapping, cache) + self.evaluateNode(node['rightChild']['nodeID'], symbolMapping, cache)
            elif(node['value'] == '*'):
                value = self.evaluateNode(node['leftChild']['nodeID'], symbolMapping, cache) * self.evaluateNode(node['rightChild']['nodeID'], symbolMapping, cache)
            else:
                raise TypeError("The operand in the node is not recoginzed")

        cache[nodeId] = value
        return value

    def addExpresionTree(self, eTree):
        "adds all nodes and edges from another expression tree to this expression tree"

        if(not isinstance(eTree, expressionTree)):
            raise ValueError("the passed object is not an expression tree")

        leafIndexes = sorted(eTree.leavesNodeMap.keys())
        for leafIndex in leafIndexes:
            nodeIds = eTree.leavesNodeMap[leafIndex]
            for nodeId in nodeIds:
                if(leafIndex == 1):
                    if(not self.hasNode(nodeId)):
                        self.addLeafNode(nodeId, eTree.getNodeValue(nodeId), eTree.idDisplayMap[nodeId])
                    else:
                        if(not self.are_identical(self.getNode(nodeId), eTree.getNode(nodeId))):
                            raise ValueError("The new tree has a leaf node id/value conflict at id: " + str(nodeId)) 
                else:
                    #handle adding nonLeaf nodes
                    if self.hasNode(nodeId):
                        if not self.are_identical(self.getNode(nodeId), eTree.getNode(nodeId)):
                            raise ValueError("The new tree has a leaf node id/value conflict at id: " + str(nodeId))
                    else:
                        node = eTree.getNode(nodeId)
                        self.addNonLeafNode(nodeId, node['value'], node['leftChild']['nodeID'], node['rightChild']['nodeID'])   

    def uniformStaging(self, disjoint_set = None):
        "perfoms staging operation on the graph under the assumption that all operators are the same and order does not matter"

        if disjoint_set != None and not isinstance(disjoint_set, ds.disjoint_set):
            raise ValueError("The disjoint set must be comprised of the class disjoint_set")

        #add a disjoint set to make staging faster
        if disjoint_set == None:
            disjoint_set = ds.disjoint_set(self.idDisplayMap.keys())
        else:
            for id in self.idDisplayMap.keys():
                disjoint_set.add(id)

        keyList = sorted(self.leavesNodeMap)
        for key in keyList:

            #get all nodes with a given dependecy
            nodeList = list(self.leavesNodeMap[key])

            #perform n^2 operation for mapping
            for i in range(0, len(nodeList)):
                leftNodeId = nodeList[i]

                if(not self.hasNode(leftNodeId)):
                    continue

                for j in range(i+1, len(nodeList)):
                    rightNodeId = nodeList[j]

                    if(not self.hasNode(rightNodeId)):
                        continue

                    equivalent = disjoint_set != None and disjoint_set.equivalent(leftNodeId, rightNodeId)
                    if(equivalent):
                        replacerId = disjoint_set.getKey(leftNodeId)
                        self.replaceNode(replacerId, rightNodeId)
                        self.replaceNode(replacerId,leftNodeId)
                    elif(self.uniformIsIsomporphic(leftNodeId, rightNodeId)):
                        if disjoint_set != None:
                            disjoint_set.add(leftNodeId)
                            disjoint_set.add(rightNodeId)
                            disjoint_set.union(leftNodeId, rightNodeId)

                            replacerId = disjoint_set.getKey(leftNodeId)
                            self.replaceNode(replacerId, rightNodeId)
                            self.replaceNode(replacerId,leftNodeId)
                        else:
                            self.replaceNode(leftNodeId,rightNodeId)

    def uniformIsIsomporphic(self, nodeAId, nodeBId):
        "determines if two nodes,nodeA and nodeB, are isomporphic. The method assumes the graph is comprised of a single operator. nodeAId - id of node A, NodeBId - id of node B"
        if(not self.hasNode(nodeAId)):
            #raise ValueError("node A does not exist in the graph")
            return False
        if(not self.hasNode(nodeBId)):
            #raise ValueError("node B does not exist in the graph")
            return False

        nodeA = self.getNode(nodeAId)
        nodeB = self.getNode(nodeBId)    
        if(not nodeA['numChildren'] == nodeB['numChildren']):
            return False
        elif(nodeA['numChildren'] == 0):
            if self.isOperator(nodeA['value']) or self.isOperator(nodeB['value']):
                return False
            return nodeA['value'] == nodeB['value']
        elif(not nodeA['numLeafDependents'] == nodeB['numLeafDependents']):
            return False
        else:        
            return nodeA['value'] == nodeB['value'] and nodeA['leafDependecies'] == nodeB['leafDependecies']
    
    def addLeafNode(self, nodeId, value, nodeDisplay = None):
        "add node to graph, nodeID, nodeDisplay display for printing, value of the leaf"

        if nodeDisplay == None:
            nodeDisplay = self.generateDisplay(nodeId, value)

        if(self.graph.has_node(nodeDisplay)):
            raise ValueError("node already exists in graph")

        
            
        #add Node to the graph
        self.idDisplayMap[nodeId] = nodeDisplay
        self.graph.add_node(nodeDisplay)        
        self.sources.add(nodeId)

        #Initialize all attributes used
        node = self.getNode(nodeId)
        self.initializeAttributes(node, nodeId, value)      
        node['numLeafDependents'] = 1
        node['leafDependecies'][node['value']] = 1
        self.addLeafMapping(node['numLeafDependents'], nodeId)

    def addNonLeafNode(self, nodeId, value, leftChildId, rightChildId, nodeDisplay = None, pushDown = False):
        "add node to graph, nodeId, nodeDisplay display for printing, value of the leaf, leftChildId, rightChildId"

        if(self.graph.has_node(nodeDisplay)):
            raise ValueError("node already exists in graph")
        if(not self.isOperator(value)):
            raise ValueError("value is not a valid operator for a nonLeafNode")
        if(not self.hasNode(leftChildId)):
            raise ValueError("left Child is not in the graph")
        if(not self.hasNode(rightChildId)):
            raise ValueError("right Child is not in the graph")

        if nodeDisplay == None:
            nodeDisplay = self.generateDisplay(nodeId, value)

        #add Node to the graph
        self.idDisplayMap[nodeId] = nodeDisplay
        self.graph.add_node(nodeDisplay) 

        #update sources set       
        self.sources.add(nodeId)
        if(leftChildId in self.sources):
            self.sources.remove(leftChildId)
        if(rightChildId in self.sources):
            self.sources.remove(rightChildId)

        #Initialize all attributes used
        node = self.getNode(nodeId)
        self.initializeAttributes(node, nodeId, value)      

        self.addChild(nodeId, leftChildId)
        self.addChild(nodeId, rightChildId)

        
        if(value == "*" and pushDown):            
            if(not nodeId in self.mult_sources):
                self.add_mult_source(nodeId)
            self.mult_push_down(nodeId)       

    def mult_push_down(self, nodeId):
        "pushes the multiplication operator all the way to the bottom or until children are either leaf children or other multiplication symbols"

        if(not self.hasNode(nodeId)):
            raise ValueError("the node does not exist")

        node = self.getNode(nodeId)
        if(node['value'] != '*'):
            raise ValueError("the node does not have a multiplication symbol.")

        leftChild = node["leftChild"]
        rightChild = node["rightChild"]

        prod_sums = leftChild['value'] == '+' and rightChild['value'] =='+'
        if(prod_sums):
            self.prod_sum_helper(nodeId, leftChild["nodeID"], rightChild["nodeID"])
        elif(leftChild['value'] == '+'):
            self.distribute_prod(nodeId, rightChild["nodeID"], leftChild["nodeID"])
        elif(rightChild['value'] == '+'):
            self.distribute_prod(nodeId, leftChild["nodeID"], rightChild["nodeID"])
        else:
            return

    def prod_sum_helper(self, nodeId, leftSumId, rightSumId):
        "method distributes a product of sums into a sum of products"

        if(not self.hasNode(nodeId)):
            raise ValueError("the node does not exist")
        if(not self.hasNode(leftSumId)):
            raise ValueError("the leftSumId does not exist")
        if(not self.hasNode(rightSumId)):
            raise ValueError("the rightSumId does not exist")

        node = self.getNode(nodeId)
        leftSumNode = self.getNode(leftSumId)
        rightSumNode = self.getNode(rightSumId)

        if(node['value'] != '*'):
            raise ValueError("the node is not a multiplication node")
        if(leftSumNode['value'] != '+'):
            raise ValueError("the leftSum is not an addition node")
        if(rightSumNode['value'] != "+"):
            raise ValueError("the rightSum is not an addition node")

        self.removeNode(nodeId, False)
        self.addNonLeafNode(nodeId, "+", leftSumId, rightSumId)

        #create new distributed multiplies
        '''format:
        1st/3rd letter: left or right
        ex: lslm_Id -> left sum left multiply
        '''
        lslm_Id = str(nodeId) + "_lslm" 
        lsrm_Id = str(nodeId) + "_lsrm"
        rslm_Id = str(nodeId) + "_rslm"
        rsrm_Id = str(nodeId) + "_rsrm"

        #get correct children
        leftSumLeftChildId = leftSumNode['leftChild']['nodeID']
        leftSumRightChildId = leftSumNode['rightChild']['nodeID']
        rightSumLeftChildId = rightSumNode['leftChild']['nodeID']
        rightSumRightChildId = rightSumNode['rightChild']['nodeID']

        #add children
        self.addNonLeafNode(lslm_Id,'*', leftSumLeftChildId, rightSumLeftChildId)
        self.addNonLeafNode(lsrm_Id,'*', leftSumLeftChildId, rightSumRightChildId)        
        self.addNonLeafNode(rslm_Id,'*', leftSumRightChildId, rightSumLeftChildId)
        self.addNonLeafNode(rsrm_Id,'*', leftSumRightChildId, rightSumRightChildId)

        #add sums
        ls_Id = str(nodeId) + "_ls" 
        rs_Id = str(nodeId) + "_rs"
        self.addNonLeafNode(ls_Id, "+", lslm_Id, lsrm_Id)
        self.addNonLeafNode(rs_Id, "+",rslm_Id, rsrm_Id)

        #update original node
        self.removeChild(nodeId,leftSumId)
        self.removeChild(nodeId,rightSumId)
        self.addChild(nodeId, ls_Id)
        self.addChild(nodeId,rs_Id)

    def add_mult_source(self, nodeId):
        "adds the nodeId to the mult_sources set"

        if(not self.hasNode(nodeId)):
            raise ValueError("the node is not in the tree " + str(nodeId))

        node = self.getNode(nodeId)

        if(node['value'] != "*"):
            raise ValueError("the node does not have a multiplication operator")

        self.mult_sources.add(nodeId)
        if(node['numChildren'] != 0):
            leftId = node['leftChild']['nodeID']
            rightId = node['rightChild']['nodeID']

            if(leftId in self.mult_sources):
                self.mult_sources.remove(leftId)
            if(rightId in self.mult_sources):
                self.mult_sources.remove(rightId)

    def remove_mult_source(self, nodeId):
        "removes the nodeID from the mult_sources set"

        if(not self.hasNode(nodeId)):
            raise ValueError("the node is not in the tree " + str(nodeId))       

        node = self.getNode(nodeId)

        if(node['value'] != "*"):
            raise ValueError("the node does not have a multiplication operator")

        self.mult_sources.remove(nodeId)
        if(node['numChildren'] != 0):
            leftNode = node['leftChild']
            rightNode = node['rightChild']

            if(len(self.getParents(leftNode['nodeID'])) == 0 and leftNode['value'] == "*" and leftNode['numChildren'] != 0):
                self.mult_sources.add(leftNode['nodeID'])
            if(len(self.getParents(rightNode['nodeID'])) == 0 and rightNode['value'] == "*" and rightNode['numChildren'] != 0):
                self.mult_sources.add(rightNode['nodeID'])

    def distribute_prod(self, nodeId, coefficientId, sumId):
        "method distributes a coefficient into a sum"

        if(not self.hasNode(nodeId)):
            raise ValueError("the node does not exist")
        if(not self.hasNode(sumId)):
            raise ValueError("the sumId does not exist")
        if(not self.hasNode(coefficientId)):
            raise ValueError("the coefficient Id does not exist")

        node = self.getNode(nodeId)
        sum_node = self.getNode(sumId)
        
        if(sum_node['value'] != '+'):
            raise ValueError("The sumId does not have an addition symbol as its value.")
        if(node['value'] != "*"):
            raise ValueError("the node does not have a multiplication symbol")
        if(self.getNode(coefficientId)['value'] == '+'):
            raise ValueError("the coefficient is actaully a sum which is not allowed")

        self.removeNode(nodeId, False)
        self.addNonLeafNode(nodeId, "+", node['leftChild']['nodeID'], node['rightChild']['nodeID'])
        

        leftSumdId = sum_node['leftChild']['nodeID']
        rightSumdId = sum_node['rightChild']['nodeID']

        #create distributed products
        h1_Id = str(nodeId) + "_dhl" 
        self.addNonLeafNode(h1_Id,'*', coefficientId, leftSumdId)
        h2_Id = str(nodeId) + "_dhr"
        self.addNonLeafNode(h2_Id, "*", coefficientId, rightSumdId) 

        #update old multiplication
        self.removeChild(nodeId, node['leftChild']['nodeID'])
        self.removeChild(nodeId, node['rightChild']['nodeID'])
        self.addChild(nodeId, h1_Id)
        self.addChild(nodeId, h2_Id)

    def generateDisplay(self, nodeId, value):
        "generate a useful display name for the node that will be used for printing the expression tree in createGraphPicture()"
        return "id " + str(nodeId) + "\nvalue " + str(value)

    def removeNode(self, nodeId, recursive = True):
        "removes a node from the graph. The node can only be removed if it does not have any parents"

        if(not self.hasNode(nodeId)):
            raise ValueError("The node does not exist in the graph")
        
        if(len(self.getParents(nodeId)) > 0):
            raise ValueError("The node cannot be removed since its parent depend on it. Use replaceNode() instead")      
        
        node = self.getNode(nodeId)
        leftChildId = None
        rightChildId = None
        if(not node['leftChild'] == None):
            leftChildId = node['leftChild']['nodeID']
        if(not node['rightChild'] == None):
            rightChildId = node['rightChild']['nodeID']

        #update sources set
        if(nodeId in self.sources):
            self.sources.remove(nodeId)  
        if(nodeId in self.mult_sources):
            self.remove_mult_source(nodeId)  

        #update leaves node map
        self.leavesNodeMap[node['numLeafDependents']].remove(nodeId)       

        self.graph.remove_node(self.idDisplayMap[nodeId])

        #update display map
        del self.idDisplayMap[nodeId]

        #remove unecessary children
        if(recursive):
            if(leftChildId != None and self.hasNode(leftChildId) and len(self.getParents(leftChildId)) == 0):
                self.removeNode(leftChildId)
            if(rightChildId != None and self.hasNode(rightChildId) and len(self.getParents(rightChildId)) == 0):
                self.removeNode(rightChildId)

    def replaceNode(self, replacerId, replaceeId):
        "replaces all instances of node with replaceeId with the corresponding node with replacerId"
        if(not self.hasNode(replaceeId)):
            raise ValueError("The replacee node does not exist in the graph")
        if(not self.hasNode(replacerId)):
            raise ValueError("The replacer node does not exist in the graph")
        if(replacerId == replaceeId):
            return

        for parentId in self.getParents(replaceeId):
            self.removeChild(parentId, replaceeId)
            self.addChild(parentId, replacerId)
        self.removeNode(replaceeId)

    def addChild(self, parentId, childId):
        parentNode = self.getNode(parentId)
        childNode = self.getNode(childId)

        #parent can have a max of 2 children
        if(parentNode['numChildren']>=2):
            raise ValueError("The parent Node already contains 2 Children")
        
        #update souces set
        if(len(self.getParents(childId)) == 0 and childId in self.sources):
            self.sources.remove(childId)

        #add an edge
        self.graph.add_edge(self.idDisplayMap[parentId], self.idDisplayMap[childId])
        parentNode['numChildren'] = parentNode['numChildren'] + 1        
            
        #update leaf dependecies
        self.updateParentDependeciesHelperAdd(parentId,childId)

        #update children
        if(parentNode['leftChild'] == None):
            parentNode['leftChild'] = childNode
        else:
            parentNode['rightChild'] = childNode

    def removeChild(self, parentId, childId):
        "remove a connection in the expression tree"

        parentNode = self.getNode(parentId)
        childNode = self.getNode(childId)

        #update leaves node map
        self.updateParentDependeciesHelperRemove(parentId,childId)    

        #remove an edge
        self.graph.remove_edge(self.idDisplayMap[parentId], self.idDisplayMap[childId])
        parentNode['numChildren'] = parentNode['numChildren'] - 1

        #update sources set
        if(len(self.getParents(childId)) == 0):
            self.sources.add(childId)

        if(parentNode['leftChild'] == childNode):
            parentNode['leftChild'] = None
        else:
            parentNode['rightChild'] = None

    def initializeAttributes(self, node, nodeId, value):
        "initializes all the attributes to be used by the node, node - graphix representation of node, nodeId, value"
        node['nodeID'] = nodeId
        node['numChildren'] = 0
        node['numLeafDependents'] = 0
        node['leftChild'] = None
        node['rightChild'] = None
        node['leafDependecies'] = {}
        node['value'] = value

    def getParents(self, nodeId):
        "returns a list of all the parentIds of the input node"

        if(not self.hasNode(nodeId)):
            raise ValueError("The node does not exist in the graph")

        nodeLabel = self.idDisplayMap[nodeId]
        parentNodes = []
        for parentLabel in self.graph.predecessors(nodeLabel):
            parentNodes.append(self.graph.nodes[parentLabel]['nodeID'])
        return parentNodes

    def get_allAncestors(self,nodeId):
        "get all ancestors until the source. This includes duplicates"

        allAncestors = []
        queue = deque()
        queue.append(nodeId)
        while len(queue) != 0:
            id = queue.popleft()
            for parentId in self.getParents(id):
                allAncestors.append(parentId)
                queue.append(parentId)
        return allAncestors

    def getMultiplyParents(self, nodeId):
        "returns a list of all parents whose value is a multiply operator"

        nodeLabel = self.idDisplayMap[nodeId]
        parentNodes = []
        for parentLabel in self.graph.predecessors(nodeLabel):
            parentNode = self.graph.nodes[parentLabel]
            if(parentNode['value'] == "*"):
                parentNodes.append(parentNode['nodeID'])
        return parentNodes

    def getNode(self, nodeId):
        if(not self.hasNode(nodeId)):
            raise ValueError("The nodeId is not in the graph.")
        else:
            nodeDisplay = self.idDisplayMap[nodeId]
            return self.graph.node[nodeDisplay]

    def hasNode(self, nodeId):
        return nodeId in self.idDisplayMap

    def getNodeValue(self, nodeId):
        if(not self.hasNode(nodeId)):
            raise ValueError("the nodeId does not exist in the graph")
        return self.getNode(nodeId)['value']

    def getNumLeafDependents(self, nodeId):
        "returns the number of individual leaves the node depends on. This does not consider repeated nodes."
        
        if(not self.hasNode(nodeId)):
            raise ValueError("the nodeId does not exist in the graph")
        return len(self.getLeafDependents(nodeId))

    def getLeafDependents(self,nodeId):
        if(not self.hasNode(nodeId)):
            raise ValueError("the nodeId does not exist in the graph")

        non_scalar_keys = []
        #note testing cases where scalars are not treated as leaf dependecies
        for key in self.getNode(nodeId)['leafDependecies'].keys():
            try:
                float(key)
            except ValueError:
                non_scalar_keys.append(key)
        return non_scalar_keys

        #return self.getNode(nodeId)['leafDependecies'].keys()

    def addLeafDependecies(self, parentNode, childNode):
        "adds all leaf dependecies from the childNode to the parentNode"

        for leafId in childNode['leafDependecies']:
            if(leafId in parentNode['leafDependecies']):
                parentNode['leafDependecies'][leafId] = parentNode['leafDependecies'][leafId] + childNode['leafDependecies'][leafId]
            else:
                parentNode['leafDependecies'][leafId] = childNode['leafDependecies'][leafId]
            parentNode['numLeafDependents'] = parentNode['numLeafDependents'] + childNode['leafDependecies'][leafId]
 
    def addLeafMapping(self, numLeaves, nodeId):
        "adds the nodeID to the leaves Mapping, numLeaves - number of leaves that the node depends on, nodeId - id of the node"
        #add leaf node to leaves mapping
        if(not numLeaves in self.leavesNodeMap):
            self.leavesNodeMap[numLeaves]= {nodeId}
        else:
            self.leavesNodeMap[numLeaves].add(nodeId)        

    def removeChildDependecies(self, parentNode, childNode):
        "removes the dependecies of the child from the parent and updates the parent's number of leaf Dependents"

        for key in childNode['leafDependecies'].keys():
            parentNode['leafDependecies'][key] = parentNode['leafDependecies'][key] - childNode['leafDependecies'][key] 
            if(parentNode['leafDependecies'][key] == 0 ):
                del parentNode['leafDependecies'][key]

            parentNode['numLeafDependents'] = parentNode['numLeafDependents'] - childNode['leafDependecies'][key]

            #this occurs after all children have been removed and the parent node becomes a leaf node
            if(parentNode['numLeafDependents'] == 0 ):
                parentNode['numLeafDependents'] = 1

    def isOperator(self, operator):
        "checks is the operator is a valid operator. Currently only supports addition"
        if(operator == "+"):
            return True
        elif(operator == "*"):
            return True
        elif(operator == '-'):
            return True
        elif(operator == '/'):
            return True
        elif(operator == 'pow'):
            return True
        elif(operator == 'sqrt'):
            return True
        else:
            return False

    def numNodes(self):
        return self.graph.number_of_nodes()

    def are_identical(self, left_node, right_node):
        "checks if two nodes are identical in terms of id, values and children if necessary. Useful when adding an additional expression tree"

        if(left_node['nodeID'] != right_node['nodeID']):
            return False
        if(left_node['value'] != right_node['value']):
            return False
        if(left_node['numChildren'] != right_node['numChildren']):
            return False
        if(left_node['numChildren']>0):
            return left_node['leftChild']['nodeID'] == right_node['leftChild']['nodeID'] and left_node['rightChild']['nodeID'] == right_node['rightChild']['nodeID']
        return True

    def createGraphPicture(self, filename):
        pyDot = nx.nx_pydot.to_pydot(self.graph)
        pyDot.write_png(str(filename) + ".png")