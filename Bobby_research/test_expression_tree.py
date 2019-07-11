import expressionTree as eTree

def test_add_single_node():
    tree = eTree.expressionTree()
    tree.addLeafNode(1, 5)
    assert tree.hasNode(1) == True

    assert tree.getNode(1)['nodeID'] == 1
    assert tree.getNode(1)['numChildren'] == 0
    assert tree.getNode(1)['numLeafDependents'] == 1
    assert tree.getNode(1)['leftChild'] == None
    assert tree.getNode(1)['rightChild'] == None
    assert tree.getNode(1)['leafDependecies'] == {5:1}
    assert tree.getNode(1)['value'] == 5

    assert tree.leavesNodeMap[1] == {1}
    assert tree.sources == {1}

def test_add_multiple_leaf_nodes():

    tree = eTree.expressionTree()

    #add leaf nodes
    tree.addLeafNode(0, 5)
    tree.addLeafNode(1, 5)
    tree.addLeafNode(2,  3)
    tree.addLeafNode(3,  4)
    tree.addLeafNode(4,  2)
    tree.addLeafNode(5,  1)
    tree.addLeafNode(6,  4)
    tree.addLeafNode(7,  2)
    tree.addLeafNode(8,  3)
    tree.addLeafNode(9,  5)
    tree.addLeafNode(10, 5)

    assert tree.hasNode(0) == True
    assert tree.getNodeValue(0) == 5
    assert tree.getNode(0)['nodeID'] == 0
    assert tree.getNode(0)['numChildren'] == 0
    assert tree.getNode(0)['numLeafDependents'] == 1
    assert tree.getNode(0)['leftChild'] == None
    assert tree.getNode(0)['rightChild'] == None
    assert tree.getNode(0)['leafDependecies'] == {5:1}
    assert tree.getNode(0)['value'] == 5

    assert tree.hasNode(1) == True
    assert tree.getNodeValue(1) == 5
    assert tree.getNode(1)['nodeID'] == 1
    assert tree.getNode(1)['numChildren'] == 0
    assert tree.getNode(1)['numLeafDependents'] == 1
    assert tree.getNode(1)['leftChild'] == None
    assert tree.getNode(1)['rightChild'] == None
    assert tree.getNode(1)['leafDependecies'] == {5:1}
    assert tree.getNode(1)['value'] == 5

    assert tree.hasNode(2) == True
    assert tree.getNodeValue(2) == 3
    assert tree.getNode(2)['nodeID'] == 2
    assert tree.getNode(2)['numChildren'] == 0
    assert tree.getNode(2)['numLeafDependents'] == 1
    assert tree.getNode(2)['leftChild'] == None
    assert tree.getNode(2)['rightChild'] == None
    assert tree.getNode(2)['leafDependecies'] == {3:1}
    assert tree.getNode(2)['value'] == 3

    assert tree.hasNode(3) == True
    assert tree.getNodeValue(3) == 4
    assert tree.getNode(3)['nodeID'] == 3
    assert tree.getNode(3)['numChildren'] == 0
    assert tree.getNode(3)['numLeafDependents'] == 1
    assert tree.getNode(3)['leftChild'] == None
    assert tree.getNode(3)['rightChild'] == None
    assert tree.getNode(3)['leafDependecies'] == {4:1}
    assert tree.getNode(3)['value'] == 4

    assert tree.hasNode(4) == True
    assert tree.getNodeValue(4) == 2
    assert tree.getNode(4)['nodeID'] == 4
    assert tree.getNode(4)['numChildren'] == 0
    assert tree.getNode(4)['numLeafDependents'] == 1
    assert tree.getNode(4)['leftChild'] == None
    assert tree.getNode(4)['rightChild'] == None
    assert tree.getNode(4)['leafDependecies'] == {2:1}
    assert tree.getNode(4)['value'] == 2

    assert tree.hasNode(5) == True
    assert tree.getNodeValue(5) == 1
    assert tree.getNode(5)['nodeID'] == 5
    assert tree.getNode(5)['numChildren'] == 0
    assert tree.getNode(5)['numLeafDependents'] == 1
    assert tree.getNode(5)['leftChild'] == None
    assert tree.getNode(5)['rightChild'] == None
    assert tree.getNode(5)['leafDependecies'] == {1:1}
    assert tree.getNode(5)['value'] == 1

    assert tree.hasNode(6) == True
    assert tree.getNodeValue(6) == 4
    assert tree.getNode(6)['nodeID'] == 6
    assert tree.getNode(6)['numChildren'] == 0
    assert tree.getNode(6)['numLeafDependents'] == 1
    assert tree.getNode(6)['leftChild'] == None
    assert tree.getNode(6)['rightChild'] == None
    assert tree.getNode(6)['leafDependecies'] == {4:1}
    assert tree.getNode(6)['value'] == 4

    assert tree.hasNode(7) == True
    assert tree.getNodeValue(7) == 2
    assert tree.getNode(7)['nodeID'] == 7
    assert tree.getNode(7)['numChildren'] == 0
    assert tree.getNode(7)['numLeafDependents'] == 1
    assert tree.getNode(7)['leftChild'] == None
    assert tree.getNode(7)['rightChild'] == None
    assert tree.getNode(7)['leafDependecies'] == {2:1}
    assert tree.getNode(7)['value'] == 2

    assert tree.hasNode(8) == True
    assert tree.getNodeValue(8) == 3
    assert tree.getNode(8)['nodeID'] == 8
    assert tree.getNode(8)['numChildren'] == 0
    assert tree.getNode(8)['numLeafDependents'] == 1
    assert tree.getNode(8)['leftChild'] == None
    assert tree.getNode(8)['rightChild'] == None
    assert tree.getNode(8)['leafDependecies'] == {3:1}
    assert tree.getNode(8)['value'] == 3

    assert tree.hasNode(9) == True
    assert tree.getNodeValue(9) == 5
    assert tree.getNode(9)['nodeID'] == 9
    assert tree.getNode(9)['numChildren'] == 0
    assert tree.getNode(9)['numLeafDependents'] == 1
    assert tree.getNode(9)['leftChild'] == None
    assert tree.getNode(9)['rightChild'] == None
    assert tree.getNode(9)['leafDependecies'] == {5:1}
    assert tree.getNode(9)['value'] == 5

    assert tree.hasNode(10) == True
    assert tree.getNodeValue(10) == 5
    assert tree.getNode(10)['nodeID'] == 10
    assert tree.getNode(10)['numChildren'] == 0
    assert tree.getNode(10)['numLeafDependents'] == 1
    assert tree.getNode(10)['leftChild'] == None
    assert tree.getNode(10)['rightChild'] == None
    assert tree.getNode(10)['leafDependecies'] == {5:1}
    assert tree.getNode(10)['value'] == 5

    assert tree.leavesNodeMap[1] == {0,1,2,3,4,5,6,7,8,9,10}
    assert tree.sources == {0,1,2,3,4,5,6,7,8,9,10}

def test_add_non_leaf_node():

    tree = eTree.expressionTree()
    tree.addLeafNode(1, 5)
    tree.addLeafNode(2,  3)

    tree.addNonLeafNode(11, '+', 1, 2)

    assert tree.hasNode(11) == True
    assert tree.getNodeValue(11) == '+'

    assert tree.getNode(11)['nodeID'] == 11
    assert tree.getNode(11)['numChildren'] == 2
    assert tree.getNode(11)['numLeafDependents'] == 2
    assert tree.getNode(11)['leftChild']["nodeID"] == 1
    assert tree.getNode(11)['rightChild']["nodeID"] == 2
    assert tree.getNode(11)['leafDependecies'] == {5:1, 3:1}
    assert tree.getNode(11)['value'] == '+'

    assert tree.leavesNodeMap[1] == {1,2}
    assert tree.leavesNodeMap[2] == {11}
    assert tree.sources == {11}
    
def test_add_multiple_non_leaf_nodes():
    tree = eTree.expressionTree()

    #add leaf nodes
    tree.addLeafNode(0, 5)
    assert tree.getNode(0)['nodeID'] == 0
    assert tree.getNode(0)['numChildren'] == 0 
    assert tree.getNode(0)['numLeafDependents'] == 1
    assert tree.getNode(0)['leftChild'] == None
    assert tree.getNode(0)['rightChild'] == None
    assert tree.getNode(0)['leafDependecies'] == {5:1}
    assert tree.getNode(0)['value'] == 5 

    tree.addLeafNode(1, 5)
    tree.addLeafNode(2,  3)
    tree.addLeafNode(3,  4)
    tree.addLeafNode(4,  2)
    tree.addLeafNode(5,  1)
    tree.addLeafNode(6,  4)
    tree.addLeafNode(7,  2)
    tree.addLeafNode(8,  3)
    tree.addLeafNode(9,  5)
    tree.addLeafNode(10, 5)
    
    #add nonleaf nodes
    tree.addNonLeafNode(11, '+', 1, 2)
    assert tree.hasNode(11) == True
    assert tree.getNodeValue(11) == '+'

    tree.addNonLeafNode(12, '+', 3, 4)
    assert tree.hasNode(12) == True
    assert tree.getNodeValue(12) == '+'

    tree.addNonLeafNode(13, '+', 6, 7)
    assert tree.hasNode(13) == True
    assert tree.getNodeValue(13) == '+'

    tree.addNonLeafNode(14, '+', 9, 10)
    assert tree.hasNode(14) == True
    assert tree.getNodeValue(14) == '+'

    tree.addNonLeafNode(15, '+', 0, 11)
    assert tree.hasNode(15) == True
    assert tree.getNodeValue(15) == '+'

    tree.addNonLeafNode(16, '+', 5, 13)
    assert tree.hasNode(16) == True
    assert tree.getNodeValue(16) == '+'

    tree.addNonLeafNode(17, '+', 8, 14)
    assert tree.hasNode(17) == True
    assert tree.getNodeValue(17) == '+'

    tree.addNonLeafNode(18, '+', 12, 15)
    assert tree.hasNode(18) == True
    assert tree.getNodeValue(18) == '+'

    tree.addNonLeafNode(19, '+', 16, 17)
    assert tree.hasNode(19) == True
    assert tree.getNodeValue(19) == '+'

    tree.addNonLeafNode(20, '+', 18, 19)
    assert tree.hasNode(20) == True
    assert tree.getNodeValue(20) == '+'

    assert tree.evaluateNode(20) == 39
    assert tree.getNode(20)

    #make sure dependecies are correct
    

    assert tree.getNode(1)['nodeID'] == 1
    assert tree.getNode(1)['numChildren'] == 0 
    assert tree.getNode(1)['numLeafDependents'] == 1
    assert tree.getNode(1)['leftChild'] == None
    assert tree.getNode(1)['rightChild'] == None
    assert tree.getNode(1)['leafDependecies'] == {5:1}
    assert tree.getNode(1)['value'] == 5 

    assert tree.getNode(2)['nodeID'] == 2
    assert tree.getNode(2)['numChildren'] == 0 
    assert tree.getNode(2)['numLeafDependents'] == 1
    assert tree.getNode(2)['leftChild'] == None
    assert tree.getNode(2)['rightChild'] == None
    assert tree.getNode(2)['leafDependecies'] == {3:1}
    assert tree.getNode(2)['value'] == 3

    assert tree.getNode(3)['nodeID'] == 3
    assert tree.getNode(3)['numChildren'] == 0 
    assert tree.getNode(3)['numLeafDependents'] == 1
    assert tree.getNode(3)['leftChild'] == None
    assert tree.getNode(3)['rightChild'] == None
    assert tree.getNode(3)['leafDependecies'] == {4:1}
    assert tree.getNode(3)['value'] == 4

    assert tree.getNode(4)['nodeID'] == 4
    assert tree.getNode(4)['numChildren'] == 0 
    assert tree.getNode(4)['numLeafDependents'] == 1
    assert tree.getNode(4)['leftChild'] == None
    assert tree.getNode(4)['rightChild'] == None
    assert tree.getNode(4)['leafDependecies'] == {2:1}
    assert tree.getNode(4)['value'] == 2

    assert tree.getNode(5)['nodeID'] == 5
    assert tree.getNode(5)['numChildren'] == 0 
    assert tree.getNode(5)['numLeafDependents'] == 1
    assert tree.getNode(5)['leftChild'] == None
    assert tree.getNode(5)['rightChild'] == None
    assert tree.getNode(5)['leafDependecies'] == {1:1}
    assert tree.getNode(5)['value'] == 1

    assert tree.getNode(6)['nodeID'] == 6
    assert tree.getNode(6)['numChildren'] == 0 
    assert tree.getNode(6)['numLeafDependents'] == 1
    assert tree.getNode(6)['leftChild'] == None
    assert tree.getNode(6)['rightChild'] == None
    assert tree.getNode(6)['leafDependecies'] == {4:1}
    assert tree.getNode(6)['value'] == 4

    assert tree.getNode(7)['nodeID'] == 7
    assert tree.getNode(7)['numChildren'] == 0 
    assert tree.getNode(7)['numLeafDependents'] == 1
    assert tree.getNode(7)['leftChild'] == None
    assert tree.getNode(7)['rightChild'] == None
    assert tree.getNode(7)['leafDependecies'] == {2:1}
    assert tree.getNode(7)['value'] == 2

    assert tree.getNode(8)['nodeID'] == 8
    assert tree.getNode(8)['numChildren'] == 0 
    assert tree.getNode(8)['numLeafDependents'] == 1
    assert tree.getNode(8)['leftChild'] == None
    assert tree.getNode(8)['rightChild'] == None
    assert tree.getNode(8)['leafDependecies'] == {3:1}
    assert tree.getNode(8)['value'] == 3 

    assert tree.getNode(9)['nodeID'] == 9
    assert tree.getNode(9)['numChildren'] == 0 
    assert tree.getNode(9)['numLeafDependents'] == 1
    assert tree.getNode(9)['leftChild'] == None
    assert tree.getNode(9)['rightChild'] == None
    assert tree.getNode(9)['leafDependecies'] == {5:1}
    assert tree.getNode(9)['value'] == 5  

    assert tree.getNode(10)['nodeID'] == 10
    assert tree.getNode(10)['numChildren'] == 0 
    assert tree.getNode(10)['numLeafDependents'] == 1
    assert tree.getNode(10)['leftChild'] == None
    assert tree.getNode(10)['rightChild'] == None
    assert tree.getNode(10)['leafDependecies'] == {5:1}
    assert tree.getNode(10)['value'] == 5 

    assert tree.getNode(11)['nodeID'] == 11
    assert tree.getNode(11)['numChildren'] == 2
    assert tree.getNode(11)['numLeafDependents'] == 2
    assert tree.getNode(11)['leftChild']["nodeID"] == 1
    assert tree.getNode(11)['rightChild']["nodeID"] == 2
    assert tree.getNode(11)['leafDependecies'] == {5:1, 3:1}
    assert tree.getNode(11)['value'] == '+'

    assert tree.getNode(12)['nodeID'] == 12
    assert tree.getNode(12)['numChildren'] == 2
    assert tree.getNode(12)['numLeafDependents'] == 2
    assert tree.getNode(12)['leftChild']["nodeID"] == 3
    assert tree.getNode(12)['rightChild']["nodeID"] == 4
    assert tree.getNode(12)['leafDependecies'] == {4:1, 2:1}
    assert tree.getNode(12)['value'] == '+' 

    assert tree.getNode(13)['nodeID'] == 13
    assert tree.getNode(13)['numChildren'] == 2
    assert tree.getNode(13)['numLeafDependents'] == 2
    assert tree.getNode(13)['leftChild']["nodeID"] == 6
    assert tree.getNode(13)['rightChild']["nodeID"] == 7
    assert tree.getNode(13)['leafDependecies'] == {2:1, 4:1}
    assert tree.getNode(13)['value'] == '+'

    assert tree.getNode(14)['nodeID'] == 14
    assert tree.getNode(14)['numChildren'] == 2
    assert tree.getNode(14)['numLeafDependents'] == 2
    assert tree.getNode(14)['leftChild']["nodeID"] == 9
    assert tree.getNode(14)['rightChild']["nodeID"] == 10
    assert tree.getNode(14)['leafDependecies'] == {5:2}
    assert tree.getNode(14)['value'] == '+'

    assert tree.getNode(15)['nodeID'] == 15
    assert tree.getNode(15)['numChildren'] == 2
    assert tree.getNode(15)['numLeafDependents'] == 3
    assert tree.getNode(15)['leftChild']["nodeID"] == 0
    assert tree.getNode(15)['rightChild']["nodeID"] == 11
    assert tree.getNode(15)['leafDependecies'] == {5:2, 3:1}
    assert tree.getNode(15)['value'] == '+'

    assert tree.getNode(16)['nodeID'] == 16
    assert tree.getNode(16)['numChildren'] == 2
    assert tree.getNode(16)['numLeafDependents'] == 3
    assert tree.getNode(16)['leftChild']["nodeID"] == 5
    assert tree.getNode(16)['rightChild']["nodeID"] == 13
    assert tree.getNode(16)['leafDependecies'] == {1:1, 4:1, 2:1}
    assert tree.getNode(16)['value'] == '+'

    assert tree.getNode(17)['nodeID'] == 17
    assert tree.getNode(17)['numChildren'] == 2
    assert tree.getNode(17)['numLeafDependents'] == 3
    assert tree.getNode(17)['leftChild']["nodeID"] == 8
    assert tree.getNode(17)['rightChild']["nodeID"] == 14
    assert tree.getNode(17)['leafDependecies'] == {3:1, 5:2}
    assert tree.getNode(17)['value'] == '+'

    assert tree.getNode(18)['nodeID'] == 18
    assert tree.getNode(18)['numChildren'] == 2
    assert tree.getNode(18)['numLeafDependents'] == 5
    assert tree.getNode(18)['leftChild']["nodeID"] == 12
    assert tree.getNode(18)['rightChild']["nodeID"] == 15
    assert tree.getNode(18)['leafDependecies'] == {5:2, 4:1, 2:1, 3:1}
    assert tree.getNode(18)['value'] == '+'

    assert tree.getNode(19)['nodeID'] == 19
    assert tree.getNode(19)['numChildren'] == 2
    assert tree.getNode(19)['numLeafDependents'] == 6
    assert tree.getNode(19)['leftChild']["nodeID"] == 16
    assert tree.getNode(19)['rightChild']["nodeID"] == 17
    assert tree.getNode(19)['leafDependecies'] == {5:2, 1:1, 4:1, 2:1, 3:1}
    assert tree.getNode(19)['value'] == '+'

    assert tree.getNode(20)['nodeID'] == 20
    assert tree.getNode(20)['numChildren'] == 2
    assert tree.getNode(20)['numLeafDependents'] == 11
    assert tree.getNode(20)['leftChild']["nodeID"] == 18
    assert tree.getNode(20)['rightChild']["nodeID"] == 19
    assert tree.getNode(20)['leafDependecies'] == {5:4, 1:1, 4:2, 2:2, 3:2}
    assert tree.getNode(20)['value'] == '+'  

    assert tree.leavesNodeMap[1] == {0,1,2,3,4,5,6,7,8,9,10}
    assert tree.leavesNodeMap[2] == {11,12,13,14}
    assert tree.leavesNodeMap[3] == {15,16,17}
    assert tree.leavesNodeMap[5] == {18}
    assert tree.leavesNodeMap[6] == {19}
    assert tree.leavesNodeMap[11] == {20}
    assert tree.sources == {20}  

def test_add_multiple_nodes():

    tree = createTree2()

    assert tree.getNode(30)['nodeID'] == 30
    assert tree.getNode(30)['numChildren'] == 0 
    assert tree.getNode(30)['numLeafDependents'] == 1
    assert tree.getNode(30)['leftChild'] == None
    assert tree.getNode(30)['rightChild'] == None
    assert tree.getNode(30)['leafDependecies'] == {5:1}
    assert tree.getNode(30)['value'] == 5 

    assert tree.getNode(31)['nodeID'] == 31
    assert tree.getNode(31)['numChildren'] == 0 
    assert tree.getNode(31)['numLeafDependents'] == 1
    assert tree.getNode(31)['leftChild'] == None
    assert tree.getNode(31)['rightChild'] == None
    assert tree.getNode(31)['leafDependecies'] == {5:1}
    assert tree.getNode(31)['value'] == 5    

    assert tree.getNode(32)['nodeID'] == 32
    assert tree.getNode(32)['numChildren'] == 0 
    assert tree.getNode(32)['numLeafDependents'] == 1
    assert tree.getNode(32)['leftChild'] == None
    assert tree.getNode(32)['rightChild'] == None
    assert tree.getNode(32)['leafDependecies'] == {3:1}
    assert tree.getNode(32)['value'] == 3    

    assert tree.getNode(33)['nodeID'] == 33
    assert tree.getNode(33)['numChildren'] == 0 
    assert tree.getNode(33)['numLeafDependents'] == 1
    assert tree.getNode(33)['leftChild'] == None
    assert tree.getNode(33)['rightChild'] == None
    assert tree.getNode(33)['leafDependecies'] == {4:1}
    assert tree.getNode(33)['value'] == 4   

    assert tree.getNode(34)['nodeID'] == 34
    assert tree.getNode(34)['numChildren'] == 0 
    assert tree.getNode(34)['numLeafDependents'] == 1
    assert tree.getNode(34)['leftChild'] == None
    assert tree.getNode(34)['rightChild'] == None
    assert tree.getNode(34)['leafDependecies'] == {2:1}
    assert tree.getNode(34)['value'] == 2    

    assert tree.getNode(35)['nodeID'] == 35
    assert tree.getNode(35)['numChildren'] == 0 
    assert tree.getNode(35)['numLeafDependents'] == 1
    assert tree.getNode(35)['leftChild'] == None
    assert tree.getNode(35)['rightChild'] == None
    assert tree.getNode(35)['leafDependecies'] == {1:1}
    assert tree.getNode(35)['value'] == 1    

    assert tree.getNode(36)['nodeID'] == 36
    assert tree.getNode(36)['numChildren'] == 0 
    assert tree.getNode(36)['numLeafDependents'] == 1
    assert tree.getNode(36)['leftChild'] == None
    assert tree.getNode(36)['rightChild'] == None
    assert tree.getNode(36)['leafDependecies'] == {4:1}
    assert tree.getNode(36)['value'] == 4    

    assert tree.getNode(37)['nodeID'] == 37
    assert tree.getNode(37)['numChildren'] == 0 
    assert tree.getNode(37)['numLeafDependents'] == 1
    assert tree.getNode(37)['leftChild'] == None
    assert tree.getNode(37)['rightChild'] == None
    assert tree.getNode(37)['leafDependecies'] == {2:1}
    assert tree.getNode(37)['value'] == 2    

    assert tree.getNode(38)['nodeID'] == 38
    assert tree.getNode(38)['numChildren'] == 0
    assert tree.getNode(38)['numLeafDependents'] == 1
    assert tree.getNode(38)['leftChild'] == None
    assert tree.getNode(38)['rightChild'] == None
    assert tree.getNode(38)['leafDependecies'] == {21:1}
    assert tree.getNode(38)['value'] == 21

    assert tree.getNode(39)['nodeID'] == 39
    assert tree.getNode(39)['numChildren'] == 0
    assert tree.getNode(39)['numLeafDependents'] == 1
    assert tree.getNode(39)['leftChild'] == None
    assert tree.getNode(39)['rightChild'] == None
    assert tree.getNode(39)['leafDependecies'] == {19:1}
    assert tree.getNode(39)['value'] == 19

    assert tree.getNode(41)['nodeID'] == 41
    assert tree.getNode(41)['numChildren'] == 2
    assert tree.getNode(41)['numLeafDependents'] == 2
    assert tree.getNode(41)['leftChild']["nodeID"] == 31
    assert tree.getNode(41)['rightChild']["nodeID"] == 32
    assert tree.getNode(41)['leafDependecies'] == {3:1, 5:1}
    assert tree.getNode(41)['value'] == '+'

    assert tree.getNode(42)['nodeID'] == 42
    assert tree.getNode(42)['numChildren'] == 2
    assert tree.getNode(42)['numLeafDependents'] == 2
    assert tree.getNode(42)['leftChild']["nodeID"] == 33
    assert tree.getNode(42)['rightChild']["nodeID"] == 34
    assert tree.getNode(42)['leafDependecies'] == {4:1, 2:1}
    assert tree.getNode(42)['value'] == '+'

    assert tree.getNode(43)['nodeID'] == 43
    assert tree.getNode(43)['numChildren'] == 2
    assert tree.getNode(43)['numLeafDependents'] == 2
    assert tree.getNode(43)['leftChild']["nodeID"] == 36
    assert tree.getNode(43)['rightChild']["nodeID"] == 37
    assert tree.getNode(43)['leafDependecies'] == {4:1, 2:1}
    assert tree.getNode(43)['value'] == '+'

    assert tree.getNode(44)['nodeID'] == 44
    assert tree.getNode(44)['numChildren'] == 2
    assert tree.getNode(44)['numLeafDependents'] == 4
    assert tree.getNode(44)['leftChild']["nodeID"] == 38
    assert tree.getNode(44)['rightChild']["nodeID"] == 46
    assert tree.getNode(44)['leafDependecies'] == {1:1, 2:1, 4:1, 21:1}
    assert tree.getNode(44)['value'] == '+'

    assert tree.getNode(45)['nodeID'] == 45
    assert tree.getNode(45)['numChildren'] == 2
    assert tree.getNode(45)['numLeafDependents'] == 3
    assert tree.getNode(45)['leftChild']["nodeID"] == 30
    assert tree.getNode(45)['rightChild']["nodeID"] == 41
    assert tree.getNode(45)['leafDependecies'] == {3:1, 5:2}
    assert tree.getNode(45)['value'] == '+'

    assert tree.getNode(46)['nodeID'] == 46
    assert tree.getNode(46)['numChildren'] == 2
    assert tree.getNode(46)['numLeafDependents'] == 3
    assert tree.getNode(46)['leftChild']["nodeID"] == 35
    assert tree.getNode(46)['rightChild']["nodeID"] == 43
    assert tree.getNode(46)['leafDependecies'] == {1:1, 2:1, 4:1}
    assert tree.getNode(46)['value'] == '+'

    assert tree.getNode(47)['nodeID'] == 47
    assert tree.getNode(47)['numChildren'] == 2
    assert tree.getNode(47)['numLeafDependents'] == 5
    assert tree.getNode(47)['leftChild']["nodeID"] == 39
    assert tree.getNode(47)['rightChild']["nodeID"] == 44
    assert tree.getNode(47)['leafDependecies'] == {1:1, 2:1, 4:1, 19:1, 21:1}
    assert tree.getNode(47)['value'] == '+'

    assert tree.getNode(48)['nodeID'] == 48
    assert tree.getNode(48)['numChildren'] == 2
    assert tree.getNode(48)['numLeafDependents'] == 5
    assert tree.getNode(48)['leftChild']["nodeID"] == 42
    assert tree.getNode(48)['rightChild']["nodeID"] == 45
    assert tree.getNode(48)['leafDependecies'] == {5:2, 3:1, 4:1, 2:1}
    assert tree.getNode(48)['value'] == '+'

    assert tree.getNode(49)['nodeID'] == 49
    assert tree.getNode(49)['numChildren'] == 2
    assert tree.getNode(49)['numLeafDependents'] == 8
    assert tree.getNode(49)['leftChild']["nodeID"] == 46
    assert tree.getNode(49)['rightChild']["nodeID"] == 47
    assert tree.getNode(49)['leafDependecies'] == {19:1, 21:1, 4:2, 2:2, 1:2}
    assert tree.getNode(49)['value'] == '+'

    assert tree.getNode(50)['nodeID'] == 50
    assert tree.getNode(50)['numChildren'] == 2
    assert tree.getNode(50)['numLeafDependents'] == 10
    assert tree.getNode(50)['leftChild']["nodeID"] == 47
    assert tree.getNode(50)['rightChild']["nodeID"] == 48
    assert tree.getNode(50)['leafDependecies'] == {19:1, 21:1, 4:2, 2:2, 1:1, 3:1, 5:2}
    assert tree.getNode(50)['value'] == '+'

    assert tree.leavesNodeMap[1] == {30, 31, 32, 33, 34, 35, 36, 37, 38, 39}
    assert tree.leavesNodeMap[2] == {41, 42, 43, }
    assert tree.leavesNodeMap[3] == {46, 45}
    assert tree.leavesNodeMap[4] == {44}
    assert tree.leavesNodeMap[5] == {48, 47}
    assert tree.leavesNodeMap[8] == {49}
    assert tree.leavesNodeMap[10] == {50}
    assert tree.sources == {49, 50} 
    assert tree.evaluateNode(49) == 54
    assert tree.evaluateNode(50) == 66

def test_uniform_staging():

    tree = createTree1()
    tree.uniformStaging()

    assert tree.evaluateNode(20) == 39
    assert tree.numNodes() == 12

    assert tree.hasNode(0) == True
    assert tree.getNodeValue(0) == 5
    assert tree.getNode(0)['nodeID'] == 0
    assert tree.getNode(0)['numChildren'] == 0
    assert tree.getNode(0)['numLeafDependents'] == 1
    assert tree.getNode(0)['leftChild'] == None
    assert tree.getNode(0)['rightChild'] == None
    assert tree.getNode(0)['leafDependecies'] == {5:1}
    assert tree.getNode(0)['value'] == 5

    assert tree.hasNode(2) == True
    assert tree.getNodeValue(2) == 3
    assert tree.getNode(2)['nodeID'] == 2
    assert tree.getNode(2)['numChildren'] == 0
    assert tree.getNode(2)['numLeafDependents'] == 1
    assert tree.getNode(2)['leftChild'] == None
    assert tree.getNode(2)['rightChild'] == None
    assert tree.getNode(2)['leafDependecies'] == {3:1}
    assert tree.getNode(2)['value'] == 3

    assert tree.hasNode(4) == True
    assert tree.getNodeValue(4) == 2
    assert tree.getNode(4)['nodeID'] == 4
    assert tree.getNode(4)['numChildren'] == 0
    assert tree.getNode(4)['numLeafDependents'] == 1
    assert tree.getNode(4)['leftChild'] == None
    assert tree.getNode(4)['rightChild'] == None
    assert tree.getNode(4)['leafDependecies'] == {2:1}
    assert tree.getNode(4)['value'] == 2

    assert tree.hasNode(3) == True
    assert tree.getNodeValue(3) == 4
    assert tree.getNode(3)['nodeID'] == 3
    assert tree.getNode(3)['numChildren'] == 0
    assert tree.getNode(3)['numLeafDependents'] == 1
    assert tree.getNode(3)['leftChild'] == None
    assert tree.getNode(3)['rightChild'] == None
    assert tree.getNode(3)['leafDependecies'] == {4:1}
    assert tree.getNode(3)['value'] == 4

    assert tree.hasNode(5) == True
    assert tree.getNodeValue(5) == 1
    assert tree.getNode(5)['nodeID'] == 5
    assert tree.getNode(5)['numChildren'] == 0
    assert tree.getNode(5)['numLeafDependents'] == 1
    assert tree.getNode(5)['leftChild'] == None
    assert tree.getNode(5)['rightChild'] == None
    assert tree.getNode(5)['leafDependecies'] == {1:1}
    assert tree.getNode(5)['value'] == 1

    assert tree.hasNode(12) == True
    assert tree.getNodeValue(12) == '+'
    assert tree.getNode(12)['nodeID'] == 12
    assert tree.getNode(12)['numChildren'] == 2
    assert tree.getNode(12)['numLeafDependents'] == 2
    assert tree.getNode(12)['leftChild']['nodeID'] == 3
    assert tree.getNode(12)['rightChild']['nodeID'] == 4
    assert tree.getNode(12)['leafDependecies'] == {4:1, 2:1}
    assert tree.getNode(12)['value'] == '+'

    assert tree.hasNode(14) == True
    assert tree.getNodeValue(14) == '+'
    assert tree.getNode(14)['nodeID'] == 14
    assert tree.getNode(14)['numChildren'] == 2
    assert tree.getNode(14)['numLeafDependents'] == 2
    assert tree.getNode(14)['leftChild']['nodeID'] == 0
    assert tree.getNode(14)['rightChild']['nodeID'] == 0
    assert tree.getNode(14)['leafDependecies'] == {5:2}
    assert tree.getNode(14)['value'] == '+'

    assert tree.hasNode(16) == True
    assert tree.getNodeValue(16) == '+'
    assert tree.getNode(16)['nodeID'] == 16
    assert tree.getNode(16)['numChildren'] == 2
    assert tree.getNode(16)['numLeafDependents'] == 3
    assert tree.getNode(16)['leftChild']['nodeID'] == 5
    assert tree.getNode(16)['rightChild']['nodeID'] == 12
    assert tree.getNode(16)['leafDependecies'] == {1:1, 2:1, 4:1}
    assert tree.getNode(16)['value'] == '+'

    assert tree.hasNode(17) == True
    assert tree.getNodeValue(17) == '+'
    assert tree.getNode(17)['nodeID'] == 17
    assert tree.getNode(17)['numChildren'] == 2
    assert tree.getNode(17)['numLeafDependents'] == 3
    assert tree.getNode(17)['leftChild']['nodeID'] == 2
    assert tree.getNode(17)['rightChild']['nodeID'] == 14
    assert tree.getNode(17)['leafDependecies'] == {5:2, 3:1}
    assert tree.getNode(17)['value'] == '+'

    assert tree.hasNode(18) == True
    assert tree.getNodeValue(18) == '+'
    assert tree.getNode(18)['nodeID'] == 18
    assert tree.getNode(18)['numChildren'] == 2
    assert tree.getNode(18)['numLeafDependents'] == 5
    assert tree.getNode(18)['leftChild']['nodeID'] == 12
    assert tree.getNode(18)['rightChild']['nodeID'] == 17
    assert tree.getNode(18)['leafDependecies'] == {5:2, 3:1, 2:1, 4:1}
    assert tree.getNode(18)['value'] == '+'

    assert tree.hasNode(19) == True
    assert tree.getNodeValue(19) == '+'
    assert tree.getNode(19)['nodeID'] == 19
    assert tree.getNode(19)['numChildren'] == 2
    assert tree.getNode(19)['numLeafDependents'] == 6
    assert tree.getNode(19)['leftChild']['nodeID'] == 16
    assert tree.getNode(19)['rightChild']['nodeID'] == 17
    assert tree.getNode(19)['leafDependecies'] == {5:2, 3:1, 1:1, 4:1, 2:1}
    assert tree.getNode(19)['value'] == '+'

    assert tree.hasNode(20) == True
    assert tree.getNodeValue(20) == '+'
    assert tree.getNode(20)['nodeID'] == 20
    assert tree.getNode(20)['numChildren'] == 2
    assert tree.getNode(20)['numLeafDependents'] == 11
    assert tree.getNode(20)['leftChild']['nodeID'] == 18
    assert tree.getNode(20)['rightChild']['nodeID'] == 19
    assert tree.getNode(20)['leafDependecies'] == {5:4, 3:2, 2:2, 4:2, 1:1}
    assert tree.getNode(20)['value'] == '+'

    assert tree.leavesNodeMap[1] == {0,2,3,4,5}
    assert tree.leavesNodeMap[2] == {12,14}
    assert tree.leavesNodeMap[3] == {16,17}
    assert tree.leavesNodeMap[5] == {18}
    assert tree.leavesNodeMap[6] == {19}
    assert tree.leavesNodeMap[11] == {20}
    assert tree.sources == {20} 

def test_uniform_staging_2():
    
    tree = createTree2()
    tree.uniformStaging()

    assert tree.getNode(30)['nodeID'] == 30
    assert tree.getNode(30)['numChildren'] == 0 
    assert tree.getNode(30)['numLeafDependents'] == 1
    assert tree.getNode(30)['leftChild'] == None
    assert tree.getNode(30)['rightChild'] == None
    assert tree.getNode(30)['leafDependecies'] == {5:1}
    assert tree.getNode(30)['value'] == 5 

    assert tree.getNode(32)['nodeID'] == 32
    assert tree.getNode(32)['numChildren'] == 0 
    assert tree.getNode(32)['numLeafDependents'] == 1
    assert tree.getNode(32)['leftChild'] == None
    assert tree.getNode(32)['rightChild'] == None
    assert tree.getNode(32)['leafDependecies'] == {3:1}
    assert tree.getNode(32)['value'] == 3   

    assert tree.getNode(33)['nodeID'] == 33
    assert tree.getNode(33)['numChildren'] == 0 
    assert tree.getNode(33)['numLeafDependents'] == 1
    assert tree.getNode(33)['leftChild'] == None
    assert tree.getNode(33)['rightChild'] == None
    assert tree.getNode(33)['leafDependecies'] == {4:1}
    assert tree.getNode(33)['value'] == 4  

    assert tree.getNode(34)['nodeID'] == 34
    assert tree.getNode(34)['numChildren'] == 0 
    assert tree.getNode(34)['numLeafDependents'] == 1
    assert tree.getNode(34)['leftChild'] == None
    assert tree.getNode(34)['rightChild'] == None
    assert tree.getNode(34)['leafDependecies'] == {2:1}
    assert tree.getNode(34)['value'] == 2    

    assert tree.getNode(35)['nodeID'] == 35
    assert tree.getNode(35)['numChildren'] == 0 
    assert tree.getNode(35)['numLeafDependents'] == 1
    assert tree.getNode(35)['leftChild'] == None
    assert tree.getNode(35)['rightChild'] == None
    assert tree.getNode(35)['leafDependecies'] == {1:1}
    assert tree.getNode(35)['value'] == 1    

    assert tree.getNode(38)['nodeID'] == 38
    assert tree.getNode(38)['numChildren'] == 0
    assert tree.getNode(38)['numLeafDependents'] == 1
    assert tree.getNode(38)['leftChild'] == None
    assert tree.getNode(38)['rightChild'] == None
    assert tree.getNode(38)['leafDependecies'] == {21:1}
    assert tree.getNode(38)['value'] == 21

    assert tree.getNode(39)['nodeID'] == 39
    assert tree.getNode(39)['numChildren'] == 0
    assert tree.getNode(39)['numLeafDependents'] == 1
    assert tree.getNode(39)['leftChild'] == None
    assert tree.getNode(39)['rightChild'] == None
    assert tree.getNode(39)['leafDependecies'] == {19:1}
    assert tree.getNode(39)['value'] == 19

    assert tree.getNode(41)['nodeID'] == 41
    assert tree.getNode(41)['numChildren'] == 2
    assert tree.getNode(41)['numLeafDependents'] == 2
    assert tree.getNode(41)['leftChild']["nodeID"] == 30
    assert tree.getNode(41)['rightChild']["nodeID"] == 32
    assert tree.getNode(41)['leafDependecies'] == {3:1, 5:1}
    assert tree.getNode(41)['value'] == '+'

    assert tree.getNode(42)['nodeID'] == 42
    assert tree.getNode(42)['numChildren'] == 2
    assert tree.getNode(42)['numLeafDependents'] == 2
    assert tree.getNode(42)['leftChild']["nodeID"] == 33
    assert tree.getNode(42)['rightChild']["nodeID"] == 34
    assert tree.getNode(42)['leafDependecies'] == {4:1, 2:1}
    assert tree.getNode(42)['value'] == '+'

    assert tree.getNode(44)['nodeID'] == 44
    assert tree.getNode(44)['numChildren'] == 2
    assert tree.getNode(44)['numLeafDependents'] == 4
    assert tree.getNode(44)['leftChild']["nodeID"] == 38
    assert tree.getNode(44)['rightChild']["nodeID"] == 46
    assert tree.getNode(44)['leafDependecies'] == {1:1, 2:1, 4:1, 21:1}
    assert tree.getNode(44)['value'] == '+'

    assert tree.getNode(45)['nodeID'] == 45
    assert tree.getNode(45)['numChildren'] == 2
    assert tree.getNode(45)['numLeafDependents'] == 3
    assert tree.getNode(45)['leftChild']["nodeID"] == 30
    assert tree.getNode(45)['rightChild']["nodeID"] == 41
    assert tree.getNode(45)['leafDependecies'] == {3:1, 5:2}
    assert tree.getNode(45)['value'] == '+'

    assert tree.getNode(46)['nodeID'] == 46
    assert tree.getNode(46)['numChildren'] == 2
    assert tree.getNode(46)['numLeafDependents'] == 3
    assert tree.getNode(46)['leftChild']["nodeID"] == 35
    assert tree.getNode(46)['rightChild']["nodeID"] == 42
    assert tree.getNode(46)['leafDependecies'] == {1:1, 2:1, 4:1}
    assert tree.getNode(46)['value'] == '+'

    assert tree.getNode(47)['nodeID'] == 47
    assert tree.getNode(47)['numChildren'] == 2
    assert tree.getNode(47)['numLeafDependents'] == 5
    assert tree.getNode(47)['leftChild']["nodeID"] == 39
    assert tree.getNode(47)['rightChild']["nodeID"] == 44
    assert tree.getNode(47)['leafDependecies'] == {1:1, 2:1, 4:1, 19:1, 21:1}
    assert tree.getNode(47)['value'] == '+'

    assert tree.getNode(48)['nodeID'] == 48
    assert tree.getNode(48)['numChildren'] == 2
    assert tree.getNode(48)['numLeafDependents'] == 5
    assert tree.getNode(48)['leftChild']["nodeID"] == 42
    assert tree.getNode(48)['rightChild']["nodeID"] == 45
    assert tree.getNode(48)['leafDependecies'] == {5:2, 3:1, 4:1, 2:1}
    assert tree.getNode(48)['value'] == '+'

    assert tree.getNode(49)['nodeID'] == 49
    assert tree.getNode(49)['numChildren'] == 2
    assert tree.getNode(49)['numLeafDependents'] == 8
    assert tree.getNode(49)['leftChild']["nodeID"] == 46
    assert tree.getNode(49)['rightChild']["nodeID"] == 47
    assert tree.getNode(49)['leafDependecies'] == {19:1, 21:1, 4:2, 2:2, 1:2}
    assert tree.getNode(49)['value'] == '+'

    assert tree.getNode(50)['nodeID'] == 50
    assert tree.getNode(50)['numChildren'] == 2
    assert tree.getNode(50)['numLeafDependents'] == 10
    assert tree.getNode(50)['leftChild']["nodeID"] == 47
    assert tree.getNode(50)['rightChild']["nodeID"] == 48
    assert tree.getNode(50)['leafDependecies'] == {19:1, 21:1, 4:2, 2:2, 1:1, 3:1, 5:2}
    assert tree.getNode(50)['value'] == '+'

    assert tree.leavesNodeMap[1] == {30, 32, 33, 34, 35, 38, 39}
    assert tree.leavesNodeMap[2] == {41, 42,}
    assert tree.leavesNodeMap[3] == {46, 45}
    assert tree.leavesNodeMap[4] == {44}
    assert tree.leavesNodeMap[5] == {48, 47}
    assert tree.leavesNodeMap[8] == {49}
    assert tree.leavesNodeMap[10] == {50}
    assert tree.sources == {49, 50} 
    assert tree.evaluateNode(49) == 54
    assert tree.evaluateNode(50) == 66

def test_merge():
    tree = createTree1()
    tree.uniformStaging()
    tree2 = createTree2()
    tree2.uniformStaging()
    tree.addExpresionTree(tree2)

    assert tree.hasNode(0) == True
    assert tree.getNodeValue(0) == 5
    assert tree.getNode(0)['nodeID'] == 0
    assert tree.getNode(0)['numChildren'] == 0
    assert tree.getNode(0)['numLeafDependents'] == 1
    assert tree.getNode(0)['leftChild'] == None
    assert tree.getNode(0)['rightChild'] == None
    assert tree.getNode(0)['leafDependecies'] == {5:1}
    assert tree.getNode(0)['value'] == 5

    assert tree.hasNode(2) == True
    assert tree.getNodeValue(2) == 3
    assert tree.getNode(2)['nodeID'] == 2
    assert tree.getNode(2)['numChildren'] == 0
    assert tree.getNode(2)['numLeafDependents'] == 1
    assert tree.getNode(2)['leftChild'] == None
    assert tree.getNode(2)['rightChild'] == None
    assert tree.getNode(2)['leafDependecies'] == {3:1}
    assert tree.getNode(2)['value'] == 3

    assert tree.hasNode(4) == True
    assert tree.getNodeValue(4) == 2
    assert tree.getNode(4)['nodeID'] == 4
    assert tree.getNode(4)['numChildren'] == 0
    assert tree.getNode(4)['numLeafDependents'] == 1
    assert tree.getNode(4)['leftChild'] == None
    assert tree.getNode(4)['rightChild'] == None
    assert tree.getNode(4)['leafDependecies'] == {2:1}
    assert tree.getNode(4)['value'] == 2

    assert tree.hasNode(3) == True
    assert tree.getNodeValue(3) == 4
    assert tree.getNode(3)['nodeID'] == 3
    assert tree.getNode(3)['numChildren'] == 0
    assert tree.getNode(3)['numLeafDependents'] == 1
    assert tree.getNode(3)['leftChild'] == None
    assert tree.getNode(3)['rightChild'] == None
    assert tree.getNode(3)['leafDependecies'] == {4:1}
    assert tree.getNode(3)['value'] == 4

    assert tree.hasNode(5) == True
    assert tree.getNodeValue(5) == 1
    assert tree.getNode(5)['nodeID'] == 5
    assert tree.getNode(5)['numChildren'] == 0
    assert tree.getNode(5)['numLeafDependents'] == 1
    assert tree.getNode(5)['leftChild'] == None
    assert tree.getNode(5)['rightChild'] == None
    assert tree.getNode(5)['leafDependecies'] == {1:1}
    assert tree.getNode(5)['value'] == 1

    assert tree.hasNode(12) == True
    assert tree.getNodeValue(12) == '+'
    assert tree.getNode(12)['nodeID'] == 12
    assert tree.getNode(12)['numChildren'] == 2
    assert tree.getNode(12)['numLeafDependents'] == 2
    assert tree.getNode(12)['leftChild']['nodeID'] == 3
    assert tree.getNode(12)['rightChild']['nodeID'] == 4
    assert tree.getNode(12)['leafDependecies'] == {4:1, 2:1}
    assert tree.getNode(12)['value'] == '+'

    assert tree.hasNode(14) == True
    assert tree.getNodeValue(14) == '+'
    assert tree.getNode(14)['nodeID'] == 14
    assert tree.getNode(14)['numChildren'] == 2
    assert tree.getNode(14)['numLeafDependents'] == 2
    assert tree.getNode(14)['leftChild']['nodeID'] == 0
    assert tree.getNode(14)['rightChild']['nodeID'] == 0
    assert tree.getNode(14)['leafDependecies'] == {5:2}
    assert tree.getNode(14)['value'] == '+'

    assert tree.hasNode(16) == True
    assert tree.getNodeValue(16) == '+'
    assert tree.getNode(16)['nodeID'] == 16
    assert tree.getNode(16)['numChildren'] == 2
    assert tree.getNode(16)['numLeafDependents'] == 3
    assert tree.getNode(16)['leftChild']['nodeID'] == 5
    assert tree.getNode(16)['rightChild']['nodeID'] == 12
    assert tree.getNode(16)['leafDependecies'] == {1:1, 2:1, 4:1}
    assert tree.getNode(16)['value'] == '+'

    assert tree.hasNode(17) == True
    assert tree.getNodeValue(17) == '+'
    assert tree.getNode(17)['nodeID'] == 17
    assert tree.getNode(17)['numChildren'] == 2
    assert tree.getNode(17)['numLeafDependents'] == 3
    assert tree.getNode(17)['leftChild']['nodeID'] == 2
    assert tree.getNode(17)['rightChild']['nodeID'] == 14
    assert tree.getNode(17)['leafDependecies'] == {5:2, 3:1}
    assert tree.getNode(17)['value'] == '+'

    assert tree.hasNode(18) == True
    assert tree.getNodeValue(18) == '+'
    assert tree.getNode(18)['nodeID'] == 18
    assert tree.getNode(18)['numChildren'] == 2
    assert tree.getNode(18)['numLeafDependents'] == 5
    assert tree.getNode(18)['leftChild']['nodeID'] == 12
    assert tree.getNode(18)['rightChild']['nodeID'] == 17
    assert tree.getNode(18)['leafDependecies'] == {5:2, 3:1, 2:1, 4:1}
    assert tree.getNode(18)['value'] == '+'

    assert tree.hasNode(19) == True
    assert tree.getNodeValue(19) == '+'
    assert tree.getNode(19)['nodeID'] == 19
    assert tree.getNode(19)['numChildren'] == 2
    assert tree.getNode(19)['numLeafDependents'] == 6
    assert tree.getNode(19)['leftChild']['nodeID'] == 16
    assert tree.getNode(19)['rightChild']['nodeID'] == 17
    assert tree.getNode(19)['leafDependecies'] == {5:2, 3:1, 1:1, 4:1, 2:1}
    assert tree.getNode(19)['value'] == '+'

    assert tree.hasNode(20) == True
    assert tree.getNodeValue(20) == '+'
    assert tree.getNode(20)['nodeID'] == 20
    assert tree.getNode(20)['numChildren'] == 2
    assert tree.getNode(20)['numLeafDependents'] == 11
    assert tree.getNode(20)['leftChild']['nodeID'] == 18
    assert tree.getNode(20)['rightChild']['nodeID'] == 19
    assert tree.getNode(20)['leafDependecies'] == {5:4, 3:2, 2:2, 4:2, 1:1}
    assert tree.getNode(20)['value'] == '+'

    assert tree.getNode(30)['nodeID'] == 30
    assert tree.getNode(30)['numChildren'] == 0 
    assert tree.getNode(30)['numLeafDependents'] == 1
    assert tree.getNode(30)['leftChild'] == None
    assert tree.getNode(30)['rightChild'] == None
    assert tree.getNode(30)['leafDependecies'] == {5:1}
    assert tree.getNode(30)['value'] == 5 

    assert tree.getNode(32)['nodeID'] == 32
    assert tree.getNode(32)['numChildren'] == 0 
    assert tree.getNode(32)['numLeafDependents'] == 1
    assert tree.getNode(32)['leftChild'] == None
    assert tree.getNode(32)['rightChild'] == None
    assert tree.getNode(32)['leafDependecies'] == {3:1}
    assert tree.getNode(32)['value'] == 3   

    assert tree.getNode(33)['nodeID'] == 33
    assert tree.getNode(33)['numChildren'] == 0 
    assert tree.getNode(33)['numLeafDependents'] == 1
    assert tree.getNode(33)['leftChild'] == None
    assert tree.getNode(33)['rightChild'] == None
    assert tree.getNode(33)['leafDependecies'] == {4:1}
    assert tree.getNode(33)['value'] == 4  

    assert tree.getNode(34)['nodeID'] == 34
    assert tree.getNode(34)['numChildren'] == 0 
    assert tree.getNode(34)['numLeafDependents'] == 1
    assert tree.getNode(34)['leftChild'] == None
    assert tree.getNode(34)['rightChild'] == None
    assert tree.getNode(34)['leafDependecies'] == {2:1}
    assert tree.getNode(34)['value'] == 2    

    assert tree.getNode(35)['nodeID'] == 35
    assert tree.getNode(35)['numChildren'] == 0 
    assert tree.getNode(35)['numLeafDependents'] == 1
    assert tree.getNode(35)['leftChild'] == None
    assert tree.getNode(35)['rightChild'] == None
    assert tree.getNode(35)['leafDependecies'] == {1:1}
    assert tree.getNode(35)['value'] == 1    

    assert tree.getNode(38)['nodeID'] == 38
    assert tree.getNode(38)['numChildren'] == 0
    assert tree.getNode(38)['numLeafDependents'] == 1
    assert tree.getNode(38)['leftChild'] == None
    assert tree.getNode(38)['rightChild'] == None
    assert tree.getNode(38)['leafDependecies'] == {21:1}
    assert tree.getNode(38)['value'] == 21

    assert tree.getNode(39)['nodeID'] == 39
    assert tree.getNode(39)['numChildren'] == 0
    assert tree.getNode(39)['numLeafDependents'] == 1
    assert tree.getNode(39)['leftChild'] == None
    assert tree.getNode(39)['rightChild'] == None
    assert tree.getNode(39)['leafDependecies'] == {19:1}
    assert tree.getNode(39)['value'] == 19

    assert tree.getNode(41)['nodeID'] == 41
    assert tree.getNode(41)['numChildren'] == 2
    assert tree.getNode(41)['numLeafDependents'] == 2
    assert tree.getNode(41)['leftChild']["nodeID"] == 30
    assert tree.getNode(41)['rightChild']["nodeID"] == 32
    assert tree.getNode(41)['leafDependecies'] == {3:1, 5:1}
    assert tree.getNode(41)['value'] == '+'

    assert tree.getNode(42)['nodeID'] == 42
    assert tree.getNode(42)['numChildren'] == 2
    assert tree.getNode(42)['numLeafDependents'] == 2
    assert tree.getNode(42)['leftChild']["nodeID"] == 33
    assert tree.getNode(42)['rightChild']["nodeID"] == 34
    assert tree.getNode(42)['leafDependecies'] == {4:1, 2:1}
    assert tree.getNode(42)['value'] == '+'

    assert tree.getNode(44)['nodeID'] == 44
    assert tree.getNode(44)['numChildren'] == 2
    assert tree.getNode(44)['numLeafDependents'] == 4
    assert tree.getNode(44)['leftChild']["nodeID"] == 38
    assert tree.getNode(44)['rightChild']["nodeID"] == 46
    assert tree.getNode(44)['leafDependecies'] == {1:1, 2:1, 4:1, 21:1}
    assert tree.getNode(44)['value'] == '+'

    assert tree.getNode(45)['nodeID'] == 45
    assert tree.getNode(45)['numChildren'] == 2
    assert tree.getNode(45)['numLeafDependents'] == 3
    assert tree.getNode(45)['leftChild']["nodeID"] == 30
    assert tree.getNode(45)['rightChild']["nodeID"] == 41
    assert tree.getNode(45)['leafDependecies'] == {3:1, 5:2}
    assert tree.getNode(45)['value'] == '+'

    assert tree.getNode(46)['nodeID'] == 46
    assert tree.getNode(46)['numChildren'] == 2
    assert tree.getNode(46)['numLeafDependents'] == 3
    assert tree.getNode(46)['leftChild']["nodeID"] == 35
    assert tree.getNode(46)['rightChild']["nodeID"] == 42
    assert tree.getNode(46)['leafDependecies'] == {1:1, 2:1, 4:1}
    assert tree.getNode(46)['value'] == '+'

    assert tree.getNode(47)['nodeID'] == 47
    assert tree.getNode(47)['numChildren'] == 2
    assert tree.getNode(47)['numLeafDependents'] == 5
    assert tree.getNode(47)['leftChild']["nodeID"] == 39
    assert tree.getNode(47)['rightChild']["nodeID"] == 44
    assert tree.getNode(47)['leafDependecies'] == {1:1, 2:1, 4:1, 19:1, 21:1}
    assert tree.getNode(47)['value'] == '+'

    assert tree.getNode(48)['nodeID'] == 48
    assert tree.getNode(48)['numChildren'] == 2
    assert tree.getNode(48)['numLeafDependents'] == 5
    assert tree.getNode(48)['leftChild']["nodeID"] == 42
    assert tree.getNode(48)['rightChild']["nodeID"] == 45
    assert tree.getNode(48)['leafDependecies'] == {5:2, 3:1, 4:1, 2:1}
    assert tree.getNode(48)['value'] == '+'

    assert tree.getNode(49)['nodeID'] == 49
    assert tree.getNode(49)['numChildren'] == 2
    assert tree.getNode(49)['numLeafDependents'] == 8
    assert tree.getNode(49)['leftChild']["nodeID"] == 46
    assert tree.getNode(49)['rightChild']["nodeID"] == 47
    assert tree.getNode(49)['leafDependecies'] == {19:1, 21:1, 4:2, 2:2, 1:2}
    assert tree.getNode(49)['value'] == '+'

    assert tree.getNode(50)['nodeID'] == 50
    assert tree.getNode(50)['numChildren'] == 2
    assert tree.getNode(50)['numLeafDependents'] == 10
    assert tree.getNode(50)['leftChild']["nodeID"] == 47
    assert tree.getNode(50)['rightChild']["nodeID"] == 48
    assert tree.getNode(50)['leafDependecies'] == {19:1, 21:1, 4:2, 2:2, 1:1, 3:1, 5:2}
    assert tree.getNode(50)['value'] == '+'

    assert tree.leavesNodeMap[1] == {0, 2, 3, 4, 5, 30, 32, 33, 34, 35, 38, 39}
    assert tree.leavesNodeMap[2] == {41, 42, 12, 14}
    assert tree.leavesNodeMap[3] == {46, 45, 16, 17}
    assert tree.leavesNodeMap[4] == {44}
    assert tree.leavesNodeMap[5] == {48, 47,18 }
    assert tree.leavesNodeMap[6] == {19}
    assert tree.leavesNodeMap[8] == {49}
    assert tree.leavesNodeMap[10] == {50}
    assert tree.leavesNodeMap[11] == {20}
    assert tree.sources == {49, 50, 20} 
    assert tree.evaluateNode(20) == 39
    assert tree.evaluateNode(49) == 54
    assert tree.evaluateNode(50) == 66

def test_merge_stage():

    tree = createTree1()
    tree.uniformStaging()
    tree2 = createTree2()
    tree2.uniformStaging()
    tree.addExpresionTree(tree2)
    tree.uniformStaging()

    assert tree.hasNode(0) == True
    assert tree.getNodeValue(0) == 5
    assert tree.getNode(0)['nodeID'] == 0
    assert tree.getNode(0)['numChildren'] == 0
    assert tree.getNode(0)['numLeafDependents'] == 1
    assert tree.getNode(0)['leftChild'] == None
    assert tree.getNode(0)['rightChild'] == None
    assert tree.getNode(0)['leafDependecies'] == {5:1}
    assert tree.getNode(0)['value'] == 5

    assert tree.hasNode(2) == True
    assert tree.getNodeValue(2) == 3
    assert tree.getNode(2)['nodeID'] == 2
    assert tree.getNode(2)['numChildren'] == 0
    assert tree.getNode(2)['numLeafDependents'] == 1
    assert tree.getNode(2)['leftChild'] == None
    assert tree.getNode(2)['rightChild'] == None
    assert tree.getNode(2)['leafDependecies'] == {3:1}
    assert tree.getNode(2)['value'] == 3

    assert tree.hasNode(4) == True
    assert tree.getNodeValue(4) == 2
    assert tree.getNode(4)['nodeID'] == 4
    assert tree.getNode(4)['numChildren'] == 0
    assert tree.getNode(4)['numLeafDependents'] == 1
    assert tree.getNode(4)['leftChild'] == None
    assert tree.getNode(4)['rightChild'] == None
    assert tree.getNode(4)['leafDependecies'] == {2:1}
    assert tree.getNode(4)['value'] == 2

    assert tree.hasNode(3) == True
    assert tree.getNodeValue(3) == 4
    assert tree.getNode(3)['nodeID'] == 3
    assert tree.getNode(3)['numChildren'] == 0
    assert tree.getNode(3)['numLeafDependents'] == 1
    assert tree.getNode(3)['leftChild'] == None
    assert tree.getNode(3)['rightChild'] == None
    assert tree.getNode(3)['leafDependecies'] == {4:1}
    assert tree.getNode(3)['value'] == 4

    assert tree.hasNode(5) == True
    assert tree.getNodeValue(5) == 1
    assert tree.getNode(5)['nodeID'] == 5
    assert tree.getNode(5)['numChildren'] == 0
    assert tree.getNode(5)['numLeafDependents'] == 1
    assert tree.getNode(5)['leftChild'] == None
    assert tree.getNode(5)['rightChild'] == None
    assert tree.getNode(5)['leafDependecies'] == {1:1}
    assert tree.getNode(5)['value'] == 1

    assert tree.hasNode(12) == True
    assert tree.getNodeValue(12) == '+'
    assert tree.getNode(12)['nodeID'] == 12
    assert tree.getNode(12)['numChildren'] == 2
    assert tree.getNode(12)['numLeafDependents'] == 2
    assert tree.getNode(12)['leftChild']['nodeID'] == 3
    assert tree.getNode(12)['rightChild']['nodeID'] == 4
    assert tree.getNode(12)['leafDependecies'] == {4:1, 2:1}
    assert tree.getNode(12)['value'] == '+'

    assert tree.hasNode(14) == True
    assert tree.getNodeValue(14) == '+'
    assert tree.getNode(14)['nodeID'] == 14
    assert tree.getNode(14)['numChildren'] == 2
    assert tree.getNode(14)['numLeafDependents'] == 2
    assert tree.getNode(14)['leftChild']['nodeID'] == 0
    assert tree.getNode(14)['rightChild']['nodeID'] == 0
    assert tree.getNode(14)['leafDependecies'] == {5:2}
    assert tree.getNode(14)['value'] == '+'

    assert tree.hasNode(16) == True
    assert tree.getNodeValue(16) == '+'
    assert tree.getNode(16)['nodeID'] == 16
    assert tree.getNode(16)['numChildren'] == 2
    assert tree.getNode(16)['numLeafDependents'] == 3
    assert tree.getNode(16)['leftChild']['nodeID'] == 5
    assert tree.getNode(16)['rightChild']['nodeID'] == 12
    assert tree.getNode(16)['leafDependecies'] == {1:1, 2:1, 4:1}
    assert tree.getNode(16)['value'] == '+'

    assert tree.hasNode(17) == True
    assert tree.getNodeValue(17) == '+'
    assert tree.getNode(17)['nodeID'] == 17
    assert tree.getNode(17)['numChildren'] == 2
    assert tree.getNode(17)['numLeafDependents'] == 3
    assert tree.getNode(17)['leftChild']['nodeID'] == 2
    assert tree.getNode(17)['rightChild']['nodeID'] == 14
    assert tree.getNode(17)['leafDependecies'] == {5:2, 3:1}
    assert tree.getNode(17)['value'] == '+'

    assert tree.hasNode(19) == True
    assert tree.getNodeValue(19) == '+'
    assert tree.getNode(19)['nodeID'] == 19
    assert tree.getNode(19)['numChildren'] == 2
    assert tree.getNode(19)['numLeafDependents'] == 6
    assert tree.getNode(19)['leftChild']['nodeID'] == 16
    assert tree.getNode(19)['rightChild']['nodeID'] == 17
    assert tree.getNode(19)['leafDependecies'] == {5:2, 3:1, 1:1, 4:1, 2:1}
    assert tree.getNode(19)['value'] == '+'

    assert tree.hasNode(20) == True
    assert tree.getNodeValue(20) == '+'
    assert tree.getNode(20)['nodeID'] == 20
    assert tree.getNode(20)['numChildren'] == 2
    assert tree.getNode(20)['numLeafDependents'] == 11
    assert tree.getNode(20)['leftChild']['nodeID'] == 48
    assert tree.getNode(20)['rightChild']['nodeID'] == 19
    assert tree.getNode(20)['leafDependecies'] == {5:4, 3:2, 2:2, 4:2, 1:1}
    assert tree.getNode(20)['value'] == '+'

    assert tree.getNode(38)['nodeID'] == 38
    assert tree.getNode(38)['numChildren'] == 0
    assert tree.getNode(38)['numLeafDependents'] == 1
    assert tree.getNode(38)['leftChild'] == None
    assert tree.getNode(38)['rightChild'] == None
    assert tree.getNode(38)['leafDependecies'] == {21:1}
    assert tree.getNode(38)['value'] == 21

    assert tree.getNode(39)['nodeID'] == 39
    assert tree.getNode(39)['numChildren'] == 0
    assert tree.getNode(39)['numLeafDependents'] == 1
    assert tree.getNode(39)['leftChild'] == None
    assert tree.getNode(39)['rightChild'] == None
    assert tree.getNode(39)['leafDependecies'] == {19:1}
    assert tree.getNode(39)['value'] == 19

    assert tree.getNode(44)['nodeID'] == 44
    assert tree.getNode(44)['numChildren'] == 2
    assert tree.getNode(44)['numLeafDependents'] == 4
    assert tree.getNode(44)['leftChild']["nodeID"] == 38
    assert tree.getNode(44)['rightChild']["nodeID"] == 16
    assert tree.getNode(44)['leafDependecies'] == {1:1, 2:1, 4:1, 21:1}
    assert tree.getNode(44)['value'] == '+'

    assert tree.getNode(47)['nodeID'] == 47
    assert tree.getNode(47)['numChildren'] == 2
    assert tree.getNode(47)['numLeafDependents'] == 5
    assert tree.getNode(47)['leftChild']["nodeID"] == 39
    assert tree.getNode(47)['rightChild']["nodeID"] == 44
    assert tree.getNode(47)['leafDependecies'] == {1:1, 2:1, 4:1, 19:1, 21:1}
    assert tree.getNode(47)['value'] == '+'

    assert tree.getNode(48)['nodeID'] == 48
    assert tree.getNode(48)['numChildren'] == 2
    assert tree.getNode(48)['numLeafDependents'] == 5
    assert tree.getNode(48)['leftChild']["nodeID"] == 12
    assert tree.getNode(48)['rightChild']["nodeID"] == 17
    assert tree.getNode(48)['leafDependecies'] == {5:2, 3:1, 4:1, 2:1}
    assert tree.getNode(48)['value'] == '+'

    assert tree.getNode(49)['nodeID'] == 49
    assert tree.getNode(49)['numChildren'] == 2
    assert tree.getNode(49)['numLeafDependents'] == 8
    assert tree.getNode(49)['leftChild']["nodeID"] == 16
    assert tree.getNode(49)['rightChild']["nodeID"] == 47
    assert tree.getNode(49)['leafDependecies'] == {19:1, 21:1, 4:2, 2:2, 1:2}
    assert tree.getNode(49)['value'] == '+'

    assert tree.getNode(50)['nodeID'] == 50
    assert tree.getNode(50)['numChildren'] == 2
    assert tree.getNode(50)['numLeafDependents'] == 10
    assert tree.getNode(50)['leftChild']["nodeID"] == 47
    assert tree.getNode(50)['rightChild']["nodeID"] == 48
    assert tree.getNode(50)['leafDependecies'] == {19:1, 21:1, 4:2, 2:2, 1:1, 3:1, 5:2}
    assert tree.getNode(50)['value'] == '+'

    assert tree.leavesNodeMap[1] == {0, 2, 3, 4, 5, 38, 39}
    assert tree.leavesNodeMap[2] == {12, 14}
    assert tree.leavesNodeMap[3] == {16, 17}
    assert tree.leavesNodeMap[4] == {44}
    assert tree.leavesNodeMap[5] == {48, 47}
    assert tree.leavesNodeMap[6] == {19}
    assert tree.leavesNodeMap[8] == {49}
    assert tree.leavesNodeMap[10] == {50}
    assert tree.leavesNodeMap[11] == {20}
    assert tree.sources == {49, 50, 20} 
    assert tree.evaluateNode(20) == 39
    assert tree.evaluateNode(49) == 54
    assert tree.evaluateNode(50) == 66

def createTree1():
    tree = eTree.expressionTree()

    #add leaf nodes
    tree.addLeafNode(0, 5)
    tree.addLeafNode(1, 5)
    tree.addLeafNode(2,  3)
    tree.addLeafNode(3,  4)
    tree.addLeafNode(4,  2)
    tree.addLeafNode(5,  1)
    tree.addLeafNode(6,  4)
    tree.addLeafNode(7,  2)
    tree.addLeafNode(8,  3)
    tree.addLeafNode(9,  5)
    tree.addLeafNode(10, 5)
    
    #add nonleaf nodes
    tree.addNonLeafNode(11, '+', 1, 2)
    tree.addNonLeafNode(12, '+', 3, 4)
    tree.addNonLeafNode(13, '+', 6, 7)
    tree.addNonLeafNode(14, '+', 9, 10)
    tree.addNonLeafNode(15, '+', 0, 11)
    tree.addNonLeafNode(16, '+', 5, 13)
    tree.addNonLeafNode(17, '+', 8, 14)
    tree.addNonLeafNode(18, '+', 12, 15)
    tree.addNonLeafNode(19, '+', 16, 17)
    tree.addNonLeafNode(20, '+', 18, 19)

    return tree

def createTree2():
    tree = eTree.expressionTree()

    #add leaf nodes
    tree.addLeafNode(30, 5)
    tree.addLeafNode(31, 5)
    tree.addLeafNode(32, 3)
    tree.addLeafNode(33, 4)
    tree.addLeafNode(34, 2)
    tree.addLeafNode(35, 1)
    tree.addLeafNode(36, 4)
    tree.addLeafNode(37, 2)
    tree.addLeafNode(38, 21)
    tree.addLeafNode(39, 19)

    #add non leaf nodes
    tree.addNonLeafNode(41, '+', 31, 32)
    tree.addNonLeafNode(42, '+', 33, 34)
    tree.addNonLeafNode(43, '+', 36, 37)
    tree.addNonLeafNode(45, '+', 30, 41)
    tree.addNonLeafNode(46, '+', 35, 43)
    tree.addNonLeafNode(44, '+', 38, 46)
    tree.addNonLeafNode(47, '+', 39, 44)
    tree.addNonLeafNode(48, '+', 42, 45)
    tree.addNonLeafNode(49, '+', 46, 47)
    tree.addNonLeafNode(50, '+', 47, 48)

    return tree