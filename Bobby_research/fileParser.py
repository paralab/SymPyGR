import expressionTree as expTree

class expressionLine:
    def __init__(self, inputString, _startingId =0):
        self.id = _startingId
        self.input = inputString.strip().replace(" ","")
        self.tokens = self.tokenize()
        self.name = self.collectName()
        #print('name ' + self.name)
    
    def tokenize(self):

        tokens =[]
        #start index of current token
        sI = 0

        #ending Index of a token
        eI = 0
        while eI < len(self.input):
            c = self.input[eI:eI+1]

            #split on Delimeter: + - * / ( , ) ;
            if self.isDelimiter(c):
                if sI != eI:
                    tokens.append(self.getToken(sI, eI))
                tokens.append(c)
                sI = eI +1

            eI = eI+1
        
        if(sI != len(self.input)):
            tokens.append(self.input[sI:])

        return tokens

    def getToken(self, sI, eI):
        '''
        Finds a returns a token from self.input.
        Performs a substring on input[sI,eI)

        sI - starting Index of the substring, inclusive
        eI - ending Index of the substring, exclusive
        '''
        if sI == eI:
            return ValueError('A token must be at least 1 character')
        return self.input[sI:eI]

    def isDelimiter(self, pD):
        '''
        determines if the passed string is a delimiter.
        returns True or False if pD is a limiter or not.

        pD - potential delimiter

        currently checks *+-/(,);= as delimiters
        '''
        if pD == '*':
            return True
        elif pD == '-':
            return True
        elif pD == '/':
            return True
        elif pD == '+':
            return True
        elif pD == '(':
            return True
        elif pD == ')':
            return True
        elif pD == ',':
            return True
        elif pD == ';':
            return True
        elif pD == '=':
            return True
        else:
            return False

    def createTree(self, eTree = None):
        '''
        creates an instance of an expression tree based on the tokens from the input.
        Note nodes that were created in the expression that do not exist in the input string are denoted with an _. ex _1542

        eTree - the expression tree to add the parsed line.
        '''

        if eTree == None:
            eTree = expTree.expressionTree()

        if len(self.tokens) == 1:
            token = self.tokens[0]
            if self.isOperator(token) or self.isFunction(token):
                raise ValueError("the only token must be a value. token: " + token)
            else:
                leafId = self.createLeaf(eTree,token, self.name)
                return eTree, leafId

        index = 0
        opStack = []
        valStack = []
        #for each token
        while index < len(self.tokens):
            token =self.tokens[index]

            #check if a function
            if self.isFunction(token):                
                #keep corresponding parens and handle the function
                if not self.isOpeningParen(self.tokens[index+1]):
                    raise ValueError('funtion must be followed by an opening paren')
                end = self.findCloseParen(index+1)
                #print('found a function: '  + str(index) +" " +str(end) + " inclusive")
                index = self.handleFunction(index, token, opStack, valStack, eTree)
                index = end

            #check if an openeing paren
            elif self.isOpeningParen(token):
                #find closing parens and call recursively. remove parens afterward                
                #print( 'found an open paren ' + str(index) + ' closing parent at : ' + str(self.findCloseParen(index)) + " inclusive")
                index = self.handleParen(index, token, opStack, valStack, eTree)

            #check if an operator
            elif self.isOperator(token):
                #print('found an operator ' + str(index))
                #determine if the operator is a negation or a subtraction
                if self.isNegation(index):
                    #print('found a negation at index ' + str(index))
                    self.handleNegateToken(token, opStack, valStack, eTree)
                elif self.isAddSub(token):
                    #print('found addition or subtraction ' + str(index))
                    self.handleAddSubToken(token, opStack, valStack, eTree)
                elif self.isMultDiv(token):
                    #print('found multiplication or division ' + str(index))
                    self.handleMultDivToken(token, opStack, valStack)

            #handle a value
            else:
                #print('found a value ' + str(index))
                self.handleValueToken(token, opStack, valStack, eTree, index)
                    
            index = index + 1
        
        #potentially an add or sub is left on the opStack
        if(len(opStack)==1):
                rightVal = valStack.pop()
                leftVal = valStack.pop()
                answer = self.handleOp(leftVal, rightVal, opStack.pop(), eTree, self.name)
                valStack.append(answer)
                                
        if(len(valStack) != 1 or len(opStack) != 0):
                raise ValueError("the stacks are misaligned")
        
        return eTree, valStack.pop()

    def handleFunction(self, index, token, opStack, valStack, eTree):
        '''
        handles a function within the expression. returns the index of the closing paren of the function
        '''

        if token == 'pow':
            return self.handlePow(index, token, opStack, valStack, eTree)
        elif token =='sqrt':
            return self.handleSqrt(index, token, opStack, valStack, eTree)
        else:
            raise ValueError('Do no recognize the function: ' + token)

    def handlePow(self, index, token, opStack, valStack, eTree):
        '''
        handles the pow function. Finds the comma and splits the parameters into two seperate expressions.
        '''

        if token != 'pow':
            raise ValueError('function must be pow: ' + token)
        if index == len(self.tokens) or self.tokens[index+1] != '(':
            raise ValueError('pow function must also contain open paren')
        openParen = index + 1
        closeParen = self.findCloseParen(openParen) 
        comma = self.findComma(openParen)
        if comma == -1:
            raise ValueError('There must be a comma for pow function')

        leftSubExp = self.createSubExpression(openParen +1, comma -1) 
        leftSubExpHandler = expressionLine(leftSubExp, self.id)
        eTree, leftId = leftSubExpHandler.createTree(eTree)
        self.id = leftSubExpHandler.id

        rightSubExp = self.createSubExpression(comma+1, closeParen-1) 
        rightSubExpHandler = expressionLine(rightSubExp, self.id)
        eTree, rightId = rightSubExpHandler.createTree(eTree)
        self.id = rightSubExpHandler.id

        if index ==0 and closeParen == len(self.tokens) -1:
            functionToken = self.handleOp(leftId, rightId, token, eTree, self.name)
            valStack.append(functionToken)
        else:
            functionToken = self.handleOp(leftId, rightId, token, eTree)
            self.handleValueToken(functionToken, opStack, valStack, eTree, closeParen)
        return closeParen

    def handleSqrt(self, index, token, opStack, valStack, eTree):
        
        if token != 'sqrt':
            raise ValueError('function must be sqrt: ' + token)
        if index == len(self.tokens) or self.tokens[index+1] != '(':
            raise ValueError('pow function must also contain open paren')
        openParen = index + 1
        closeParen = self.findCloseParen(openParen) 

        subExpression = self.createSubExpression(openParen + 1, closeParen -1)
        subExpHandler = expressionLine(subExpression, self.id)

        eTree, returnId = subExpHandler.createTree(eTree)
        self.id = subExpHandler.id

        

        if index ==0 and closeParen == len(self.tokens) -1:
            returnId = self.handleOp(returnId, 1.0, 'sqrt', eTree, self.name)
            valStack.append(returnId)
        else:
            returnId = self.handleOp(returnId, 1.0, token, eTree)
            self.handleValueToken(returnId, opStack, valStack, eTree, closeParen)
        return closeParen

    def handleParen(self, index, token, opStack, valStack, eTree):
        '''
        handles a parenthesis within the expression. returns the index of the closing paren
        '''
        endIndex = self.findCloseParen(index)
        subExpression = self.createSubExpression(index + 1, endIndex -1)
        subExpHandler = expressionLine(subExpression, self.id)

        eTree, returnId = subExpHandler.createTree(eTree)
        self.id = subExpHandler.id

        self.handleValueToken(returnId, opStack, valStack, eTree, endIndex)
        return endIndex

    def createSubExpression(self, leftIndex, rightIndex):
        '''
        isolates the subexpression from leftIndex to rightIndex, inclusive.
        The subexpression is returned as a string

        leftIndex = index to start the subexpression.
        rightIndex = index to end the subexpression.

        Note both idicies are inclusive
        '''

        if leftIndex > rightIndex:
            raise ValueError("there is no subspace since right index is too large leftIndex: " + str(leftIndex) + " rightIndex: " + str(rightIndex))

        subExp = ""
        while leftIndex <= rightIndex:
            subExp = subExp + str(self.tokens[leftIndex])
            leftIndex = leftIndex + 1
        return subExp

    def handleNegateToken(self, token, opStack, valStack, eTree):
        
        if token != '-':
            raise ValueError('Negate token must be a negative sign. token: ' + str(token))
        self.handleValueToken('Negate', opStack, valStack, eTree)
        self.handleMultDivToken('*', opStack, valStack)

    def handleValueToken(self, token, opStack, valStack, eTree, index = -1):    
        #remove multiply or divide operator if possble
        if len(opStack) > 0 and self.isMultDiv(opStack[-1]):
            #print('handle Muliplication or division after adding: ' + token)
            if index == len(self.tokens) -1  and len(opStack) == 1:
                token = self.handleOp(valStack.pop(), token, opStack.pop(), eTree, self.name)
            else:
                token = self.handleOp(valStack.pop(), token, opStack.pop(), eTree)
        valStack.append(token)

    def handleAddSubToken(self, token, opStack, valStack, eTree):

        #remove previous addition or subtraction token if possible
        if len(opStack) > 0 and self.isAddSub(opStack[-1]):
            rightVal = valStack.pop()
            leftVal = valStack.pop()
            oldToken = self.handleOp(leftVal, rightVal, opStack.pop(), eTree)
            valStack.append(oldToken)
        opStack.append(token)

    def handleMultDivToken(self, token, opStack, valStack):
        if not self.isMultDiv(token):
            raise ValueError("cannot handle multplity or divide operation on token: " + str(token))
        opStack.append(token)

    def handleOp(self, leftTok, rightTok, opToken, tree, opId = None):
        '''
        Takes an operation and 2 values and adds the operation to the expression tree.
        Method checks if each id already exists or assigns a new one if necessary.

        leftTok - the left hand operand of the expression. If not represented in the tree then it will be added to the tree.
        rightTok - the left hand operand of the expression. If not represented in the tree then it will be added to the tree.
        opToken - the operation to be performed on the tokens.
        tree - the expression tree that operation will be added to.
        opId - an optional parameter to give the operation a specific id
        '''

        if not tree.hasNode(leftTok):
            leftTok = self.createLeaf(tree, leftTok)
        if not tree.hasNode(rightTok):
            rightTok = self.createLeaf(tree, rightTok)
        if opId == None:
            opId = "_" + str(self.getUniqueId())
        else:
            if tree.hasNode(opId):
                raise ValueError("the opId already exists in the tree: " + opId)

        tree.addNonLeafNode(opId, opToken, leftTok, rightTok)
        return opId

    def createLeaf(self, tree, value, id = None):
        '''
        Adds a value to the expression tree to the expression tree.

        value - leaf node value to add to the expression tree
        tree - the expression tree that operation will be added to.
        '''

        if id == None:
            id = "_" + str(self.getUniqueId())
        if not tree.hasNode(id):
            tree.addLeafNode(id, value)
        '''    
        else:
            if tree.getNode(id)['value'] != value:
                raise ValueError('Houston We have a problem id: ' + id +'\n value mismatch value: ' + value + " node[value]: " + tree.getNode(id)['value'])
        '''
        return id

    def findCloseParen(self, index):

        if index > len(self.tokens):
            raise ValueError('index out of bounds for CloseParen')

        token = self.tokens[index]
        if not self.isOpeningParen(token):
            raise ValueError('start of the search must begin with an open Paren: ' + token)
        index = index+1
        token = self.tokens[index]

        #represent how many new opening parens have been found
        openCount = 0
         
        while openCount > 0 or not self.isCloseParen(token):
            if self.isOpeningParen(token):
                openCount = openCount+1
            if self.isCloseParen(token):
                openCount = openCount -1
            index = index+1
            token = self.tokens[index]
        
        if self.tokens[index] != ')':
            raise ValueError("There is a Paren mismatch.")
        return index

    def findComma(self, index):
        '''
        starts at index and looks for a comma token. If no comma can be found then -1 returned.
        Meant to only be used by the pow function.
        '''

        if index > len(self.tokens):
            raise ValueError('index out of bounds for CloseParen')

        nestedPowCount = 0
        while index < len(self.tokens):
            if self.tokens[index] == ',':
                if nestedPowCount ==0:
                    return index
                else:
                    nestedPowCount = nestedPowCount -1
            if self.tokens[index] == 'pow':
                nestedPowCount = nestedPowCount + 1
            index = index + 1
        return -1

    def isFunction(self, token):
        return token == 'sqrt' or token == 'pow'
    
    def isOpeningParen(self, token):
        return token == '('
    
    def isCloseParen(self, token):
        return token == ')'

    def isOperator(self, token):
        return self.isAddSub(token) or self.isMultDiv(token)
    
    def isNegation(self, index):
        if self.tokens[index] != '-':
            return False
        if index == 0:
            return True
        else:
            prevTok = self.tokens[index-1]
            return self.isOperator(prevTok) or self.isOpeningParen(prevTok)

    def isAddSub(self, token):
        return token == '+' or token == "-"

    def isMultDiv(self, token):
        return token == '*' or token == '/'

    def collectName(self):
        if self.tokens[len(self.tokens) - 1] == ';':
            self.tokens.pop()
        if len(self.tokens) >= 3:
            # check for assignment
            if self.tokens[1]=='=':
                name = self.tokens[0]
                if name.startswith('double'):
                        name = name[6:]
                #remove name and = from the tokens
                self.tokens.pop(0)
                self.tokens.pop(0)
                return name            
        if len(self.tokens) == 1:
            return self.tokens[0]
        return '_' + str(self.getUniqueId())
    
    def getUniqueId(self):
        self.id = self.id+1
        return self.id