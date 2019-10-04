import expressionTree as eTree
import disjoint_set as djs
import fileParser as fP
import sympy as sympy
import sympy.parsing.sympy_parser as sympy_parser
import sys

def main1():
        #pow(-pow(chi[pp], BSSN_ETA_POWER[0]) + 1, -BSSN_ETA_POWER[1])
        
        id =0
        global_vars_dictionary ={}
        for line in open('bssn.cpp'):
                if not line.startswith('//') and line != '\n':
                        tree = eTree.expressionTree()
                        parse = fP.expressionLine(line,id)
                        
                        updated_tokens = parse.replace_pp_tokens()
                        global_vars_dictionary = merge_dictionaries(global_vars_dictionary, updated_tokens)
                        #print(parse.tokens)
                        tree, junk = parse.createTree(tree)
                        id = parse.id

                        for source in tree.sources:
                                print(source)
                                line = tree.createCodeOutput(source)

                                generated_code = tree.createCodeOutput(source, doubleStar=True)
                                expression = sympy.parsing.sympy_parser.parse_expr(generated_code, evaluate=False)
                                answer = sympy.cse(expression)
                                #print(ccode(Assignment(' DENDRO_972', answer[1])))
                                code = sympy.printing.ccode(answer[1])
                                code = code[1:len(code) -1]
                                print(code)
                                updated_code = replace_tokens(code, updated_tokens)
                                
                                print(updated_code)
                                #print(answer[0][0][0])

def replace_tokens(string, tokens):
        for key in tokens.keys():
                string = string.replace(key, tokens[key])
        return string

def merge_dictionaries(d1, d2):
        for key in d2.keys():
                if key in d1:
                        if d1[key] != d2[key]:
                                raise ValueError('key mismatch')
                else:
                        d1[key] = d2[key]
        return d1


def main():
        #new main file

        '''
        line =  fP.expressionLine("double DENDRO_972 = BSSN_ETA_R0*sqrt(DENDRO_41*(DENDRO_106*DENDRO_971 - DENDRO_274*DENDRO_55 - DENDRO_278*DENDRO_35 - DENDRO_274*DENDRO_55 - DENDRO_357*DENDRO_971 + 2*DENDRO_274*DENDRO_55))*pow(-pow(chi[pp], BSSN_ETA_POWER[0]) + 1, -BSSN_ETA_POWER[1]);")
        tree, id = line.createTree()
        tree.createGraphPicture('pictures/pow')
        '''

        '''
        tree1 = createTree1()
        tree1.createGraphPicture('original_tree')
        subtrees_1 = tree1.cacheAdaptTrees(4)
        count = 0
        for subtree in subtrees_1:
                subtree.createGraphPicture('subtree' + str(count))
                count = count + 1
        '''
        
        tree = eTree.expressionTree()
        id =0
        global_vars_dictionary ={}
        for line in open('bssn.cpp'):
                if not line.startswith('//') and line != '\n':
                        parse = fP.expressionLine(line,id)
                        updated_tokens = parse.replace_pp_tokens()
                        global_vars_dictionary = merge_dictionaries(global_vars_dictionary, updated_tokens)
                        #print(parse.tokens)
                        tree, junk = parse.createTree(tree)
                        id = parse.id
        #tree.createGraphPicture('pictures/basicAdd')
        
        #cache = int(sys.argv[1])
        #registers = int(sys.argv[2])
        cache = 100
        #registers = 100
        passSize = 1

        #print(cache)

              
        counter = 0
        #while cache <=95:
        subtrees = tree.cacheAdaptTrees(cache)
        

        

        # prints all the similarities between nodes
        '''
        sources = []
        source_tree_map = {}
        for subtree in subtrees:
                for source in subtree.sources:
                        if subtree.getNumLeafDependents(source) >=2 :
                                source_tree_map[source] = subtree
                                sources.append(source)
        '''           
        '''
        print("total sources: " + str(len(sources)))
        for source in sources:
                print(source + " " +str(source_tree_map[source].getNumLeafDependents(source)))
        '''
        '''
        for i in range(len(sources)):
                left = sources[i]
                left_dependecies = source_tree_map[left].getNumLeafDependents(left)
                for j in range(i+1, len(sources)):
                        right = sources[j]
                        right_dependecies = source_tree_map[right].getNumLeafDependents(right)
                        similarity_percent = source_tree_map[left].jaccard_similarity(left, source_tree_map[right], right)
                        similarity_count = source_tree_map[left].jaccard_similiarity_count(left,  source_tree_map[right], right)

                        print(str(left) + " " + str(right) + " subtree similarity analysis")
                        print(str(left) + " dependents " + str(left_dependecies))
                        print(str(right) + " dependents " + str(right_dependecies))
                        print('similarity% ' + str(similarity_percent))
                        print('similarity count ' + str(similarity_count))

                #subtree.createGraphPicture('pictures/subT'+str(counter)) 
                #print(len(subtree.sources))  
                 
                if registers > 1:                                        
                        register_trees = subtree.registerAdaptTrees(registers)
                        for register_subtree in register_trees:
                                for register_source in register_subtree.sources: 
                                        output = output + '            ' + 'double ' + register_source + ' = ' + register_subtree.createCodeOutput(register_source, globals=subtree_sources)+ ';\n'
                        
        '''
        
        

        
        
        
        #to autogenerate code
        
        print("num subtrees " + str(len(subtrees)))
        file = open("staged_codes/" + "cache" + str(cache) + "/Code"+str(cache)+".cpp", "w")
        alloc_file = open("staged_codes/" + "cache" + str(cache) + "/allocate"+str(cache)+".cpp", "w")
        dealloc_file = open("staged_codes/" + "cache" + str(cache) + "/deallocate"+str(cache)+".cpp", "w")
        subtree_sources=set()
        for subtree in subtrees:   
                for source in subtree.sources: 
                        if subtree.getNumLeafDependents(source) >=2 : 
                                
                                output = ''
                                var_end = ''
                                end = ''
                                if not source.endswith('[pp]'):
                                        #subtree_sources.add(source)
                                        #output = 'double '
                                        var_end = '[pp]'     

                                
                                end = '        }\n'
                                end = end + '    }\n'
                                end = end + '}\n'
                                end = end + "bssn::timer::t_rhs.stop();"
                                

                                output = "\n\nbssn::timer::t_rhs.start();"                                     
                                output = output + '\nfor (unsigned int k = 3; k < nz-3; k++) {\n'
                                output = output + '    for (unsigned int j = 3; j < ny-3; j++) {\n'
                                output = output + '        for (unsigned int i = 3; i < nx-3; i++) {\n'
                                output = output + '            double x = pmin[0] + i*hx;\n'
                                output = output + '            double y = pmin[1] + j*hy;\n'
                                output = output + '            double z = pmin[2] + k*hz;\n'
                                output = output + '            double r_coord = sqrt(x*x + y*y + z*z);\n'
                                output = output + '            double eta=ETA_CONST;\n'
                                output = output + '            if (r_coord >= ETA_R0) {\n'
                                output = output + '                eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);\n'
                                output = output + '            }\n'
                                output = output + '            int pp = i + nx*(j + ny*k);\n'                     
                                                              
                                generated_code = tree.createCodeOutput(source, doubleStar=True)
                                expression = sympy.parsing.sympy_parser.parse_expr(generated_code, evaluate=False)
                                answer = sympy.cse(expression)
                                #print(ccode(Assignment(' DENDRO_972', answer[1])))

                                print(len(answer[0]))
                                for temp_var in answer[0]:
                                        name = str(temp_var[0])
                                        expression = temp_var[1]
                                        expression = sympy.printing.ccode(temp_var[1])
                                        #expression = expression[1:len(expression) -1]
                                        expression = replace_tokens(expression, global_vars_dictionary)
                                        output = output + '            double ' + name + ' = ' + str(expression)+ ';\n'
                                print('finish temps')

                                code = sympy.printing.ccode(answer[1])
                                code = code[1:len(code) -1]
                                updated_code = replace_tokens(code, global_vars_dictionary)

                                output = output + '            ' + source + var_end + ' = ' + updated_code + ';\n' 
                                if source.startswith('_') or source.startswith("DENDRO"):
                                        
                                        alloc_str = "double* " + source + "= (double*)malloc(sizeof(double)*n);\n"
                                        alloc_file.write(alloc_str)
                                        
                                        dealloc_str = 'free(' + source + ');\n'
                                        dealloc_file.write(dealloc_str)
                                #generate code
                                counter = counter + 1
                                

                                
                output = output + end
                file.write(output)
                #counter = counter + 1
        print('number stages ' + str(counter))
                #cache = cache +5
                   
        

def calculate_rhs():

        a = 1
        K = 2
        b = [3,4,5]
        l1 = 6
        e_i = [0,1,2]
        ad = 7,8,9

        lie =  sum([b[i] * ad[i] for i in e_i]) 
        value = l1*lie - 2*a*K 

        symbols = {}
        symbols['a'] = a
        symbols['K'] = K
        symbols['l1'] = l1
        symbols['b0'] = b[0]
        symbols['b1'] = b[1]
        symbols['b2'] = b[2]
        symbols['ad(0,a)'] = ad[0]
        symbols['ad(1,a)'] = ad[1]
        symbols['ad(2,a)'] = ad[2]

        return value, symbols

def create_arhs():
        tree = eTree.expressionTree()
        
        #create left side term
        tree.addLeafNode(1,"-2")
        tree.addLeafNode(2,"a")
        tree.addLeafNode(3, "K")
        tree.addNonLeafNode(4,"*",2,3)
        tree.addNonLeafNode(5,"*",1,4)

        #implement tree for dendro.lie
        tree.addLeafNode(6,"b0")
        tree.addLeafNode(7,"b1")
        tree.addLeafNode(8,"b2")

        tree.addLeafNode(9,"ad(0,a)")
        tree.addLeafNode(10,"ad(1,a)")
        tree.addLeafNode(11,"ad(2,a)")

        tree.addNonLeafNode(12,"*",6,9)
        tree.addNonLeafNode(13,"*",7,10)
        tree.addNonLeafNode(14,"*",8,11)
        tree.addNonLeafNode(15,"+",12,13)
        tree.addNonLeafNode(16,"+",15,14)

        #create right side term
        tree.addLeafNode(0,"l1")
        tree.addNonLeafNode(17,"*",0,16)

        #merfe right side and left side
        tree.addNonLeafNode(18,"+", 17, 5)

        return tree

def run_add_tree_example():

    tree1 = createTree1()
    print(tree1.evaluateNode(20))
    tree1.createGraphPicture('kTree1')
    tree1.uniformStaging()
    print(tree1.evaluateNode(20))
    tree1.createGraphPicture('kTree1Stage')
    
    


    tree2 = createTree2()
    print(str(tree2.evaluateNode(49)) + " " + str(tree2.evaluateNode(50)))
    tree2.createGraphPicture('kTree2')
    tree2.uniformStaging()
    print(str(tree2.evaluateNode(49)) + " " + str(tree2.evaluateNode(50)))
    tree2.createGraphPicture('kTree2Stage')
    

    
    tree1.addExpresionTree(tree2)
    tree1.createGraphPicture('kMerge')
    tree1.uniformStaging()
    tree1.createGraphPicture('kMerge')
    

    subtrees = tree1.cacheAdaptTrees(5)

    cache = {}
    count = 0
    for subtree in subtrees:
        subtree.createGraphPicture("kMergeSubtree " + str(count))
        count = count + 1

        for source in subtree.sources:
            sourceValue = subtree.evaluateNode(source,cache)
            cache[source] = sourceValue
    print(str(cache[20]) + " " + str(cache[49]) + " " + str(cache[50]))

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

def create_mult_tree():
    
    tree = eTree.expressionTree()
    tree.addLeafNode(0,2)
    tree.addLeafNode(1,5)
    tree.addLeafNode(2,2)

    
    tree.addLeafNode(3,4)
    tree.addLeafNode(4,5)
    
    tree.addLeafNode(5,2)
    tree.addLeafNode(6,4)
    tree.addLeafNode(7,5)
    tree.addLeafNode(8,2)
    
    tree.addLeafNode(9,5)
    tree.addLeafNode(10,4)
    

    tree.addNonLeafNode(11,"+",0,1)
    tree.addNonLeafNode(12,"*",2,11)
    
    tree.addNonLeafNode(13,"*",4,3)
    tree.addNonLeafNode(14,"+",13,12)
    
    tree.addNonLeafNode(15,"+",6,5)
    tree.addNonLeafNode(16,"+",8,7)
    tree.addNonLeafNode(17,"*",16,15)
    tree.addNonLeafNode(18,"+",17,14)
    
    tree.addNonLeafNode(19,"*",10,9)
    
    tree.addNonLeafNode(20,"*",19,18)
    

    return tree

def disjoint_set():
    dj = djs.disjoint_set()
    dj.add(5)
    dj.add(6)
    dj.add(7)
    dj.add(8)
    dj.add(9)
    dj.add(10)

    dj.union(5,6)
    dj.union(7,8)
    dj.union(9,10)
    dj.union(10,6)

    print(dj.find(10))
    print(dj.find(9))
    print(dj.find(8))
    print(dj.find(7))
    print(dj.find(6))
    print(dj.find(5))


if __name__ == '__main__':
    main()