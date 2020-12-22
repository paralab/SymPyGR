"""
@author: Milinda Fernando
@brief: graph utils for optimized code generation. 
"""

import matplotlib.pyplot as plt
import sympy



"""
Computes the nodes and rank them based on their in node degree and out node degree threshold values. 
"""
def sorted_nodes_by_in_degree(expr_g,thresh_in,thresh_out=0):
    G=expr_g._G_
    n_list=list()
    for  n in G.nodes():
        if(G.in_degree(n)>=thresh_in  and G.out_degree(n) >thresh_out ):
            n_list.append((G.in_degree(n),n))
    
    n_list.sort(key=lambda  x : x[0],reverse=True)
    return n_list

"""
expr = (x+y + x*y)**2 + x*y
print(gutils.expr_term_rewrite(expr,x*y,Symbol('DENDRO_0')))

expr: sympy expression 
sub_expr: sub expression to rewrite
replace_expr: replace expr sub_expr will get replaced by replace_expr and expr will get regenerated. 
"""
def expr_term_rewrite(expr,sub_expr,replace_expr):

    def _preorder_search(expr,arg_list,sub_expr,replace_expr):
        #print("befor : ", arg_list)

        for i,arg in enumerate(arg_list):
            if(arg == sub_expr):
                arg_list[i] = replace_expr
            else:
                aa_list = list(arg.args)
                if len(aa_list) > 0:
                    _preorder_search(arg,aa_list,sub_expr,replace_expr)
                    aa = arg.func(* tuple(aa_list))
                    arg_list[i] = aa

        #print("after: ",arg_list)
    
    arg_list = list(expr.args)
    _preorder_search(expr,arg_list,sub_expr,replace_expr)
    return expr.func(*tuple(arg_list))
    

"""
returns True if sub_expr is contained in expr
"""
def is_sub_expr(expr,sub_expr):
    for e in sympy.preorder_traversal(expr):
        if e == sub_expr:
            return True
    
    return False


"""
subexpression elimination based on the high in node value. 
"""
def high_node_cse(expr_g, reuse_thresh, rename_prefix="DENDRO_", dep_thresh=0):
    
    G = expr_g._G_
    expr_dict = dict(expr_g._sympy_expr)
    
    n_list=list()
    for  n in G.nodes():
        if(G.in_degree(n)>=reuse_thresh  and G.out_degree(n) >dep_thresh ):
            n_list.append((G.in_degree(n),n))
    
    n_list.sort(key=lambda  x : x[0],reverse=True)
    cse_list=list()

    for node in n_list:
        expr = sympy.parse_expr(str(node[1]))
        cse_list.append(expr)

    # original list of expressions that we will eliminate
    cse_renamed_list = list(cse_list)
    
    # tuple list containing replacement (j,i) expression_i is a subset in expression_j
    replace_idx=list()
    for i,expr_i in enumerate(cse_list):
        for j,expr_j in enumerate(cse_list):
            if (i < j  and is_sub_expr(expr_j,expr_i)):
                replace_idx.append((j,i))
    
    # elimination of the CSE in the selected sub expressions. 
    for (expr_i,sub_i) in replace_idx:
        cse_renamed_list[expr_i] = expr_term_rewrite(cse_renamed_list[expr_i],cse_renamed_list[sub_i],sympy.Symbol(rename_prefix + str(sub_i)))


    # now replace sub expressions in the actual main expressions. 
    for (k,expr) in expr_dict.items():
        print("expr renaming for ", k)
        for i,sub_expr in enumerate(cse_list):
            if is_sub_expr(expr,sub_expr):
                expr=expr_term_rewrite(expr,sub_expr,sympy.Symbol(rename_prefix + str(i)))
            
        expr_dict[k]=expr


    #for (k,expr) in expr_dict.items():
    #    print("expr name : %s  expression %s " %(k,expr))
    cse_dict=dict()
    for i,expr in enumerate(cse_renamed_list):
        cse_dict[rename_prefix+str(i)]=expr

    return [cse_dict,expr_dict]


"""
Generate C sequential code based on high degree node cse
expr_g : expression graph
replace_dic : variable renaming to perform during code generation. 
custom_rename_func: function to perform andy renames. 
cse_mode (=None): if none will use the in node degree to cse metric
reuse_thresh: reuse threshold 
"""
def generate_cpu_c_code(expr_g, reuse_thresh, fname, replace_dict={}, custom_rename_func=None, user_func={}, cse_mode=None):
    cse_dict  = dict()
    exp_dict  = dict()
    fo = open(fname, "w")
    if cse_mode==None:
        [cse_dict,exp_dict] = high_node_cse(expr_g,reuse_thresh,rename_prefix="DENDRO_",dep_thresh=0)
        print('// Dendro: printing temp variables',file=fo)
        for (l,r) in cse_dict.items():
            print("const double %s = " %l, end="",file=fo)
            c_code=sympy.ccode(r,user_functions=user_func)
            
            for (f,r) in replace_dict.items():
                c_code=c_code.replace(str(f),str(r))
            
            if custom_rename_func!= None:
                c_code = custom_rename_func(c_code)

            print(c_code+";",file=fo)

        print('// Dendro: printing RHS variables\n\n\n',file=fo)
        
        for (l,r) in exp_dict.items():
            c_code=sympy.ccode(r,user_functions=user_func)
            for (f,r) in replace_dict.items():
                c_code=c_code.replace(str(f),str(r))
            
            if custom_rename_func!= None:
                c_code = custom_rename_func(c_code)
            
            print(str(l) + " = " + c_code + ";",file=fo)
            print("",file=fo)
    
    fo.close()



"""
Computes the number of minimum registors needed to evaluate an expression. 
"""
# def min_registers(expr):

#     def __pre_traversal(self,expr,reg_count):
#         if isinstance(expr.func, sympy.core.function.UndefinedFunction):
#             sym_name=str(expr.func)
#             for a in expr.args:
#                 sym_name = sym_name + '_' + str(a)
            
#             reg_count.append(1)

#         else if isinstance(expr.func, sympy.core.)
#         else:
#             node_list.append(expr)
        
#         for arg in expr.args:
#             if isinstance(arg.func, sympy.core.function.UndefinedFunction):
#                 f=arg.func
#                 sym_name=str(f)
#                 for a in arg.args:
#                     sym_name = sym_name + '_' + str(a)
                
#                 node_list.append(sympy.Symbol(sym_name))
#                 edge_list.append((expr,sympy.Symbol(sym_name)))
#             else:
#                 edge_list.append((expr,arg))
#                 self.__pre_traversal(arg,node_list,edge_list)
    
    

    


            