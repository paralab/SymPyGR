"""codegen.py

This file contains the functions that generate the C++ or CUDA code
that can be used by Dendro to run the programmed simulations.
"""

import re as regex
from typing import List, Tuple, Union

import sympy as sym
from dendrosym import nr


def construct_cse(ex: Union[list, sym.Matrix, sym.Expr], vnames: List[str],
                  idx: str) -> Tuple[list, int]:
    """Construct the common sub-expression ellimination tree

    TODO: detailed explanation

    Parameters
    ----------
    ex : list, sympy.Matrix, sympy.Expr
        The expression to parse for the corresponding varibles in the next
        parameter. e.g. (as Sympy expressions) [x + 1, y + x**2]
    vnames : list
        The variable names corresponding to each expression, e.g.
        ['alpha_rhs', 'beta_rhs']
    idx : str
        The string used for indexing. e.g. '[pp]'

    Returns
    -------


    """

    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']

    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i] + repr(j) + idx)
        elif type(e) == sym.Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i] + midx[j] + idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i] + idx)

    # ''.join(random.choice(string.ascii_uppercase) for _ in range(5))
    ee_name = 'DENDRO_'
    ee_syms = sym.numbered_symbols(prefix=ee_name)
    _v = sym.cse(lexp, symbols=ee_syms, optimizations='basic')

    return _v, sym.count_ops(lexp)


def generate_cpu(ex: Union[list, sym.Matrix, sym.Expr], vnames: List[str],
                 idx: str) -> Tuple[list, int]:
    """Generate the CPU C++ code by simplifying the expressions

    TODO: expand the documentation

    Parameters
    ----------
    ex : list, sympy.Matrix, sympy.Expr
        The expression to parse for the corresponding varibles in the next
        parameter. e.g. (as Sympy expressions) [x + 1, y + x**2]
    vnames : list
        The variable names corresponding to each expression, e.g.
        ['alpha_rhs', 'beta_rhs']
    idx : str
        The string used for indexing. e.g. '[pp]'
    """
    # print(ex)

    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']

    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i] + repr(j) + idx)
        elif type(e) == sym.Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i] + midx[j] + idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i] + idx)

    cse = construct_cse(ex, vnames, idx)
    _v = cse[0]

    print("// Dendro: {{{ ")
    print("// Dendro: original ops: %d " % (cse[1]))

    ee_name = 'DENDRO_'
    ee_syms = sym.utilities.numbered_symbols(prefix=ee_name)

    custom_functions = {
        'grad': 'grad',
        'grad2': 'grad2',
        'agrad': 'agrad',
        'kograd': 'kograd'
    }

    rops = 0
    print('// Dendro: printing temp variables')
    for (v1, v2) in _v[0]:
        print('const double ', end='')
        print(
            change_deriv_names(
                sym.ccode(v2, assign_to=v1, user_functions=custom_functions)))
        rops = rops + sym.count_ops(v2)

    print()
    print('// Dendro: printing variables')
    for i, e in enumerate(_v[1]):
        print("//--")
        print(
            change_deriv_names(
                sym.ccode(e,
                          assign_to=lname[i],
                          user_functions=custom_functions)))
        rops = rops + sym.count_ops(e)

    print('// Dendro: reduced ops: %d' % (rops))
    print('// Dendro: }}} ')


def change_deriv_names(in_str: str) -> str:
    """Change derivative names within a string

    This code is used to replace all instances of plain
    derivative names with their advanced counterparts. This is
    a fallback in case users define their derivatives as 'd' or
    'd2' or something similar.

    Parameters
    ----------
    in_str : str
        The input string that has the text to be modified

    Returns
    -------
    str
        The output that has the replaced text
    """

    # TODO: this may need to be changed, since it's replacing things
    c_str = in_str
    derivs = ['agrad', 'grad', 'kograd']
    for deriv in derivs:
        key = deriv + r'\(\d, \w+\[pp\]\)'
        slist = regex.findall(key, c_str)
        for s in slist:
            # print(s)
            w1 = s.split('(')
            w2 = w1[1].split(')')[0].split(',')
            # print(w1[0]+'_'+w2[0].strip()+'_'+w2[1].strip()+';')
            rep = w1[0]
            for v in w2:
                rep = rep + '_' + v.strip()
            # rep=rep+';'
            c_str = c_str.replace(s, rep)

    derivs2 = ['grad2']
    for deriv in derivs2:
        key = deriv + r'\(\d, \d, \w+\[pp\]\)'
        slist = regex.findall(key, c_str)
        for s in slist:
            # print(s)
            w1 = s.split('(')
            w2 = w1[1].split(')')[0].split(',')
            # print(w1[0]+'_'+w2[0].strip()+'_'+w2[1].strip()+';')
            rep = w1[0]
            for v in w2:
                rep = rep + '_' + v.strip()
            # rep=rep+';'
            c_str = c_str.replace(s, rep)
    return c_str


def generate_fpcore(ex, vnames, idx):
    """Gennerate FPCore code

    FPCore is the formate used in FPBench benchmarks. It is a basic
    programming language with conditionals and simple loops. FPBench
    is a good benchmark for floating-point computation to understand
    how effective the code is.
    """

    # TODO: this may need to be updated to match the power of CPU gen

    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']

    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i] + repr(j) + idx)
        elif type(e) == sym.Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i] + midx[j] + idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i] + idx)

    cse = construct_cse(ex, vnames, idx)
    _v = cse[0]

    # print("// Dendro: {{{ ")
    # print("// Dendro: original ops: %d " %(cse[1]))

    ee_name = 'DENDRO_'
    ee_syms = sym.utilities.numbered_symbols(prefix=ee_name)

    custom_functions = {
        'grad': 'grad',
        'grad2': 'grad2',
        'agrad': 'agrad',
        'kograd': 'kograd'
    }
    rops = 0

    # re_symbol=regex.compile(r"Symbol\('[a-z,A-Z,_]+[0-9,\[pp\],\[0-9\]]*'\)")
    re_symbol = regex.compile(r"Symbol\('([a-z,A-Z,0-9,_,\[\]]*)'\)")
    re_integer = regex.compile(r"Integer\(([\-,0-9]+)\)")
    re_float = regex.compile(r"Float\('([\-,0-9]*\.[0-9]*)'\s prec=([0-9]+)\)")
    re_grad = regex.compile(
        r"Function\('([a-z]+[0-9]*)'\)\(Integer\(([0-9]+)\)"
        r",\s*Symbol\('([a-z,A-Z]+[0-9]*\[pp\])'\)\)")

    subs_functions = {
        "Add(": "(+ ",
        "Integer(-1)": "-1 ",
        "Mul(": "(* ",
        "Div(": "(/ ",
        "Pow(": "(pow ",
        "Rational(": "(/ "
    }

    # print('// Dendro: printing temp variables')
    tmp_vars = list()
    for (v1, v2) in _v[0]:
        tmp_vars.append(str(v1))
        sym_sub = dict()
        srep = sym.srepr(v2)
        # print(srep)

        res = re_grad.findall(srep)
        for g in res:
            s = "Function('%s')(Integer(%s), Symbol('%s'))" % (g[0], g[1],
                                                               g[2])
            # print(s)
            ss = "Symbol('%s')" % (g[0] + "_" + g[1] + "_" + g[2])
            srep = srep.replace(s, ss)

        srep = srep.replace(",", " ")
        # print(srep)

        res = re_symbol.findall(srep)
        inp_params = list()
        # print(res)
        for s in res:
            ss = s.replace("[pp]", "")
            for index in range(0, 6):
                ss = ss.replace("[" + str(index) + "]", str(index))
            inp_params.append(ss)
            tmp_vars.append(ss)
            sym_sub["Symbol(\'%s\')" % (s)] = ss

        int_sub = dict()
        res = re_integer.findall(srep)
        for s in res:
            int_sub["Integer(%s)" % (s)] = s

        float_sub = dict()
        res = re_float.findall(srep)

        for s in res:
            float_sub["Float('%s'  prec=%s)" % (s[0], s[1])] = s[0]

        for key, val in sym_sub.items():
            # print("{%s: %s} "%(key,val))
            srep = srep.replace(key, val)

        for key, val in int_sub.items():
            # print("{%s: %s} "%(key,val))
            srep = srep.replace(key, val)

        for key, val in float_sub.items():
            # print("{%s: %s} "%(key,val))
            srep = srep.replace(key, val)

        for key, val in subs_functions.items():
            srep = srep.replace(key, val)

        print("(FPCore (%s)" % (" ".join(inp_params)))
        print("\t%s" % (srep))
        print(")\n")

    # print(tmp_vars)
    tmp_vars.clear()
    tmp_vars = list()
    for i, e in enumerate(_v[1]):
        srep = sym.srepr(e)

        res = re_grad.findall(srep)
        for g in res:
            s = "Function('%s')(Integer(%s), Symbol('%s'))" % (g[0], g[1],
                                                               g[2])
            # print(s)
            ss = "Symbol('%s')" % (g[0] + "_" + g[1] + "_" + g[2])
            srep = srep.replace(s, ss)

        srep = srep.replace(",", " ")

        res = re_symbol.findall(srep)
        inp_params = list()
        # print(res)
        for s in res:
            ss = s.replace("[pp]", "")
            for index in range(0, 6):
                ss = ss.replace("[" + str(index) + "]", str(index))
            inp_params.append(ss)
            tmp_vars.append(ss)
            sym_sub["Symbol(\'%s\')" % (s)] = ss

        int_sub = dict()
        res = re_integer.findall(srep)
        for s in res:
            int_sub["Integer(%s)" % (s)] = s

        float_sub = dict()
        res = re_float.findall(srep)

        for s in res:
            float_sub["Float('%s'  prec=%s)" % (s[0], s[1])] = s[0]

        for key, val in sym_sub.items():
            # print("{%s: %s} "%(key,val))
            srep = srep.replace(key, val)

        for key, val in int_sub.items():
            # print("{%s: %s} "%(key,val))
            srep = srep.replace(key, val)

        for key, val in float_sub.items():
            # print("{%s: %s} "%(key,val))
            srep = srep.replace(key, val)

        for key, val in subs_functions.items():
            srep = srep.replace(key, val)

        tmp_vars = list(set(tmp_vars))
        print("(FPCore (%s)" % (" ".join(tmp_vars)))
        print("\t%s" % (srep))
        print(")")
        # print(")")
        # print(")")
        # print(change_deriv_names(ccode(e, assign_to=lname[i],
        #       user_functions=custom_functions)))


def generate_avx(ex, vnames, idx):
    """Generate AVX C++ code

    AVX is the 'advanced vector extensions' library that you can
    use in C++.

    Notes
    -----
    I'm unsure if this is still useful to the project, but has been
    kept for future reference and backwards compatibility.
    """

    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']

    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i] + repr(j) + idx)
        elif type(e) == sym.Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i] + midx[j] + idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i] + idx)

    cse = construct_cse(ex, vnames, idx)
    _v = cse[0]

    print('// Dendro: {{{ ')
    print("// Dendro: original ops: %d " % (cse[1]))

    ee_name = 'DENDRO_'
    ee_syms = sym.utilities.numbered_symbols(prefix=ee_name)

    print('// Dendro vectorized code: {{{')
    oper = {'mul': 'dmul', 'add': 'dadd', 'load': '*'}
    prevdefvars = set()
    for (v1, v2) in _v[0]:
        vv = sym.utilities.numbered_symbols('v')
        vlist = []
        gen_vector_code(v2, vv, vlist, oper, prevdefvars, idx)
        print('  double ' + repr(v1) + ' = ' + repr(vlist[0]) + ';')
    for i, e in enumerate(_v[1]):
        print("//--")
        vv = sym.utilities.numbered_symbols('v')
        vlist = []
        gen_vector_code(e, vv, vlist, oper, prevdefvars, idx)
        # st = '  ' + repr(lname[i]) + '[idx] = ' + repr(vlist[0]) + ';'
        st = '  ' + repr(lname[i]) + " = " + repr(vlist[0]) + ';'
        print(st.replace("'", ""))

    print('// Dendro vectorized code: }}} ')


def generate_separate(ex, vnames, idx, prefix=""):
    """Generate 'separate' C++ code after simplification

    Note
    ----
    I'm not sure what this 'separate' means, but it also includes
    hard coded references to the `bssn` namespace, so this may
    need to be modified in the future.
    """
    # print(ex)
    if len(ex) != 1:
        print('pass each variable separately ', end='\n')
        return

    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']

    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                lexp.append(ev)
                lname.append(vnames[i] + repr(j) + idx)
        elif type(e) == sym.Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                lexp.append(e[k])
                lname.append(vnames[i] + midx[j] + idx)
        else:
            num_e = num_e + 1
            lexp.append(e)
            lname.append(vnames[i] + idx)

    # print(num_e)
    # print(len(lname))
    c_file = open(prefix + vnames[0] + '.cpp', 'w')
    print('generating code for ' + vnames[0])
    print('    bssn::timer::t_rhs.start();', file=c_file)
    print('for (unsigned int k = 3; k < nz-3; k++) { ', file=c_file)
    print('    z = pmin[2] + k*hz;', file=c_file)

    print('for (unsigned int j = 3; j < ny-3; j++) { ', file=c_file)
    print('    y = pmin[1] + j*hy; ', file=c_file)

    print('for (unsigned int i = 3; i < nx-3; i++) {', file=c_file)
    print('    x = pmin[0] + i*hx;', file=c_file)
    print('    pp = i + nx*(j + ny*k);', file=c_file)
    print('    r_coord = sqrt(x*x + y*y + z*z);', file=c_file)
    print('    eta=ETA_CONST;', file=c_file)
    print('    if (r_coord >= ETA_R0) {', file=c_file)
    print('    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);', file=c_file)
    print('    }', file=c_file)

    print('// Dendro: {{{ ', file=c_file)
    print('// Dendro: original ops: ', sym.count_ops(lexp), file=c_file)

    # print("--------------------------------------------------------")
    # print("Now trying Common Subexpression Detection and Collection")
    # print("--------------------------------------------------------")

    # Common Subexpression Detection and Collection
    # for i in range(len(ex)):
    #     # print("--------------------------------------------------------")
    #     # print(ex[i])
    #     # print("--------------------------------------------------------")
    #     ee_name = ''.join(
    #         random.choice(string.ascii_uppercase) for _ in range(5))
    #     ee_syms = sym.utilities.numbered_symbols(prefix=ee_name)
    #     _v = cse(ex[i],symbols=ee_syms)
    #     # print(type(_v))
    #     for (v1,v2) in _v[0]:
    #         print("double %s = %s;" % (v1, v2))
    #     print("%s = %s" % (vnames[i], _v[1][0]))

    # mex = Matrix(ex)
    ee_name = 'DENDRO_'
    # (ABOVE) ''.join(random.choice(string.ascii_uppercase) for _ in range(5))
    ee_syms = sym.utilities.numbered_symbols(prefix=ee_name)
    _v = construct_cse(lexp, symbols=ee_syms, optimizations='basic')

    custom_functions = {
        'grad': 'grad',
        'grad2': 'grad2',
        'agrad': 'agrad',
        'kograd': 'kograd'
    }

    rops = 0
    print('// Dendro: printing temp variables', file=c_file)
    for (v1, v2) in _v[0]:
        # print("double %s = %s;" % (v1, v2)) # replace_pow(v2)))
        print('double ', end='', file=c_file)
        print(change_deriv_names(
            sym.ccode(v2, assign_to=v1, user_functions=custom_functions)),
              file=c_file)
        rops = rops + sym.count_ops(v2)

    print('// Dendro: printing variables', file=c_file)
    for i, e in enumerate(_v[1]):
        print("//--", file=c_file)
        # print("%s = %s;" % (lname[i], e)) # replace_pow(e)))
        f = open(str(vnames[0]) + '.gv', 'w')
        print(sym.printing.dot.dotprint(e), file=f)
        f.close()
        print(change_deriv_names(
            sym.ccode(e, assign_to=lname[i], user_functions=custom_functions)),
              file=c_file)
        # c_file.write('\n')
        rops = rops + sym.count_ops(e)

    print('// Dendro: reduced ops: ', rops, file=c_file)
    print('// Dendro: }}} ', file=c_file)

    print('  }', file=c_file)
    print(' }', file=c_file)
    print('}', file=c_file)
    print('     bssn::timer::t_rhs.stop();', file=c_file)
    c_file.close()
    print('generating code for ' + vnames[0] + ' completed')


def replace_pow(exp_in):
    """Convert integer powers to multiplications

    This function finds all instances of integer power in expressions
    and converts them within the expression to multiplication.
    This is done because the C++ implementation of `pow` is much
    slower than pure multiplication of expressions.

    Parameters
    ----------
    exp_in : sympy.Expression
        The full expression that should be analyzed

    Returns
    -------
    sympy.Expression
        The new output expression that has replaced the "pow" with
        multiplication.
    """
    pows = list(exp_in.atoms(sym.Pow))
    if any(not e.is_Integer for b, e in (i.as_base_exp() for i in pows)):
        raise ValueError("Dendro: Non integer power encountered.")
    repl = zip(pows, (sym.Mul(*[b] * e, evaluate=False)
                      for b, e in (i.as_base_exp() for i in pows)))
    return exp_in.xreplace(dict(repl))


def generate_debug(ex, vnames):
    """Debug version of generating code

    I believe this is depreciated and never used, since generate_cpu
    and other functions have been declared and fleshed out. Kept
    for potential use.
    """
    # print(ex)

    mi = [0, 1, 2, 4, 5, 8]
    midx = ['00', '01', '02', '11', '12', '22']

    # total number of expressions
    # print("--------------------------------------------------------")
    num_e = 0
    lexp = []
    lname = []
    print('// Dendro: {{{ ')
    for i, e in enumerate(ex):
        if type(e) == list:
            num_e = num_e + len(e)
            for j, ev in enumerate(e):
                # lexp.append(ev)
                print(vnames[i] + repr(j), end='')
                print(' = ', end='')
                print(replace_pow(ev), ';')
        elif type(e) == sym.Matrix:
            num_e = num_e + len(e)
            for j, k in enumerate(mi):
                # lexp.append(e[k])
                print(vnames[i] + midx[j], end='')
                print(' = ', end='')
                print(replace_pow(e[k]), ';')
        else:
            num_e = num_e + 1
            # lexp.append(e)
            print(vnames[i], end='')
            print(' = ', end='')
            print(replace_pow(e), ';')

    print('// Dendro: }}} ')


def vec_print_str(tv, pdvars):
    """Generate vector string

    This returns a string that will be used to print a line of code. If the
    variable 'tv' has not yet been used before, then the declaration of this
    variable must be included in the string.

    Parameters
    ----------
    tv : str
        The new temporary variable to print
    pdvars : list
        List of previously declared variables
    """
    st = '  '
    if tv not in pdvars:
        st += 'double '
        pdvars.add(tv)
    return st


def gen_vector_code(ex, vsym, vlist, oper, prevdefvars, idx):
    """Generate vectorized code from an expression

    This function takes the expressions and generates vector
    code to be used in the C++ implementation.

    Parameters
    ----------
    ex : sympy.Function, sympy.Pow, sympy.Expression
        The expression to create the code for. Note that "expression"
        here just means SymPy generated values.
    vsym : list
        Numbered symbols for use
    vlist : list
        An empty list that will be used to process the tree on return
    oper : dict
        Dictionary used to map the '+' and '*' operators
    prevdefvars : dict
        An empty set used to identify previously defined temporary
        variables.
    idx : str
        The name of the index used for accessing arrays. '[pp]'
    """

    one = sym.symbols('one')
    negone = sym.symbols('negone')
    # print (vlist)
    if isinstance(ex, sym.Function):
        # check to see if we are processing a derivative
        if isinstance(ex, nr.ad) or isinstance(ex, nr.d) or isinstance(
                ex, nr.kod) or isinstance(ex, nr.d2s):
            # print('...ex and args: ',ex,ex.func,ex.args)
            tv = next(vsym)
            vlist.append(tv)
            st = vec_print_str(tv, prevdefvars)
            str_args = [repr(a) for a in ex.args]
            o1 = oper['load']
            o1s = repr(o1).replace("'", "")
            idxn = idx.replace("[", "")
            idxn = idxn.replace("]", "")
            st += repr(tv) + ' = ' + o1s + '(' + repr(
                ex.func) + '_' + '_'.join(str_args) + '+' + idxn + ' );'
            # st += repr(tv) + ' = ' + repr(ex) + ';'
            print(st.replace(idx, ""))
            return

    if isinstance(ex, sym.Pow):
        # check to see if we are processing a simple pow
        a1, a2 = ex.args
        # print('processing pow...',ex,a1,a2)
        if isinstance(a1, sym.Symbol) and isinstance(a2, sym.Number):
            # This is a simple Pow function. Process it here and return
            tv = next(vsym)
            vlist.append(tv)
            st = vec_print_str(tv, prevdefvars)
            if (a2 == -1):
                st += repr(tv) + ' = 1.0 / ' + repr(a1) + ';'
            elif (a2 == 2):
                st += repr(tv) + ' = ' + repr(a1) + ' * ' + repr(a1) + ';'
            else:
                st += repr(tv) + ' = pow( ' + repr(a1) + ', ' + repr(a2) + ');'
            print(st)
            return

    # recursively process the arguments of the function or operator
    for arg in ex.args:
        gen_vector_code(arg, vsym, vlist, oper, prevdefvars, idx)

    if isinstance(ex, sym.Number):
        if isinstance(ex, sym.Integer) and ex == 1:
            vlist.append(one)
        elif isinstance(ex, sym.Number) and ex == -1:
            vlist.append(negone)
        else:
            tv = next(vsym)
            vlist.append(tv)
            st = vec_print_str(tv, prevdefvars)
            if isinstance(ex, sym.Rational):
                st += repr(tv) + ' = ' + repr(float(ex)) + ';'
            else:
                st += repr(tv) + ' = ' + repr(ex) + ';'
            print(st)

    elif isinstance(ex, sym.Symbol):
        tv = next(vsym)
        vlist.append(tv)
        st = vec_print_str(tv, prevdefvars)
        st += repr(tv) + ' = ' + repr(ex) + ';'
        print(st)

    elif isinstance(ex, sym.Mul):
        nargs = len(ex.args)
        # print('mul..',len(vlist))
        for i in range(nargs - 1):
            tv = next(vsym)
            st = vec_print_str(tv, prevdefvars)
            st += repr(tv) + ' = '
            v1 = vlist.pop()
            v2 = vlist.pop()
            # st += repr(v1) + ' * ' + repr(v2) + ';'
            o1 = oper['mul']
            st += repr(o1) + '(' + repr(v1) + ', ' + repr(v2) + ');'
            print(st.replace("'", ""))
            vlist.append(tv)

    elif isinstance(ex, sym.Add):
        nargs = len(ex.args)
        # print('add..',len(vlist))
        for i in range(nargs - 1):
            tv = next(vsym)
            st = vec_print_str(tv, prevdefvars)
            st += repr(tv) + ' = '
            v1 = vlist.pop()
            v2 = vlist.pop()
            o1 = oper['add']
            st += repr(o1) + '(' + repr(v1) + ', ' + repr(v2) + ');'
            print(st.replace("'", ""))
            vlist.append(tv)

    elif isinstance(ex, sym.Pow):
        tv = next(vsym)
        qexp = vlist.pop()
        qman = vlist.pop()
        a1, a2 = ex.args
        o1 = oper['mul']
        if isinstance(a2, sym.Integer):
            if (a2 == -1):
                st = vec_print_str(tv, prevdefvars)
                st += repr(tv) + ' =  1.0 / ' + repr(qman) + ';'
            elif (a2 == 2):
                st = vec_print_str(tv, prevdefvars)
                st += repr(tv) + ' = ' + repr(o1) + '(' + repr(
                    qman) + ', ' + repr(qman) + ');'
            elif (a2 == -2):
                v1 = next(vsym)
                st = vec_print_str(v1, prevdefvars)
                st += repr(v1) + ' = ' + repr(o1) + '(' + repr(
                    qman) + ', ' + repr(qman) + ');'
                print(st.replace("'", ""))
                st = vec_print_str(tv, prevdefvars)
                st += repr(tv) + ' = 1.0 / ' + repr(v1) + ';'
            elif (a2 > 2 and a2 < 8):
                v1 = next(vsym)
                st = vec_print_str(v1, prevdefvars)
                st += repr(v1) + ' = ' + repr(o1) + '(' + repr(
                    qman) + ', ' + repr(qman) + ');'
                print(st.replace("'", ""))
                for i in range(a2 - 3):
                    v2 = next(vsym)
                    st = vec_print_str(v2, prevdefvars)
                    st += repr(v2) + ' = ' + repr(o1) + '(' + repr(
                        v1) + ', ' + repr(qman) + ');'
                    print(st.replace("'", ""))
                    v1 = v2
                st = vec_print_str(tv, prevdefvars)
                st += repr(tv) + ' = ' + repr(o1) + '(' + repr(
                    v1) + ', ' + repr(qman) + ');'
            else:
                st = vec_print_str(tv, prevdefvars)
                st += repr(tv) + ' = pow(' + repr(qman) + ',' + repr(
                    qexp) + ');'
        else:
            st = vec_print_str(tv, prevdefvars)
            st = repr(tv) + ' = pow(' + repr(qman) + ',' + repr(qexp) + ');'

        print(st.replace("'", ""))
        vlist.append(tv)
