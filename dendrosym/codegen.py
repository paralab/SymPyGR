"""codegen.py

This file contains the functions that generate the C++ or CUDA code
that can be used by Dendro to run the programmed simulations.
"""

import re as regex
from typing import List, Tuple, Union

import sympy as sym


# def construct_cse(ex: List[Union(List[sym.Expr], sym.Matrix, sym.Expr)], vnames: List[str], idx: str) -> Tuple[list, int]:
def construct_cse(ex, vnames, idx) -> Tuple[list, int]:
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


# def generate_cpu(ex: List[Union(List[sym.Expr], sym.Matrix, sym.Expr)], vnames: List[str], idx: str):
def generate_cpu(ex, vnames, idx):
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


def change_deriv_names(in_str):
    c_str = in_str
    derivs = ['agrad', 'grad', 'kograd']
    for deriv in derivs:
        key = deriv + '\(\d, \w+\[pp\]\)'
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
        key = deriv + '\(\d, \d, \w+\[pp\]\)'
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
