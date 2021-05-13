'''
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use,copy, modify, merge, publish, distribute, sublicense,and/or sell copies
of the Software,and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
'''

##########################################################################
import sympy as sym
import re as regex
import numpy as np
# from import abc import ABC, abstractmethod

# class DType(enum.Enum):
#     SCALAR   = 1
#     VECTOR   = 2
#     MAT      = 3
#     SYM_MAT  = 4
#     ASYM_MAT = 5

# '''
# basic abstract class for data types supported.
# '''
# class DSymbol():
#     @abstractmethod
#     def get_sym_var(self):
#         pass

# '''
# Scalar variable
# '''
# class DScalar(DSymbol):
#     def __init__(self,name):
#         self._name = name
#         self._type = DType.SCALAR

#     def get_sym_var(self):
#         return self._sym

# '''
# Vector variable
# '''
# class DVector(DSymbol):
#     def __init__(self,name,idx_range):
#         self._name = name
#         self._type = DType.VECTOR
#         self._sym  = list()
#         for i in idx_range:
#             nameR = self._name + repr(i)
#             self._sym.append(nameR)

#     def get_sym_var(self):
#         return self._sym

##########################################################################
# variable initialization functions
##########################################################################


def get_name_suffix(name):
    """Used to extract out the name and suffix of a variable so that they can
    include indexing information

    """
    # NOTE: if name has our indexing information, we need to rip it out
    if name.endswith("]"):
        re_match = regex.search(r"\[[^\]]*\]", name)
        suffix = re_match.group(0)
        name = name[0:re_match.start()]
    else:
        suffix = ""
    return name, suffix


def scalar(name):
    """
    Create a scalar variable with the corresponding name.
    The 'name' will be during code generation, so should match the
    variable name used in the C++ code.
    """
    tname = name
    return sym.symbols(tname)


def vec(name, n):
    """
    Create a nD vector variable with the corresponding name.
    The 'name' will be during code generation, so should match
    the variable name used in the C++ code. The returned
    variable can be indexed(0,1,2), i.e.,

    b = dendro.vec("beta")
    b[1] = x^2
    """

    name, suffix = get_name_suffix(name)

    vname = list()

    for i in range(0, n):
        nameR = name + repr(i) + suffix
        vname.append(nameR)

    return sym.Matrix([sym.symbols(vname[i]) for i in range(0, n)])


def vec3(name):
    """Simple function wrapper of vec for 3D

    Simply calls vec and passes through dimensionality of three
    """

    return vec(name, 3)


def mat(name, n, m):
    """
    Creates a symbolic matrix of size nxm
    """
    vname = list()

    for i in range(0, n):
        nameR = name + repr(i)
        nameC = ' '.join([nameR + repr(j) for j in range(0, m)])
        vname.append(nameC)

    return sym.Matrix([sym.symbols(vname[i]) for i in range(0, n)])


def matI(n):
    """
    Creates a Idendity matrix of size nxn
    """

    return sym.eye(n)


def sym_mat(name, n):
    """
    Creates a symbolic symmetric matrix of size nxn
    """

    name, suffix = get_name_suffix(name)

    vname = list()

    for i in range(0, n):
        nameR = name + repr(i)
        nameC = ' '.join([nameR + repr(j) + suffix for j in range(i, n)])
        vname.append(nameC.split(' '))

    mat_out = sym.zeros(n, n)

    # then iterate through our values to fill the sides
    for i in range(0, n):
        mat_out[i, i:n] = sym.symbols([vname[i]])

    # add the transpose to "diagonalize"
    mat_out = mat_out + mat_out.T
    # convert to numpy for a minute just to make things easier
    mat_out = np.array(mat_out)
    # then divide the diagonal by two
    mat_out[range(n), range(n)] /= 2

    # return the symbolic part of the matrix
    return sym.Matrix(mat_out)


def sym_3x3(name):
    """Creates a 3x3 symbolic symmetric matrix

    Simple wrapper of sym_mat passing through 3 as the `n` parameter
    """
    return sym_mat(name, 3)


def antisym_mat(name, n):
    """
    Creates a symbolic anti-symmetric matrix of size nxn
    """

    vname = list()

    for i in range(0, n):
        nameR = name + repr(i)
        nameC = ' '.join([nameR + repr(j) for j in range(i, n)])
        vname.append(nameC)

    return sym.Matrix()
