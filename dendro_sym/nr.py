'''
@author Hari Sundar
@author Eric Hirschmann
@author David Neilsen
@author David Van Komen
@author Milinda Fernando

@note: The code is take from the original dendro.py for the SymPyGR refactoring. 
@brief: all the numerical relativity operators should go here.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use,copy, modify, merge, publish, distribute, sublicense,and/or sell copies
of the Software,and to permit persons to whom the Software is furnished to do so,
subject to the following conditions: 

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.

'''

from sympy import *
from sympy.tensor.array import *
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.utilities import numbered_symbols
from sympy.printing import print_ccode
from sympy.printing.dot import dotprint

undef = symbols('undefined')

metric = undef
inv_metric = undef
C1 = undef
C2 = undef
# C2_spatial
C3 = undef

# first derivative
d = undef
# second derivative
d2s = undef
# advective derivative
ad = undef

# Kreiss-Oliger dissipation operator
kod = undef

one = symbols('one_')
negone = symbols('negone_')

e_i = [0, 1, 2]
e_ij = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]

Ricci = undef


def d2(i, j, a):
    global d2s
    if (i>j):
        return d2s(j,i,a)
    else:
        return d2s(i,j,a)

def set_first_derivative(g):
    """
    Set how the stencil for the first derivative will be called. Here g is a string

    Typically,

    d_i u =  g(i, u)
    """
    global d
    d = Function(g)
    return d


def set_second_derivative(g):
    """
    Set how the stencil for the second derivative will be called. Here g is a string

    Typically,

    d_ij u =  g(i, j, u)
    """
    global d2s
    d2s = Function(g)
    return d2s

def set_advective_derivative(g):
    """
    Set how the stencil for the second derivative will be called. Here g is a string

    Typically,

    ad_i u =  g(i, u)
    """
    global ad
    ad = Function(g)
    return ad

def set_kreiss_oliger_dissipation(g):
    """
    Set how the stencil for Kreiss-Oliger dissipation will be called. Here g is a string.

    Typically,

    kod_i u = g(i, u)
    """
    global kod
    kod = Function(g)
    return kod

# Covariant Derivatives
def DiDj(a):
    """
    Defines the covariant derivative for a scalar a with respect to the full metric.
    [ewh] Actually this defines two covariant derivatives acting on a scalar.
    The derivative in this case is built from the full (non-conformal) metric.
    Thus C3 is built from the full metric.  This object is symmetric in both
    indices.
    """
    global d, C3

    m = Matrix([d2(i, j, a) - sum([C3[l, i, j] * d(l, a) for l in e_i]) for i, j in e_ij])
    return m.reshape(3, 3)


def _Di_Dj(a):
    """
    Defines the covariant derivative.
    [ewh]  Actually, this defines two covariant derivatives acting on a scalar.
    The use of C2 below, however, suggests that this derivative is built
    from the conformal metric.  Such an operator and term shows up in the
    definition of the Ricci scalar which, in turn shows up in the trace-free
    term in the At evolution equation.  As with DiDj, this object is symmetric
    in both indices when acting on a scalar.
    """
    #[ewh] shouldn't this be C2 instead of C3, i.e.:
    global d, C2
    #global d, d2, C3

    m = Matrix([d2(i, j, a) - sum([C2[l, i, j] * d(l, a) for l in e_i]) for i, j in e_ij])
    return m.reshape(3, 3)


# Index Raising
def up_up(A):
    """
    raises both the indices of A, i.e., A_{ij} --> A^{ij}
    """
    global inv_metric

    m = Matrix([sum([inv_metric[i, k]*inv_metric[j, l]*A[k, l] for k, l in e_ij]) for i, j in e_ij])
    return m.reshape(3, 3)

# One index rasing
def up_down(A):
    """
    raises one index of A, i.e., A_{ij} --> A^i_j
    """
    global inv_metric

    m = Matrix([sum([inv_metric[i, k]*A[k, j] for k in e_i]) for i, j in e_ij])
    return m.reshape(3, 3)


def lie(b, a, weight=0):
    """
    Computes the Lie derivative of a field, a, along the vector b.  Assumes
    the metric has been set.  An optional weight for the field can be
    specified.

    b must be of type dendro.vec3
    a can be scalar, vec3 or sym_3x3

    Computes L_b(v)
    """
    global d, ad

    # e_ij = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]

    if type(b) != tuple:
        raise ValueError('Dendro: The field wrt which the Lie derivative is calculated needs to be vec3.')

    if type(a) == Symbol:
        return sum([b[i] * ad(i, a) for i in e_i]) + weight*a*sum([d(i, b[i]) for i in e_i])
    elif type(a) == tuple:
        return [sum([b[j] * ad(j, a[i]) - a[j] * d(j, b[i]) + weight*a[i]*d(j, b[j]) for j in e_i]) for i in e_i]
    elif type(a) == Matrix:
        m = Matrix([sum([b[k]*ad(k, a[i, j]) + a[i, k]*d(j, b[k]) + a[k, j]*d(i, b[k]) + weight*a[i, j]*d(k, b[k]) for k in e_i]) for i, j in e_ij])
        return m.reshape(3, 3)
    else:
        raise ValueError('Dendro: Unknown type for input field to compute Lie derivative for.')

def kodiss(a):
    """
    Kreiss-Oliger dissipation operator
    """
    global kod

    if type(a) == Symbol:
        return sum( [ kod(i, a) for i in e_i ] )
    elif type(a) == tuple:
        return [ sum ( [ kod(i, a[j]) for i in e_i ] ) for j in e_i ]
    elif type(a) == Matrix:
        return Matrix( [ sum( [ kod(k, a[i, j]) for k in e_i ] ) for i, j in e_ij ]).reshape(3, 3)
    else:
        raise ValueError('Dendro: Unknown type for input to computer kodiss.')


def laplacian(a, chi):
    """
    Computes the laplacian of a scalar function with respect to the 3D metric
    gamma_ij.  Assumes that the conformally rescaled metric (called gt in
    various places) and the conformal factor (chi) is set.  Note that C3 is
    built from the same 3D metric.  The only place that this laplacian is
    used in the bssn equations is in the evolution equation for K and is
    the laplacian of alpha (the lapse).
    """
    global d, metric, C3

    full_metric = metric/chi
    inv_full_metric = simplify(full_metric.inv('ADJ'))

    #return sum([(inv_full_metric[i, j] * d2(i, j, a) - sum([C3[l, i, j] * d(l, a) for l in e_i])) for i, j in e_ij])
    return sum([ inv_full_metric[i, j] * ( d2(i, j, a) - sum([C3[l, i, j] * d(l, a) for l in e_i]) ) for i, j in e_ij])


def laplacian_conformal(a):
    """
    Computes the (conformal) laplacian of a scalar function with respect
    to the tilded or conformally rescaled metric (called gt in various
    places).  We assume the rescaled metric is set as well the conformal
    factor, chi.  Note that C2 is built from the conformally rescaled
    metrci.  This (conformal) laplacian is only used in the definition of
    Ricci that shows up in the evolution equation for At (under the trace
    free operation), and even then only in the part that multiplies the
    metric and which will drop out on taking the trace free part.  So, in
    fact, the code could be written to completely ignore this operation
    in the evolution equations themselves.  However, if the constraints
    are included or the full Ricci is needed for another reason, this
    would be needed.
    """
    global d, inv_metric, C2

    if inv_metric == undef:
        inv_metric = get_inverse_metric()

    #ewh3    return sum([(inv_metric[i, j] * d2(i, j, a) - sum([C2[l, i, j] * d(l, a) for l in e_i])) for i, j in e_ij])
    return sum([ inv_metric[i, j] * (d2(i, j, a) - sum([C2[l, i, j] * d(l, a) for l in e_i])) for i, j in e_ij])


def sqr(a):
    """
    Computes the square of the matrix. Assumes metric is set.
    """
    global inv_metric

    if inv_metric == undef:
        inv_metric = get_inverse_metric()

    return sum([a[i, j]*sum([inv_metric[i, k] * inv_metric[j, l] * a[k, l] for k in e_i for l in e_i]) for i, j in e_ij])


def trace_free(x):
    """
    makes the operator trace-free
    """
    global metric, inv_metric

    if inv_metric == undef:
        inv_metric = get_inverse_metric()

    trace = sum([ inv_metric[i, j] * x[i, j] for i, j in e_ij])

    # X_{ab} - 1/3 gt_{ab} X.
    # tf = Matrix([x[i, j] - 1/3*metric[i,j]*trace for i, j in e_ij])
    tf = Matrix([x[i, j] - metric[i,j]*trace/3 for i, j in e_ij])

    return tf.reshape(3, 3)

def vec_j_del_j(b, a):
    """
    expands to  $\beta^i\partial_i \alpha$
    """
    return sum([b[i]*d(i, a) for i in e_i])


#[ewh] Adding this as this term needs to be in the beta equation as an
#      advective derivative ... and not as a regular (partial) derivative.
def vec_j_ad_j(b, f):
    """
    expands to  $\beta^i\partial_i f$
    """
    return sum([b[i]*ad(i, f) for i in e_i])

    #vec_k_del_k = vec_j_del_j

##########################################################################
# metric related functions
##########################################################################


def set_metric(g):
    """
    sets the metric variable, so that dendro knows how to compute the derived variables. This should be done fairly
    early on. e.g.,

    gt = dendro.sym_3x3("gt")
    dendro.set_metric(gt)
    """
    global metric

    metric = g


def get_inverse_metric():
    """
    Computes and returns the inverse metric. The variables need for be defined in advance. e.g.,

    gt = dendro.sym_3x3("gt")
    dendro.set_metric(gt)
    igt = dendro.get_inverse_metric()
    """
    global metric, inv_metric, undef

    if metric == undef:
        raise ValueError('Dendro: Metric not defined.')

    if inv_metric == undef:
        # method : ('GE', 'LU', or 'ADJ')
        inv_metric = simplify(metric.inv('ADJ'))

    return inv_metric


def get_first_christoffel():
    """
    Computes and returns the first Christoffel Symbols. Assumes the metric has been set. e.g.,

    dendro.set_metric(gt);

    C1 = dendro.get_first_christoffel();
    """
    global metric, inv_metric, undef, C1, d

    if inv_metric == undef:
        get_inverse_metric()

    if C1 == undef:
        C1 = MutableDenseNDimArray(range(27), (3, 3, 3))

        for k in e_i:
            for j in e_i:
                for i in e_i:
                    #C1[k, i, j] = 1 / 2 * (d(j, metric[k, i]) + d(i, metric[k, j]) - d(k, metric[i, j]))
                    C1[k, i, j] = 0.5 * (d(j, metric[k, i]) + d(i, metric[k, j]) - d(k, metric[i, j]))

    return C1


def get_second_christoffel():
    """
    Computes and returns the second Christoffel Symbols. Assumes the metric has been set. Will compute the first
    Christoffel if not already computed. e.g.,

    dendro.set_metric(gt);

    C2 = dendro.get_second_christoffel();
    """
    global C2, C1, inv_metric

    if C2 == undef:
        if C1 == undef:
            get_first_christoffel()

        igt_t = Array(inv_metric, (3, 3))
        C2 = tensorcontraction(tensorproduct(igt_t, C1), (1, 2))

    return C2


def get_complete_christoffel(chi):
    """
    Computes and returns the second Christoffel Symbols. Assumes the metric has been set. Will compute the first/second
    Christoffel if not already computed. e.g.,

    dendro.set_metric(gt);

    C2_spatial = dendro.get_complete_christoffel();
    """
    global metric, inv_metric, undef, C1, C2, C3, d

    if C3 == undef:
        C3 = MutableDenseNDimArray(range(27), (3, 3, 3))

        if C2 == undef:
            get_second_christoffel()

        for k in e_i:
            for j in e_i:
                for i in e_i:
                    #C3[i, j, k] = C2[i, j, k] - 1/(2*chi)*(KroneckerDelta(i, j) * d(k, chi) +
                    C3[i, j, k] = C2[i, j, k] - 0.5/(chi)*(KroneckerDelta(i, j) * d(k, chi) +
                                                           KroneckerDelta(i, k) * d(j, chi) -
                                                           metric[j, k]*sum([inv_metric[i, m]*d(m, chi) for m in e_i])
                                                           )

    return C3


def compute_ricci(Gt, chi):
    """
    Computes the Ricci tensor. e.g.,

    dendro.set_metric(gt)

    R = dendro.compute_ricci(Gt, chi)

    or

    dendro.compute_ricci(Gt, chi)

    and use

    dendro.ricci

    The conformal connection coefficient and the conformal variable needs to be supplied.
    """
    global metric, inv_metric, C1, C2

    Lchi = laplacian_conformal(chi)

    #print(type(Lchi))

    #print('Done with Lphi') #simplify(Lchi))


    #ewh4 DKchiDkchi = Matrix([4*metric[i, j]*sum([sum([inv_metric[k, l]*d(l, chi) for l in e_i])*d(k, chi) for k in e_i]) for i, j in e_ij])
    DKchiDkchi = Matrix([0.25/chi/chi*metric[i, j]*sum([sum([inv_metric[k, l]*d(l, chi) for l in e_i])*d(k, chi) for k in e_i]) for i, j in e_ij])

    #print('done with DKchi') # simplify(DKchiDkchi))

    CalGt = [sum(inv_metric[k,l]*C2[i,k,l] for k, l in e_ij) for i in e_i]

    Rt = Matrix([-0.5*sum([inv_metric[l, m]*d2(l, m, metric[i, j]) for l, m in e_ij]) +
              0.5*sum([metric[k,i]*d(j, Gt[k]) + metric[k,j]*d(i, Gt[k]) for k in e_i]) +
              0.5*sum([CalGt[k]*(C1[i,j,k] + C1[j,i,k]) for k in e_i]) +
              sum([inv_metric[l,m]*(C2[k,l,i]*C1[j,k,m] + C2[k,l,j]*C1[i,k,m] + C2[k,i,m]*C1[k,l,j])
                   for k in e_i for l,m in e_ij]) for i,j in e_ij])

    #print('done with Rt') #simplify(Rt))

    #ewh5    Rphi_tmp = Matrix([2*metric[i, j]*Lchi - 4*d(i, chi)*d(j, chi) for i, j in e_ij])
    #dwn    Rphi_tmp = Matrix([ 0.5*metric[i, j]*Lchi/chi - 0.25*d(i, chi)*d(j, chi)/chi/chi for i, j in e_ij])

    #print(simplify(Rphi_tmp))

    #ewh6    Rphi = -2*_Di_Dj(chi) - Rphi_tmp.reshape(3, 3) - DKchiDkchi.reshape(3, 3)
    #dwn    Rphi = -0.5*_Di_Dj(chi)/chi - Rphi_tmp.reshape(3, 3) - DKchiDkchi.reshape(3, 3)
    xRphi = Matrix( [ 1/(2*chi)*(d2(i,j,chi) -
          sum(C2[k,j,i]*d(k,chi) for k in e_i)) -
          1/(4*chi*chi)*d(i,chi)*d(j,chi) for i, j in e_ij]).reshape(3,3)

    Rphi = xRphi + Matrix( [
           1/(2*chi)*metric[i,j] * ( sum(inv_metric[k,l]*(d2(k,l,chi) -
           3/(2*chi)*d(k,chi)*d(l,chi))  for k, l in e_ij) -
           sum(CalGt[m]*d(m,chi) for m in e_i))
           for i, j in e_ij ] ).reshape(3,3)

    return [Rt.reshape(3, 3) + Rphi, Rt.reshape(3,3), Rphi, CalGt]