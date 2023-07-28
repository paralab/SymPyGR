"""
@author Hari Sundar
@author Eric Hirschmann
@author David Neilsen
@author David Van Komen
@author Milinda Fernando

@note: The code is take from the original dendro.py for the
       SymPyGR refactoring.
@brief: all the numerical relativity operators should go here.

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
"""

import sympy as sym
from sympy.functions.special.tensor_functions import KroneckerDelta

undef = sym.symbols("undefined")

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

# OLD VERSIONS OF THE DERIVATIVE CODE
# first derivative
d_old = sym.Function("grad")
# second derivative
d2s_old = sym.Function("grad2")
# advective derivative
ad_old = sym.Function("agrad")

# Kreiss-Oliger dissipation operator
kod = undef

one = sym.symbols("one_")
negone = sym.symbols("negone_")

e_i = [0, 1, 2]
e_ij = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]

Ricci = undef


def d2(i, j, a):
    """Apply the second derivative stencil

    This is more for the C++ implementation, where the
    second derivative stencil is a function that takes a
    first index, a second index, and a variable.

    Typically, this is the same as inputting

    .. math:: d_{ij}u = g(i, j, u)

    Parameters
    ----------
    i : int
        The first index for the derivative
    j : int
        The second index for the derivative
    a : sympy.Symbol
        The variable over which the derivative should be calculated

    Returns
    -------
    sympy.Symbol
        The symbolic representation of the function being applied
    """
    global d2s
    if i > j:
        return d2s(j, i, a)
    else:
        return d2s(i, j, a)


def set_first_derivative(g):
    """Set the stencil function for the first derivative

    This is more for the C++ implementation, where the
    first derivative stencil is a function that takes an
    index, and a variable. This function takes a string which
    is the function name in the C++ implementation that will
    be applied.

    This is the same as inputting

    .. math:: d_{i} = g(i, u)

    This function returns the Sympy Function that can be used
    in your implementation and also sets the DendroSym global
    for the first derivative.

    Parameters
    ----------
    g : str
        The name of the function that will be applied as
        the stencil for the first derivative

    Returns
    -------
    sympy.Function
        The Sympy function that can be used directly later
    """

    global d

    if type(g) is sym.Symbol:
        d = sym.Function(g)
    else:
        d = g

    return d


def set_second_derivative(g):
    """Set the stencil function for the second derivative

    This is more for the C++ implementation, where the
    second derivative stencil is a function that takes an
    index, a second index, and a variable. This function takes
    a string which is the function name in the C++ implementation
    that will be applied.

    This is the same as inputting

    .. math:: d_{ij} u =  g(i, j, u)

    This function returns the Sympy Function that can be used
    in your implementation and also sets the DendroSym global
    for the second derivative.

    Parameters
    ----------
    g : str
        The name of the function that will be applied as
        the stencil for the first derivative.

    Returns
    -------
    sympy.Function
        The Sympy function that can be used directly later.
    """

    global d2s

    if type(g) is sym.Symbol:
        d2s = sym.Function(g)
    else:
        d2s = g

    return d2s


def set_advective_derivative(g):
    """Set the stencil function for the advective derivative

    This is more for the C++ implementation, where the
    advective derivative stencil is a function that takes an
    index and a variable. This function takes
    a string which is the function name in the C++ implementation
    that will be applied.

    This is the same as inputting

    .. math:: d_{i} u =  g(i, u)

    This function returns the Sympy Function that can be used
    in your implementation and also sets the DendroSym global
    for the advective derivative.

    Parameters
    ----------
    g : str
        The name of the function that will be applied as
        the stencil for the first derivative

    Returns
    -------
    sympy.Function
        The Sympy function that can be used directly later
    """

    global ad
    ad = sym.Function(g)
    return ad


def set_kreiss_oliger_dissipation(g):
    """Set the stencil function for the Kriess-Oliger dissipation

    This is more for the C++ implementation, where the
    advective derivative stencil is a function that takes an
    index and a variable. This function takes
    a string which is the function name in the C++ implementation
    that will be applied.

    This is the same as inputting

    .. math:: kod_{i} u =  g(i, u)

    This function returns the Sympy Function that can be used
    in your implementation and also sets the DendroSym global
    for the Kriess-Oliger dissipation.

    Parameters
    ----------
    g : str
        The name of the function that will be applied as
        the stencil for the first derivative

    Returns
    -------
    sympy.Function
        The Sympy function that can be used directly later
    """

    global kod
    kod = sym.Function(g)
    return kod


# Covariant Derivatives


def covariant_divergence(T):
    """Determine the covariant divergence of a 3 vector

    This function will calculate the covariant divergence
    of a given 3 vector. Please note that it requires the
    "metric" to be set via the `dendrosym.gr.set_metric`
    function as the covariant divergence calculation requires
    a metric.

    Parameters
    ----------
    T : sympy.Matrix, tuple
        The 3 dimensional input vector

    Returns
    -------
    sympy.Expression
        The output expression that would compute the
        covariant divergence
    """

    global metric, e_i

    if metric == undef:
        raise ValueError("You need to set the metric first")

    # NOTE: norm already has a square root in it, so I think that covers the
    # formulation as in the link below, will need to further research that
    metric_norm = metric.norm()
    one_over_norm = 1 / metric_norm

    # then from
    # https://ned.ipac.caltech.edu/level5/March01/Carroll3/Carroll3.html
    # equation
    return one_over_norm * sum([d(ii, metric_norm * T[ii]) for ii in e_i])


def DiDj(a):
    """Defines two covariant derivatives acting on a scalar

    Calculates the covariant derivative of a scalar. This is
    an important operation in Numerical Relativity. It requires the
    first derivative and the Christoffel symbols to be calculated
    and stored.

    Parameters
    ----------
    a : sympy.Symbol
        The scalar over which the covariant derivatives should
        be calculated

    Returns
    -------
    sympy.Expression
        A SymPy expression calculated that is the covariant derivative

    Note
    ----
    [ewh] Actually this defines two covariant derivatives acting on a scalar.
    The derivative in this case is built from the full (non-conformal) metric.
    Thus C3 is built from the full metric.  This object is symmetric in both
    indices.
    """

    global d, C3

    if d == undef:
        raise ValueError("First derivative was not defined, please define.")
    if C3 == undef:
        raise ValueError("Third Christoffel symbols were not defined, please define.")

    m = sym.Matrix(
        [
            d2(ii, jj, a) - sum([C3[ll, ii, jj] * d(ll, a) for ll in e_i])
            for ii, jj in e_ij
        ]
    )
    return m.reshape(3, 3)


def _Di_Dj(a):
    """Defines two covariant derivatives acting on a scalar

    Calculates the covariant derivative of a scalar. This is
    an important operation in Numerical Relativity. It requires the
    first derivative and the Christoffel symbols to be calculated
    and stored.

    Parameters
    ----------
    a : sympy.Symbol
        The scalar over which the covariant derivatives should
        be calculated

    Returns
    -------
    sympy.Expression
        A SymPy expression calculated that is the covariant derivative

    Note
    ----
    [ewh] Actually, this defines two covariant derivatives acting on a scalar.
    The use of C2 below, however, suggests that this derivative is built
    from the conformal metric.  Such an operator and term shows up in the
    definition of the Ricci scalar which, in turn shows up in the trace-free
    term in the At evolution equation.  As with DiDj, this object is symmetric
    in both indices when acting on a scalar.
    """

    # [ewh] shouldn't this be C2 instead of C3, i.e.:
    global d, C2

    if d == undef:
        raise ValueError("First derivative was not defined, please define.")

    if C2 == undef:
        raise ValueError("Second Christoffel symbols were not defined, please define.")
    # global d, d2, C3

    m = sym.Matrix(
        [
            d2(ii, jj, a) - sum([C2[ll, ii, jj] * d(ll, a) for ll in e_i])
            for ii, jj in e_ij
        ]
    )

    return m.reshape(3, 3)


def up_up(A):
    """Raises both indices of a matrix A

    This takes a 3x3 matrix and raises both of the indices:

    .. math:: A_{ij} \rightarrow A^{ij}

    Parameters
    ----------
    A : sympy.Matrix
        The input 3x3 matrix to raise the indices of

    Returns
    -------
    sympy.Matrix
        The matrix of expressions to raise both indices of A
    """
    # Index Raising

    global inv_metric

    if inv_metric == undef:
        # if the inverse metric is still undefined, we have to get it
        inv_metric = get_inverse_metric()

    if type(A) is not sym.Matrix:
        raise ValueError("Input variable to UP UP needs to be a matrix.")

    if A.shape != (3, 3):
        raise ValueError("Input matrix into UP UP needs to be 3x3 in shape.")

    m = sym.Matrix(
        [
            sum(
                [inv_metric[ii, kk] * inv_metric[jj, ll] * A[kk, ll] for kk, ll in e_ij]
            )
            for ii, jj in e_ij
        ]
    )

    return m.reshape(3, 3)


def up_down(A):
    """Raises one indices of a matrix A

    This takes a 3x3 matrix and raises one of the indices:

    .. math:: A_{ij} \rightarrow A^i_j

    Parameters
    ----------
    A : sympy.Matrix
        The input 3x3 matrix to transform

    Returns
    -------
    sympy.Matrix
        The matrix of expressions to raise one index of A
    """

    # One index rasing
    global inv_metric

    m = sym.Matrix(
        [sum([inv_metric[ii, kk] * A[kk, jj] for kk in e_i]) for ii, jj in e_ij]
    )

    return m.reshape(3, 3)


def lie(b, a, weight=0, use_advective=False):
    """Compute the Lie derivative of a field along a vector

    Takes a field, a, and computes the Lie derivative along the vector b.
    This function assumes that the metric has already been set.

    An optional weight for the field is also allowed.

    Parameters
    ----------
    b : tuple, sympy.Matrix
        The vector b the derivative will be along
    a : sympy.Symbol, tuple, sympy.Matrix
        The field a to take the derivative of
    weight : number, optional
        A number (could be a sympy rational) for the weighting
        of the Lie derivative
    use_advective : bool, optional
        A flag for enabling advective derivative calculation.
        In some cases, using the advective derivative may be
        beneficial, but after discussion the stability provided
        often is not worth the computational expense. The option
        remains for those that wish to use it.

    Returns
    -------
    sympy.Expression, list, sympy.Matrix
        The expressions that compute the derivative
    """

    global d, ad

    # assign if we're using the advective derivative function or not
    ad_use = ad if use_advective else d

    if type(b) == sym.Matrix:
        # NOTE: a matrix *will* work, but it must be 3, or (3x1)
        # temporarily convert b to a 1D array
        assert b.shape[0] == 3, "The field wrt is too large!"

    elif type(b) != tuple:
        raise ValueError(
            "Dendro: The field wrt which the Lie derivative is calculated"
            " needs to be vec3."
        )

    if type(a) == sym.Symbol:
        return sum([b[ii] * ad_use(ii, a) for ii in e_i]) + weight * a * sum(
            [d(ii, b[ii]) for ii in e_i]
        )

    # if it's a 3-vector, this is what we calculate
    elif type(a) == tuple or (type(a) == sym.Matrix and a.shape == (3, 1)):
        return [
            sum(
                [
                    b[jj] * ad_use(jj, a[ii])
                    - a[jj] * d(jj, b[ii])
                    + weight * a[ii] * d(jj, b[jj])
                    for jj in e_i
                ]
            )
            for ii in e_i
        ]

    elif type(a) == sym.Matrix:
        m = sym.Matrix(
            [
                sum(
                    [
                        b[kk] * ad_use(kk, a[ii, jj])
                        + a[ii, kk] * d(jj, b[kk])
                        + a[kk, jj] * d(ii, b[kk])
                        + weight * a[ii, jj] * d(kk, b[kk])
                        for kk in e_i
                    ]
                )
                for ii, jj in e_ij
            ]
        )
        return m.reshape(3, 3)

    else:
        raise ValueError(
            "Dendro: Unknown type for input field to compute Lie derivative."
        )


def kodiss(a):
    """Apply Kreiss-Oliger dissipation to a variable

    This takes the variable and automatically applies the Kreiss-Oliger
    dissipation. It assumes that it was already set through the
    `dendrosym.nr.set_kreiss_oliger_dissipation`

    Parameter
    ---------
    a : sympy.Symbol, tuple, sympy.Matrix
        The variable to apply Kreiss-Oliger dissipation to

    Returns
    -------
    sympy.Expression, list, sympy.Matrix
        The expression that applies the Kreiss-Oliger dissipation.
        The output type depends on what was input to the function.
    """

    global kod

    if type(a) == sym.Symbol:
        return sum([kod(ii, a) for ii in e_i])

    elif type(a) == tuple:
        return [sum([kod(ii, a[jj]) for ii in e_i]) for jj in e_i]

    elif type(a) == sym.Matrix:
        return sym.Matrix(
            [sum([kod(kk, a[ii, jj]) for kk in e_i]) for ii, jj in e_ij]
        ).reshape(3, 3)

    else:
        raise ValueError("Dendro: Unknown type for input to compute kodiss.")


def laplacian(a, chi):
    """Compute the Laplacian of a scalar with respect to the metric

    Assumes that the conformally rescaled metric (called gt in
    various places) and the conformal factor (chi) is set.  Note that C3 is
    built from the same 3D metric.  The only place that this Laplacian is
    used in the BSSN equations is in the evolution equation for K and is
    the Laplacian of alpha (the lapse).

    Parameters
    ----------
    a : sympy.Symbol
        The scalar that will have the Laplacian calculated of
    chi : sympy.Symbol
        The conformal factor required in the BSSN equation

    Returns
    -------
    sympy.Expression
        The final expression that calculates the Laplacian of the inputs
    """
    global d, metric, C3

    full_metric = metric / chi
    inv_full_metric = sym.simplify(full_metric.inv("ADJ"))

    # return sum([(inv_full_metric[i, j] * d2(i, j, a) -
    #            sum([C3[l, i, j] * d(l, a) for l in e_i])) for i, j in e_ij])
    return sum(
        [
            inv_full_metric[ii, jj]
            * (d2(ii, jj, a) - sum([C3[ll, ii, jj] * d(ll, a) for ll in e_i]))
            for ii, jj in e_ij
        ]
    )


def laplacian_conformal(a):
    """Compute the conformal Laplacian of a scalar with respect to the metric

    We assume the rescaled metric is set as well the conformal
    factor, chi.  Note that C2 is built from the conformally rescaled
    metrci.  This (conformal) laplacian is only used in the definition of
    Ricci that shows up in the evolution equation for At (under the trace
    free operation), and even then only in the part that multiplies the
    metric and which will drop out on taking the trace free part.  So, in
    fact, the code could be written to completely ignore this operation
    in the evolution equations themselves.  However, if the constraints
    are included or the full Ricci is needed for another reason, this
    would be needed.

    Parameters
    ----------
    a : sympy.Symbol
        The scalar that will have the Laplacian calculated of

    Returns
    -------
    sympy.Expression
        The final expression that calculates the Laplacian of the inputs
    """

    global d, inv_metric, C2

    if inv_metric == undef:
        inv_metric = get_inverse_metric()

    # ewh3
    # return sum([(inv_metric[i, j] * d2(i, j, a) - sum([C2[l, i, j]\
    #     * d(l, a) for l in e_i])) for i, j in e_ij])
    return sum(
        [
            inv_metric[ii, jj]
            * (d2(ii, jj, a) - sum([C2[ll, ii, jj] * d(ll, a) for ll in e_i]))
            for ii, jj in e_ij
        ]
    )


def sqr(A):
    """Computes the square of the matrix.

    This will take a 3x3 matrix and will compute its square.
    Assumes the metric is set.

    Parameters
    ----------
    A : sympy.Matrix
        The input matrix to compute the square of

    Returns
    -------
    sympy.Expression
        The expression that computes the matrix square
    """

    global inv_metric

    if inv_metric == undef:
        inv_metric = get_inverse_metric()

    return sum(
        [
            A[ii, jj]
            * sum(
                [
                    inv_metric[ii, kk] * inv_metric[jj, ll] * A[kk, ll]
                    for kk in e_i
                    for ll in e_i
                ]
            )
            for ii, jj in e_ij
        ]
    )


def trace(A):
    """Calculate the trace of a Matrix

    Assumes the metric is set, will get the inverse metric.

    Parameters
    ----------
    A : sympy.Matrix
        The input matrix

    Returns
    -------
    sympy.Expression
        The expression that computes the trace of the matrix.
    """

    global inv_metric

    if inv_metric == undef:
        inv_metric = get_inverse_metric()

    trace = sum([inv_metric[ii, jj] * A[ii, jj] for ii, jj in e_ij])

    return trace


def trace_free(A):
    """Transforms an input matrix to be trace free

    If an equation requires a matrix to be trace free, then this
    function is necessary. Assumes the metric is set.

    Parameters
    ----------
    A : sympy.Matrix
        The input matrix that needs to be converted to have no
        trace

    Returns
    -------
    sympy.Matrix
        The now-trace-free matrix
    """
    global metric, inv_metric

    if inv_metric == undef:
        inv_metric = get_inverse_metric()

    trace = sum([inv_metric[ii, jj] * A[ii, jj] for ii, jj in e_ij])

    # X_{ab} - 1/3 gt_{ab} X.
    # tf = sym.Matrix([x[i, j] - 1/3*metric[i,j]*trace for i, j in e_ij])
    tf = sym.Matrix([A[ii, jj] - metric[ii, jj] * trace / 3 for ii, jj in e_ij])

    return tf.reshape(3, 3)


def vec_j_del_j(b, a):
    """Computes the regular derivative of a scalar over all indices with a vector

    This function calculates the following:

    .. math:: \\beta^i \\partial_i \\alpha

    Parameters
    ----------
    b : tuple, sympy.Matrix
        The 3 vector "beta"
    a : sympy.Symbol
        The scalar "alpha"

    Returns
    -------
    sympy.Expression
        The output expression
    """
    return sum([b[ii] * d(ii, a) for ii in e_i])


def vec_j_ad_j(b, f):
    """Computes the advective derivative of a scalar over all indices with a vector

    Note that this function calculates the advective derivative when
    using the following equation:

    .. math:: \\beta^i \\partial_i \\alpha

    Parameters
    ----------
    b : tuple, sympy.Matrix
        The 3 vector "beta"
    a : sympy.Symbol
        The scalar "alpha"

    Returns
    -------
    sympy.Expression
        The output expression

    Note
    ----
    [ewh] Adding this as this term needs to be in the beta equation as an
    advective derivative and not as a regular (partial) derivative.
    """

    return sum([b[ii] * ad(ii, f) for ii in e_i])


def calc_symmetric_part_rank2(M):
    """Compute the symmetric part of a rank 2 tensor

    Parameters
    ----------
    M : sympy.Matrix
        The rank 2 tensor (3x3 matrix) to compute the symmetric
        part of.

    Returns
    -------
    sympy.Matrix
        The newly-symmetric tensor (3x3 matrix)
    """

    return (1 / 2) * sym.Matrix([M[ii, jj] + M[jj, ii] for ii, jj in e_ij]).reshape(
        3, 3
    )


##########################################################################
# metric related functions
##########################################################################


def set_metric(g):
    """Set the metric variable to be used in the NR calculations

    Sets the metric variable, so that DendroSym knows how to compute
    various derived variables. This should be done early on in the
    generating script.

    Parameters
    ----------
    g : sympy.Matrix
        The 3x3 matrix that will be used as the metric

    Example
    -------
    >>> gt = dendrosym.dtypes.sym_3x3("gt")
    >>> dendrosym.nr.set_metric(gt)
    """

    global metric

    metric = g


def get_inverse_metric():
    """Computes and returns the inverse metric

    Sets the inverse metric variable, so that DendroSym knows how to compute
    various derived variables. This should be done early on in the
    generating script. It requires the metric to already be defined.

    Returns
    -------
    sympy.Matrix
        The 3x3 matrix returned is the inverse metric

    Example
    -------
    >>> gt = dendrosym.dtypes.sym_3x3("gt")
    >>> dendrosym.nr.set_metric(gt)
    >>> igt = dendrosym.nr.get_inverse_metric()
    """

    global metric, inv_metric, undef

    if metric == undef:
        raise ValueError("Dendro: Metric not defined.")

    if inv_metric == undef:
        # method : ('GE', 'LU', or 'ADJ')
        inv_metric = sym.simplify(metric.inv("ADJ"))

    return inv_metric


def get_first_christoffel():
    """Computes and returns the first Christoffel Symbols

    This function will take the metric, inverse metric,
    first derivative, and then calculate the first Christoffel
    Symbols. It will store these symbols internally for
    other calculations as well.

    If these first Christoffel Symbols were already computed,
    it will just return the already-stored Symbols.

    Returns
    -------
    sympy.Matrix
        The 3x3 matrix returned is the first Christoffel Symbols

    Example
    -------
    >>> C1 = dendrosym.nr.get_first_christoffel()
    """
    global metric, inv_metric, undef, C1, d

    if inv_metric == undef:
        get_inverse_metric()

    if C1 == undef:
        C1 = sym.tensor.array.MutableDenseNDimArray(range(27), (3, 3, 3))

        for kk in e_i:
            for jj in e_i:
                for ii in e_i:
                    # C1[k, i,
                    #    j] = 1/2 * (d(j, metric[k, i]) + d(i, metric[k, j]) -
                    #                  d(k, metric[i, j]))
                    C1[kk, ii, jj] = 0.5 * (
                        d(jj, metric[kk, ii])
                        + d(ii, metric[kk, jj])
                        - d(kk, metric[ii, jj])
                    )

    return C1


def get_second_christoffel():
    """Computes and returns the second Christoffel Symbols

    This function will take the metric, inverse metric,
    first derivative, and then calculate the second Christoffel
    Symbols. It will store these symbols internally for
    other calculations as well.

    If the first Christoffel Symbols were not already computed,
    it will compute them and store them internally.

    If these second Christoffel Symbols were already computed,
    it will just return the already-stored Symbols.

    Returns
    -------
    sympy.Matrix
        The 3x3 matrix returned is the second Christoffel Symbols

    Example
    -------
    >>> C2 = dendrosym.nr.get_second_christoffel()
    """

    global C2, C1, inv_metric

    if C2 == undef:
        if C1 == undef:
            get_first_christoffel()

        igt_t = sym.Array(inv_metric, (3, 3))
        C2 = sym.tensor.array.tensorcontraction(
            sym.tensor.array.tensorproduct(igt_t, C1), (1, 2)
        )

    return C2


def get_complete_christoffel(chi):
    """Computes and returns the complete Christoffel Symbols

    This function will take the metric, inverse metric,
    first derivative, and then calculate the complete Christoffel
    Symbols. It will store these symbols internally for
    other calculations as well.

    If the second Christoffel Symbols were not already computed,
    it will compute them and store them internally.

    If these complete Christoffel Symbols were already computed,
    it will just return the already-stored Symbols.

    Parameters
    ----------
    chi : sympy.Scalar
        The input scalar used in the calculation

    Returns
    -------
    sympy.Matrix
        The 3x3 matrix returned is the complete Christoffel Symbols

    Example
    -------
    >>> C3 = dendrosym.nr.get_complete_christoffel()
    """

    global metric, inv_metric, undef, C1, C2, C3, d

    if C3 == undef:
        C3 = sym.tensor.array.MutableDenseNDimArray(range(27), (3, 3, 3))

        if C2 == undef:
            get_second_christoffel()

        for kk in e_i:
            for jj in e_i:
                for ii in e_i:
                    # C3[i, j, k] = C2[i, j, k] - \
                    # 1/(2*chi)*(KroneckerDelta(i, j) * d(k, chi) +
                    C3[ii, jj, kk] = C2[ii, jj, kk] - 0.5 / (chi) * (
                        KroneckerDelta(ii, jj) * d(kk, chi)
                        + KroneckerDelta(ii, kk) * d(jj, chi)
                        - metric[jj, kk]
                        * sum([inv_metric[ii, mm] * d(mm, chi) for mm in e_i])
                    )

    return C3


def compute_ricci(Gt, chi):
    """Computes the Ricci tensor

    The conformal connection coefficient and the conformal variable
    need to be supplied.

    The metric needs to already have been set,as well as the first
    and second Christoffel symbols.

    Parameters
    ----------
    Gt : tuple, sympy.Matrix
        The 3-vector input Gt values
    chi : sympy.Scalar
        The input scalar used in the calculations

    Returns
    -------
    R : sympy.Matrix
        The Ricci tensor (includes the additional addon term of
        Rphi)
    Rt : sympy.Matrix
        The Ricci tensor (without Rphi added)
    Rphi : sympy.Matrix
        Just the Rphi portion
    CalGt : sympy.Matrix
        The calculated Gt matrix

    Example
    -------
    >>> dendrosym.nr.compute_ricci(Gt, chi)
    """
    global metric, inv_metric, C1, C2

    # TODO: is this necessary? They are never referenced or returned
    # NOTE: DFVK has commented this piece out because it's never used!!!
    # Lchi = laplacian_conformal(chi)

    # print(type(Lchi))

    # print('Done with Lphi') #simplify(Lchi))

    # ewh4 DKchiDkchi = sym.Matrix([
    #     4 * metric[i, j] * sum([
    #         sum([inv_metric[k, l] * d(l, chi) for l in e_i]) * d(k, chi)
    #         for k in e_i
    #     ]) for i, j in e_ij
    # ])
    # TODO: is this necessary? THey are never referenced or returned
    # NOTE: DFVK has commented this piece out because it's never used!!!!!
    # DKchiDkchi = sym.Matrix([
    #     0.25 / chi / chi * metric[i, j] * sum([
    #         sum([inv_metric[k, l] * d(l, chi) for l in e_i]) * d(k, chi)
    #         for k in e_i
    #     ]) for i, j in e_ij
    # ])

    # print('done with DKchi') # simplify(DKchiDkchi))

    CalGt = [sum(inv_metric[kk, ll] * C2[ii, kk, ll] for kk, ll in e_ij) for ii in e_i]

    Rt = sym.Matrix(
        [
            -0.5
            * sum([inv_metric[ll, mm] * d2(ll, mm, metric[ii, jj]) for ll, mm in e_ij])
            + 0.5
            * sum(
                [
                    metric[kk, ii] * d(jj, Gt[kk]) + metric[kk, jj] * d(ii, Gt[kk])
                    for kk in e_i
                ]
            )
            + 0.5 * sum([CalGt[kk] * (C1[ii, jj, kk] + C1[jj, ii, kk]) for kk in e_i])
            + sum(
                [
                    inv_metric[ll, mm]
                    * (
                        C2[kk, ll, ii] * C1[jj, kk, mm]
                        + C2[kk, ll, jj] * C1[ii, kk, mm]
                        + C2[kk, ii, mm] * C1[kk, ll, jj]
                    )
                    for kk in e_i
                    for ll, mm in e_ij
                ]
            )
            for ii, jj in e_ij
        ]
    )

    # print('done with Rt') #simplify(Rt))

    # ewh5  Rphi_tmp = sym.Matrix(
    #     [2 * metric[i, j] * Lchi - \
    #        4 * d(i, chi) * d(j, chi) for i, j in e_ij])
    # dwn Rphi_tmp = sym.Matrix([
    #     0.5 * metric[i, j] * Lchi / chi -
    #     0.25 * d(i, chi) * d(j, chi) / chi / chi for i, j in e_ij
    # ])

    # print(simplify(Rphi_tmp))

    # ewh6    Rphi = -2*_Di_Dj(chi) - Rphi_tmp.reshape(3, 3) - \
    #         DKchiDkchi.reshape(3, 3)
    # dwn    Rphi = -0.5*_Di_Dj(chi)/chi - Rphi_tmp.reshape(3, 3) - \
    #        DKchiDkchi.reshape(3, 3)
    xRphi = sym.Matrix(
        [
            1
            / (2 * chi)
            * (d2(ii, jj, chi) - sum(C2[kk, jj, ii] * d(kk, chi) for kk in e_i))
            - 1 / (4 * chi * chi) * d(ii, chi) * d(jj, chi)
            for ii, jj in e_ij
        ]
    ).reshape(3, 3)

    Rphi = xRphi + sym.Matrix(
        [
            1
            / (2 * chi)
            * metric[ii, jj]
            * (
                sum(
                    inv_metric[kk, ll]
                    * (d2(kk, ll, chi) - 3 / (2 * chi) * d(kk, chi) * d(ll, chi))
                    for kk, ll in e_ij
                )
                - sum(CalGt[mm] * d(mm, chi) for mm in e_i)
            )
            for ii, jj in e_ij
        ]
    ).reshape(3, 3)

    return Rt.reshape(3, 3) + Rphi, Rt.reshape(3, 3), Rphi, CalGt
