'''
@author : Milinda Fernando
School of Computing, University of Utah. 
@date: 22 Jan. 2019

@package: Symbolic rerferance element class, for FEM computation. 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense,and/or sell copies
of the Software,and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software. 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.

'''

import dendrosym.dtypes as dtypes


class RefEelement:
    def __init__(self, dim, eleOrder):

        # dimention of the reference element
        self._dim = dtypes.scalar("dim")

        # nodal locations of the element
        self._r = dtypes.vec("r", eleOrder + 1)

        # weights for the gll quadrature
        self._wgll = dtypes.vec("wgll", eleOrder + 1)

        # gauss points for quadrature
        self._g = dtypes.vec("g", eleOrder + 1)

        # weights for the gauss quadrature
        self._wg = dtypes.vec("w", eleOrder + 1)

        # polynomial Vandermonde matrix evaluated at the reference points
        self._Vr = dtypes.mat("Vr", eleOrder + 1, eleOrder + 1)

        # gradient of the polynomial Vandermonde matrix evaluated at the reference points
        self._gradVr = dtypes.mat("gradVr", eleOrder + 1, eleOrder + 1)

        # polynomial Vandermonde matrix evaluated at the gauss points
        self._Vg = dtypes.mat("Vg", eleOrder + 1, eleOrder + 1)

        # polynomial Vandermonde matrix evaluated at the gauss points
        self._gradVg = dtypes.mat("gradVg", eleOrder + 1, eleOrder + 1)

        # Gradient interpolation to the reference points from gauss points
        self._Dr = dtypes.mat("gradDr", eleOrder + 1, eleOrder + 1)

        # Gradient interpolation to the gauss points from reference points
        self._Dg = dtypes.mat("gradDg", eleOrder + 1, eleOrder + 1)

        # matrix Identity
        self._matI = dtypes.matI(eleOrder + 1)

        # parent to child (one child) interpolation
        self._Vph = dtypes.mat("Vph", eleOrder + 1, eleOrder + 1)

        # parent to all children interpolation
        self._Vpp = dtypes.mat("Vpp", eleOrder + 1, eleOrder + 1)
