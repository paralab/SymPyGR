import math

import quaternion
import spherical_functions as sf
from sympy import pprint

import dendrosym


class GWExtract:
    """Object to define Gravitational Wave Extraction

    Originally written by Milinda Fernando
    """

    def __init__(self, N, qprec, lmodes, spin, varNames) -> None:

        # number of different radii for psi4 polynomial approx
        self.N = N

        # precision for the lebedev quadrature
        self.qprec = qprec

        # l modes for SWSH
        self.lmodes = lmodes

        self.numLModes = len(self.lmodes)

        # spin values for SWSH (2)
        self.spin = spin

        # computing the SWSH decomposition
        # self.r=dtypes.scalar(varNames[0])
        # self.psi4_s2=dtypes.scalar(varNames[1])
        # self.psi4_lm=dtypes.scalar(varNames[2])

        # part of generating symbolic expression for the polynomial fit.
        # self.V=dtypes.mat(varNames[3],N,N)
        self.vR = dendrosym.dtypes.vec("r", N)  # radius of the extraction
        self.psi4 = dendrosym.dtypes.vec(
            "psi4_", N
        )  # different lm values at the extraction radius
        self.argpsi4 = dendrosym.dtypes.vec(
            "apsi4_", N
        )  # different lm values at the extraction radius
        self.A = dendrosym.dtypes.vec("a", N)  # coefficients for the polynomial fit.
        self.B = dendrosym.dtypes.vec(
            "b", N
        )  # coefficennts for the phase polynomical fit.

        # fname=sys.path[0]+"/Lebedev/lebedev_" +"{:03d}".format(self.qprec)+".txt"
        fname = "Lebedev/lebedev_" + "{:03d}".format(self.qprec) + ".txt"
        self.lebedev_theata = []
        self.lebedev_phi = []
        self.lebedev_w = []

        deg_radians = math.pi / 180.0
        # read lebedev points
        with open(fname, "r") as f:
            for row in f:
                row = row.strip().split()
                self.lebedev_theata.append(deg_radians * float(row[1]))
                self.lebedev_phi.append(deg_radians * (float(row[0]) + 180.0))
                self.lebedev_w.append(float(row[2]))

    def swsh(self, s, l, m, theta, phi):
        # get the rotor related for Wigner matrix
        # Note: We need to specify the (-theta,-phi direction to match the convenction used in HAD)
        # 10/12/20: David has pointed out this should be theta, phi to match the SWSH used in the NR community. so changed back to theta phi.
        # this is the reason that we have a sign difference between the LazEv code.
        R_tp = quaternion.from_spherical_coords(theta, phi)
        W = sf.Wigner_D_element(R_tp, l, m, -s)
        # print(W)
        return ((-1) ** (s)) * math.sqrt((2 * l + 1) / (4 * math.pi)) * W

    def initVars(self, varnames):

        outputstr = ""

        outputstr += f"#define {varnames[0]} {len(self.lebedev_theata)} \n"

        # allocate Lebedev quadrature points
        outputstr += f"static const double {varnames[1]} [] = {{ {','.join(['{:.20f}'.format(x) for x in self.lebedev_theata])} }};"
        outputstr += f"static const double {varnames[2]} [] = {{ {','.join(['{:.20f}'.format(x) for x in self.lebedev_phi])} }};"
        outputstr += f"static const double {varnames[3]} [] = {{ {','.join(['{:.20f}'.format(x) for x in self.lebedev_w])} }};"

        for l in self.lmodes:
            for m in range(-l, l + 1):
                lm_value = []
                for q in range(0, len(self.lebedev_theata)):
                    lm_value.append(
                        self.swsh(-2, l, m, self.lebedev_theata[q], self.lebedev_phi[q])
                    )

                if m < 0:
                    outputstr += f"static const DendroComplex {varnames[4]}_sp_{2}_l_{l}_m_minus_{abs(m)} [] = {{ {','.join(['DendroComplex({:.20f},{:.20f})'.format(x.real,x.imag) for x in lm_value])} }};"
                else:
                    outputstr += f"static const DendroComplex {varnames[4]}_sp_{2}_l_{l}_m_plus_{abs(m)} [] = {{ {','.join(['DendroComplex({:.20f},{:.20f})'.format(x.real,x.imag) for x in lm_value])} }};"

        ptr_str = []
        for i in range(0, len(self.lmodes)):
            for m in range(-self.lmodes[i], self.lmodes[i] + 1):
                if m < 0:
                    ptr_str.append(
                        "%s_sp_%d_l_%d_m_minus_%d"
                        % (varnames[4], 2, self.lmodes[i], abs(m))
                    )
                else:
                    ptr_str.append(
                        "%s_sp_%d_l_%d_m_plus_%d"
                        % (varnames[4], 2, self.lmodes[i], abs(m))
                    )

        outputstr += (
            f"static const DendroComplex* {varnames[4]}[]={{ {','.join(ptr_str)} }};"
        )

        return outputstr

    def psi4PolyFit(self):
        """Compute the polynomial fit for the psi4 extraction at the far radius"""

        for r in range(0, self.N):
            for j in range(0, self.N):
                self.V[r, j] = 1 / (self.vR[r] ** (j + 1))

        self.Vinv = self.V ** (-1)
        self.A = self.Vinv * self.psi4
        self.B = self.Vinv * self.argpsi4

        pprint(self.A)
        pprint(self.B)

    def propTimeCompute(self):
        """Compute the t* propagation time computation"""
        print("not implemented yet")
