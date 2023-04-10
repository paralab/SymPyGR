# # Energetics of the Einstein Rosen Spacetime by Herrera et al
# ##  Geoff Cope
# ##  Univeristy of Utah
# ##  January 2, 2022

# https://arxiv.org/abs/gr-qc/0606052

# In[1]:

import warnings

import matplotlib.cbook
import sympy

import einsteinpy
from einsteinpy.symbolic import *

k, m, n = sympy.symbols('k m n', integer=True)
f, g, h = sympy.symbols('f g h', cls=sympy.Function)

# init_session(use_latex=True)

# In[2]:


warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

# In[3]:

t, r, phi, z, dt, dr, dphi, dz = sympy.symbols('t r phi z dt dr dphi dz')

# In[4]:

psi = sympy.Function('psi')(r, t)
psi

# In[5]:

gamma = sympy.Function('gamma')(r, t)
gamma

# In[6]:

variables = sympy.Matrix([t, r, phi, z])
variables

# In[7]:

differentials = sympy.Matrix([dt, dr, dphi, dz])
differentials

# In[8]:

lineElement = sympy.expand(-sympy.exp(2 * (gamma - psi)) * (dt**2 - dr**2) +
                     sympy.exp(2 * psi) * dz**2 + r**2 * sympy.exp(-2 * psi) * dphi**2)
lineElement

# In[9]:

g = sympy.zeros(4)

for i in range(4):
    for j in range(4):
        if i == j:
            g[i, j] = lineElement.coeff(differentials[i], 2)
        else:
            g[i, j] = sympy.Rational(1, 2) * lineElement.coeff(
                differentials[i] * differentials[j], 1)

g

# In[10]:

# In[11]:

m = sympy.Array(g)
m

# In[12]:



# In[13]:

syms = sympy.symbols("t r phi z")
t, r, phi, z = syms

# In[14]:

metric = MetricTensor(m, syms)

# In[15]:

ch = ChristoffelSymbols.from_metric(metric)
sympy.simplify(ch.tensor())

# In[16]:

Ric = RicciTensor.from_metric(metric)
sympy.simplify(Ric.tensor())

# In[17]:

R = RicciScalar.from_riccitensor(Ric)
R.simplify()
R.expr

# In[18]:

einst = EinsteinTensor.from_metric(metric)
einst.tensor()

# In[19]:

#  rm1 = RiemannCurvatureTensor.from_christoffels(ch)
#  rm1.tensor()

# In[20]:

#  weyl = WeylTensor.from_metric(metric)
#  weyl.tensor()

# In[21]:

einsteinSimplifed = sympy.simplify(einst.tensor())
einsteinSimplifed

# In[22]:

equation2 = sympy.Eq(sympy.expand((-1 / r) * sympy.simplify(Ric.tensor())[2, 2].args[1]), 0)
equation2

# In[23]:

equation3 = sympy.Eq(sympy.diff(gamma, t),
               sympy.solve(einsteinSimplifed[1, 0], sympy.diff(gamma, t))[0])
equation3

# In[24]:

equation4 = sympy.Eq(sympy.diff(gamma, r),
               sympy.solve(einsteinSimplifed[1, 1], sympy.diff(gamma, r))[0])
equation4

# In[34]:

vacuumFieldEquations = [equation2, equation3, equation4]
vacuumFieldEquations

# In[27]:

X, Y = map(sympy.Function, 'XY')

# In[35]:

eq = vacuumFieldEquations[0]
eq

# In[29]:

xODE = sympy.Eq(sympy.pde_separate(eq, psi, [X(r), Y(t)])[0], k**2)
xODE

# In[30]:

xSolution = sympy.dsolve(xODE, X(r))
xSolution

# In[31]:

yODE = sympy.Eq(sympy.pde_separate(eq, psi, [X(r), Y(t)])[1], -k**2)
print("yODE:", yODE)

# In[32]:

ySolution = sympy.dsolve(yODE, Y(t))
print("ySolution:", ySolution)

# In[ ]:
