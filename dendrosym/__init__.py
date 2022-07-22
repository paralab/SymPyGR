"""__init__.py

This is the initialization file for the DendroSym Python package.
It helps us make sure that everything is ready to roll
"""

# the datatypes that includes all of the sympy pieces
from . import dtypes

# numerical relativity functions
from . import nr

# code generation
from . import codegen

# parameter information
from . import params

# ==== MISC. FUNCTIONS AND OPERATIONS ====
# the basic memory manager
from . import memoryManager

# nxgraph generation
from . import nxgraph

# reference element information
from . import refEl

# sympy cache simulations
# from . import sympy_cachesim
# TODO: cachesim currently does not work on some machines

from . import utils

from . import params

# ====== DENDRO CONFIGURATION ======
# base configuration class
from .general_configs import DendroConfiguration

# and then the numerical relativity class
from .nr_configs import NRConfig

# TODO: the package that is used with gw is not currently supported!
# from . import gw


# TODO: quaternion and cachesym are not working on MaryLou currently
