"""__init__.py

This is the initialization file for the DendroSym Python package.
It helps us make sure that everything is ready to roll
"""

# the datatypes that includes all of the sympy pieces
import dendrosym.dtypes
# numerical relativity functions
import dendrosym.nr

# code generation
import dendrosym.codegen
# parameter information
import dendrosym.params

# ==== MISC. FUNCTIONS AND OPERATIONS ====
# the basic memory manager
import dendrosym.memoryManager
# nxgraph generation
import dendrosym.nxgraph
# reference element information
import dendrosym.refEl
# sympy cache simulations
import dendrosym.sympy_cachesim
