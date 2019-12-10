import math
from enum import Enum
import warnings

import sys    
import os
import datetime
import ntpath
from sympy import *
from sympy.printing.cxxcode import cxxcode
from sympy.codegen.ast import String

from collections import OrderedDict, namedtuple

#holds dict of parameters, contains mostly output functions
class Parameters:
	def __init__(self, paramDict=OrderedDict(), initialData=OrderedDict(), varsDict = OrderedDict()):
		self.paramDict = paramDict
		self.initialData = initialData
		self.varsDict = varsDict

		self.indent = "\t"

		# if set, will be used as the default category for any newly added parameters
		self.category = None
		self.generationNote = "This file was generated at {0}, through the parameters object defined in {1}".format(
								datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"),
								os.path.basename(sys.argv[0]))

	def __getitem__(self, key):
		return self.paramDict[key]

	#add is probably the simpler thing to use generally. But this gives direct access
	def __setitem__(self, key, value):
		self.paramDict[key] = value

	def add(self, id, value, description=None, category=None, cppType=None, arraySize = None):
		#set default value if it makes sense
		if category is None and self.category is not None:
			category = self.category

		if isinstance(id, PresetParam):
			#get pre-made parameter object from enum value
			param = id.value
			param.category = category
			param.setValue(value, cppType, arraySize)
			#note the dict key will be the PresetParam enum, not a string
			self.paramDict[id] = param
		else:
			param = Parameter(id, value, description, category, cppType, arraySize)
			self.paramDict[param.id] = param

	def addInitialData(self, id, expression, category = None):
		if category is None:
			category = self.category
		if category is None:
			category = "Initial Data"

		if isinstance(expression, Expr):

			evalExpression = expression
			for initialData in self.initialData.values():
				evalExpression = evalExpression.subs(initialData.symbol,initialData.value)
			
			self.initialData[id] = ExpressionParameter(id,evalExpression.evalf(), expression, category)
		else:
			self.initialData[id] = ExpressionParameter(id,expression, category=category)

		return self.initialData[id].symbol
	def addVar(self, id, initialValue, rhsExpression):
		if isinstance(initialValue, Expr):

			evalExpression = initialValue
			for var in list(self.varsDict.values()) + list(self.initialData.values()):
				evalExpression = evalExpression.subs(var.symbol,var.value)
			
			self.varsDict[id] = Var(id, evalExpression.evalf(), valueExpression=initialValue, rhsExpression=rhsExpression)
		else:
			self.varsDict[id] = Var(id, initialValue, rhsExpression = rhsExpression)

		return self.varsDict[id].symbol

	#func is a PrecomputeFunc enum
	#varSymbol is sympy Symbol object associated with a var, returned by addVar
	#in c++, an array will be created and computed ahead of time, so it can then be used to compute var rhs
	def addRhsPrecompute(self, func, varSymbol):
		self.varsDict[varSymbol.name].addPrecompute(func)
		#returned string can be used in sympy expressions for RHS
		return String(func.getArrayName(self.varsDict[varSymbol.name]) + "[pp]")

	def setCategory(self, category):
		self.category = category

	def clearCategory(self):
		self.category = None

	def values(self):
		return self.paramDict.values()
	def initialDataValues(self):
		return self.initialData.values()
	def varsValues(self):
		return self.varsDict.values()

	#generate parameters.cpp and parameters.h equivalent files
	def writeCpp(self, hPath, cppPath, namespace):

		hName = ntpath.basename(hPath)
		cppName = ntpath.basename(cppPath)

		# Initialize output and first lines for .h file
		os.makedirs(os.path.dirname(hPath), exist_ok=True)
		outH = open(hPath, "w")
		outH.write("//" + self.generationNote + "\n")
		outH.write("#ifndef SFCSORTBENCH_PARAMETERS_H\n")
		outH.write("#define SFCSORTBENCH_PARAMETERS_H\n\n")

		outH.write("#include <string.h>\n")
		outH.write("#include <iostream>\n\n")

		outH.write("namespace {0}\n{{\n\n".format(namespace))

		# Initialize output and first lines for .cpp file.
		os.makedirs(os.path.dirname(cppPath), exist_ok=True)
		outCpp = open(cppPath, "w")
		outCpp.write("//" + self.generationNote + "\n")

		outCpp.write("#include \"{0}\"\n".format(hName))
		outCpp.write("namespace {0}\n{{\n\n".format(namespace))

		for parameter in list(self.paramDict.values()) + list(self.initialData.values()):
			if parameter.isComment():
				continue
			outH.write(parameter.toStringH(1))
			outCpp.write(parameter.toStringCpp(1))

		outH.write("\tstatic const unsigned int NUM_VARS={0};\n".format(len(self.varsDict)))

		outH.write("\tenum VAR {")
		i = 0
		for var in self.varsValues():
			if i > 0:
				outH.write(", ")
			outH.write("{0}={1}".format(var.id, i))
			i+=1

		outH.write("};\n")

		outH.write("\tstatic const char * VAR_NAMES[] = {")
		i = 0
		for var in self.varsValues():
			if i > 0:
				outH.write(", ")
					
			outH.write("\"" + var.id + "\"")
			i+=1

		outH.write("};\n")

		outH.write("\tstatic const unsigned int RK3_STAGES=3;\n")
		outH.write("\tstatic const unsigned int RK4_STAGES=4;\n")
		outH.write("\tstatic const unsigned int RK45_STAGES=6;\n")

		outH.write("}\n")
		outH.write("#endif //SFCSORTBENCH_PARAMETERS_H")
		outH.close()

		outCpp.write("}")
		outCpp.close()

	def writeJson(self, filepath):
		categories = OrderedDict()
		for parameter in list(self.paramDict.values()) + list(self.initialData.values()):
			if parameter.category not in categories:
				categories[parameter.category] = ""
			categories[parameter.category] += parameter.toStringJson(2)

		os.makedirs(os.path.dirname(filepath), exist_ok=True)

		out = open(filepath, "w")
		out.write("{\n")
		out.write('"__comment__" : "{0}",\n'.format(self.generationNote))

		first = True
		for key in categories:
			if not first:
				out.write(",\n")
				first = False

			keyLength = len(key) if key is not None else 0
			if keyLength > 0:
				out.write('{0}"__comment__" : "{1} {2} {3}",\n'.format(self.indent, int(math.floor((50 - keyLength) / 2) - 1) * "=", key,
																	int(math.ceil((50 - keyLength) / 2) - 1) * "="))
			out.write(categories[key])
			if keyLength > 0:
				out.write('{0}"__comment__" : "{1}"\n\n'.format(self.indent, "=" * 50))

		out.write("\n}")
		out.close()

	def writeCMakeLists(self, filepath, projectName):
		out = open(filepath, "w")


class Parameter:
	def __init__(self, id, value, description=None, category=None, cppType=None, arraySize=None):
		self.id = id
		self.description = description
		self.category = category
		self.indent = "\t"
		self.setValue(value, cppType, arraySize)

	#cppType and arraySize may need to be inferred from the value, so handle that all in one function
	def setValue(self, value, cppType=None, arraySize=None):
		self.value = value

		if isinstance(arraySize, int):
			self.arraySize=arraySize
		# if value is a list but array size isn't explicitly set, get it from values
		elif isinstance(value, list) and len(value) > 0:
			self.arraySize = len(value)
			# set temp value to first element in list, so we can get the array type below
			value = value[0]
		else:
			self.arraySize = None

		if isinstance(cppType, CppType):
			self.cppType = cppType
		else:
			if cppType is not None:
				warnings.warn("Invalid value for cppType arg when creating parameter. Should use CppType enum.")

			# infer data type if not set
			if isinstance(value, int):
				self.cppType = CppType.unsignedInt
			elif isinstance(value, float) or isinstance(value, Expr):
				self.cppType = CppType.double
			#default to string if nothing else matches
			else:
				self.cppType = CppType.string

	def isComment(self):
		return self.id.startswith("__comment__", 0, 11)

	def toStringJson(self, indentCount=0):
		output = ""
		if (self.description is not None):
			output += '{0}"           " : "{1}",\n'.format(self.indent * indentCount, self.description)
		output += '{0}"{1}" : {2},\n'.format(self.indent * indentCount, self.id,
											 '"'+self.value +'"' if isinstance(self.value, str) else self.value)

		return output

	#get the value of the parameter, enclosed in quotes if necessary, for use in c++
	def getCppValue(self):
		valueQuote = '"' if self.cppType is CppType.string else ""
		return "{0}{1}{0}".format(valueQuote,self.value)

	#get the c++ string definition for this parameter to be used in parameters.h equivalent
	def toStringH(self, indentCount = 0):
		output = ""
		if self.description is not None:
			output += "{0}//{1}\n".format(self.indent * indentCount, self.description)

		output += "{0}extern {1} {2}{3};\n\n".format("\t" * indentCount, self.cppType.value, self.id,
											  "[" + str(self.arraySize) +"]" if self.arraySize is not None else "")
		return output

	#get the c++ string definition for this parameter to be used in parameters.cpp equivalent
	def toStringCpp(self, indentCount = 0):
		output = ""
		if self.description is not None:
			output += "{0}//{1}\n".format(self.indent * indentCount, self.description)

		if self.arraySize is None:
			output += "{0}{1} {2}={3};\n\n".format(self.indent * indentCount, self.cppType.value, self.id, self.getCppValue())
		else:
			# surround value with quotes if it's a string, otherwise nothing
			valueQuote = '"' if self.cppType is CppType.string else ""
			arrayValue = valueQuote + "{0},{0}".format(valueQuote).join(str(v) for v in self.value) + valueQuote
			output += "{0}{1} {2}[{3}]={{{4}}};\n\n".format(self.indent * indentCount, self.cppType.value, self.id, self.arraySize, arrayValue)
			
		return output

#parameter capable of handling sympy expression as its value
class ExpressionParameter(Parameter):
	def __init__(self, id, value, expression=None, category = None):
		Parameter.__init__(self, id, value, category=category)
		self.symbol = symbols(id)
		self.expression = expression

	def getCppExpression(self, expression):
		return cxxcode(expression, standard="C++11")

	def toStringJson(self, indentCount=0):
		#params defined with an expression will have their value calculated
		#in the cpp directly. Do not want to include them in json
		if self.expression is not None:
			return ""
		return Parameter.toStringJson(self, indentCount)

	def toStringCpp(self, indentCount=0):
		if self.expression is not None:
			return ""
		return Parameter.toStringCpp(self, indentCount)

	def toStringH(self, indentCount=0):
		if self.expression is not None:
			return ""
		return Parameter.toStringH(self, indentCount)

class Var(ExpressionParameter):
	def __init__(self, id, value, valueExpression=None, rhsExpression=None, category = None):
		Parameter.__init__(self, id, value, category=category)
		self.symbol = symbols(id)
		self.expression = valueExpression
		self.rhs = rhsExpression
		self.precomputeDict = OrderedDict()

	#func should be a PrecomputeFunc enum
	def addPrecompute(self, func):
		#treating this like an ordered set, so just going to use the keys and ignore values
		self.precomputeDict[func] = None

class CppType(Enum):
	unsignedInt = "unsigned int"
	double = "double"
	string = "std::string"

class PrecomputeFunc(Enum):
	deriv_xx = "deriv_xx"
	deriv_yy = "deriv_yy"
	deriv_zz = "deriv_zz"
	def getArrayName(self, var):
		return self.value + "_" + var.id

	def getFunctionCallString(self, var):
		if self is PrecomputeFunc.deriv_xx:
			return "{0}({1}, {2}, hx, sz, bflag);".format(self.value, self.getArrayName(var), var.id)
		elif self is PrecomputeFunc.deriv_yy:
			return "{0}({1}, {2}, hy, sz, bflag);".format(self.value, self.getArrayName(var), var.id)
		elif self is PrecomputeFunc.deriv_zz:
			return "{0}({1}, {2}, hz, sz, bflag);".format(self.value, self.getArrayName(var), var.id)
		else:
			return ""


class PresetParam(Enum):
	#note values start out None, will be filled in before added to paramDict
	DIM = Parameter("DIM", None, description="dimentionality of the octree, (meshing is supported only for 3D)", cppType = CppType.unsignedInt)
	MAXDEPTH = Parameter("MAXDEPTH", None, description="maximum level of refinement of the mesh", cppType = CppType.unsignedInt)	
	ASYNC_COMM_K = Parameter("ASYNC_COMM_K", None, "variable group size for the asynchronous unzip operation", cppType = CppType.unsignedInt)
	CFL_FACTOR = Parameter("CFL_FACTOR", None, description="CFL Factor", cppType = CppType.double)
	GRID_MIN_X = Parameter("GRID_MIN_X", None, cppType = CppType.double)
	GRID_MIN_Y = Parameter("GRID_MIN_Y", None, cppType = CppType.double)
	GRID_MIN_Z = Parameter("GRID_MIN_Z", None, cppType = CppType.double)
	GRID_MAX_X = Parameter("GRID_MAX_X", None, cppType = CppType.double)
	GRID_MAX_Y = Parameter("GRID_MAX_Y", None, cppType = CppType.double)
	GRID_MAX_Z = Parameter("GRID_MAX_Z", None, cppType = CppType.double)
	COMPD_MIN = Parameter("COMPD_MIN", None, description="Should be an array defined from values of GRID_MIN", cppType = CppType.double)
	COMPD_MAX = Parameter("COMPD_MAX", None, description="Should be an array defined from values of GRID_MAX", cppType = CppType.double)
	OCTREE_MIN = Parameter("OCTREE_MIN", None, description="", cppType = CppType.double)
	OCTREE_MAX = Parameter("OCTREE_MAX", None, description="", cppType = CppType.double)
	RK45_TIME_STEP_SIZE = Parameter("RK45_TIME_STEP_SIZE", None, description="value should be calculated using the CFL Factor", cppType = CppType.double)
	RK45_DESIRED_TOL = Parameter("RK45_DESIRED_TOL", None, description="used in adaptive time stepping (not currently used)", cppType = CppType.double)
	KO_DISS_SIGMA = Parameter("KO_DISS_SIGMA", None, description="Kreiss-Oliger dissipation", cppType = CppType.double)
	ENABLE_BLOCK_ADAPTIVITY = Parameter("ENABLE_BLOCK_ADAPTIVITY", None, description="Set to 1 disable AMR and use block adaptivity (not recomended).", cppType = CppType.unsignedInt)
	BLK_MIN_X = Parameter("BLK_MIN_X", None, cppType = CppType.double)
	BLK_MIN_Y = Parameter("BLK_MIN_Y", None, cppType = CppType.double)
	BLK_MIN_Z = Parameter("BLK_MIN_Z", None, cppType = CppType.double)
	BLK_MAX_X = Parameter("BLK_MAX_X", None, cppType = CppType.double)
	BLK_MAX_Y = Parameter("BLK_MAX_Y", None, cppType = CppType.double)
	BLK_MAX_Z = Parameter("BLK_MAX_Z", None, cppType = CppType.double)
	NUM_REFINE_VARS = Parameter("NUM_REFINE_VARS", None, description="number of refinement variables", cppType = CppType.unsignedInt)
	REFINE_VARIABLE_INDICES = Parameter("REFINE_VARIABLE_INDICES", None, description="refinement variable IDs", cppType = CppType.unsignedInt)
	RK_MIN_TOL = Parameter("RK_MIN_TOL", None, description="RK solver tolerance", cppType = CppType.double)
	ELE_ORDER = Parameter("ELE_ORDER", None, description="element order", cppType = CppType.unsignedInt)
	DENDRO_GRAIN_SZ = Parameter("DENDRO_GRAIN_SZ", None, description="grain size N/p , Where N number of total octants, p number of active cores", cppType = CppType.unsignedInt)
	LOAD_IMB_TOL = Parameter("LOAD_IMB_TOL", None, description="dendro load imbalance tolerance for flexible partitioning", cppType = CppType.double)
	SPLIT_FIX = Parameter("SPLIT_FIX", None, description="Splitter fix value", cppType = CppType.unsignedInt)
	RK45_TIME_BEGIN = Parameter("RK45_TIME_BEGIN", None, description="simulation time begin", cppType = CppType.double)
	RK45_TIME_END = Parameter("RK45_TIME_END", None, description="simulation time end", cppType = CppType.double)
	RESTORE_SOLVER = Parameter("RESTORE_SOLVER", None, description="Set to 1 to restore solver from a checkpoint. 0 otherwise", cppType = CppType.unsignedInt)
	CHKPT_FILE_PREFIX = Parameter("CHKPT_FILE_PREFIX", None, description="file prefix for the checkpoint files", cppType = CppType.string)
	VTU_FILE_PREFIX = Parameter("VTU_FILE_PREFIX", None, description="file prefix for the vtu files", cppType = CppType.string)
	PROFILE_FILE_PREFIX = Parameter("PROFILE_FILE_PREFIX", None, description="file prefix for the intermediate profile files", cppType = CppType.string)
	IO_OUTPUT_FREQ = Parameter("IO_OUTPUT_FREQ", None, description="frequency for VTU output", cppType = CppType.unsignedInt)
	TIME_STEP_OUTPUT_FREQ = Parameter("TIME_STEP_OUTPUT_FREQ", None, description="", cppType = CppType.unsignedInt)
	REMESH_TEST_FREQ = Parameter("REMESH_TEST_FREQ", None, description="frequency for remeshing test based on wavelets", cppType = CppType.unsignedInt)
	CHECKPT_FREQ = Parameter("CHECKPT_FREQ", None, description="frequency for checkpoint output", cppType = CppType.unsignedInt)
	IO_OUTPUT_GAP = Parameter("IO_OUTPUT_GAP", None, description="VTU file output gap. (Not currently used. Might be usefull in adaptive timestepping)", cppType = CppType.unsignedInt)
	DENDRO_AMR_FAC = Parameter("DENDRO_AMR_FAC", None, description="dendro coarsening factor, corsent if computed wavelet tol < DENDRO_AMR_FAC*WAVELET_TOL ", cppType = CppType.double)
	NUM_EVOL_VARS_VTU_OUTPUT = Parameter("NUM_EVOL_VARS_VTU_OUTPUT", None, description="number of variables (evolution) to output in vtu files", cppType = CppType.unsignedInt)
	VTU_OUTPUT_EVOL_INDICES = Parameter("VTU_OUTPUT_EVOL_INDICES", None, description="evolution variable ids", cppType = CppType.unsignedInt)
