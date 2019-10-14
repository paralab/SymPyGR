import math
from enum import Enum
import warnings

import sys    
import os
import datetime
import ntpath

from collections import OrderedDict

#holds dict of parameters, contains mostly output functions
class Parameters:
	def __init__(self, paramDict=None):
		if paramDict is None:
			paramDict = OrderedDict()
		self.paramDict = paramDict
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
			param.setValue(value, cppType, arraySize)
			#note the dict key will be the PresetParam enum, not a string
			self.paramDict[id] = param
		else:
			param = Parameter(id, value, description, category, cppType, arraySize)
			self.paramDict[param.id] = param

	def setCategory(self, category):
		self.category = category

	def clearCategory(self):
		self.category = None

	def values(self):
		return self.paramDict.values()

	def writeCpp(self, hPath, cppPath, namespace):

		hName = ntpath.basename(hPath)
		cppName = ntpath.basename(cppPath)

		# Initialize output and first lines for .h file
		outH = open(hPath, "w")
		outH.write("//" + self.generationNote + "\n")
		outH.write("#ifndef SFCSORTBENCH_PARAMETERS_H\n")
		outH.write("#define SFCSORTBENCH_PARAMETERS_H\n\n")

		outH.write("#include <string.h>\n")
		outH.write("#include <iostream>\n\n")

		outH.write("namespace {0}\n{{\n\n".format(namespace))

		# Initialize output and first lines for .cpp file.
		outCpp = open(cppPath, "w")
		outCpp.write("//" + self.generationNote + "\n")

		outCpp.write("#include \"{0}\"\n".format(hName))
		outCpp.write("namespace {0}\n{{\n\n".format(namespace))

		for parameter in self.paramDict.values():
			if parameter.isComment():
				continue
			outH.write(parameter.toStringH(1))
			outCpp.write(parameter.toStringCpp(1))

		outH.write("}\n")
		outH.write("#endif //SFCSORTBENCH_PARAMETERS_H")
		outH.close()

		outCpp.write("}")
		outCpp.close()

	def writeJson(self, filepath):
		categories = OrderedDict()
		for parameter in self.paramDict.values():
			if parameter.category not in categories:
				categories[parameter.category] = ""
			categories[parameter.category] += parameter.toStringJson(2)

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
			elif isinstance(value, float):
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
	def toStringH(self, indentCount = 0):
		output = ""
		if self.description is not None:
			output += "{0}//{1}\n".format(self.indent * indentCount, self.description)

		output += "{0}extern {1} {2}{3};\n\n".format("\t" * indentCount, self.cppType.value, self.id,
											  "[" + str(self.arraySize) +"]" if self.arraySize is not None else "")
		return output

	def toStringCpp(self, indentCount = 0):
		output = ""
		if self.description is not None:
			output += "{0}//{1}\n".format(self.indent * indentCount, self.description)

		# surround value with quotes if it's a string, otherwise nothing
		valueQuote = '"' if self.cppType is CppType.string else ""
		if self.arraySize is None:
			output += "{0}{1} {2}={4}{3}{4};\n\n".format(self.indent * indentCount, self.cppType.value, self.id, self.value, valueQuote)
		else:
			arrayValue = valueQuote + "{0},{0}".format(valueQuote).join(str(v) for v in self.value) + valueQuote
			output += "{0}{1} {2}[{3}]={{{4}}};\n\n".format(self.indent * indentCount, self.cppType.value, self.id, self.arraySize, arrayValue)
			
		return output

class CppType(Enum):
	unsignedInt = "unsigned int"
	double = "double"
	string = "std::string"

class PresetParam(Enum):
	#note values start out None, will be filled in before added to paramDict
	NUM_VARS = Parameter("NUM_VARS", None, cppType = CppType.unsignedInt)
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
	WAVELET_TOL = Parameter("WAVELET_TOL", None, description="wavelet tolerance", cppType = CppType.double)
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
	REMESH_TEST_FREQ = Parameter("REMESH_TEST_FREQ", None, description="frequency for remeshing test based on wavelets", cppType = CppType.unsignedInt)
	CHECKPT_FREQ = Parameter("CHECKPT_FREQ", None, description="frequency for checkpoint output", cppType = CppType.unsignedInt)
	IO_OUTPUT_GAP = Parameter("IO_OUTPUT_GAP", None, description="VTU file output gap. (Not currently used. Might be usefull in adaptive timestepping)", cppType = CppType.unsignedInt)
	DENDRO_AMR_FAC = Parameter("DENDRO_AMR_FAC", None, description="dendro coarsening factor, corsent if computed wavelet tol < DENDRO_AMR_FAC*WAVELET_TOL ", cppType = CppType.double)
	NUM_EVOL_VARS_VTU_OUTPUT = Parameter("NUM_EVOL_VARS_VTU_OUTPUT", None, description="number of variables (evolution) to output in vtu files", cppType = CppType.unsignedInt)
	VTU_OUTPUT_EVOL_INDICES = Parameter("VTU_OUTPUT_EVOL_INDICES", None, description="evolution variable ids", cppType = CppType.unsignedInt)
