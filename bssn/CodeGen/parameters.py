import math
from collections import OrderedDict

class Parameters:
	def __init__(self, paramDict=None):
		if paramDict is None:
			paramDict = OrderedDict()
		self.paramDict = paramDict
		self.indent = "\t"
	
	def add(self, id, value, description=None, category=None):
		param = Parameter(id, value, description, category)
		self.paramDict[param.id] = param
		
	def write(self, filepath):
		categories = OrderedDict()
		for parameter in self.paramDict.values():
			if parameter.category not in categories:
				categories[parameter.category] = ""
			categories[parameter.category] += parameter.toString(2)
		
		out = open(filepath, "w")
		out.write("{\n")
		
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
				out.write('{0}"__comment__" : "{1}"'.format(self.indent, "=" * 50))
		
		out.write("\n}")
		out.close()

class Parameter:
	def __init__(self, id, value, description=None, category=None):
		self.id = id
		self.value = value
		self.description = description
		self.category = category

	def toString(self, indentCount=0):
		indent = "\t"
		output = ""
		if (self.description is not None):
			output += '{0}"           " : "{1}",\n'.format(indent * indentCount, self.description)
		output += '{0}"{1}" : {2},\n'.format(indent * indentCount, self.id,
											 '"' + self.value + '"' if isinstance(self.value, str) else self.value)

		return output
