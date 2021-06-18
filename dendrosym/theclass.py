"""the_class.py

Temporary name of the file for the class that can be used
to interface with the code generation stuff.
"""

import sympy as sym
from sympy.core.symbol import var

import dendrosym


class DendroGRCodeGen:
    """Class for generating code files for Dendro computations


    """
    def __init__(self, project_short: str):

        # project short is for namespace info, basically like an abbreviation
        # for the entire project
        self.project_short = project_short.lower()
        self.upper_project_short = project_short.upper()

        self.var_types = ["VAR"]

        self.constraint_vars = []
        self.constraint_var_names = []
        self.other_vars = []
        self.other_var_names = []

        self.functions = {}

        self.advective_der_var = None

        self.parameter_vars = []

        pass

    def set_advective_der_var(self, in_var):
        self.advective_der_var = str(in_var).split("[")[0]

    def add_parameter_variables(self, in_vars: list):
        # TODO: clean up the parameters

        self.parameter_vars += in_vars

    def generate_code(self, arc_type="cpu"):

        if arc_type == "cpu":
            # generate for the CPU

            return_str = dendrosym.codegen.generate_cpu(
                *self.functions["other"](), '[pp]')

            pass

        elif arc_type == "gpu":
            # generate for the GPU

            pass

        return return_str

    def generate_var_enum(self):

        # start with the constraints
        return_str = ""

        return_str += dendrosym.codegen.gen_enum_info(
            self.constraint_var_names, "VAR_CONSTRAINT", "C") + "\n"

        return_str += dendrosym.codegen.gen_enum_info(self.other_var_names,
                                                      "VAR", "U")

        return return_str

    def generate_var_name_array(self):
        """Method that generates the var name array

        """

        return_str = ""

        return_str += dendrosym.codegen.gen_var_name_array(
            self.constraint_var_names, self.project_short, "C",
            "CONSTRAINT_VAR")

        return_str += dendrosym.codegen.gen_var_name_array(
            self.other_var_names, self.project_short, "U", "VAR")

        return return_str

    def generate_memory_alloc(self, var_type: str):
        """This generates the memory allocation code

        Please note that there needs to be an assigned function
        relating to this variable type for this to work
        """

        func = self.functions.get(var_type, None)

        if func is None:
            raise NotImplementedError(
                "Function for this var type not assigned")

        # get all of the possible gradient variables as uniques
        # TODO: I'm thinking we just save the text file of results
        # and then read that in to find the gradient vars
        # TODO: fix this!
        with open("test1.txt", 'r') as fh:
            in_text = fh.read()

        grad_list = dendrosym.utils.find_derivtype_in_text(in_text, "grad_")
        grad2_list = dendrosym.utils.find_derivtype_in_text(in_text, "grad2_")
        agrad_list = dendrosym.utils.find_derivtype_in_text(in_text, "agrad_")

        return_str = ""

        return_str += dendrosym.codegen.generate_deriv_comp(grad_list)
        return_str += dendrosym.codegen.generate_deriv_comp(grad2_list)
        return_str += dendrosym.codegen.generate_deriv_comp(
            agrad_list, self.advective_der_var)

        return return_str

    def assign_rhs_function(self, func, var_type: str = "other"):

        self.functions[var_type] = func

    def generate_rhs_vars(self,
                          var_type: str = "other",
                          zip_var_name="unzipVarsRHS"):

        return_str = ""

        func = self.functions.get(var_type, None)

        if func is None:
            raise ValueError("Sorry, but that function hasn't been assigned")

        if var_type == "other":
            expre, var_names = func()
            _, exp_names, _ = dendrosym.codegen.construct_expression_list(
                expre, var_names, "")

            # now I gotta match up the lists to be sure
            other_var_names = self.other_var_names.copy()
            other_var_names.sort()
            other_var_names = list(map(lambda x: x.upper(), other_var_names))
            exp_names.sort()

            return_str = dendrosym.codegen.gen_var_info(
                exp_names,
                zip_var_name=zip_var_name,
                offset_name="offset",
                enum_name="VAR",
                enum_prefix="U",
                use_const=False,
                enum_var_names=other_var_names)
        else:
            raise NotImplementedError("Fix this")

        return return_str

    def generate_variable_extraction(self,
                                     var_type="other",
                                     use_const=False,
                                     zip_var_name="uZipVars"):
        return_str = ""
        if var_type == "other":
            return_str += dendrosym.codegen.gen_var_info(
                self.other_var_names,
                zip_var_name=zip_var_name,
                enum_prefix="U",
                use_const=use_const)
        else:
            raise NotImplementedError("Not yet implemented")

        return return_str

    def generate_phys_constraint_vars(self):

        return_str = ""

        # start with the physical constraint vars
        return_str += dendrosym.codegen.gen_var_info(
            self.constraint_var_names,
            zip_var_name="uZipConVars",
            enum_name="VAR_CONSTRAINT",
            enum_prefix="C",
            use_const=False)

        return return_str

    def generate_parameter_code(self):

        return_str = ""

        for param in self.parameter_vars:
            return_str += param.generate_cpp_line(
                global_param_prefix=self.upper_project_short)
            return_str += "\n"

        return return_str

    def generate_ko_calculations(self,
                                 var_type="other",
                                 ko_sigma_name="sigma",
                                 idx="pp"):

        return_str = ""

        if var_type == "other":
            for var_name in self.other_var_names:
                temp_str = f"{var_name}_rhs[{idx}]"
                temp_str += f" += {ko_sigma_name} * ("
                temp_str += " + ".join(f"grad_{ii}_{var_name}[{idx}]"
                                       for ii in range(3))
                temp_str += ");\n"

                return_str += temp_str

        return return_str

    def add_var_type(self, var_prefix):
        """Add a variable type to Dendro computations

        This method allows you to add another subset of variable type
        to be used in the computations. Meaning you can split variables
        into different portions based on constraints or something similar.

        Note that it will force this to be all caps because this will be
        used in global variable definitions.
        """

        self.var_types.append(var_prefix)

    def add_constraint_var(self, the_var, add_multi=False):
        """Add a single variable to the collection of constraints

        This includes variables such as Ham, Mom, Psi4, etc.
        """

        if not add_multi:
            # make it a list if not incoming
            the_var = [the_var]

        # iterate through the list
        for va in the_var:
            if type(va) not in [sym.Symbol, sym.Matrix, tuple]:
                raise NotImplementedError("Invalid variable type")

            self.constraint_vars.append(va)

        self._create_raw_strings_for_vars()

    def add_var(self, the_var, add_multi=False):
        """Add a single variable to the collection of variables

        This is for variables not considered constraints in the
        dendro GR code.

        This includes variables such as alpha, chi, k, Gt, beta, Beta,
        and others.
        """

        if not add_multi:
            the_var = [the_var]

        for va in the_var:
            if type(va) not in [sym.Symbol, sym.Matrix, tuple]:
                raise NotImplementedError("Invalid variable type")

            self.other_vars.append(va)

        self._create_raw_strings_for_vars()

    def _create_raw_strings_for_vars(self):

        self.constraint_var_names = self.clean_var_names(self.constraint_vars)
        self.other_var_names = self.clean_var_names(self.other_vars)

    @staticmethod
    def clean_var_names(in_vars):
        """Create list of strings for the input variables

        """
        # first find all strings that correspond to the vars
        all_var_names = []
        for the_var in in_vars:
            if type(the_var) is sym.Matrix:
                unique_vars = []
                for ii in range(the_var.shape[0]):
                    for jj in range(the_var.shape[1]):
                        if str(the_var[ii, jj]) not in unique_vars:
                            unique_vars.append(str(the_var[ii, jj]))
                all_var_names += unique_vars
            elif type(the_var) is tuple or type(the_var) is list:
                all_var_names += [str(the_name) for the_name in the_var]
            elif type(the_var) is sym.Symbol:
                all_var_names.append(str(the_var))
            else:
                raise NotImplementedError(
                    "That variable type isn't implemented yet")

        # now we clean away potential [idx] information
        for ii in range(len(all_var_names)):
            curr_var_name = all_var_names[ii]
            if "[" in curr_var_name:
                curr_var_name = curr_var_name.split("[")[0]
                all_var_names[ii] = curr_var_name

        return all_var_names
