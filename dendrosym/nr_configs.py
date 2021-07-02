"""nr_configs.py

This file contains the configuration class to generate
and store the necessary pieces of information regarding
the set-up and configuration for a numerical relativity
project using the Dendro framework.
"""

import re

from numpy import single
import dendrosym
import sympy as sym


class ImproperInitalization(Exception):
    pass


class AdvectiveDerivativeNotSet(Exception):
    pass


class NRConfig:
    """A class to store numerical relativity configurations

    In particular, this class can store infomration about
    the numerical relativity project such as RHS equations,
    store variables, and then generate the C++ code
    from the symbolic equations.
    """
    def __init__(self, project_name: str):
        self.project_name = project_name
        self.project_upper = project_name.upper()

        # make all vars a dictionary in case there are
        # other types to store, but by default we will have
        # "constraint", "evolution", and "parameter" which are initialized
        # by default
        self.all_vars = {
            "constraint": [],
            "evolution": [],
            "parameter": {
                "constraint": [],
                "evolution": []
            }
        }

        # also create the same type of dictionary but
        # just for string representations
        self.all_var_names = {
            "constraint": [],
            "evolution": [],
            "parameter": {
                "constraint": [],
                "evolution": []
            }
        }

        self.enum_prefixes = {"constraint": "C", "evolution": "U"}

        # then also store the functions that return the rhs
        # portions
        self.all_rhs_functions = {"constraint": None, "evolution": None}

        # initialize the advective derivative variable as well
        self.advective_der_var = None

        self.idx_str = ""

        self.bcs_info = {"constraint": {}, "evolution": {}}

        self.metric_var = None

        self.evolution_constraint_info = {"trace_zero": [], "pos_floor": []}

    def add_parameter_variables(self,
                                in_vars: list,
                                eqn_type: str = "evolution"):
        """Add a parameter variable to the list

        Use this function when there is a constant parameter variable
        that you have in your equations. Be sure to input
        the type of equation that it belongs to. Currently
        accepts "evolution" or "constraint"
        """

        if type(in_vars) is not list:
            in_vars = [in_vars]

        if eqn_type not in ["evolution", "constraint"]:
            raise ValueError("Parameter variables can only be assigned" +
                             " to evolution and constraint")

        for one_var in in_vars:
            # get the variable name from the object
            var_names_clean = one_var.var_name

            if not self.check_repeat_items_ignore_case(
                    var_names_clean,
                    self.all_var_names["parameter"][eqn_type]):
                self.all_var_names["parameter"][eqn_type].append(
                    var_names_clean)
            else:
                raise Exception("That parameter has already been added")
            # then append the paramter to the list
            self.all_vars["parameter"][eqn_type].append(one_var)

    def add_evolution_variables(self, in_vars: list):

        if type(in_vars) is not list:
            in_vars = [in_vars]

        for one_var in in_vars:
            # get the clean var name
            var_names_clean = self.clean_var_names([one_var])

            for vname_clean in var_names_clean:
                if not self.check_repeat_items_ignore_case(
                        vname_clean, self.all_var_names["evolution"]):
                    self.all_var_names["evolution"].append(vname_clean)
                else:
                    raise Exception(
                        f"'{vname_clean}' has already been added or is" +
                        " too similar to another existing one")
            # then append the paramter to the list
            self.all_vars["evolution"].append(one_var)

    def add_constraint_variables(self, in_vars: list):

        if type(in_vars) is not list:
            in_vars = [in_vars]

        for one_var in in_vars:
            # get the clean var name
            var_names_clean = self.clean_var_names([one_var])

            for vname_clean in var_names_clean:
                if not self.check_repeat_items_ignore_case(
                        vname_clean, self.all_var_names["constraint"]):
                    self.all_var_names["constraint"].append(vname_clean)
                else:
                    raise Exception(
                        f"'{vname_clean}' has already been added or is" +
                        " too similar to another existing one")
            # then append the paramter to the list
            self.all_vars["constraint"].append(one_var)

    def set_advective_derivative_var(self, in_var):
        self.advective_der_var = str(in_var).split("[")[0]

    def set_idx_str(self, idx_str):
        self.idx_str = idx_str

    def set_metric(self, in_var):
        assert in_var in self.all_vars.get(
            "evolution",
            []), "Incoming metric variable is not in evolution variables"
        self.metric_var = in_var
        # also set the metric in the NR code for use there
        dendrosym.nr.set_metric(in_var)

    def get_rhs_var_names(self, var_type: str):
        """Creates the RHS variable names for a variable type

        This does not generate C++ code, just creates the list
        of rhs variable names.

        This also does not guarantee that all RHS equations
        will line up, it just generates the names of the RHS
        variables associated with the variables provided during
        configuration of the project.
        """

        if var_type not in self.all_var_names.keys():
            raise ValueError("Invalid variable type in RHS Var Name Creation")

        if var_type == "parameter":
            raise ValueError("RHS Variables Cannot Exist for Parameters")

        vars_use = self.all_var_names.get(var_type)
        rhs_vars = []

        for the_var in vars_use:
            rhs_vars.append(the_var + "_rhs")

        return rhs_vars

    def generate_rhs_code(self, var_type: str, arc_type="cpu"):

        if var_type not in ["evolution", "constraint"]:
            raise ValueError(f"Unfortunately {var_type} doesn't work yet")

        # so now we can get started by getting the rhs information
        all_exp, all_rhs_names, orig_n_exp = self._extract_rhs_expressions(
            var_type)

        # count the number of original operations on the expressions
        orig_n_ops = sym.count_ops(all_exp)

        # construct cse from the list
        cse_exp = dendrosym.codegen.construct_cse_from_list(all_exp)

        if arc_type == "cpu":
            # then we're to send this information to the cpu gen function
            return dendrosym.codegen.generate_cpu_preextracted(
                cse_exp, all_rhs_names, self.idx_str, orig_n_ops)

        else:
            raise NotImplementedError(
                "That archetecture generation code isn't ready yet")

    def set_bhs_falloff_and_asymptotic(self, var_type, var_list, var_info):
        """Sets the falloff and asymptotic values to the individual variables

        To use this function, say what type of variables you're providing information
        for, and then pass through a dictionary of values for each of the varibles
        you are assigning the data to. That dictionary should be of structure

        var_information = {beta: {"falloff": 1.0, "asymptotic": 2.0}}

        and include all variables through the dictionary construction. This
        function will throw an error if not all variables are accounted for.
        """

        if var_type == "parameter":
            raise ValueError("Cannot set BHS data to parameters")

        temp_vars = self.all_vars.get(var_type, []).copy()

        if len(temp_vars) == 0:
            raise ImproperInitalization(
                f"{var_type} does not have variables assigned to it.")

        # let's make sure that they're all there
        for ii, var_ in enumerate(var_list):

            found_match = False

            # then we have to iterate through our temp_vars list and find a match
            for temp_var in temp_vars:
                if temp_var == var_:
                    # say we found a match, and then remove it
                    found_match = True
                    temp_vars.remove(temp_var)

                    break

            if not found_match:
                raise ImproperInitalization(
                    f"Couldn't find identical variable {var_} in stored list")

        # now make sure that we've cleared both lists
        if len(temp_vars) != 0:
            raise ImproperInitalization(
                "Not all assigned variables to '" + var_type +
                "' were included in RHS function. Remaining variables were:" +
                repr(temp_vars))

        self.bcs_info.update({var_type: {"vars": var_list, "info": var_info}})

    def generate_bcs_calculations(self,
                                  var_type="evolution",
                                  pmin="pmin",
                                  pmax="pmax",
                                  sz="sz",
                                  bflag="bflag"):
        """Generate the BCS calculation functions

        This one is tricky as different variables require using
        """

        # just iterate through the data we have stored
        all_var_info = self.bcs_info.get(var_type, None)

        if all_var_info is None:
            raise ImproperInitalization(
                f"BCS information was not initialized for {var_type}.")

        if not all_var_info:
            raise ImproperInitalization(
                f"BCS information was not initialized for {var_type}.")

        return_str = ""

        for ii, the_var in enumerate(all_var_info["vars"]):
            var_info = all_var_info["info"][ii]
            cleaned_var_name = self.clean_var_names([the_var])

            if len(cleaned_var_name) > 1:

                for clean_var in cleaned_var_name:
                    rhs_var = self.generate_rhs_var_names([clean_var])

                    grad_vars = self.create_grad_var_names([clean_var], "grad",
                                                           3)

                    if len(var_info) == 2 and (type(var_info[0]) is int
                                               or type(var_info[0]) is float):
                        falloff = var_info[0]
                        asymptotic = var_info[1]
                    else:
                        idxs = self.get_indices_from_var_name(clean_var)

                        if len(idxs) == 1:
                            falloff = var_info[idxs[0]][0]
                            asymptotic = var_info[idxs[0]][1]
                        else:
                            falloff = var_info[idxs[0]][idxs[1]][0]
                            asymptotic = var_info[idxs[0]][idxs[1]][1]

                    return_str += dendrosym.codegen.generate_bcs_function_call(
                        rhs_var[0], clean_var, falloff, asymptotic, grad_vars,
                        self.project_name, pmin, pmax, sz, bflag)

                pass

            else:
                rhs_var = self.generate_rhs_var_names(cleaned_var_name)

                grad_vars = self.create_grad_var_names(cleaned_var_name,
                                                       "grad", 3)
                return_str += dendrosym.codegen.generate_bcs_function_call(
                    rhs_var[0], cleaned_var_name[0], var_info[0], var_info[1],
                    grad_vars, self.project_name, pmin, pmax, sz, bflag)

        return return_str

    @staticmethod
    def get_indices_from_var_name(var_name):

        var_digits = re.findall('[0-9]+', var_name)
        var_end = var_digits[-1]

        if len(var_end) == 1:
            return (int(var_end), )

        # take the length of this string and divide it by two for how many
        # characters each index is
        assert len(var_end) % 2 == 0, "incoming string not divisible by two"
        len_idx = int(len(var_end) / 2)

        # we're only dealing with 2D for now
        return (int(var_end[:len_idx]), int(var_end[len_idx:]))

    def generate_variable_extraction(self,
                                     var_type="evolution",
                                     use_const=True,
                                     zip_var_name="uZipVars",
                                     dtype="double"):
        """Generate the C++ Code for Variable Extraction

        The variable extraction refers to the zipping and
        unzipping that happens with a variable like 'uZipVars'.
        This function can generate all of the code necessary
        to create pointers to where that data is stored.
        """

        if var_type == "parameter":
            raise ValueError("Parameters are generated with another function")

        named_enums = self.get_enum_var_names(var_type)

        if var_type == "evolution":
            enum_name = "VAR"
        elif var_type == "constraint":
            enum_name = "VAR_CONSTRAINT"
        else:
            raise NotImplementedError("Not Yet Implemented")

        return_str = dendrosym.codegen.gen_var_info(self.all_var_names.get(
            var_type, []),
                                                    zip_var_name=zip_var_name,
                                                    use_const=use_const,
                                                    enum_name=enum_name,
                                                    enum_var_names=named_enums,
                                                    dtype=dtype)

        return return_str

    def generate_rhs_var_extraction(self,
                                    var_type="evolution",
                                    use_const=False,
                                    zip_var_name="unZipVarsRHS",
                                    dtype="double"):
        """Generate the C++ Code for Variable Extraction

        The variable extraction refers to the zipping and
        unzipping that happens with a variable like 'uZipVars'.
        This function can generate all of the code necessary
        to create pointers to where that data is stored.
        """

        if var_type == "parameters":
            raise ValueError("Parameters are generated with another function")

        named_enums = self.get_enum_var_names(var_type)
        vars_use = self.get_rhs_var_names(var_type)

        return_str = dendrosym.codegen.gen_var_info(vars_use,
                                                    zip_var_name=zip_var_name,
                                                    use_const=use_const,
                                                    dtype=dtype,
                                                    enum_var_names=named_enums)

        return return_str

    def gen_enum_code(self,
                      var_type: str,
                      enum_start_idx: int = 0,
                      enum_name: str = "VAR"):

        if var_type == "parameter":
            raise NotImplementedError("Not yet implemented")

        # get the enum names first
        enum_names = self.get_enum_var_names(var_type)

        enum_text = f"enum {enum_name}\n{{\n"

        for ii, enum_name in enumerate(enum_names):
            enum_line = f"    {enum_name}"
            enum_line += f" = {enum_start_idx}" if ii == 0 else ""
            enum_line += ",\n" if ii != len(enum_names) - 1 else "\n"
            enum_text += enum_line

        return enum_text + "};"

    def gen_enum_names(self, var_type: str, enum_name: str = "VAR"):

        if var_type == "parameter":
            raise ValueError("This function can't handle parameter enum names")

        enum_names = self.get_enum_var_names(var_type)

        return_str = dendrosym.codegen.gen_var_name_array(
            enum_names, self.project_upper, enum_name)

        return return_str

    def gen_parameter_code(self, var_type="evolution"):

        return_str = ""

        for param in self.all_vars["parameter"].get(var_type, []):
            return_str += param.generate_cpp_line(
                global_param_prefix=self.project_upper)
            return_str += "\n"

        return return_str

    def get_enum_var_names(self, var_type: str):

        if var_type == "parameter":
            raise NotImplementedError("Not yet ready")

        orig_vars = self.all_var_names.get(var_type, [])

        out_enum_names = []
        enum_prefix = self.enum_prefixes.get(var_type, "")

        for orig_var in orig_vars:
            out_enum_names.append(f"{enum_prefix}_{orig_var.upper()}")

        return out_enum_names

    def set_rhs_equation_function(self, var_type: str, rhs_func):

        if var_type == "parameter":
            raise ValueError("Cannot set RHS function to parameters")

        # first we need to check that the function gives us what we want
        # so we need to evaluate it first
        rhs_list, var_list = rhs_func()

        if len(rhs_list) != len(var_list):
            raise ImproperInitalization(
                "The RHS function does not have the same number" +
                " of expressions as variables")

        # so, now we make a copy of our variables for this variable type
        temp_vars = self.all_vars.get(var_type, []).copy()

        if len(temp_vars) == 0:
            raise ImproperInitalization(
                f"{var_type} does not have variables assigned to it.")

        # we can't change the order, but we need to make sure they're all there
        # so we traverse the variable list
        for ii, var_ in enumerate(var_list):

            found_match = False

            # then we have to iterate through our temp_vars list and find a match
            for temp_var in temp_vars:
                if temp_var == var_:
                    # say we found a match, and then remove it
                    found_match = True
                    temp_vars.remove(temp_var)
                    break

            if not found_match:
                raise ImproperInitalization(
                    f"Couldn't find identical variable {var_} in stored list")

        # now make sure that we've cleared both lists
        if len(temp_vars) != 0:
            raise ImproperInitalization(
                "Not all assigned variables to '" + var_type +
                "' were included in RHS function. Remaining variables were:" +
                repr(temp_vars))

        # if that's all good, then we're good to add the function to our list
        self.all_rhs_functions.update({var_type: rhs_func})

    def _extract_rhs_expressions(self, var_type: str):
        """An internal function that extracts the expressions for a variable type

        It does this to match up expressions with their corresponding RHS variables
        and put the list in order. This is because the class manages the
        creation of RHS and grad variables.
        """

        if var_type == "parameter":
            raise ValueError("Cannot extract expressions from parameters")

        rhs_func = self.all_rhs_functions.get(var_type, None)

        if rhs_func is None:
            raise ImproperInitalization(
                f"A RHS function was not assigned to {var_type}")

        # now we need the variables, and the good news is that we know that the
        # function we extracted can't exist without having gone through
        # our filtering function to make sure all variables are present
        # and accounted for.
        # TODO: one problem would be adding more variables after setting the
        # function... need to fix that

        # evaluate the function so we can get the pieces
        rhs_list, var_list = rhs_func()

        # now we iterate through it and put together all expressions and var names
        # in order
        all_expressions = []
        all_rhs_var_names = []
        original_number_expressions = 0

        for ii in range(len(rhs_list)):
            expression = rhs_list[ii]
            the_var = var_list[ii]

            list_expressions, num_e = dendrosym.codegen.extract_expression(
                expression)

            # note: if we have a sym.Matrix as our variable, we need to then
            # keep in mind the indexing so we can build the RHS variables
            rhs_vars = self.generate_rhs_var_names(
                self.clean_var_names([the_var]))

            all_expressions += list_expressions
            all_rhs_var_names += rhs_vars
            original_number_expressions += num_e

        # for ii in range(len(all_expressions)):
        #     print(all_rhs_var_names[ii], "=", all_expressions[ii])

        # okay, now that we have them we can return them
        return all_expressions, all_rhs_var_names, original_number_expressions

    def gen_grad_memory_alloc(self, var_type: str, grad_type: str = "grad"):
        """Generates memory allocation code

        This method takes the internally-stored set of variables as
        given by the user and generates the C++ code to allocate
        memory for their calculation. However, do note that it will
        generate allocation code for every single variable in the list,
        which may not be necessary depending on the equations that they
        will be used and called from.

        See `gen_grad_memory_alloc_from_code` for an alternative
        method to generating this memory allocation code.
        """

        if var_type == "parameter":
            raise ValueError("Parameters cannot use gradients")

        orig_vars = self.all_var_names.get(var_type, [])

        grad_vars = self.create_grad_var_names(orig_vars, grad_type, 3)

        return_text = dendrosym.codegen.generate_memory_alloc(
            grad_vars, "double", grad_type == "grad")

        return return_text

    @staticmethod
    def gen_grad_memory_alloc_from_code(filename: str,
                                        grad_type: str = "grad",
                                        include_byte_declaration: bool = False,
                                        dtype: str = "double"):
        """Alternative method to generating grad memory alloc

        Instead of generating the variable list from internally stored
        data, this method will take a pre-generated file for a set of
        equations and find all instances of the derivatives that start
        with the 'grad_type' prefix included.

        This might be preferable to the method that uses internally stored
        configurations of variables, simply because it will only allocate
        the variables that are absolutely needed. It also sorts them
        by alphabetical order.
        """

        with open(filename, "r") as fh:
            in_text = fh.read()

        # first get the grad list from the "find deriv type in text" function
        grad_list = dendrosym.utils.find_derivtype_in_text(
            in_text, grad_type + "_")

        # sort it for clarity
        grad_list.sort()

        if len(grad_list) == 0:
            # TODO: make sure to actually call some kind of logger
            # or warning system
            print("WARNING: no variables found with " +
                  f"'{grad_type}' string in file {filename}")

        # then we can create the grad variable names directly
        return_text = dendrosym.codegen.generate_memory_alloc(
            grad_list,
            var_type=dtype,
            include_byte_declaration=include_byte_declaration)

        return return_text

    def gen_grad_memory_dealloc(self, var_type: str, grad_type: str = "grad"):

        if var_type == "parameter":
            raise ValueError("Parameters cannot use gradients")

        orig_vars = self.all_var_names.get(var_type, [])

        grad_vars = self.create_grad_var_names(orig_vars, grad_type, 3)

        return_text = dendrosym.codegen.generate_memory_dealloc(grad_vars)

        return return_text

    @staticmethod
    def gen_grad_memory_dealloc_from_code(filename: str,
                                          grad_type: str = "grad"):

        with open(filename, "r") as fh:
            in_text = fh.read()

        # first get the grad list from the "find deriv type in text" function
        grad_list = dendrosym.utils.find_derivtype_in_text(
            in_text, grad_type + "_")

        # sort it for clarity
        grad_list.sort()

        if len(grad_list) == 0:
            # TODO: make sure to actually call some kind of logger
            # or warning system
            print("WARNING: no variables found with " +
                  f"'{grad_type}' string in file {filename}")

        # then we can create the grad variable names directly
        return_text = dendrosym.codegen.generate_memory_dealloc(grad_list)

        return return_text

    def gen_grad_calculations(self, var_type: str, grad_type: str = "grad"):

        orig_vars = self.all_var_names.get(var_type, [])

        grad_vars = self.create_grad_var_names(orig_vars, grad_type, 3)

        return_text = dendrosym.codegen.generate_deriv_comp(grad_vars)

        return return_text

    def generate_ko_derivs(self,
                           var_type: str,
                           sz="sz",
                           bflag="bflag",
                           ko_func_name="ko_deriv"):

        if var_type == "parameter":
            raise ValueError("Cannot generate KO Calculations with Parameters")

        return_text = ""
        grad_dirs = ["x", "y", "z"]

        # get the variables
        var_list = self.all_var_names.get(var_type, [])

        for ii, var_name in enumerate(var_list):
            grad_vars = self.create_grad_var_names([var_name])

            for jj, grad_var in enumerate(grad_vars):
                g_dir = grad_dirs[jj]
                return_text += f"{ko_func_name}_{g_dir}("
                return_text += f"{grad_var}, {var_name}, h{g_dir}"
                return_text += f", {sz}, {bflag});\n"

        return return_text

    def generate_ko_calculations(self, var_type: str, ko_sigma_name="sigma"):

        return_str = ""

        if var_type == "parameter":
            raise ValueError("Cannot generate KO Calculations with Parameters")

        var_list = self.all_var_names.get(var_type, [])
        rhs_var_list = self.get_rhs_var_names(var_type)

        for ii, var_name in enumerate(var_list):
            temp_str = f"{rhs_var_list[ii]}{self.idx_str}"
            temp_str += f" += {ko_sigma_name} * ("
            temp_str += " + ".join(f"grad_{ii}_{var_name}{self.idx_str}"
                                   for ii in range(3))
            temp_str += ");\n"

            return_str += temp_str

        return return_str

    def add_evolution_constraint(self, var_adjust, constraint_info):
        """Adds an evolution constraint of some kind

        Takes in the input variable that needs to have some kind of constraint
        and then information about the constraint. The constraint info
        currently supports 'trace_zero', 'pos_floor'

        Trace of zero requires the metric to have been set, and the
        generation code will ensure that it has been properly extracted
        and updated due to the metric's determinant being 1.

        Positive threshold only works for single values and requires
        a parameter be set *elsewhere* (in parameters.cpp for example)
        with the name "{VAR_NAME}_FLOOR" where 'VAR_NAME' is the same
        name that was assigned to the variable but in all capitals.
        """

        if constraint_info == "trace_zero":
            assert var_adjust in self.all_vars.get(
                "evolution", []), "Variable not in set evolution variables"
            
            if type(var_adjust) is not sym.Matrix:
                raise ValueError("For zero trace, the incoming variable needs to be a matrix")
            
            self.evolution_constraint_info["trace_zero"].append(var_adjust)
        elif constraint_info == "pos_floor":
            assert var_adjust in self.all_vars.get(
                "evolution", []), "Variable not in set evolution variables"
            
            if type(var_adjust) is not sym.Symbol:
                raise ValueError("For positive floor, the variable needs to not have dimensions")
            
            self.evolution_constraint_info["pos_floor"].append(var_adjust)
        else:
            raise NotImplementedError("This type of constraint has not been implemented")


    def generate_evolution_constraints(self):

        if self.metric_var is None:
            raise ImproperInitalization(
                "Metric was not set, cannot generate evolution constraints")

        return_str = ""

        metric_var_names = self.clean_var_names([self.metric_var])

        metric_name = metric_var_names[0][:-2]

        # then we also need to get the metric enum information
        enum_prefix = self.enum_prefixes.get("evolution", "")
        metric_enums = [
            f"VAR::{enum_prefix}_{mvar.upper()}" for mvar in metric_var_names
        ]

        # first things first, the metric determinant needs to be zero,
        # so we generate that code first
        return_str += "////\n//// CODE TO REQUIRE METRIC DETERMINANT TO BE ONE\n"
        return_str += dendrosym.codegen.generate_update_sym_mat_extract(
            metric_name, metric_enums)

        # then we need to generate the determinant code
        return_str += dendrosym.codegen.generate_force_sym_matrix_det_to_one(
            metric_name, metric_enums)

        # now we can start adding our other constraints
        for var_use in self.evolution_constraint_info.get("trace_zero", []):

            return_str += "////\n////APPLYING TRACE ZERO TO "
            # for trace of zero, we need to get all of the variables
            var_names = self.clean_var_names([var_use])
            # then get the name of the variable
            single_var_name = var_names[0][:-2]
            return_str += single_var_name + "\n"
            # then generate the enums
            var_enums = [
                f"VAR::{enum_prefix}_{ivar.upper()}" for ivar in var_names
            ]

            # then add the code to extract the sym matrix
            return_str += dendrosym.codegen.generate_update_sym_mat_extract(
                single_var_name, var_enums)

            # then the code to force the traceless symmat
            return_str += dendrosym.codegen.generate_force_symmat_traceless(
                single_var_name, metric_name
            )

            # then the code to actually update the matrix
            return_str += dendrosym.codegen.generate_update_sym_mat_code(
                single_var_name, var_enums, include_up=False)
            return_str += "//// END APPLICATION OF TRACE ZERO\n////\n\n"
        
        for var_use in self.evolution_constraint_info.get("pos_floor"):

            return_str += "////\n////APPLYING POSITIVE FLOOR TO "
            single_var_name = self.clean_var_names([var_use])[0]
            return_str += single_var_name + "\n"

            var_enum = f"VAR::{enum_prefix}_{single_var_name.upper()}"

            return_str += dendrosym.codegen.generate_variable_always_positive(
                var_enum, floor_var=f"{single_var_name.upper()}_FLOOR"
            )
            return_str += "//// END APPLICATION OF POSITIVE FLOOR\n////\n\n"


        # now generate the metric update code
        return_str += "//// NOW UPDATING THE METRIC VALUES\n"
        return_str += dendrosym.codegen.generate_update_sym_mat_code(
            metric_name, metric_enums)

        return return_str

    @staticmethod
    def create_grad_var_names(in_vars: list, grad_type="grad", ndim=3):
        """Generate the gradient variable names

        This also takes in the gradient type, which is specifically for
        generating the full list of gradient variable names.

        Notes
        -----
        TODO: this might not be entirely necessary, there is a backup
        method that scans the output code for the gradient types, but
        at least this way we can generate the actual gradient names
        based on the full list.
        """

        grad_vars = []

        # 1d gradient
        if grad_type == "grad":
            for curr_var in in_vars:
                for ii in range(ndim):
                    grad_vars.append(f"grad_{ii}_{curr_var}")
        # 2d gradient
        elif grad_type == "grad2":
            for curr_var in in_vars:
                for ii in range(ndim):
                    for jj in range(ndim):
                        grad_vars.append(f"grad2_{ii}_{jj}_{curr_var}")
        # advective gradient variables
        elif grad_type == "agrad":
            for curr_var in in_vars:
                for ii in range(ndim):
                    grad_vars.append(f"agrad_{ii}_{curr_var}")

        return grad_vars

    @staticmethod
    def generate_rhs_var_names(in_vars: list):

        new_var_list = []
        for the_var in in_vars:
            new_var_list.append(the_var + "_rhs")

        return new_var_list

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
            elif type(the_var) is dendrosym.dtypes.ParameterVariable:
                all_var_names
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

    @staticmethod
    def check_repeat_items_ignore_case(new_item, existing_items):

        for existing in existing_items:
            if new_item.lower() == existing.lower():
                return True

        return False

    def __repr__(self):
        return f"<NRConfigs for '{self.project_name}'>"
