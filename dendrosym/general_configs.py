"""general.py

This file contains the configuration class to generate
and store the necessary pieces of information regarding
the set-up and configuration for a project using the
Dendro framework. The numerical relativity class inherits
many of the classes and methods of those in this file.
"""

import re
import sys
import hashlib
from pathlib import Path
import dill as pickle

import sympy as sym
import numpy as np

import dendrosym

# pickle.detect.trace(True)


class ImproperInitalization(Exception):
    pass


class DendroConfiguration:
    """Store and use configurations for Dendro projects

    This class is used to define and store all of the pieces necessary
    for the Dendro code generation tool. Dendro is an adaptive mesh
    framework written in C/C++ that can divide a problem across
    multiple nodes efficiently. The purpose of this class is to generate
    some of the Dendro-ready C++ code from symbolic Python code to make
    it easier to go from equations to project.

    More information on how to use this class will be forthcoming.

    This class is the base configuration class. For use with general
    relativity projects, use the NRConfigs child class within `nr_configs`.


    """

    def __init__(self, project_name: str, project_description: str = ""):
        self.project_name = project_name
        self.project_upper = project_name.upper()

        self.project_description = project_description

        self.all_vars = {"general": [], "parameter": {"general": []}}
        self.all_var_names = {"general": [], "parameter": {"general": []}}

        self.enum_prefixes = {"general": "G"}

        self.all_rhs_functions = {"general": None}

        self.all_initial_data_functions = {"general": []}

        self.idx_str = "[pp]"

        self.bcs_info = {"general": {}}

        # NOTE: stored_rhs_function is the only thing that is restored by the pickle!
        self.stored_rhs_function = {}

        self.stored_staged_exprs = {}

        # by default we want to replace and or expand the derivatives
        self.replace_and_expand_derivatives = True

        # by default we also want to pull the derivatives from the derivative workspace
        self.use_deriv_workspace = True

        # attempt to restore from pickle
        self.load_all_rhs_from_pickle()

    def create_project_config(self, config_file=""):
        output_dict = {}

        output_dict["project_name"] = self.project_name
        output_dict["project_description"] = self.project_description

        output_dict["all_var_names"] = self.all_var_names

        output_dict["var_counts"] = {}
        for key, item in self.all_var_names.items():
            if key == "parameter":
                continue
            output_dict["var_counts"][key] = len(item)

        output_dict["dendrosym_version"] = dendrosym.__version__

        if config_file != "":
            config_file_hash = dendrosym.utils.calculate_sha256(config_file)
        else:
            config_file_hash = ""

        output_dict["config_file_hash"] = config_file_hash

        output_dict["unique_deriv_counts"] = {}
        largest_num_derivs = 0
        for var_type in self.all_var_names.keys():
            if var_type == "parameter":
                continue
            print(
                f"Now calculating the number of found unique derivatives for var type {var_type}..."
            )
            uq_deriv_count = self.get_num_unique_derivs(var_type)
            print("  ...done")

            if uq_deriv_count > largest_num_derivs:
                largest_num_derivs = uq_deriv_count
            output_dict["unique_deriv_counts"][var_type] = uq_deriv_count

        output_dict["num_derivs"] = largest_num_derivs

        return output_dict

    def create_parameter_variable_config(self):
        # gather all parameters
        pars = self.all_vars["parameter"]

        all_params = {}

        for typ, additional_pars in pars.items():
            # print(typ, additional_pars)
            for par in additional_pars:
                if par.num_params > 1:
                    dtype = par.dtype + "[]"
                else:
                    dtype = par.dtype

                par_name = self.project_upper + "_" + par.var_name.upper()

                all_params[par_name] = {
                    "class": "variant",
                    "dtype": dtype,
                    "desc": par.description,
                    "default": par.default,
                }

                if par.num_params > 1:
                    all_params[par_name]["size"] = par.num_params

        # now we add in the constraint info
        # TODO: might want to update the naming here
        for key, items in self.evolution_constraint_info.items():
            if key == "trace_zero":
                continue
            elif key == "pos_floor":
                # generate params for each of these
                for symb in items:
                    var_name = str(symb)
                    if var_name.endswith(self.idx_str):
                        var_name = var_name.replace(self.idx_str, "")
                    par_name = self.project_upper + "_" + var_name.upper() + "_FLOOR"
                    all_params[par_name] = {
                        "class": "variant",
                        "dtype": "double",
                        "desc": "Positive floor value for variable " + var_name,
                        "default": 1e-5,
                    }

        dsolve_additional_params = {}

        # print(all_params)
        all_params = {self.project_name + "-NMSPC": all_params}

        # then we also need the base dendro additional parameters for

        all_params["requirements"] = []

        for var_type in self.all_var_names.keys():
            if len(self.all_var_names[var_type]) == 0:
                continue
            if var_type == "parameter":
                continue

            num_vars = len(self.all_var_names[var_type])

            dsolve_additional_params["DSOLVER_" + var_type.upper() + "_NUM_VARS"] = {
                "class": "invariant",
                "dtype": "unsigned int",
                "default": num_vars,
                "desc": "Number of " + var_type + " variables.",
            }

            dsolve_additional_params[
                "DSOLVER_VTU_OUTPUT_NUM_" + var_type.upper() + "_VARS"
            ] = {
                "class": "semivariant",
                "dtype": "unsigned int",
                "default": num_vars,
                "desc": "The number of "
                + var_type
                + " variables to be written to the VTU output. Corresponds to the output indices parameter.",
            }

            dsolve_additional_params[
                "DSOLVER_VTU_OUTPUT_" + var_type.upper() + "_INDICES"
            ] = {
                "class": "semivariant",
                "dtype": "unsigned int[]",
                "size": num_vars,
                "default": [ii for ii in range(num_vars)],
                "desc": "The indices of the variables you wish to have as output. The numbers are: "
                + ", ".join(
                    [
                        f"{ii} -> {var}"
                        for ii, var in enumerate(self.all_var_names[var_type])
                    ]
                ),
            }

            all_params["requirements"].append(
                {
                    "cond": "DSOLVER_VTU_OUTPUT_NUM_"
                    + var_type.upper()
                    + "_VARS <= DSOLVER_NUM_"
                    + var_type.upper()
                    + "_NUM_VARS",
                    "message": "Number of VTU "
                    + var_type
                    + " output variables must be less than or equal to the number of variables!",
                }
            )

            if var_type == "evolution" or var_type == "general":
                all_params["requirements"].append(
                    {
                        "cond": "DSOLVER_NUM_REFINE_VARS <= DSOLVER_"
                        + var_type.upper()
                        + "_NUM_VARS",
                        "message": "Number of refinement variables must be less than or equal to the number of variables!",
                    }
                )

                dsolve_additional_params["DSOLVER_NUM_REFINE_VARS"] = {
                    "class": "semivariant",
                    "dtype": "unsigned int",
                    "default": num_vars,
                    "desc": "The number of "
                    + var_type
                    + " variables to used during refinement. Corresponds to the refine variable indices parameter.",
                }

                dsolve_additional_params["DSOLVER_REFINE_VARIABLE_INDICES"] = {
                    "class": "semivariant",
                    "dtype": "unsigned int[]",
                    "size": num_vars,
                    "default": [ii for ii in range(num_vars)],
                    "desc": "The indices of the variables for use in refinement. Be sure to update DSOLVER_NUM_REFINE_VARS. The numbers are: "
                    + ", ".join(
                        [
                            f"{ii} -> {var}"
                            for ii, var in enumerate(self.all_var_names[var_type])
                        ]
                    ),
                }

        # if constraint is in our variables...
        if "constraint" in self.all_var_names.keys():
            num_constraint = len(self.all_var_names["constraint"])
            dsolve_additional_params["DSOLVER_CONSTRAINT_NUM_VARS"] = {
                "class": "invariant",
                "dtype": "unsigned int",
                "default": num_constraint,
                "desc": "Number of constraint variables",
            }
            dsolve_additional_params["DSOLVER_NUM_CONST_VARS_VTU_OUTPUT"] = {
                "class": "semivariant",
                "dtype": "unsigned int",
                "default": num_constraint,
                "desc": "The number of constraint variables to be written to the VTU output. Corresponds to the output indices parameter.",
            }

        all_params["dsolve-NMSPC"] = dsolve_additional_params

        return all_params

    def save_rhs_to_pickle(self, var_type: str):
        """Saves a RHS set of functions to a Pickle

        Parameters
        ----------
        var_type : str
            _description_
        """
        (
            all_exp,
            all_rhs_names,
            staged_exp,
            staged_names,
            orig_n_exp,
        ) = self.get_rhs_functions_all(var_type, include_rhs_in_name=False)

        collection = {
            "all_exp": all_exp,
            "all_rhs_names": all_rhs_names,
            "staged_exp": staged_exp,
            "staged_names": staged_names,
            "orig_n_exp": orig_n_exp,
        }

        # get the dendro folder from working directory
        dendroFolder = Path.cwd() / ".dendro"
        # make the directory, if it doesn't exist
        dendroFolder.mkdir(exist_ok=True)

        rhsPickleFile = dendroFolder / f"{var_type}_rhs.pkl"

        # print(rhsPickleFile)

        with open(rhsPickleFile, "wb") as f:
            pickle.dump(collection, f, recurse=True)

    def load_rhs_from_pickle(self, var_type: str):
        # get the dendro folder from working directory
        dendroFolder = Path.cwd() / ".dendro"

        if not dendroFolder.exists():
            # raise a directory error, then that can be handled in a try-except block
            raise NotADirectoryError(
                "The dendro folder does not exist. Cannot load RHS from pickle"
            )

        rhsPickleFile = dendroFolder / f"{var_type}_rhs.pkl"

        print("loading", rhsPickleFile)

        if not rhsPickleFile.exists():
            print(f"WARNING: Could not find RHS pickle file for '{var_type}'")
            return

        # then load it if it works
        with open(rhsPickleFile, "rb") as f:
            collection = pickle.load(f)

        print("loaded", rhsPickleFile)

        temp_dict = self.stored_rhs_function.get(var_type, {})

        # load it in then
        temp_dict["exprs"] = collection["all_exp"]
        temp_dict["all_rhs_names"] = collection["all_rhs_names"]
        temp_dict["staged_exprs"] = collection["staged_exp"]
        temp_dict["staged_exprs_names"] = collection["staged_names"]
        temp_dict["orig_n_exp"] = collection["orig_n_exp"]

        self.stored_rhs_function[var_type] = temp_dict
        # and done!

    def load_all_rhs_from_pickle(self):
        # get the dendro folder from working directory
        dendroFolder = Path.cwd() / ".dendro"

        if not dendroFolder.exists():
            # if the folder doesn't exist, there's nothing to restore :(
            return

        print("Found the Dendro folder, attempting to load saved pickles")
        # find all of the pickle files that end in _rhs
        for file in dendroFolder.glob("*_rhs.pkl"):
            print("   attempting ", file)
            var_name = file.name[:-8]

            self.load_rhs_from_pickle(var_name)

            print("    Loaded - ", var_name)

        # and done!

    def save_all_rhs_to_pickle(self):
        for var_name in self.all_rhs_functions.keys():
            if var_name == "parameters" or self.all_rhs_functions[var_name] is None:
                continue

            print("Saving", var_name, "to a pickle file to save work...")

            self.save_rhs_to_pickle(var_name)

        # and done!

    def set_idx_str(self, idx_str):
        """Store the string used for indexing into the variables

        For the Dendro projects, this should be "[pp]", but the flexibility
        is provided.

        TODO: it might be convenient to make this "update idx str"
        for if someone wants to use something *other* than "[pp]".
        """
        self.idx_str = idx_str

        dendrosym.derivs.idx_str = idx_str

    def override_derivative_expansion(self):
        """This function enables the override for replacing and expanding derivatives

        This should **not** be called unless you know for a fact that you have no nested
        derivatives that should be expanded or computed before hand. In the BSSN case,
        we know this to be true so we can skip these calculations. But other theories
        or even ways to express the theory can be different.
        """
        self.replace_and_expand_derivatives = False

    def add_parameter_variables(self, in_vars: list, eqn_type: str = "general"):
        """Add a parameter variable to the list

        Use this function when there is a constant parameter variable
        that you have in your equations. Be sure to input
        the type of equation that it belongs to. Currently
        accepts "evolution" or "constraint"
        """

        if eqn_type not in self.all_vars.keys():
            raise ImproperInitalization(
                f"Equation type {eqn_type} has not been declared"
            )

        if type(in_vars) is not list:
            in_vars = [in_vars]

        for one_var in in_vars:
            # get the variable name from the object
            var_names_clean = one_var.var_name

            if not self.check_repeat_items_ignore_case(
                var_names_clean, self.all_var_names["parameter"][eqn_type]
            ):
                self.all_var_names["parameter"][eqn_type].append(var_names_clean)
            else:
                raise Exception(f"Parameter '{one_var}' has already been added")
            # then append the paramter to the list
            self.all_vars["parameter"][eqn_type].append(one_var)

    def add_variable(self, in_vars: list, eqn_type: str = "general"):
        """Add a symbolic variable to a particular equation type

        This is necessary so the code generation can allocate memory,
        keep track of, and store variables throughout the problem's evolution.

        Parameters
        ----------
        in_vars:
        """

        if eqn_type not in self.all_vars.keys():
            raise ImproperInitalization(
                f"Equation type {eqn_type} has not been declared"
            )

        if type(in_vars) is not list:
            in_vars = [in_vars]

        for one_var in in_vars:
            # get the clean var name
            var_names_clean = self.clean_var_names([one_var])

            for vname_clean in var_names_clean:
                # if not vname_clean.endswith(self.idx_str):
                #     vname_clean += self.idx_str
                if not self.check_repeat_items_ignore_case(
                    vname_clean, self.all_var_names[eqn_type]
                ):
                    self.all_var_names[eqn_type].append(vname_clean)
                else:
                    raise Exception(
                        f"'{vname_clean}' has already been added or is"
                        + " too similar to another existing one"
                    )
            # then append the paramter to the list
            self.all_vars[eqn_type].append(one_var)

        # then update the global variables names in derivatives for good measure
        flattened_list = []
        for tt, vars in self.all_var_names.items():
            flattened_list += vars

        dendrosym.derivs.variable_strs = flattened_list

    def get_rhs_var_names(self, var_type: str = "general"):
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

    def generate_rhs_code(
        self, var_type: str, arc_type="cpu", include_rhs_in_name=True
    ):
        """Generate the code that calculates the 'RHS' variables

        This method takes the SymPy equations that were stored during configuration
        and then generates the C++ code that calculates them. It is referred to
        RHS even if every equation type is the RHS of a partial differential
        system of equations.

        For example, if "evolution" is a variable type, that can be sent in
        and the class will then generate the code.

        As an additional bonus, if used in conjunction with other methods,
        it will take into account other "optimizations" used to either simplify,
        precalculate, or likewise the equations automatically.

        Please note that while `arc_type` is a possible parameter, it currently
        only supports CPU code generation. GPU is planned.
        """

        if var_type not in self.all_var_names.keys():
            raise ValueError(f"Unfortunately {var_type} doesn't work yet")

        if self.stored_rhs_function.get(var_type, None) is not None:
            temp_funcs = self.stored_rhs_function[var_type]
            all_exp = temp_funcs["exprs"]
            all_rhs_names = temp_funcs["all_rhs_names"]
            if not include_rhs_in_name:
                all_rhs_names_tmp = []
                for rhs_name in all_rhs_names:
                    # replace/remove the _rhs side
                    all_rhs_names_tmp.append(rhs_name.replace("_rhs", ""))
                all_rhs_names = all_rhs_names_tmp
            orig_n_exp = temp_funcs["orig_n_exp"]

            staged_exp = temp_funcs["staged_exprs"]
            staged_exprs_names = temp_funcs["staged_exprs_names"]
            print("Found stored RHS info", file=sys.stderr)
        else:
            # so now we can get started by getting the rhs information
            (
                all_exp,
                all_rhs_names,
                staged_exp,
                staged_names,
                orig_n_exp,
            ) = self._extract_rhs_expressions(
                var_type, append_rhs_to_var=include_rhs_in_name
            )

        # now to restore the 'symbols' inside the derivitaves to their original functionality
        if (
            type(dendrosym.nr.d) is not sym.core.function.UndefinedFunction
            or type(dendrosym.nr.d2s) is not sym.core.function.UndefinedFunction
        ):
            print(".... Now restoring some of the symbols from derivatives...")
            temp_exprs = []
            for ii, expr in enumerate(staged_exp):
                expr_temp = dendrosym.derivs.restore_only_symbols(expr)
                temp_exprs.append(expr_temp)

            staged_exp = temp_exprs

            temp_exprs = []

            for ii, expr in enumerate(all_exp):
                print(f"     on expr... {ii + 1}/{len(all_exp)}")
                expr_temp = dendrosym.derivs.restore_only_symbols(expr)
                temp_exprs.append(expr_temp)

            all_exp = temp_exprs
            temp_exprs = None

            print("... finished replacing the derivatives")

        # start with the staged expressions
        if len(staged_exp) > 0:
            orig_n_ops = sym.count_ops(staged_exp)
            cse_exp = dendrosym.codegen.construct_cse_from_list(
                staged_exp, temp_var_prefix="DENDRO_STAGED_VAR_"
            )

            # NO INDEX STRING HERE, WE AREN'T INDEXING INTO THEM
            output_str = dendrosym.codegen.generate_cpu_preextracted(
                cse_exp, staged_exprs_names, "", orig_n_ops
            )
        else:
            orig_n_ops = 0
            output_str = ""

        # count the number of original operations on the expressions
        orig_n_ops += sym.count_ops(all_exp)

        # construct cse from the list
        cse_exp = dendrosym.codegen.construct_cse_from_list(all_exp)

        if arc_type == "cpu":
            # then we're to send this information to the cpu gen function

            main_output_str = dendrosym.codegen.generate_cpu_preextracted(
                cse_exp, all_rhs_names, self.idx_str, orig_n_ops
            )

            output_str += main_output_str
            return output_str

        else:
            raise NotImplementedError(
                "That archetecture generation code isn't ready yet"
            )

    def generate_rhs_separate(self, var_type: str, arc_type="cpu"):
        if var_type not in self.all_var_names.keys():
            raise ValueError(f"'{var_type}' is not a valid RHS type")

        # get the expressions from the rhs info
        (
            all_exp,
            all_rhs_names,
            staged_exprs,
            staged_vars,
            original_number_expressions,
        ) = self._extract_rhs_expressions(var_type)

        # now we iterate through all of them individually to get their new code

        if arc_type == "cpu":
            return dendrosym.codegen.generate_separate_cpu(
                all_exp, all_rhs_names, self.idx_str, orig_n_exp
            )

    def set_bhs_falloff_and_asymptotic(self, var_type, var_list, var_info):
        """Sets the falloff and asymptotic values to the individual variables

        To use this function, say what type of variables you're providing
        information for, and then pass through a dictionary of values for
        each of the varibles you are assigning the data to. That dictionary
        should be of structure:

        var_information = {beta: {"falloff": 1.0, "asymptotic": 2.0}}

        and include all variables through the dictionary construction. This
        function will throw an error if not all variables are accounted for.
        """

        if var_type == "parameter":
            raise ValueError("Cannot set BHS data to parameters")

        temp_vars = self.all_vars.get(var_type, []).copy()

        if len(temp_vars) == 0:
            raise ImproperInitalization(
                f"{var_type} does not have variables assigned to it."
            )

        # let's make sure that they're all there
        for ii, var_ in enumerate(var_list):
            found_match = False

            # iterate through our temp_vars list and find a match
            for temp_var in temp_vars:
                if temp_var == var_:
                    # say we found a match, and then remove it
                    found_match = True
                    temp_vars.remove(temp_var)

                    break

            if not found_match:
                raise ImproperInitalization(
                    f"Couldn't find identical variable {var_} in stored list"
                )

        # now make sure that we've cleared both lists
        if len(temp_vars) != 0:
            raise ImproperInitalization(
                "Not all assigned variables to '"
                + var_type
                + "' were included in RHS function. Remaining variables were:"
                + repr(temp_vars)
            )

        self.bcs_info.update({var_type: {"vars": var_list, "info": var_info}})

    def generate_bcs_calculations(
        self,
        var_type="general",
        pmin="pmin",
        pmax="pmax",
        sz="sz",
        bflag="bflag",
    ):
        """Generate the BCS calculation functions

        This one is tricky as different variables require using
        """

        # just iterate through the data we have stored
        all_var_info = self.bcs_info.get(var_type, None)

        if all_var_info is None:
            raise ImproperInitalization(
                f"BCS information was not initialized for {var_type}."
            )

        if not all_var_info:
            raise ImproperInitalization(
                f"BCS information was not initialized for {var_type}."
            )

        return_str = ""

        for ii, the_var in enumerate(all_var_info["vars"]):
            var_info = all_var_info["info"][ii]
            cleaned_var_name = self.clean_var_names([the_var])

            if len(cleaned_var_name) > 1:
                for clean_var in cleaned_var_name:
                    rhs_var = self.generate_rhs_var_names([clean_var])

                    grad_vars = self.create_grad_var_names([clean_var], "grad", 3)

                    if len(var_info) == 2 and (
                        type(var_info[0]) is int or type(var_info[0]) is float
                    ):
                        falloff = var_info[0]
                        asymptotic = var_info[1]
                    else:
                        idxs = self.get_indices_from_var_name(clean_var)

                        if len(idxs) == 1:
                            falloff = var_info[idxs[0]][0]
                            asymptotic = var_info[idxs[0]][1]
                        else:
                            # print(clean_var, file=sys.stderr)
                            try:
                                falloff = var_info[idxs[0]][idxs[1]][0]
                                asymptotic = var_info[idxs[0]][idxs[1]][1]
                            except:
                                print(f"FAILURE {var_info}")
                                raise Exception

                    return_str += dendrosym.codegen.generate_bcs_function_call(
                        rhs_var[0],
                        clean_var,
                        falloff,
                        asymptotic,
                        grad_vars,
                        self.project_name,
                        pmin,
                        pmax,
                        sz,
                        bflag,
                    )

                pass

            else:
                rhs_var = self.generate_rhs_var_names(cleaned_var_name)

                grad_vars = self.create_grad_var_names(cleaned_var_name, "grad", 3)
                return_str += dendrosym.codegen.generate_bcs_function_call(
                    rhs_var[0],
                    cleaned_var_name[0],
                    var_info[0],
                    var_info[1],
                    grad_vars,
                    self.project_name,
                    pmin,
                    pmax,
                    sz,
                    bflag,
                )

        return return_str

    def generate_variable_extraction(
        self,
        var_type="general",
        use_const=True,
        zip_var_name="uZipVars",
        dtype="double",
    ):
        """Generate the C++ Code for Variable Extraction

        The variable extraction refers to the zipping and
        unzipping that happens with a variable like 'uZipVars'.
        This function can generate all of the code necessary
        to create pointers to where that data is stored.
        """

        if var_type == "parameter":
            raise ValueError("Parameters are generated with another function")

        named_enums = self.get_enum_var_names(var_type)

        if var_type in ["evolution", "general"]:
            enum_name = "VAR"
        elif var_type == "constraint":
            enum_name = "VAR_CONSTRAINT"
        else:
            raise NotImplementedError("Not Yet Implemented")

        return_str = dendrosym.codegen.gen_var_info(
            self.all_var_names.get(var_type, []),
            zip_var_name=zip_var_name,
            use_const=use_const,
            enum_name=enum_name,
            enum_var_names=named_enums,
            dtype=dtype,
        )

        return return_str

    def generate_rhs_var_extraction(
        self,
        var_type="general",
        use_const=False,
        zip_var_name="unZipVarsRHS",
        dtype="double",
    ):
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

        return_str = dendrosym.codegen.gen_var_info(
            vars_use,
            zip_var_name=zip_var_name,
            use_const=use_const,
            dtype=dtype,
            enum_var_names=named_enums,
        )

        return return_str

    def gen_enum_code(
        self, var_type: str, enum_start_idx: int = 0, enum_name: str = "VAR"
    ):
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
            enum_names, self.project_upper, enum_name
        )

        return return_str

    def gen_enum_iterable_list(self, var_type: str, enum_name: str = "VAR"):
        if var_type == "parameter":
            raise ValueError("This function can't handle parameter enum names")

        enum_names = self.get_enum_var_names(var_type)

        return_str = dendrosym.codegen.gen_var_iterable_list(
            enum_names, self.project_upper, enum_name
        )

        return return_str

    def gen_parameter_code(self, var_type="evolution"):
        return_str = ""

        for param in self.all_vars["parameter"].get(var_type, []):
            return_str += param.generate_cpp_line(
                global_param_prefix=self.project_upper
            )
            return_str += "\n"

        return return_str

    def get_enum_var_names(self, var_type: str):
        if var_type == "parameter":
            raise NotImplementedError("Not yet ready")

        orig_vars = self.all_var_names.get(var_type, [])

        enum_prefix = self.enum_prefixes.get(var_type, "")

        return self._generate_enum_var_names(orig_vars, enum_prefix)

    def _generate_enum_var_names(self, list_vars, enum_prefix):
        out_enum_names = []
        for the_var in list_vars:
            out_enum_names.append(f"{enum_prefix}_{the_var.upper()}")

        return out_enum_names

    def _gen_enum_var_names_from_sympy_vars(self, list_vars, enum_prefix):
        # first need to conver them
        list_vars_clean = []
        for one_var in list_vars:
            tmp = self.clean_var_names([one_var])

            list_vars_clean += tmp

        return self._generate_enum_var_names(list_vars_clean, enum_prefix)

    def set_rhs_equation_function(
        self, var_type: str, rhs_func, override_checks: bool = False
    ):
        if var_type == "parameter":
            raise ValueError("Cannot set RHS function to parameters")

        if not override_checks:
            # first we need to check that the function gives us what we want
            # so we need to evaluate it first
            rhs_list, var_list = rhs_func()

            if len(rhs_list) != len(var_list):
                raise ImproperInitalization(
                    "The RHS function does not have the same number"
                    + " of expressions as variables"
                )

            # so, now we make a copy of our variables for this variable type
            temp_vars = self.all_vars.get(var_type, []).copy()

            if len(temp_vars) == 0:
                raise ImproperInitalization(
                    f"{var_type} does not have variables assigned to it."
                )

            # we can't change the order, but we need to make sure they're all there
            # so we traverse the variable list
            for ii, var_ in enumerate(var_list):
                found_match = False

                # iterate through our temp_vars list and find a match
                for temp_var in temp_vars:
                    if temp_var == var_:
                        # say we found a match, and then remove it
                        found_match = True
                        temp_vars.remove(temp_var)
                        break

                if not found_match:
                    raise ImproperInitalization(
                        f"Couldn't find identical variable {var_} in stored list"
                    )

            # now make sure that we've cleared both lists
            if len(temp_vars) != 0:
                raise ImproperInitalization(
                    "Not all assigned variables to '"
                    + var_type
                    + "' were included in RHS function. Remaining variables were:"
                    + repr(temp_vars)
                )

        # if that's all good, then we're good to add the function to our list
        self.all_rhs_functions.update({var_type: rhs_func})

    def get_rhs_eqns(self, var_type: str):
        """Get the plain output from the RHS equation"""
        rhs_func = self.all_rhs_functions.get(var_type, None)

        rhs_list, var_list = rhs_func()

        return rhs_list, var_list

    def get_rhs_eqns_flat(self, var_type: str):
        """This function gets the RHS equations and makes them "flat"

        Useful for evaluating each possible variable.
        """

        rhs_func = self.all_rhs_functions.get(var_type, None)

        rhs_list, var_list = rhs_func()

        mi = [0, 1, 2, 4, 5, 8]

        lexp = []
        lname = []

        num_e = 0

        # flatten the expressions to lists with corresponding names
        for i, e in enumerate(rhs_list):
            if type(e) == list:
                num_e = num_e + len(e)
                for j, ev in enumerate(e):
                    lexp.append(ev)
                    lname.append(str(var_list[i][j]))
            elif type(e) == sym.Matrix:
                if len(e.shape) == 2:
                    first_dim, second_dim = e.shape
                    print(var_list[i])
                    for j in range(first_dim):
                        for k in range(second_dim):
                            if var_list[i][j, k] not in lname:
                                lexp.append(e[j, k])
                                lname.append(var_list[i][j, k])
                                num_e += 1
                else:
                    num_e = num_e + np.prod(e.shape)
                    for j, k in enumerate(mi):
                        lexp.append(e[k])
                        lname.append(str(var_list[i][j]))
            else:
                num_e = num_e + 1
                lexp.append(e)
                lname.append(str(var_list[i]))

        return lexp, lname

    def get_rhs_functions_all(self, var_type: str, include_rhs_in_name: bool = True):
        """A more advanced version of getting RHS functions, returns 5 things"""

        if self.stored_rhs_function.get(var_type, None) is not None:
            temp_funcs = self.stored_rhs_function[var_type]
            all_exp = temp_funcs["exprs"]
            all_rhs_names = temp_funcs["all_rhs_names"]
            if not include_rhs_in_name:
                all_rhs_names_tmp = []
                for rhs_name in all_rhs_names:
                    # replace/remove the _rhs side
                    all_rhs_names_tmp.append(rhs_name.replace("_rhs", ""))
                all_rhs_names = all_rhs_names_tmp
            orig_n_exp = temp_funcs["orig_n_exp"]

            staged_exp = temp_funcs["staged_exprs"]
            staged_names = temp_funcs["staged_exprs_names"]
            print("Found stored RHS info", file=sys.stderr)
        else:
            # so now we can get started by getting the rhs information
            (
                all_exp,
                all_rhs_names,
                staged_exp,
                staged_names,
                orig_n_exp,
            ) = self._extract_rhs_expressions(
                var_type, append_rhs_to_var=include_rhs_in_name
            )

        return all_exp, all_rhs_names, staged_exp, staged_names, orig_n_exp

    def _extract_rhs_expressions(self, var_type: str, append_rhs_to_var=True):
        """An internal function that extracts the expressions for a variable type

        It does this to match up expressions with their corresponding RHS
        variables and put the list in order. This is because the class
        manages the creation of RHS and grad variables.
        """

        if var_type == "parameter":
            raise ValueError("Cannot extract expressions from parameters")

        rhs_func = self.all_rhs_functions.get(var_type, None)

        if rhs_func is None:
            raise ImproperInitalization(
                f"A RHS function was not assigned to {var_type}"
            )

        # now we need the variables, and the good news is that we know that the
        # function we extracted can't exist without having gone through
        # our filtering function to make sure all variables are present
        # and accounted for.
        # TODO: one problem would be adding more variables after setting the
        # function... need to fix that

        # evaluate the function so we can get the pieces
        rhs_list, var_list = rhs_func()

        # iterate through it and put together all expressions and var names
        # in order
        all_expressions = []
        all_rhs_var_names = []
        original_number_expressions = 0

        for ii in range(len(rhs_list)):
            expression = rhs_list[ii]
            the_var = var_list[ii]

            # check if var is symmetric
            is_sym = False
            if type(the_var) == sym.Matrix:
                if the_var.shape == (3, 3):
                    # check if sym 3x3
                    if (
                        the_var[0, 1] == the_var[1, 0]
                        and the_var[0, 2] == the_var[2, 0]
                        and the_var[1, 2] == the_var[2, 1]
                    ):
                        is_sym = True

            list_expressions, num_e = dendrosym.codegen.extract_expression(
                expression, is_symmetric_matrix=is_sym
            )

            # note: if we have a sym.Matrix as our variable, we need to then
            # keep in mind the indexing so we can build the RHS variables
            if append_rhs_to_var:
                rhs_vars = self.generate_rhs_var_names(self.clean_var_names([the_var]))
            else:
                rhs_vars = self.clean_var_names([the_var])

            all_expressions += list_expressions
            all_rhs_var_names += rhs_vars
            original_number_expressions += num_e

        # for ii in range(len(all_expressions)):
        #     print(all_rhs_var_names[ii], "=", all_expressions[ii])

        # now add the staged variables, which are used by the user to simplify things by hand
        if self.stored_staged_exprs.get(var_type, None) is None:
            staged_vars = []
            staged_exprs = []
        else:
            # run the function, and assume it's the only one
            staged_exprs, staged_vars = self.stored_staged_exprs.get(var_type)()

            # now we can make sure they have the same number
            if len(staged_exprs) != len(staged_vars):
                raise ImproperInitalization(
                    f"The staged expression function for {var_type} has an invalid number of outputs"
                )

            # quick check on the staged expressions, make sure they're not in a sympy matrix
            for expr in staged_exprs:
                if isinstance(expr, sym.Matrix):
                    raise ImproperInitalization(
                        "Staged expressions must not be vectors or matrices, ensure that each staged variable is its own single expression"
                    )

            # otherwise we're good to continue

        # okay, now that we have them we can return them
        return (
            all_expressions,
            all_rhs_var_names,
            staged_exprs,
            staged_vars,
            original_number_expressions,
        )

    def gen_grad_memory_alloc(
        self,
        var_type: str,
        grad_type: str = "grad",
        include_byte_declaration=False,
    ):
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

        return_text, current_index = dendrosym.codegen.generate_memory_alloc(
            grad_vars,
            "double",
            not self.use_deriv_workspace,
            include_byte_declaration,
            start_id=self.current_index,
        )

        self.current_index = current_index

        return return_text

    @staticmethod
    def gen_grad_memory_alloc_from_code(
        filename: str,
        grad_type: str = "grad",
        include_byte_declaration: bool = False,
        dtype: str = "double",
    ):
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
        grad_list = dendrosym.utils.find_derivtype_in_text(in_text, grad_type + "_")

        # sort it for clarity
        grad_list.sort()

        if len(grad_list) == 0:
            # TODO: make sure to actually call some kind of logger
            # or warning system
            print(
                "WARNING: no variables found with "
                + f"'{grad_type}' string in file {filename}"
            )

        # then we can create the grad variable names directly
        return_text = dendrosym.codegen.generate_memory_alloc(
            grad_list,
            var_type=dtype,
            include_byte_declaration=include_byte_declaration,
        )

        return return_text

    def gen_grad_memory_dealloc(self, var_type: str, grad_type: str = "grad"):
        if var_type == "parameter":
            raise ValueError("Parameters cannot use gradients")

        orig_vars = self.all_var_names.get(var_type, [])

        grad_vars = self.create_grad_var_names(orig_vars, grad_type, 3)

        return_text = dendrosym.codegen.generate_memory_dealloc(grad_vars)

        return return_text

    @staticmethod
    def gen_grad_memory_dealloc_from_code(filename: str, grad_type: str = "grad"):
        with open(filename, "r") as fh:
            in_text = fh.read()

        # first get the grad list from the "find deriv type in text" function
        grad_list = dendrosym.utils.find_derivtype_in_text(in_text, grad_type + "_")

        # sort it for clarity
        grad_list.sort()

        if len(grad_list) == 0:
            # TODO: make sure to actually call some kind of logger
            # or warning system
            print(
                "WARNING: no variables found with "
                + f"'{grad_type}' string in file {filename}"
            )

        # then we can create the grad variable names directly
        return_text = dendrosym.codegen.generate_memory_dealloc(grad_list)

        return return_text

    def gen_grad_calculations(
        self, var_type: str, grad_type: str = "grad", use_eqns=False
    ):
        """"""
        if use_eqns:
            pass

        else:
            orig_vars = self.all_var_names.get(var_type, [])

            grad_vars = self.create_grad_var_names(orig_vars, grad_type, 3)

            grad_vars.sort()

            return_text = dendrosym.codegen.generate_deriv_comp(grad_vars)

        return return_text

    def generate_ko_derivs(
        self, var_type: str, sz="sz", bflag="bflag", ko_func_name="ko_deriv"
    ):
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
            temp_str += " + ".join(
                f"grad_{ii}_{var_name}{self.idx_str}" for ii in range(3)
            )
            temp_str += ");\n"

            return_str += temp_str

        return return_str

    @staticmethod
    def get_indices_from_var_name(var_name):
        var_digits = re.findall("[0-9]+", var_name)
        var_end = var_digits[-1]

        if len(var_end) == 1:
            return (int(var_end),)

        # take the length of this string and divide it by two for how many
        # characters each index is
        assert len(var_end) % 2 == 0, "incoming string not divisible by two"
        len_idx = int(len(var_end) / 2)

        # we're only dealing with 2D for now
        return (int(var_end[:len_idx]), int(var_end[len_idx:]))

    @staticmethod
    def create_grad_var_names(
        in_vars: list, grad_type="grad", ndim=3, assume_symmetry=True
    ):
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
                collected_dirs = []
                for ii in range(ndim):
                    for jj in range(ndim):
                        if assume_symmetry:
                            if (jj, ii) in collected_dirs:
                                continue

                        grad_vars.append(f"grad2_{ii}_{jj}_{curr_var}")

                        collected_dirs.append((ii, jj))

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
        """Create list of strings for the input variables"""
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
                all_var_names.append(the_var.var_name)
            else:
                raise NotImplementedError("That variable type isn't implemented yet")

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
        return f"<DendroConfigs for '{self.project_name}'>"

    def add_initial_data(
        self,
        in_func,
        func_name,
        var_type="general",
        initial_data_id=None,
        has_explicit_solution=False,
    ):
        if initial_data_id is None:
            initial_data_id = 0
            initial_data_funcs = self.all_initial_data_functions.get(var_type, [])

            for _ in initial_data_funcs:
                initial_data_id += 1

        # then we need to assign it and make sure the sizes align like the
        # regular functions
        exprs, vars_list = in_func()

        # now get a copy of our variables for this type
        temp_vars = self.all_vars.get(var_type, []).copy()

        if len(temp_vars) == 0:
            raise ImproperInitalization(
                f"{var_type} does not have variables assigned to it."
            )

        for ii, var_ in enumerate(vars_list):
            found_match = False

            # then we can iterate through our temp vars list to find a match
            for temp_var in temp_vars:
                if temp_var == var_:
                    # say we found the match, and then remove it from the temp
                    found_match = True
                    temp_vars.remove(temp_var)
                    break

            if not found_match:
                raise ImproperInitalization(
                    f"Couldn't find identical variable {var_} in stored list"
                )

        if len(temp_vars) != 0:
            raise ImproperInitalization(
                "Not all assigned variables to '" + var_type + "' were included"
                " in Initial Data function. Remaining variables were:" + repr(temp_vars)
            )

        self.all_initial_data_functions[var_type].append(
            (initial_data_id, func_name, in_func, has_explicit_solution)
        )

    def generate_initial_data_declaration(self, var_type="general", dtype="double"):
        # first we need to get the functions
        all_init_funcs = self.all_initial_data_functions.get(var_type, [])

        if len(all_init_funcs) == 0:
            # if there's nothing we need to just return an empty string, nothing to generate
            return ""

        outstr = ""

        for init_data_id, func_name, init_func in all_init_funcs:
            outstr += f"/**\n * @brief {func_name} initialization function\n"
            outstr += " *\n * @param x : x coord (octree coord)\n"
            outstr += " * @param y : y coord (octree coord)\n"
            outstr += " * @param z : z coord (octree coord)\n"
            outstr += (
                " * @param var : initialized dsolve variable pointer for grid points\n"
            )
            outstr += " */\n"
            outstr += f"void {func_name}(const {dtype} x, const {dtype} y, const {dtype} z, {dtype} *var);\n\n"

        return outstr

    def generate_initial_data_code(self, var_type="general", dtype="double"):
        # first we need to get the functions
        all_init_funcs = self.all_initial_data_functions.get(var_type, [])

        if len(all_init_funcs) == 0:
            # if there's nothing we need to just return an empty string, nothing to generate
            return ""

        outstr = ""

        named_enums = self.get_enum_var_names(var_type)

        if var_type in ["evolution", "general"]:
            enum_name = "VAR"
        elif var_type == "constraint":
            enum_name = "VAR_CONSTRAINT"
        else:
            raise NotImplementedError("Not yet implemented")

        enum_extraction_str = dendrosym.codegen.gen_var_info(
            self.all_var_names.get(var_type, []),
            zip_var_name="var",
            use_const=False,
            enum_name=enum_name,
            enum_var_names=named_enums,
            dtype=dtype,
            offset_name="",
            num_spaces=4,
        )

        # now we iterate through all functions
        for init_data_id, func_name, init_func, _ in all_init_funcs:
            outstr += f"void {func_name}(const {dtype} x, const {dtype} y, const {dtype} z, {dtype} *var){{\n"
            outstr += f"    const unsigned int {self.idx_str[1:-1]} = 0;\n"
            outstr += enum_extraction_str

            par_code = self.gen_parameter_code(var_type)
            outstr += par_code

            # then we have to generate the code
            expressions, variables = init_func()
            # # need to sympify the expressions in the list first:
            # # just in case they are constants!
            # expressions = [sym.sympify(expr) for expr in expressions]
            if (
                type(dendrosym.nr.d) is not sym.core.function.UndefinedFunction
                or type(dendrosym.nr.d2s) is not sym.core.function.UndefinedFunction
            ):
                print(".... Now restoring some of the symbols from derivatives...")
                temp_exprs = []
                for ii, expr in enumerate(expressions):
                    expr_temp = dendrosym.derivs.restore_only_symbols(expr)
                    temp_exprs.append(expr_temp)

                expressions = temp_exprs

                temp_exprs = None

                print("... finished replacing the derivatives")

            # convert the variables to clean variables
            tmp_vars = []
            for the_var in variables:
                tmp_vars += self.clean_var_names([the_var])
            variables = tmp_vars

            # then convert the expressions to be extracted
            tmp_expressions = []
            original_number_expressions = 0
            for expr in expressions:
                list_expressions, num_e = dendrosym.codegen.extract_expression(expr)
                tmp_expressions += list_expressions
                original_number_expressions += num_e

            expressions = tmp_expressions

            # then construct the cse from the list
            cse_exp = dendrosym.codegen.construct_cse_from_list(expressions)

            # TODO: modify this call (or use a separate function?)
            # NOTE: this needs to be redone because we're generating full extraction vars but var is already at the coords we want
            outstr += "    "
            outstr += "    ".join(
                dendrosym.codegen.generate_cpu_preextracted(
                    cse_exp, variables, "[0]", original_number_expressions
                ).splitlines(True)
            )

            outstr += "}\n\n"

        return outstr

    def generate_explicit_solution_code(self, var_type="general", dtype="double"):
        # first we need to get the functions
        all_init_funcs = self.all_initial_data_functions.get(var_type, [])

        if len(all_init_funcs) == 0:
            # if there's nothing we need to just return an empty string, nothing to generate
            return ""

        outstr = ""

        named_enums = self.get_enum_var_names(var_type)

        if var_type in ["evolution", "general"]:
            enum_name = "VAR"
        elif var_type == "constraint":
            enum_name = "VAR_CONSTRAINT"
        else:
            raise NotImplementedError("Not yet implemented")

        enum_extraction_str = dendrosym.codegen.gen_var_info(
            self.all_var_names.get(var_type, []),
            zip_var_name="var",
            use_const=False,
            enum_name=enum_name,
            enum_var_names=named_enums,
            dtype=dtype,
            offset_name="",
            num_spaces=4,
        )

        # now we iterate through all functions
        for init_data_id, func_name, init_func, has_explicit_solution in all_init_funcs:
            if not has_explicit_solution:
                continue

            outstr += f"void {func_name}Analytical(const {dtype} xx, const {dtype} yy, const {dtype} zz, {dtype} t, {dtype} *var, bool varsAreGrid = true){{\n"
            outstr += f"    const unsigned int {self.idx_str[1:-1]} = 0;\n"
            outstr += f"    {dtype} x, y, z;\n"
            outstr += """    if (varsAreGrid) {
        x = GRIDX_TO_X(xx);
        y = GRIDY_TO_Y(yy);
        z = GRIDZ_TO_Z(zz);
    } else {
        x = xx;
        y = yy;
        z = zz;
    }
\n\n"""
            outstr += enum_extraction_str

            par_code = self.gen_parameter_code(var_type)
            outstr += par_code

            # then we have to generate the code
            expressions, variables = init_func(t=sym.Symbol("t"))
            # # need to sympify the expressions in the list first:
            # # just in case they are constants!
            # expressions = [sym.sympify(expr) for expr in expressions]
            if (
                type(dendrosym.nr.d) is not sym.core.function.UndefinedFunction
                or type(dendrosym.nr.d2s) is not sym.core.function.UndefinedFunction
            ):
                print(".... Now restoring some of the symbols from derivatives...")
                temp_exprs = []
                for ii, expr in enumerate(expressions):
                    expr_temp = dendrosym.derivs.restore_only_symbols(expr)
                    temp_exprs.append(expr_temp)

                expressions = temp_exprs

                temp_exprs = None

                print("... finished replacing the derivatives")

            # convert the variables to clean variables
            tmp_vars = []
            for the_var in variables:
                tmp_vars += self.clean_var_names([the_var])
            variables = tmp_vars

            # then convert the expressions to be extracted
            tmp_expressions = []
            original_number_expressions = 0
            for expr in expressions:
                list_expressions, num_e = dendrosym.codegen.extract_expression(expr)
                tmp_expressions += list_expressions
                original_number_expressions += num_e

            expressions = tmp_expressions

            # then construct the cse from the list
            cse_exp = dendrosym.codegen.construct_cse_from_list(expressions)

            # TODO: modify this call (or use a separate function?)
            # NOTE: this needs to be redone because we're generating full extraction vars but var is already at the coords we want
            outstr += "    "
            outstr += "    ".join(
                dendrosym.codegen.generate_cpu_preextracted(
                    cse_exp, variables, "[0]", original_number_expressions
                ).splitlines(True)
            )

            outstr += "}\n\n"

        return outstr

    def find_derivatives(self, var_type, do_extra_check=False):
        """This method finds other derivatives

        This is particularly useful for when a user uses an
        equation that takes the derivative of something beyond
        just the variables used in the problem.

        For a simple example, d(x * y) is something that may happen,
        but would need to be calculated beforehand.
        """

        # start by checking if we have them stored to avoid recalculating
        # things for the pre-extraction
        if self.stored_rhs_function.get(var_type, None) is not None:
            temp_funcs = self.stored_rhs_function[var_type]
            return (
                temp_funcs["exprs"],
                temp_funcs["all_rhs_names"],
                temp_funcs["found_derivatives"],
                temp_funcs["orig_n_exp"],
                temp_funcs["staged_exprs"],
                temp_funcs["staged_exprs_names"],
            )

        # then check the stored data for the functions
        if var_type not in self.all_var_names.keys():
            raise ValueError(f"Unfortunately {var_type} doesn't work yet")

        # so now we can get started by getting the rhs information
        (
            all_exp,
            all_rhs_names,
            staged_exp,
            staged_names,
            orig_n_exp,
        ) = self._extract_rhs_expressions(var_type)

        # collection of "new" (modified) expressions for each of the variables
        new_exprs = []

        if self.replace_and_expand_derivatives:
            print(
                "... Now finding derivatives to attempt chain rule expansion or staging..."
            )

            for ii, expr in enumerate(all_exp):
                print(
                    "    "
                    + str(all_rhs_names[ii])
                    + f"  - prog. {ii + 1}/{len(all_exp)} : {(ii + 1) / len(all_exp):.2%}",
                    file=sys.stderr,
                )

                # call the find and replace complicated derivatives function
                # this will update and modify the collection of derivatives
                new_expr = self.find_and_replace_complex_ders(
                    expr, self.every_var_name, 0, self.idx_str
                )

                # add these new expressions to our list
                new_exprs.append(new_expr)

            # iterate through the staged ones as well
            new_staged_exprs = []
            print("    Finding derivitaves in staged expressions if there are any")
            for ii, expr in enumerate(staged_exp):
                new_expr = self.find_and_replace_complex_ders(
                    expr, self.every_var_name, 0, self.idx_str
                )

                # these new expressions are then "replaced"
                new_staged_exprs.append(new_expr)

            if do_extra_check:
                # go again with the original method
                # collection of "new" (modified) expressions for each of the variables
                new_exprs_again = []
                # the list of all derivatives that we've found that will need to be
                # precalculated
                found_derivatives = []

                for ii, expr in enumerate(new_exprs):
                    print(
                        str(all_rhs_names[ii])
                        + f" {ii + 1}/{len(all_exp)} : {(ii + 1) / len(all_exp):.2%}",
                        file=sys.stderr,
                    )

                    # call the find and replace complicated derivatives function
                    # this will update and modify the collection of derivatives
                    (
                        new_expr,
                        found_derivatives,
                    ) = self.find_and_replace_complex_ders_staged(
                        expr, found_derivatives, 0, self.idx_str
                    )

                    # add these new expressions to our list
                    new_exprs_again.append(new_expr)

                # TODO: implement extra check for staged expressions

            else:
                new_exprs_again = new_exprs
                found_derivatives = []

        else:
            print("NOTE: Override for finding and replacing derivatives was set!")
            new_exprs_again = all_exp
            found_derivatives = []
            new_staged_exprs = staged_exp

        # store them internally so we don't lose them
        self.stored_rhs_function[var_type] = {
            "exprs": new_exprs_again,
            "all_rhs_names": all_rhs_names,
            "found_derivatives": found_derivatives,
            "orig_n_exp": orig_n_exp,
            "staged_exprs": new_staged_exprs,
            "staged_exprs_names": staged_names,
        }

        # return the new expressions with the replacements, the names of the
        # expressions in order, and the found derivatives as well
        # as number of operations
        return (
            new_exprs_again,
            all_rhs_names,
            found_derivatives,
            orig_n_exp,
            new_staged_exprs,
            staged_names,
        )

    @staticmethod
    def get_only_derivatives(expr, var_names, idx_str="[pp]"):
        funcs_to_find = [dendrosym.nr.d, dendrosym.nr.d2s, dendrosym.nr.ad]

        # start by finding all expressions using the atoms funcion
        all_funcs = DendroConfiguration.find_and_sort_atoms(expr, funcs_to_find)

        pass

    def replace_derivatives_with_stencil(self, var_type, stencil_order=6):
        # make sure we've done the derivative mess to expand things out if possible
        (
            exprs,
            rhs_names,
            found_complex_deriv_exprs,
            orig_n_exp,
            staged_exprs,
            staged_exprs_names,
        ) = self.find_derivatives(var_type)

        # then go ahead and go through them and replace them

        new_exprs = []

        for ii, expr in enumerate(exprs):
            print("On expression ", ii)

            new_expr = self.find_and_replace_derivatives_with_stencil(
                expr,
                self.every_var_name,
                stencil_order,
                self.idx_str,
            )

            new_exprs.append(new_expr)

        # update stored_rhs_function's exprs
        self.stored_rhs_function[var_type]["exprs"] = new_exprs

    @staticmethod
    def find_and_replace_derivatives_with_stencil(
        expr, var_names, stencil_order=6, idx_str="[pp]"
    ):
        # lets build up the var names and create all of the stencils and just replace that way
        deriv_funcs_find = []
        deriv_funcs_replace = []
        for xx in var_names:
            if xx.endswith(idx_str):
                xx = xx[: -len(idx_str)]
                # if it ends with the index string, then we want to add it and the non idx string version

            sym_no_idx_str = sym.Symbol(xx)
            sym_w_idx_str = sym.Symbol(xx + idx_str)

            # first order derivatives
            for ii in range(3):
                to_replace = dendrosym.derivs.build_first_deriv_stencil_by_dim_idx(
                    xx, stencil_order, ii, idx_str=idx_str
                )

                deriv_funcs_find.append(dendrosym.nr.d(ii, sym_no_idx_str))
                deriv_funcs_replace.append(to_replace)
                deriv_funcs_find.append(dendrosym.nr.d(ii, sym_w_idx_str))
                deriv_funcs_replace.append(to_replace)

            # then second order derivatives
            for ii in range(3):
                for jj in range(3):
                    to_replace = dendrosym.derivs.build_second_deriv_stencil_by_dim_idx(
                        xx, stencil_order, ii, jj, idx_str=idx_str
                    )

                    deriv_funcs_find.append(dendrosym.nr.d2(ii, jj, sym_no_idx_str))
                    deriv_funcs_replace.append(to_replace)
                    deriv_funcs_find.append(dendrosym.nr.d2(ii, jj, sym_w_idx_str))
                    deriv_funcs_replace.append(to_replace)

        all_repl = dict(zip(deriv_funcs_find, deriv_funcs_replace))
        reverse_repl = dict(zip(deriv_funcs_replace, deriv_funcs_find))

        expr = expr.subs(all_repl)

        return expr

    @staticmethod
    def find_and_replace_complex_ders_staged(
        expr, found_derivatives, depth=0, idx_str="[pp]"
    ):
        funcs_to_find = [dendrosym.nr.d, dendrosym.nr.d2s, dendrosym.nr.ad]

        # start by finding all expressions using the atoms funcion
        all_funcs = DendroConfiguration.find_and_sort_atoms(expr, funcs_to_find)

        while len(all_funcs) > 0:
            # get the first element in the functions list
            func = all_funcs.pop(0)

            term_differentiate = func.args[-1]

            if isinstance(term_differentiate, sym.Symbol):
                if not term_differentiate.name.startswith("DENDRO_STAGED_"):
                    continue

            elif (
                isinstance(term_differentiate, sym.Function)
                or len(term_differentiate.atoms(*funcs_to_find)) > 0
            ):
                # if there are any sub expressions in this one, then we need
                # otherwise, we need to store the operation here and move on
                needs_to_go_deeper = False
                if len(term_differentiate.atoms(*funcs_to_find)) > 0:
                    needs_to_go_deeper = True

                if needs_to_go_deeper:
                    (
                        mini_expr,
                        found_derivatives,
                    ) = DendroConfiguration.find_and_replace_complex_ders_staged(
                        term_differentiate,
                        found_derivatives,
                        depth + 1,
                        idx_str=idx_str,
                    )

                    # replace it within the expression
                    expr = expr.xreplace({term_differentiate: mini_expr})

                    # then reupdate the all functions and continue
                    all_funcs = DendroConfiguration.find_and_sort_atoms(
                        expr, funcs_to_find
                    )

                    continue

            # otherwise we're good to check if this already exists in our list
            if func.name == "grad":
                # generate a staged gradient name for it
                temp_var_name = sym.Symbol(
                    "DENDRO_STAGED_GRAD_" + f"{len(found_derivatives):03d}" + idx_str
                )

                func_args = func.args
                term_to_differentiate = func_args[1]
                index_order = func_args[0]

            elif func.name == "grad2":
                # generate a staged gradient name for it
                temp_var_name = sym.Symbol(
                    "DENDRO_STAGED_GRAD2_" + f"{len(found_derivatives):03d}" + idx_str
                )

                func_args = func.args
                term_to_differentiate = func_args[2]
                index_order = (func_args[0], func_args[1])

            elif func.name == "agrad":
                # generate a staged gradient name for it
                temp_var_name = sym.Symbol(
                    "DENDRO_STAGED_AGRAD_" + f"{len(found_derivatives):03d}" + idx_str
                )

                func_args = func.args
                term_to_differentiate = func_args[1]
                index_order = func_args[0]

            else:
                # anything else we'll continue
                continue

            # check to see if the term already exists in our collection
            (
                already_exists,
                found_var_name,
            ) = DendroConfiguration.find_repeat_derivative_terms(
                found_derivatives,
                term_to_differentiate,
                index_order,
                func.name,
            )

            if not already_exists:
                # if it didn't already exist
                found_derivatives.append(
                    {
                        "operation": func.name,
                        "orig_exp": term_to_differentiate,
                        "index_order": index_order,
                        "temp_var_name": temp_var_name,
                        "depth": depth,
                    }
                )
            else:
                # if it *did* exist, then we replace the temp_var_name
                # with the one that was found
                temp_var_name = found_var_name

            # go ahead and replace the found function with the piece that we found
            expr = expr.xreplace({func: temp_var_name})

            # refresh the list of functions we're interested in finding
            all_funcs = DendroConfiguration.find_and_sort_atoms(expr, funcs_to_find)

        return expr, found_derivatives

    @staticmethod
    def find_and_replace_complex_ders(expr, var_names, depth=0, idx_str="[pp]"):
        funcs_to_find = [dendrosym.nr.d, dendrosym.nr.d2s, dendrosym.nr.ad]

        # start by finding all expressions using the atoms funcion
        all_funcs = DendroConfiguration.find_and_sort_atoms(expr, funcs_to_find)

        if len(all_funcs) == 0:
            print(f"        no derivatives to assess here!")
            return expr

        # print("found", all_funcs)
        print(f"        found {len(all_funcs)} to assess...")

        something_changed = False

        while len(all_funcs) > 0:
            # get the first element in the functions list
            func = all_funcs.pop(0)

            term_differentiate = func.args[-1]

            if isinstance(term_differentiate, sym.Symbol):
                if term_differentiate.name.startswith("DENDRO_STAGED_"):
                    continue

            elif (
                isinstance(term_differentiate, sym.Function)
                or len(term_differentiate.atoms(*funcs_to_find)) > 0
            ):
                # if there are any sub expressions in this one, then we need
                # otherwise, we need to store the operation here and move on
                needs_to_go_deeper = False
                if len(term_differentiate.atoms(*funcs_to_find)) > 0:
                    needs_to_go_deeper = True

                if needs_to_go_deeper:
                    mini_expr = DendroConfiguration.find_and_replace_complex_ders(
                        term_differentiate,
                        var_names,
                        depth + 1,
                        idx_str=idx_str,
                    )

                    # replace it within the expression
                    expr = expr.xreplace({term_differentiate: mini_expr})

                    # then reupdate the all functions and continue
                    all_funcs = DendroConfiguration.find_and_sort_atoms(
                        expr, funcs_to_find
                    )

                    continue

            # symbols we need and will replace
            x, y, z = sym.symbols("xtem ytem ztem")
            d_order = (x, y, z)  # the index goes in this order

            # any sort of symbol that we want to find should also check if it has strings
            symbols_find = []
            symbols_replace = []
            for xx in var_names:
                if xx.endswith(idx_str):
                    xx = xx[: -len(idx_str)]
                    # if it ends with the index string, then we want to add it and the non idx string version

                symbols_find.append(sym.Symbol(xx))
                symbols_find.append(sym.Symbol(xx + idx_str))
                # replace it with:
                to_replace = sym.Function(sym.Symbol(xx))(x, y, z)
                # add twice for the lineup
                symbols_replace.append(to_replace)
                symbols_replace.append(to_replace)

            # symbols_find = [sym.Symbol(xx + idx_str) for xx in var_names]
            # symbols_replace = [
            #     sym.Function(sym.Symbol(xx))(x, y, z) for xx in var_names
            # ]
            all_repl = dict(zip(symbols_find, symbols_replace))
            reverse_repl = dict(zip(symbols_replace, symbols_find))

            # otherwise we're good to check if this already exists in our list
            if func.name == "grad":
                # FIRST ORDER DERIVATIVE
                # remember: grad(x, y) is what's incoming here...
                # x is the "dimension" over which it is being differentiated
                # y is the actual term to differentiate

                func_args = func.args
                term_to_differentiate = func_args[1]
                index_order = func_args[0]

                # also need to find pre-existing derivatives that have been replaced
                term_to_differentiate = term_to_differentiate.subs(all_repl)

                # then take the derivate of that term over the index_order
                diff_term = sym.diff(term_to_differentiate, d_order[index_order])

            elif func.name == "grad2":
                # generate a staged gradient name for it
                func_args = func.args
                term_to_differentiate = func_args[2]
                index_order = (func_args[0], func_args[1])

                term_to_differentiate = term_to_differentiate.subs(all_repl)

                # then take the second order derivative in that direction
                diff_term = sym.diff(
                    term_to_differentiate,
                    d_order[index_order[0]],
                    d_order[index_order[1]],
                )

            elif func.name == "agrad":
                print(
                    "    WARNING: AN EXPRESSION USING AGRAD WAS FOUND "
                    + "AND CANNOT BE SIMPLIFIED, PLEASE CONSIDER REWRITING"
                )
                continue

            else:
                # anything else we'll continue
                continue

            # go ahead and replace the found function with the piece that we found
            expr = expr.xreplace({func: diff_term})

            something_changed = True

            # refresh the list of functions we're interested in finding
            all_funcs = DendroConfiguration.find_and_sort_atoms(expr, funcs_to_find)

        if depth == 0:
            # now need to replace the sympy derivatives with the grad and grad2 functions again for code gen
            expr = DendroConfiguration.replace_sympy_derivatives(
                expr, var_names, idx_str
            )

            if something_changed:
                # make sure everything else is replaced in the expression
                expr = expr.xreplace(reverse_repl)

        return expr

    @staticmethod
    def replace_sympy_derivatives(expr, var_names, idx_str):
        print("  Now attempting to replace derivatives...")
        # symbols we need and will replace
        x, y, z = sym.symbols("xtem ytem ztem")
        d_order = (x, y, z)  # the index goes in this order

        symbols_find = [sym.Symbol(xx + idx_str) for xx in var_names]
        symbols_replace = [sym.Function(sym.Symbol(xx))(x, y, z) for xx in var_names]
        all_repl = dict(zip(symbols_find, symbols_replace))
        reverse_repl = dict(zip(symbols_replace, symbols_find))

        # we need to replace all derivative terms with variables
        all_derivatives = expr.atoms(sym.Derivative)

        while len(all_derivatives) > 0:
            # print("    Length of Remaining Derivatives: ", len(all_derivatives))

            curr_deriv = all_derivatives.pop()

            if len(curr_deriv.args) == 2:
                # this is only in one single direction, so what will matter is the "diff_dir"
                diff_term, diff_dir = curr_deriv.args

                if diff_dir[1] > 2:
                    raise NotImplemented(
                        "A 3rd or higher order derivative was found after expanding!"
                    )
                idx_use = -1

                for ii, val in enumerate(d_order):
                    if diff_dir[0] == val:
                        idx_use = ii

                if diff_dir[1] == 1:
                    new_expr = dendrosym.nr.d(idx_use, diff_term)
                else:
                    new_expr = dendrosym.nr.d2s(idx_use, idx_use, diff_term)

            elif len(curr_deriv.args) == 3:
                # this is a multi-directional derivative, so we use grad2 in replacement
                diff_term, diff_dir1, diff_dir2 = curr_deriv.args

                if diff_dir1[1] > 1 or diff_dir2[1] > 1:
                    raise NotImplemented(
                        "A 3rd or higher order derivative was found after expanding!"
                    )

                idx_1 = -1
                for ii, val in enumerate(d_order):
                    if diff_dir1[0] == val:
                        idx_1 = ii

                idx_2 = -1
                for ii, val in enumerate(d_order):
                    if diff_dir2[0] == val:
                        idx_2 = ii

                # swap them if they're out of order, idx_1 should be smaller for efficiency
                idxtmp = 0
                if idx_1 > idx_2:
                    print("Swapping the indices of derivatives:", idx_1, idx_2)
                    idxtmp = idx_1
                    idx_1 = idx_2
                    idx_2 = idxtmp

                new_expr = dendrosym.nr.d2s(idx_1, idx_2, diff_term)

            else:
                raise NotImplemented(
                    "A more complicated derivative was found requiring 3rd or higher order!"
                )

            new_expr = new_expr.subs(reverse_repl)

            expr = expr.xreplace({curr_deriv: new_expr})

            all_derivatives = expr.atoms(sym.Derivative)

        return expr

    @staticmethod
    def find_and_sort_atoms(expr, objs: list):
        # TODO: fix this if expression isn't a sympy object

        # find the atoms
        atoms = np.array(list(expr.atoms(*objs)))

        num_sub_exprs = []

        # find the total number of sub expressions
        for ii, at in enumerate(atoms):
            num_sub_exprs.append(len(at.atoms(*objs)))

        num_sub_exprs = np.array(num_sub_exprs)
        # then we arg sort and flip so it's largest to shortest
        sorted_idxs = np.flip(np.argsort(num_sub_exprs))

        # then we return the sorted list
        return atoms[sorted_idxs].tolist()

    @staticmethod
    def find_repeat_derivative_terms(
        found_derivatives, term_to_differentiate, index_order, operation
    ):
        # iterate through the found derivatives
        for der_info in found_derivatives:
            # if the operation doesn't match, no need to continue
            if der_info["operation"] != operation:
                continue

            if (
                sym.simplify(der_info["orig_exp"] - term_to_differentiate) == 0
                and der_info["index_order"] == index_order
            ):
                return True, der_info["temp_var_name"]

        return False, None

    def generate_pre_necessary_derivatives(
        self, var_type, dtype="double", include_byte_declaration=False
    ):
        (
            exprs,
            all_rhs_names,
            found_derivatives,
            orig_n_exp,
            new_staged_exprs,
            staged_names,
        ) = self.find_derivatives(var_type)

        if include_byte_declaration:
            # get the number of bytes we need for allocation
            outstr = f"const unsigned int bytes = n * sizeof({dtype});\n\n"
        else:
            outstr = ""

        if len(found_derivatives) == 0:
            return outstr + "// NO INTERMEDIATE DERIVATIVES FOUND\n", ""

        # add a declaration for computing the boundaries based on the boundary flag
        outstr += "// initializing start and end points for the grids\n"
        outstr += "unsigned int kstart = 0, jstart = 0, istart = 0;\n"
        outstr += "unsigned int kend = nz, jend = ny, iend = nx;\n\n"

        outstr += "if (bflag & (1u << OCT_DIR_LEFT)) {\n    "
        outstr += "istart += PW;\n}\n\n"
        outstr += "if (bflag & (1u << OCT_DIR_RIGHT)) {\n    "
        outstr += "iend -= PW;\n}\n\n"
        outstr += "if (bflag & (1u << OCT_DIR_DOWN)) {\n    "
        outstr += "jstart += PW;\n}\n\n"
        outstr += "if (bflag & (1u << OCT_DIR_UP)) {\n    "
        outstr += "jend -= PW;\n}\n\n"
        outstr += "if (bflag & (1u << OCT_DIR_BACK)) {\n    "
        outstr += "kstart += PW;\n}\n\n"
        outstr += "if (bflag & (1u << OCT_DIR_FRONT)) {\n    "
        outstr += "kend -= PW;\n}\n\n"

        out_dealloc = ""

        allocate_str = "// Allocating memory for STAGED variables\n"
        # allocation for temporary intermediate steps that need to be
        # calculated before calculating the derivatives
        allocate_tmp_str = "// " + "Allocating memory for intermediate calculations\n"

        deallocate_str = "// Deallocating memory for STAGED variables\n"
        # TODO: do we really want to deallocate right after calculating?
        deallocate_tmp_str = (
            "// " + "Deallocating memory for intermediate calculations\n"
        )

        t = "    "

        calculation_str = "/**\n * CALCULATING INTERMEDIATE EXPRESSIONS\n" + " */\n"
        calculation_str += "for (unsigned int k = kstart; k < kend; k++)\n"
        calculation_str += t + "{\n"
        calculation_str += t + "for (unsigned int j = jstart; j < jend; j++)\n"
        calculation_str += t + "{\n"
        calculation_str += t * 2 + "for (unsigned int i = istart; i < iend; i++)\n"
        calculation_str += t * 2 + "{\n"
        calculation_str += t * 3 + f"const {dtype} x = pmin[0] + i * hx;\n"
        calculation_str += t * 3 + f"const {dtype} y = pmin[1] + j * hy;\n"
        calculation_str += t * 3 + f"const {dtype} z = pmin[2] + k * hz;\n\n"
        calculation_str += t * 3 + "const unsigned int pp = i + nx * (j + ny * k);\n\n"
        # TODO: add eta check????

        # then create the string for calculating the derivatives
        deriv_str = "// Staged gradient calculation\n"

        # sort the found derivatives somewhat via bubble sort based on *type*
        # of derivative calculation. We want to do grad before grad2, for example
        found_derivatives = self.sort_found_derivatives(found_derivatives)

        # now we need to iterate through the found derivatives and start with largest
        # "depth" field
        max_depth = self.find_max_depth_value(found_derivatives)
        # DEBUG:
        # print("Max Depth is:" + str(max_depth), file=sys.stderr)

        already_calculated_inters = []
        already_added_derivs = []

        final_intermediate_dealloc = ""

        for curr_depth in reversed(range(max_depth + 1)):
            curr_all = "" + allocate_str
            curr_temp_all = "" + allocate_tmp_str
            curr_deall = "" + deallocate_str
            curr_temp_deall = "" + deallocate_tmp_str
            curr_calc = "" + calculation_str
            curr_deriv = "" + deriv_str

            ders_remove = []

            for ii, deriv_info in enumerate(found_derivatives):
                if deriv_info["depth"] != curr_depth:
                    continue

                (
                    already_included,
                    inter_var_name,
                ) = self.find_if_inter_processed(
                    already_calculated_inters, deriv_info, self.idx_str
                )

                # otherwise, get the information from the function
                (
                    all_str,
                    all_tmp_str,
                    calc_str,
                    deriv_str,
                    deall_str,
                    deall_tmp_str,
                ) = self.gen_individual_der_strs(
                    deriv_info,
                    inter_var_name,
                    already_calculated_inters=already_added_derivs,
                )

                curr_all += all_str
                curr_temp_all += all_tmp_str
                curr_deall += deall_str
                curr_temp_deall += deall_tmp_str
                curr_calc += calc_str
                curr_deriv += deriv_str

                if not already_included:
                    already_calculated_inters.append(deriv_info)
                already_added_derivs.append(deriv_info)

                ders_remove.append(ii)

            # then we cap off the calculation string
            curr_calc += t * 2 + "}\n" + t * 1 + "}\n" + "}\n"

            # then append them all together
            if curr_depth > 0:
                outstr += (
                    "//\n// CALCULATING EXPRESSIONS THAT INCLUDE"
                    + " DERIVATIVES WITHIN\n"
                    + f"// CURRENT DEPTH: {curr_depth}\n"
                )
            outstr += curr_all + "\n"
            outstr += curr_temp_all + "\n"
            outstr += curr_calc + "\n"
            outstr += curr_deriv + "\n"
            final_intermediate_dealloc += curr_temp_deall + "\n"

            out_dealloc += curr_deall + "\n"

            # and remove the set of derivatives
            if max_depth > 0:
                for idx in sorted(ders_remove, reverse=True):
                    del found_derivatives[idx]

        # then stitch on the intermediate deallocation string at the end
        outstr += final_intermediate_dealloc

        return outstr, out_dealloc

    @staticmethod
    def find_max_depth_value(deriv_info):
        max_depth = 0

        for der in deriv_info:
            if der["depth"] > max_depth:
                max_depth = der["depth"]

        return max_depth

    def gen_individual_der_strs(
        self,
        deriv_info,
        inter_var_name=None,
        already_calculated_inters=[],
        dtype="double",
        t="    ",
    ):
        # little string used to map direction index to real direction
        idx_to_dir = {0: "x", 1: "y", 2: "z"}

        var_name = str(deriv_info["temp_var_name"])
        # remove the indexing string
        var_name_no_idx = var_name.replace(self.idx_str, "")

        # allocate memory string addition
        allocate_str = f"{dtype} *{var_name_no_idx}" + f" = ({dtype} *)malloc(bytes);\n"

        # if intervar name isn't set, then we need to include it for
        # calculations
        if inter_var_name is None:
            inter_var_name = var_name_no_idx + "_intermediate"
            # intermediate allocation string addition
            allocate_tmp_str = (
                f"{dtype} *{inter_var_name}" + f" = ({dtype} *)malloc(bytes);\n"
            )

            # then the calculation string addition
            calculation_str = (
                t * 3
                + inter_var_name
                + self.idx_str
                + " = "
                + sym.ccode(deriv_info["orig_exp"])
                + ";\n"
            )

            # then deallocate the temporary var
            deallocate_tmp_str = f"free({var_name_no_idx}_intermediate);\n"
        else:
            allocate_tmp_str = ""
            calculation_str = ""
            deallocate_tmp_str = ""

        # then we do the derivative string based on the type of derivative
        if deriv_info["operation"] == "grad":
            func = f"deriv_{idx_to_dir[deriv_info['index_order']]}("
            deriv_str = func + var_name_no_idx
            deriv_str += ", " + inter_var_name
            deriv_str += ", h" + idx_to_dir[deriv_info["index_order"]]
            deriv_str += ", sz, bflag);\n"
        elif deriv_info["operation"] == "grad2":
            # now we calculate the derivative order
            if deriv_info["index_order"][0] == deriv_info["index_order"][1]:
                func = f"deriv_{idx_to_dir[deriv_info['index_order'][0]] * 2}("
                deriv_str = func + var_name_no_idx
                deriv_str += ", " + inter_var_name
                deriv_str += ", h" + idx_to_dir[deriv_info["index_order"][0]]
                deriv_str += ", sz, bflag);\n"
            else:
                # so we have to start by checking if the first direction was already calculated
                dir1, dir2 = deriv_info["index_order"]

                # inter_var_name
                found_first_dir = False
                # base string so that the compiler throws and error
                intermediate_value_name = "INVALID_VARIABLE_NAME_DO_NOT_USE"
                for completed_deriv_info in already_calculated_inters:
                    if type(completed_deriv_info["index_order"]) == tuple:
                        # we have a tuple for the index order
                        completed_dir1 = completed_deriv_info["index_order"][0]
                    else:
                        completed_dir1 = completed_deriv_info["index_order"]

                    # now if it's the same direction
                    if completed_dir1 == dir1:
                        # then we check the operation is grad 1
                        if completed_deriv_info["operation"] == "grad":
                            # then we check if the
                            if (
                                sym.simplify(
                                    completed_deriv_info["orig_exp"]
                                    - deriv_info["orig_exp"]
                                )
                                == 0
                            ):
                                found_first_dir = True
                                # get the intermediate value name, and then we're golden
                                intermediate_value_name = str(
                                    completed_deriv_info["temp_var_name"]
                                ).split(self.idx_str)[0]
                                break

                if found_first_dir:
                    deriv_str = f"deriv_{idx_to_dir[dir2]}("
                    deriv_str += var_name_no_idx
                    deriv_str += ", " + intermediate_value_name
                    deriv_str += ", h" + idx_to_dir[dir2]
                    deriv_str += ", sz, bflag);\n"

                else:
                    # TODO: handle edge case where we have some weird derivative that wasn't already calculated
                    deriv_str = f"TODO: {deriv_info['index_order']}, {inter_var_name}\n"

        elif deriv_info["operation"] == "agrad":
            adv_der_var_use = getattr(self, "advective_der_var", "UNDEF") + str(
                deriv_info["index_order"]
            )

            func = f"adv_deriv_{idx_to_dir[deriv_info['index_order']]}("
            deriv_str = func + var_name_no_idx
            deriv_str += ", " + inter_var_name
            deriv_str += ", h" + idx_to_dir[deriv_info["index_order"]]
            deriv_str += f", sz, {adv_der_var_use}, bflag);\n"

        # then we add the deallocation
        deallocate_str = f"free({var_name_no_idx});\n"

        return (
            allocate_str,
            allocate_tmp_str,
            calculation_str,
            deriv_str,
            deallocate_str,
            deallocate_tmp_str,
        )

    @staticmethod
    def find_if_inter_processed(curr_list, incoming_der_info, idx_str):
        for curr_test in curr_list:
            if sym.simplify(curr_test["orig_exp"] - incoming_der_info["orig_exp"]) == 0:
                inter_var_name = str(curr_test["temp_var_name"])
                inter_var_name = inter_var_name.replace(idx_str, "")
                return True, inter_var_name + "_intermediate"

        return False, None

    @staticmethod
    def sort_found_derivatives(found_derivatives):
        # quick bubble sort to put all grad2 operations after grad1's
        for mx in range(len(found_derivatives) - 1, -1, -1):
            swapped = False
            for i in range(mx):
                if (
                    found_derivatives[i]["operation"] == "grad2"
                    and found_derivatives[i + 1]["operation"] == "grad"
                ):
                    found_derivatives[i], found_derivatives[i + 1] = (
                        found_derivatives[i + 1],
                        found_derivatives[i],
                    )
                    swapped = True
            if not swapped:
                break

        return found_derivatives

    def generate_deriv_allocation_and_calc(
        self,
        var_type="evolution",
        include_byte_declaration=False,
        use_old_method=False,
    ):
        """Generates all of the C++ code for allocation and calculation of derivatives"""

        # get the RHS stuff
        temp_funcs = self.stored_rhs_function[var_type]
        grad_list, grad2_list, agrad_list = self.find_all_unique_ders(
            temp_funcs["exprs"], temp_funcs["all_rhs_names"]
        )

        self.current_index = 0

        # gradient level 1
        if var_type == "evolution":
            # NOTE: for evolution, we need every derivitive **no matter what**
            # this is because we need to calculate KO diss.
            grad_alloc = self.gen_grad_memory_alloc(
                var_type,
                "grad",
                include_byte_declaration=include_byte_declaration,
            )
            grad_calc = self.gen_grad_calculations(var_type, "grad")
            grad_dealloc = self.gen_grad_memory_dealloc(var_type, "grad")

            orig_vars = self.all_var_names.get(var_type, [])
            grad_grad_names = self.create_grad_var_names(orig_vars, "grad", 3)
        else:
            (
                grad_alloc,
                grad_dealloc,
                grad_calc,
                grad_grad_names,
                start_index,
            ) = self.create_func_list_code(
                grad_list, self.idx_str, start_index=self.current_index
            )
            self.current_index = start_index

        # grad 2
        (
            grad2_alloc,
            grad2_dealloc,
            grad2_calc,
            grad2_grad_names,
            start_index,
        ) = self.create_func_list_code(
            grad2_list,
            self.idx_str,
            grad1_list=grad_grad_names,
            start_index=self.current_index,
        )
        self.current_index = start_index

        # agrad
        (
            agrad_alloc,
            agrad_dealloc,
            agrad_calc,
            agrad_grad_names,
            start_index,
        ) = self.create_func_list_code(
            agrad_list,
            self.idx_str,
            agrad_var=getattr(self, "advective_der_var", "NONE"),
            start_index=self.current_index,
        )
        self.current_index = start_index

        # stitch them together and return the strings
        all_alloc = grad_alloc + grad2_alloc + agrad_alloc
        all_calc = grad_calc + grad2_calc + agrad_calc
        all_dealloc = grad_dealloc + grad2_dealloc + agrad_dealloc

        if not use_old_method:
            all_dealloc = "// NO DEALLOCATION REQUIRED DUE TO DERIVATIVE BASE"

        return all_alloc, all_calc, all_dealloc

    @staticmethod
    def create_func_list_code(
        input_funcs,
        idx_str,
        dtype="double",
        grad1_list=[],
        agrad_var="beta",
        use_old_method=False,
        start_index=0,
    ):
        # if it's not an input
        if not isinstance(input_funcs, set) and not isinstance(input_funcs, list):
            input_funcs = [input_funcs]

        out_alloc = ""
        out_dealloc = ""
        out_calc = ""
        all_grad_names = []

        dir_dict = {0: "x", 1: "y", 2: "z"}

        for fun in input_funcs:
            if fun.name == "grad":
                direction, var = fun.args

                var_name = str(var).split(idx_str)[0]

                grad_name = f"grad_{direction}_{var_name}"

                calc_str = (
                    f"deriv_{dir_dict[direction]}("
                    + grad_name
                    + ", "
                    + var_name
                    + ", "
                    + "h"
                    + dir_dict[direction]
                    + ", sz, bflag);\n"
                )

                all_grad_names.append(grad_name)

            elif fun.name == "grad2":
                dir1, dir2, var = fun.args

                var_name = str(var).split(idx_str)[0]

                grad_name = f"grad2_{dir1}_{dir2}_{var_name}"

                # now the calculation is a bit different
                if dir1 == dir2:
                    calc_str = "deriv_" + 2 * dir_dict[dir1]
                    calc_str += "(" + grad_name + ", " + var_name
                    calc_str += ", h" + dir_dict[dir1] + ", sz, bflag);\n"

                # if they are not equal, we need to do a few things
                else:
                    first_dir_name = f"grad_{dir1}_{var_name}"
                    # first check if the first direction is in our list of grad1
                    if first_dir_name in grad1_list:
                        # then we can just add the second direction
                        calc_str = f"deriv_{dir_dict[dir2]}("
                        calc_str += (
                            grad_name
                            + ", "
                            + f"grad_{dir1}_{var_name}"
                            + ", "
                            + "h"
                            + dir_dict[dir2]
                            + ", sz, bflag);\n"
                        )

                        all_grad_names.append(grad_name)
                    else:
                        # then we need to add to the allocation string
                        # a temporary variable
                        # NOTE: this is a fallback in case it wasn't actually made
                        if use_old_method:
                            out_alloc += f"{dtype} *{first_dir_name} = ({dtype} *)malloc(bytes);\n"
                            out_dealloc += f"free({first_dir_name});\n"
                        else:
                            out_alloc += f"{dtype} *{first_dir_name} = deriv_base + {start_index} * BLK_SZ;\n"
                            start_index += 1

                        # then the calculation string
                        calc_str = (
                            f"deriv_{dir_dict[dir1]}("
                            + first_dir_name
                            + ", "
                            + var_name
                            + ", "
                            + "h"
                            + dir_dict[dir1]
                            + ", sz, bflag);\n"
                        )

                        calc_str += (
                            f"deriv_{dir_dict[dir2]}("
                            + grad_name
                            + ", "
                            + first_dir_name
                            + ", "
                            + "h"
                            + dir_dict[dir2]
                            + ", sz, bflag);\n"
                        )

            elif fun.name == "agrad":
                direction, var = fun.args

                var_name = str(var).split(idx_str)[0]

                grad_name = f"agrad_{direction}_{var_name}"

                calc_str = (
                    "adv_deriv_"
                    + dir_dict[direction]
                    + "("
                    + grad_name
                    + ", "
                    + var_name
                    + ", h"
                    + dir_dict[direction]
                    + ", sz, "
                    + agrad_var
                    + str(direction)
                    + ", bflag);\n"
                )

            # grad allocation
            if use_old_method:
                alloc_str = f"{dtype} *{grad_name} = ({dtype} *)malloc(bytes);\n"

                # grad deallocation
                dealloc_str = f"free({grad_name});\n"
            else:
                alloc_str = (
                    f"{dtype} *{grad_name} = deriv_base + {start_index} * BLK_SZ;\n"
                )
                start_index += 1
                dealloc_str = ""

            out_alloc += alloc_str
            out_dealloc += dealloc_str
            out_calc += calc_str

        return out_alloc, out_dealloc, out_calc, all_grad_names, start_index

    def get_unique_ders_from_var_type(self, var_type):
        (
            all_exp,
            all_rhs_names,
            staged_exp,
            staged_names,
            orig_n_exp,
        ) = self.get_rhs_functions_all(var_type, True)

        # now that we have them, we find all unique derivatives below

        grad_list, grad2_list, agrad_list = self.find_all_unique_ders(
            all_exp, all_rhs_names
        )

        (
            grad_list_staged,
            grad2_list_staged,
            agrad_list_staged,
        ) = self.find_all_unique_ders(staged_exp, all_rhs_names)

        # count all of the unique ones
        uq_grads = set(grad_list)
        uq_grads.update(grad_list_staged)

        uq_grad2s = set(grad2_list)
        uq_grad2s.update(grad2_list_staged)

        uq_agrads = set(agrad_list)
        uq_agrads.update(agrad_list_staged)

        return uq_grads, uq_grad2s, uq_agrads

    def get_num_unique_derivs(self, var_type):
        uq_grads, uq_grad2s, uq_agrads = self.get_unique_ders_from_var_type(var_type)

        return len(uq_grads) + len(uq_grad2s) + len(uq_agrads)

    # TODO: potentially remove this from static and make it real
    @staticmethod
    def find_all_unique_ders(rhs_funcs, rhs_names, sort=True):
        # NOTE: this is assumed to be after the complicated ones are
        # gathered
        # this will then help us keep track of all of the different derivatives
        # we have to calculate

        funcs_to_find = [dendrosym.nr.d, dendrosym.nr.d2s, dendrosym.nr.ad]
        uses_orig_deriv = True
        uses_orig_deriv2 = True

        if type(dendrosym.nr.d) is not sym.core.function.UndefinedFunction:
            # if we're using the built-in derivative stuff, we need to check for that
            funcs_to_find[0] = sym.Derivative
            uses_orig_deriv = False
        else:
            # if it's not a derivative, we want to check for advective derivative
            # and derivative being the same
            if funcs_to_find[0].name == funcs_to_find[2].name:
                # setting it to a random string that we will certainly never find!
                funcs_to_find[2] = sym.Function("7OobSHkNgRuNhk7M01yu")

        # also need to replace d2s
        if type(dendrosym.nr.d2s) is not sym.core.function.UndefinedFunction:
            # if we're using the built-in derivative stuff, we need to check for that
            funcs_to_find[1] = sym.Derivative
            uses_orig_deriv2 = False

        grad_list = set()
        grad2_list = set()
        agrad_list = set()

        for ii, eqn in enumerate(rhs_funcs):
            # this creates sets for all of the gradient atoms
            all_current_grads = eqn.atoms(funcs_to_find[0])
            all_current_grad2s = eqn.atoms(funcs_to_find[1])
            all_current_agrads = eqn.atoms(funcs_to_find[2])

            # print(all_current_grad2s)

            if not uses_orig_deriv and not uses_orig_deriv2:
                updated_grads = set()
                updated_grad2s = set()

                for the_grad in all_current_grads:
                    temp = dendrosym.derivs.restore_original_derivatives(the_grad)

                    if temp.name == "grad":
                        updated_grads.update(
                            {
                                temp,
                            }
                        )
                    elif temp.name == "grad2":
                        idx1, idx2, symbol = temp.args
                        if idx1 > idx2:
                            idxtmp = temp.args[1]
                            temp.args[1] = temp.args[0]
                            temp.args[0] = idxtmp
                        updated_grad2s.update(
                            {
                                temp,
                            }
                        )
                    else:
                        raise ValueError("Invalid type of derivative coming in!")

                grad_list.update(updated_grads)
                grad2_list.update(updated_grad2s)
                agrad_list.update(all_current_agrads)

            elif uses_orig_deriv and not uses_orig_deriv2:
                raise NotImplementedError(
                    "Mixed derivative types is not currently allowed,"
                    + " make sure both are either sym.Derivative or 'original' function"
                )
            elif not uses_orig_deriv and uses_orig_deriv2:
                raise NotImplementedError(
                    "Mixed derivative types is not currently allowed,"
                    + " make sure both are either sym.Derivative or 'original' function"
                )
            else:
                # for grad2s, we need to make sure the args are in the right order
                updated_grad2s = set()
                for the_grad2 in all_current_grad2s:
                    idx1, idx2, symbol = the_grad2.args
                    if idx1 > idx2:
                        idxtmp = idx1
                        idx2 = idx1
                        idx1 = idxtmp
                    updated_grad2s.update(
                        {
                            dendrosym.nr.d2s(idx1, idx2, symbol),
                        }
                    )

                # print(updated_grad2s)

                # then add them to our set which forces uniqueness, but not order
                grad_list.update(all_current_grads)
                grad2_list.update(updated_grad2s)
                agrad_list.update(all_current_agrads)

        # then let's sort them
        if sort:
            grad_list = list(grad_list)
            grad_list.sort(key=lambda x: x.args[1].name)

            grad2_list = list(grad2_list)
            grad2_list.sort(key=lambda x: x.args[2].name)

            agrad_list = list(agrad_list)
            agrad_list.sort(key=lambda x: x.args[1].name)

        return grad_list, grad2_list, agrad_list

    @property
    def every_var_name(self):
        all_vars = []
        for key, the_vars in self.all_var_names.items():
            # ignore parameters
            if key == "parameter":
                pass
            all_vars += the_vars
        return all_vars

    # def add_initial_data(init_func, name="", desc=""):
    #     # NOTE: just like the RHS function part, we need the
    #     # init function to handle returning a list of expressions
    #     # for each of the variables

    #     # the name is for when we save a file and create a function
    #     # for the initial data. The description is to give some
    #     # helpful info about how things

    #     # NOTE: I also need to decide if we're going to be using
    #     # x or xx. The general ctx files are using xx, so I'll just
    #     # consider that

    #     pass

    def add_staged_function(self, var_type, expr_func):
        self.stored_staged_exprs.update({var_type: expr_func})
