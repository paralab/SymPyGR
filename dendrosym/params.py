"""Parameter file generation functions

This file contains the functions that can generate entire parameter
files based on information provided by the user. This simplifies adding
a parameter to a project because it can place it exactly where it needs
to go while also providing a routine for reading the configurable values
during runtime.

Please see each function for more information, but the most important
ones for preparing this generation are `generate_all_parameter_text`
and `generate_sample_config_file_text`. The first generates all of the
header and source file code while the second creates a sample configuration
file to be used during runtime.
"""

import os
import tomlkit as toml

TAB = "    "

BHOLE_VARS = [
    "MASS", "X", "Y", "Z", "V_X", "V_Y", "V_Z", "SPIN", "SPIN_THETA",
    "SPIN_PHI"
]


def get_toml_data(filename):
    """Function that can read a TOML file and create a table

    The main purpose of this function is to make it easier to
    load a TOML file in Python. Instead of having to write the
    `open` block of code, you can just pass through a filename
    and it will read the entire file.

    Parameters
    ----------
    filename : str
        The input file name for the TOML data to load.

    Returns
    -------
    tomlkit.Table
        A "wrapped" dict object that contains the key and
        item pairs found inside the TOML file.
    """
    with open(filename, "r") as in_file:
        toml_data = toml.parse(in_file.read())

    return toml_data


def get_rank_npes(n_tabs=2, tab_char="    "):
    """Generates the code for getting MPI rank and size

    This is to keep the block clean later on, but can also
    handle having a different number of indents with custom
    tab characters.

    Parameters
    ----------
    n_tabs : int, optional
        The indentation size that you would like for this block.
        Defaults to two.
    tab_char : str, optional
        The character to use for the indentation. Defaults to
        four spaces.

    Returns
    -------
    tomlkit.Table
        A "wrapped" dict object that contains the key and
        item pairs found inside the TOML file.
    """
    t = tab_char * n_tabs
    outstr = f"{t}int rank, npes;\n"
    outstr += f"{t}MPI_Comm_rank(comm, &rank);\n"
    outstr += f"{t}MPI_Comm_size(comm, &npes);\n\n"
    return outstr


def get_variant_text(namespaces: list,
                     table_name: str,
                     vname: str,
                     vinfo: toml.table,
                     indent_ph: str,
                     indent_pc: str,
                     curr_namespace_list_ph,
                     base_t: int = 3):
    """Generate text pieces for a variant parameter

    This function will create the different pieces of text
    that are required for a variant parameter. A variant
    parameter is one that is declared as an external variable
    in the header file, initialilzed in the source file, but
    required to be read in at run time.

    This returns three different output strings that define
    the behavior for each of those different files.

    Parameters
    ----------
    namespaces : list
        The list of namespaces that belong to this parameter.
    table_name : str
        The "table" name for this parameter. This is typically
        when there's a block of parameters that all use the same
        prefix. I.e. bssn::BH1_MASS and bssn::BH2_MASS can be defined
        in the template file in one table alone.
    vname : str
        The name of the variable (the key for the corresponding table)
    vinfo : tomlkit.Table or dict
        The information for the current parameter
    indent_ph : str
        The current indentation for the parameter header
    indent_pc : str
        The current indentation for the parameter source
    curr_namespace_list_ph : list of str
        The current list of namespaces captured inside the
        header file.
    base_t : int
        The number of indents that are currently used in this block.

    Returns
    -------
    paramh_str : str
        The parameters header string
    paramc_str : str
        The parameters source string
    param_read : str
        The parameter read string
    """

    # In parameters.h: declare extern
    # In parameters.cpp: declare but don't define
    # In grUtils.cpp: define from user input; throw
    # error if value input not given

    paramh_str = indent_pc
    if vinfo.get("desc", "") != "":
        paramh_str += "/** @brief: " + vinfo.get(
            "desc", "No description given") + " */\n"
        paramh_str += indent_ph
    paramh_str += "extern std::" if vinfo["dtype"][0] == 's' else "extern "
    paramh_str += vinfo["dtype"][:-2] if vinfo["dtype"][-2] == '[' else vinfo[
        "dtype"]
    paramh_str += " "
    paramh_str += table_name
    paramh_str += vname

    # this writes the size declaration for arrays
    if vinfo["dtype"][-2] == '[':
        paramh_str += "[" + str(vinfo["size"]) + "]"

    paramh_str += ';\n\n'

    # add to the c version
    paramc_str = indent_pc
    if vinfo["dtype"][0] == 's':
        paramc_str += "std::"  # strings need "std::"
    # arrays: chop "[]"
    paramc_str += vinfo["dtype"][:-2] if vinfo["dtype"][-2] == '[' else vinfo[
        "dtype"]
    paramc_str += " "
    paramc_str += table_name
    paramc_str += vname

    # this writes the size declaration for arrays
    if vinfo["dtype"][-2] == '[':
        paramc_str += "[" + str(vinfo["size"]) + "]"

    paramc_str += ";\n"

    # now we update the param_read string
    param_read = f"{TAB * base_t}if(file.contains(\""
    param_read += get_full_vname(namespaces, table_name, vname)
    param_read += f"\"))\n{TAB * base_t}{{\n{TAB * (base_t + 1)}"

    if (vinfo["dtype"][0] == 'i'
            or vinfo["dtype"][0] == 'd') and vinfo["dtype"][-2] != '[':
        param_read += get_bound_text(namespaces, table_name, vname, vinfo)

    # only for arrays: write a loop to assign values
    if vinfo["dtype"][-2] == '[':
        param_read += get_array_assignment_start(vinfo, base_t + 1)

    param_read += get_full_vname(namespaces, table_name, vname)

    # add the read definition
    param_read += "[i]" if vinfo["dtype"][-2] == '[' else ""
    param_read += " = file[\""
    param_read += get_full_vname(namespaces, table_name, vname)
    param_read += "\"][i]" if vinfo["dtype"][-2] == '[' else "\"]"
    param_read += ".as_"
    if vinfo["dtype"][0] == 'd':
        param_read += "floating"
    elif vinfo["dtype"][0] == 'i' or vinfo["dtype"][0] == 'u':
        param_read += "integer"
    elif vinfo["dtype"][0] == 's':
        param_read += "string"
    else:
        param_read += "boolean"

    if vinfo["dtype"][-2] == '[':
        param_read += f"();\n{TAB * (base_t + 1)}}}\n{TAB * base_t}}}\n\n"
    else:
        param_read += f"();\n{TAB * base_t}}}\n\n"

    param_read += f"{TAB * base_t}else\n{TAB * base_t}{{\n{TAB * (base_t + 1)}"
    param_read += "std::cerr << R\"(No value for \""
    param_read += table_name + vname
    param_read += "\"; \""
    param_read += table_name + vname
    param_read += "\" must be given a value)\" " + "<< std::endl;"
    param_read += f"\n{TAB * (base_t + 1)}exit(-1);\n{TAB * base_t}}}\n\n"

    return paramh_str, paramc_str, param_read


def get_semivariant_text(namespaces: list,
                         table_name: str,
                         vname: str,
                         vinfo: toml.table,
                         indent_ph: str,
                         indent_pc: str,
                         curr_namespace_list_ph: list,
                         base_t: int = 3):
    """Generate text pieces for a semivariant parameter

    This function will create the different pieces of text
    that are required for a semivariant parameter. A semivariant
    parameter is one that is declared as an external variable
    in the header file, initialilzed in the source file, and
    able to be modified through an input parameter file at runtime.

    This returns three different output strings that define
    the behavior for each of those different files.

    Parameters
    ----------
    namespaces : list
        The list of namespaces that belong to this parameter.
    table_name : str
        The "table" name for this parameter. This is typically
        when there's a block of parameters that all use the same
        prefix. I.e. bssn::BH1_MASS and bssn::BH2_MASS can be defined
        in the template file in one table alone.
    vname : str
        The name of the variable (the key for the corresponding table)
    vinfo : tomlkit.Table or dict
        The information for the current parameter
    indent_ph : str
        The current indentation for the parameter header
    indent_pc : str
        The current indentation for the parameter source
    curr_namespace_list_ph : list of str
        The current list of namespaces captured inside the
        header file.
    base_t : int
        The number of indents that are currently used in this block.

    Returns
    -------
    paramh_str : str
        The parameters header string
    paramc_str : str
        The parameters source string
    param_read : str
        The parameter read string
    """

    paramh_str = ""
    paramc_str = ""
    param_read = ""

    # header needs extern declaration
    # source needs declaration and define
    # grUtils.cpp define from user input
    paramh_str += indent_ph
    if vinfo.get("desc", "") != "":
        paramh_str += "/** @brief: " + vinfo.get(
            "desc", "No description given") + " */\n"
        paramh_str += indent_ph
    paramh_str += "extern std::" if vinfo["dtype"][0] == 's' else "extern "
    paramh_str += vinfo["dtype"][:-2] if vinfo["dtype"][-2] == '[' else vinfo[
        "dtype"]
    paramh_str += " " + table_name + vname

    # size declaration for arrays
    if vinfo["dtype"][-2] == '[':
        paramh_str += "[" + str(vinfo["size"]) + "]"
    paramh_str += ';\n\n'

    # now for the declarations in the cpp file
    paramc_str += indent_pc
    if vinfo["dtype"][0] == 's':
        paramc_str += "std::"
    paramc_str += vinfo["dtype"][:-2] if vinfo["dtype"][-2] == '[' else vinfo[
        "dtype"]
    paramc_str += " " + table_name + vname

    # if we have an array
    if vinfo["dtype"][-2] == '[':
        paramc_str += "[" + str(vinfo["size"]) + "]"
        paramc_str += " = {" + str(vinfo["default"])[1:-1] + "}"
    # if it isn't an array
    else:
        paramc_str += " = \"" if vinfo["dtype"][0] == 's' else " = "
        paramc_str += str(
            vinfo["default"]).lower() if vinfo["dtype"][0] == 'b' else str(
                vinfo["default"])

    paramc_str += "\";\n" if vinfo["dtype"][0] == 's' else ";\n"

    # then the parameter reading chunk
    param_read += f"{TAB*3}if(file.contains(\""
    param_read += get_full_vname(namespaces, table_name, vname)
    param_read += f"\"))\n{TAB*base_t}{{\n{TAB*(base_t+1)}"

    if (vinfo["dtype"][0] == 'i'
            or vinfo["dtype"][0] == 'd') and vinfo["dtype"][-2] != '[':
        param_read += get_bound_text(namespaces, table_name, vname, vinfo)

    # only for arrays, we need to write a loop to assign values
    if vinfo["dtype"][-2] == '[':
        param_read += get_array_assignment_start(vinfo)

    param_read += get_full_vname(curr_namespace_list_ph, table_name, vname)

    # add the read definition
    param_read += "[i]" if vinfo["dtype"][-2] == '[' else ""
    param_read += " = file[\""
    param_read += get_full_vname(namespaces, table_name, vname)
    param_read += "\"][i]" if vinfo["dtype"][-2] == '[' else "\"]"
    param_read += ".as_"
    if vinfo["dtype"][0] == 'd':
        param_read += "floating"
    elif vinfo["dtype"][0] == 'i' or vinfo["dtype"][0] == 'u':
        param_read += "integer"
    elif vinfo["dtype"][0] == 's':
        param_read += "string"
    else:
        param_read += "boolean"

    if vinfo["dtype"][-2] == "[":
        param_read += f"();\n{TAB * (base_t + 1)}}}\n{TAB * base_t}}}\n\n"
    else:
        param_read += f"();\n{TAB * base_t}}}\n\n"

    return paramh_str, paramc_str, param_read


def get_invariant_text(table_name: str, vname: str, vinfo: toml.table,
                       indent_ph: str):
    """Generate text pieces for an invariant parameter

    This function will create the different pieces of text
    that are required for an invariant parameter. An invariant
    parameter is one that is declared as a constant variable
    in the header file. It cannot be modified unless the program
    is recompiled.

    This returns the output string that define
    the behavior for the header file.

    Parameters
    ----------
    table_name : str
        The "table" name for this parameter. This is typically
        when there's a block of parameters that all use the same
        prefix. I.e. bssn::BH1_MASS and bssn::BH2_MASS can be defined
        in the template file in one table alone.
    vname : str
        The name of the variable (the key for the corresponding table)
    vinfo : tomlkit.Table or dict
        The information for the current parameter
    indent_ph : str
        The current indentation for the parameter header

    Returns
    -------
    paramh_str : str
        The parameters header string
    """

    paramh_str = indent_ph
    if vinfo.get("desc", "") != "":
        paramh_str += "/** @brief: " + vinfo.get(
            "desc", "No description given") + " */\n"
        paramh_str += indent_ph
    paramh_str += "static const "
    paramh_str += "std::" if vinfo["dtype"][0] == 's' else ""
    paramh_str += vinfo["dtype"][:-2] if vinfo["dtype"][-2] == '[' else vinfo[
        "dtype"]
    paramh_str += " "

    paramh_str += table_name + vname

    if vinfo["dtype"][-2] == "[":
        paramh_str += "[" + str(vinfo["size"]) + "[ = {"
        paramh_str += str(vinfo["default"])[1:-1]
        paramh_str += "};\n"
    else:
        paramh_str += " = \"" if vinfo["dtype"] == 's' else " = "
        paramh_str += str(vinfo["default"])
        paramh_str += "\";\n" if vinfo["dtype"] == 's' else ";\n"

    return paramh_str + "\n"


def get_hyperinvariant_text(table_name: str, vname: str, vinfo: toml.table,
                            indent_ph: str):
    """Generate text pieces for a hyperinvariant parameter

    This function will create the different pieces of text
    that are required for an invariant parameter. An invariant
    parameter is one that is declared as a constant variable
    in the header file with options from the compiler. It cannot be
    modified unless the program is recompiled.

    This returns the output string that define
    the behavior for the header file.

    Parameters
    ----------
    table_name : str
        The "table" name for this parameter. This is typically
        when there's a block of parameters that all use the same
        prefix. I.e. bssn::BH1_MASS and bssn::BH2_MASS can be defined
        in the template file in one table alone.
    vname : str
        The name of the variable (the key for the corresponding table)
    vinfo : tomlkit.Table or dict
        The information for the current parameter
    indent_ph : str
        The current indentation for the parameter header

    Returns
    -------
    paramh_str : str
        The parameters header string
    """

    paramh_str = indent_ph

    # start with the brief comment
    if vinfo.get("desc", "") != "":
        paramh_str += "/** @brief: " + vinfo.get(
            "desc", "No description given") + " */\n"

    for ii, curr_item in enumerate(vinfo["compiler_options"]):

        # initial string for the compiler option
        paramh_str += "#ifdef " + str(curr_item) + "\n"

        paramh_str += indent_ph
        paramh_str += "static const "
        paramh_str += "std::" if vinfo["dtype"][0] == 's' else ""
        paramh_str += vinfo["dtype"][:-2] if vinfo["dtype"][
            -2] == '[' else vinfo["dtype"]
        paramh_str += " "

        paramh_str += table_name + vname

        if vinfo["dtype"][-2] == "[":
            paramh_str += "[" + str(vinfo["size"]) + "[ = {"
            paramh_str += str(vinfo["default"][ii])[1:-1]
            paramh_str += "};\n"
        else:
            paramh_str += " = \"" if vinfo["dtype"] == 's' else " = "
            paramh_str += str(vinfo["default"][ii])
            paramh_str += "\";\n" if vinfo["dtype"] == 's' else ";\n"

    if len(vinfo["compiler_options"]) < len(vinfo["default"]):
        default_option = len(vinfo["compiler_options"])
        paramh_str += "#else\n"

        paramh_str += indent_ph
        paramh_str += "static const "
        paramh_str += "std::" if vinfo["dtype"][0] == 's' else ""
        paramh_str += vinfo["dtype"][:-2] if vinfo["dtype"][
            -2] == '[' else vinfo["dtype"]
        paramh_str += " "

        paramh_str += table_name + vname

        if vinfo["dtype"][-2] == "[":
            paramh_str += "[" + str(vinfo["size"]) + "[ = {"
            paramh_str += str(vinfo["default"][default_option])[1:-1]
            paramh_str += "};\n"
        else:
            paramh_str += " = \"" if vinfo["dtype"] == 's' else " = "
            paramh_str += str(vinfo["default"][default_option])
            paramh_str += "\";\n" if vinfo["dtype"] == 's' else ";\n"

    paramh_str += "#endif\n"

    return paramh_str + "\n"


def get_dependent_text(namespaces: list,
                       table_name: str,
                       vname: str,
                       vinfo: toml.table,
                       indent_ph: str,
                       curr_namespace_list_ph: list,
                       base_t: int = 3):
    """Generate text pieces for a dependent parameter

    This function will create the different pieces of text
    that are required for an invariant parameter. A dependent
    parameter is one that requires information from a previously
    read one. It is declared in the header file, and then
    updated in the broadcast text (just a simple place to add it
    as broadcasting is done after reading but before continuing).

    This returns the output string that define
    the behavior for the header file and C++ file.

    Parameters
    ----------
    table_name : str
        The "table" name for this parameter. This is typically
        when there's a block of parameters that all use the same
        prefix. I.e. bssn::BH1_MASS and bssn::BH2_MASS can be defined
        in the template file in one table alone.
    vname : str
        The name of the variable (the key for the corresponding table)
    vinfo : tomlkit.Table or dict
        The information for the current parameter
    indent_ph : str
        The current indentation for the parameter header

    Returns
    -------
    paramh_str : str
        The parameters header string
    """

    paramh_str = ""

    # header needs extern declaration
    # source needs declaration and define
    # grUtils.cpp define from user input
    paramh_str += indent_ph
    if vinfo.get("desc", "") != "":
        paramh_str += "/** @brief: " + vinfo.get(
            "desc", "No description given") + " */\n"
        paramh_str += indent_ph
    paramh_str += "extern std::" if vinfo["dtype"][0] == 's' else "extern "
    paramh_str += vinfo["dtype"][:-2] if vinfo["dtype"][-2] == '[' else vinfo[
        "dtype"]
    paramh_str += " " + table_name + vname

    # size declaration for arrays
    if vinfo["dtype"][-2] == '[':
        paramh_str += "[" + str(vinfo["size"]) + "]"
    paramh_str += ';\n\n'

    # now the "broadcast" parameter text to add
    param_bcast = f"\n{TAB*3}/* NOTE: this is not a broadcast," + \
        " but an assignment due to previous dependencies! */\n"
    param_bcast += f"{TAB*3}"

    # only for arrays, we need to write a loop to assign values
    if vinfo["dtype"][-2] == '[':
        raise NotImplementedError("Dependent arrays are not yet supported")

    param_bcast += get_full_vname(curr_namespace_list_ph, table_name, vname)

    # add the read definition
    param_bcast += " = "
    param_bcast += vinfo["default"]

    param_bcast += ";\n"

    return paramh_str, param_bcast


def get_broadcast(table_name: str,
                  vname: str,
                  vinfo: toml.table,
                  curr_namespace_list_ph: str,
                  base_t: int = 3):
    """Generate broadcasting text for MPI implementations

    This function will create the code that broadcasts the parameters
    through MPI to the other processes. This is necessary for each
    individual process (across multiple nodes) as it is the only
    way for the parameters to be consistent.

    This returns the output string that defines the behavior
    of broadcasting the parameters.

    Parameters
    ----------
    table_name : str
        The "table" name for this parameter. This is typically
        when there's a block of parameters that all use the same
        prefix. I.e. bssn::BH1_MASS and bssn::BH2_MASS can be defined
        in the template file in one table alone.
    vname : str
        The name of the variable (the key for the corresponding table)
    vinfo : tomlkit.Table or dict
        The information for the current parameter
    curr_namespace_list_ph : list of str
        The current list of namespaces captured inside the
        header file.
    base_t : int
        The number of indents that are currently used in this block.

    Returns
    -------
    bcasts : str
        The broadcasting source string
    """
    bcasts = f"{TAB * base_t}"

    # if we have an array
    if vinfo["dtype"][-1] == ']':
        bcasts += "MPI_Bcast(&("
        bcasts += get_full_vname(curr_namespace_list_ph, table_name, vname)
        bcasts += "), " + str(vinfo["size"]) + ", "
        bcasts += "MPI_DOUBLE" if vinfo["dtype"][0] == "d" else "MPI_INT"
        bcasts += ", 0, comm);"

    # if it's a string
    elif vinfo["dtype"][0] == "s":
        bcasts += "MPI_Bcast(const_cast<char*>("
        bcasts += get_full_vname(curr_namespace_list_ph, table_name, vname)
        bcasts += ".c_str()), "
        bcasts += get_full_vname(curr_namespace_list_ph, table_name, vname)
        bcasts += ".size() + 1, MPI_CHAR, 0, comm);"

    # if it's anything else, we have to use the par
    # implementation of mpi_bcast in dendro
    else:
        bcasts += "par::Mpi_Bcast(&("
        bcasts += get_full_vname(curr_namespace_list_ph, table_name, vname)
        bcasts += "), 1, 0, comm);"

    return bcasts


def get_full_vname(namespaces: list, table_name: str, vname: str):
    """Generates the string that is the full variable name

    This full variable name is what is used in the C++ code across
    all of the different namespaces and if there's a table prefix.
    Putting it all here reduces repeated code.

    Parameters
    ----------
    namespaces : list
        The list of namespaces that belong to this parameter.
    table_name : str
        The "table" name for this parameter. This is typically
        when there's a block of parameters that all use the same
        prefix. I.e. bssn::BH1_MASS and bssn::BH2_MASS can be defined
        in the template file in one table alone.
    vname : str
        The name of the variable (the key for the corresponding table)

    Returns
    -------
    str
        The output string of the full variable name
    """
    out_str = "::".join(namespaces)
    out_str += "::" + table_name + vname
    return out_str


def get_array_assignment_start(vinfo: toml.table, base_t: int = 4):
    """Generates the text for starting an array's assignment

    This simple function creates the C++ code for populating
    a parameter that's an array with the data from the input
    file. The rest of the block needs to be handled separately.

    Parameters
    ----------
    vinfo : tomlkit.Table or dict
        The information for the current parameter
    base_t : int
        The number of indents that are currently used in this block.

    Returns
    -------
    str
        The output string for array population
    """

    out_str = "for (int i = 0; i < "
    out_str += str(vinfo["size"])
    out_str += f"; ++i)\n{TAB*base_t}{{\n{TAB*(base_t+1)}"
    return out_str


def get_bound_text(namespaces, table_name, vname, vinfo, base_t=4):
    """Creates the code for checking the bounds of a parameter

    A parameter has the option of having a minimum and a maximum
    set in the template. This function creates the code that
    will then check if the input value by the incoming parameter
    file is within those bounds.

    Parameters
    ----------
    namespaces : list
        The list of namespaces that belong to this parameter.
    table_name : str
        The "table" name for this parameter. This is typically
        when there's a block of parameters that all use the same
        prefix. I.e. bssn::BH1_MASS and bssn::BH2_MASS can be defined
        in the template file in one table alone.
    vname : str
        The name of the variable (the key for the corresponding table)
    vinfo : tomlkit.Table or dict
        The information for the current parameter
    base_t : int
        The number of indents that are currently used in this block.

    Returns
    -------
    str
        The output code that checks the input against the bounds
    """
    out_str = "if (" + str(vinfo["min"]) + " > file[\""
    out_str += get_full_vname(namespaces, table_name, vname)
    out_str += "\"].as_"

    out_str += "floating" if vinfo["dtype"][0] == 'd' else "integer"

    out_str += "() || " + str(vinfo["max"]) + " < file[\""

    out_str += get_full_vname(namespaces, table_name, vname)

    out_str += "\"].as_"
    out_str += "floating" if vinfo["dtype"][0] == 'd' else "integer"

    # now for the text that declares it's an invalid value
    out_str += f"())\n{TAB * base_t}{{\n"
    out_str += f"{TAB * (base_t + 1)}std::cerr << R\"(Invalid value for \""

    out_str += get_full_vname(namespaces, table_name, vname)

    out_str += f"\")\" << std::endl;\n{TAB * (base_t + 1)}"
    out_str += f"exit(-1);\n{TAB * base_t}}}\n\n"

    return out_str + f"{TAB * base_t}"


def generate_all_parameter_text(project_short: str, filename: str):
    """Generate all of the C++ code for parameters

    This function takes a project "short", which is just a short
    name like "bssn" or "ccz4". It also takes in a template file
    in TOML format that defines all of the behavior of *all* parameters
    corresponding to the project.

    Please note that the output strings of source code *do not* include
    includes or preprocessor definitions. Those need to be added by
    hand! This was done so that the parameter files could have different
    names than "parameters.cpp" or "parameters.h" if necessary.

    Parameters
    ----------
    project_short : str
        A short string for prefixing various things. "bssn" or
        "ccz4" for example.
    filename : str
        The input filename that has the template for all parameters.

    Returns
    -------
    paramh_str : str
        The header file string. This should then be placed in the header
        file for parameters.
    paramc_str : str
        The source file string. This should then be placed in the source
        file for parameters.
    """
    namespace_tag = "-NMSPC"

    ps = project_short.lower()
    ps_u = project_short.upper()

    setup_dict = get_toml_data(filename)
    # tab delimiter
    t = TAB

    # start with the namespace block
    paramh_str = f"namespace {ps}\n{{\n"
    paramh_str += f"{t}void readParamFile(const char* inFile, " + \
        "MPI_Comm comm);\n"
    paramh_str += f"{t}void dumpParamFile(std::ostream& sout, int" + \
        " root, MPI_Comm comm);\n\n"
    paramh_str += f"{t}extern mem::memory_pool<double> {ps_u}_MEM_POOL;\n"
    paramh_str += f"{t}extern BH BH1;\n{t}extern BH BH2;\n"
    paramh_str += f"{t}extern Point {ps_u}_BH_LOC[2];\n"
    paramh_str += f"{t}//extern RefinementMode {ps_u}_REFINEMENT_MODE;\n}}\n"

    paramc_str = f"namespace {ps}\n{{\n"
    paramc_str += f"{t}mem::memory_pool<double> {ps_u}_MEM_POOL = " + \
        "mem::memory_pool<double>(0,16);\n"
    paramc_str += f"{t}BH BH1;\n{t}BH BH2;\n"
    paramc_str += f"{t}//RefinementMode {ps_u}_REFINEMENT_MODE = " + \
        "RefinementMode::WAMR;\n"
    paramc_str += f"{t}Point {ps_u}_BH_LOC[2];\n}}\n"

    # now we define the param read definition
    param_read = f"namespace {ps}\n{{\n"
    param_read += f"{t}void readParamFile(const char* inFile, " + \
        f"MPI_Comm comm)\n{t}{{\n\n"
    param_read += get_rank_npes()
    param_read += f"{t*2}auto file = toml::parse(inFile);\n\n"
    param_read += f"{t*2}if(!rank)\n{t*2}{{\n"

    # and now we start the dump param file information
    param_dump = f"{t}void dumpParamFile(std::ostream& sout, " + \
        "int root, MPI_Comm comm)\n"
    param_dump += f"{t}{{\n\n"
    param_dump += get_rank_npes()
    param_dump += f"{t*2}if(rank==root)\n{t*2}{{\n"
    param_dump += f'{t*3}sout << "Parameters read: " << std::endl;\n'

    # now we iteratively pass through the contents of setup_dict, while
    # extracting the desired information and saving it to all of the different
    # variables that represent the header file, source file, and the read and
    # dump functions.

    # what follows here are helper variables and their definitions
    stack = []  # used much as the call stack is in a recursive algorithm
    table_names = ""  # store names of the tables containing params
    namespace_list = []  # stores the list of namespaces the parameter is in
    curr_namespace_list_ph = []  # namespace for params_h
    curr_namespace_list_pc = []  # namespace for params_c
    indent_ph = ""
    indent_pc = ""
    bcasts = []  # code for mpi broadcasts for the parameter
    first_ph = True  # indicates if we're in the first namespace for params_h
    first_pc = True  # indicates if we're in the first namespace for params_c
    PNIP = False  # previous non-invariant parameter in curr namespace
    CNS_ph = False  # see if we've changed namespaces in params_h
    CNS_pc = False  # see if we've chanced namespaces in params_c

    # start by putting the setup dict on the stack
    # with the empty names and namespace list
    stack.append((setup_dict, table_names, namespace_list))

    # as long as we aren't empty, we iterate
    while stack:
        # pop the top item off of the stack
        (dikt, table_names, namespace_list) = stack.pop()

        # NOTE: both "tables" and namespaces correspond to toml tables, it's
        # just that some tables correspond to namespaces in parameters.h
        # and/or parameters.cpp while others are just used for organizational
        # purposes in the toml file. To indicate that a table designates a
        # namespace, the tag "-NMSP" is appended to its name in the toml file.
        # The names of purly organizational tables are prepended to that of
        # the parameter (when it is declared/defined in parameters.h
        # and parameters.cpp). So, if there is a parameter whose name
        # (in the toml file) is "param1" that is in "table2", which is inside
        # "table1" (both organizational tables) then in parameters.h and
        # parameters.cpp it is given the name table1_table2_param1; if
        # it is furthermore inside the table "nam2-NMSPC" which is itself
        # inside the table "nam1-NMSPC", it would be defined/declared
        # inside the namespace nam2, which would be defined within
        # the namespace nam1.

        # start by checking what kind of stuff is inside dikt (tables and
        # namespaces get put on the stack; parameters have their info
        # written to parameters.h, parameters.cpp param_read, and param_dump)
        for k, v in dikt.items():
            # v is a dict, k is a table/namespace
            if isinstance(v, dict):
                for kk, vv in v.items():
                    # if vv is a dict, then k has nested tables/namespaces
                    if isinstance(vv, dict):
                        # if k has the "-NMSPC" tab, it is a namespace
                        if k[-(len(namespace_tag)):] == namespace_tag:
                            namespace_list.append(k.split(namespace_tag)[0])
                            stack.append((v, table_names, namespace_list[:]))
                            namespace_list.pop()

                        # if we don't have the namespace tag, then
                        # k is a normal table
                        else:
                            stack.append((v, "".join([table_names, k,
                                                      "_"]), namespace_list))
                        # break for the iterative testing
                        break

                    # if k isn't nested tables or namespaces, we have a param
                    else:
                        # for parameters first add the code for output
                        # parameter name
                        tmp_param_dump = ""
                        tmp_param_dump += f'{t*3}sout << "\\t'
                        tmp_param_dump += get_full_vname(
                            namespace_list, table_names, k)

                        if v['dtype'][-2] == "[":  # arrays
                            tmp_param_dump += ": [\";"
                            tmp_param_dump += f"\n{t*3}for (unsigned" + \
                                " int i = 0; i < "
                            tmp_param_dump += str(v["size"])
                            tmp_param_dump += f"; ++i)\n{t*3}{{\n{t*4}sout << "
                            tmp_param_dump += get_full_vname(
                                namespace_list, table_names, k)
                            tmp_param_dump += "[i] << (i<"
                            tmp_param_dump += str(v["size"])
                            tmp_param_dump += "-1"
                            tmp_param_dump += "?',':']');\n"
                            tmp_param_dump += f"{t*3}}}\n{t*3}sout"
                        else:  # non-arrays
                            tmp_param_dump += ": \" << "
                            tmp_param_dump += get_full_vname(
                                namespace_list, table_names, k)
                        tmp_param_dump += " << std::endl;\n"

                        param_dump += tmp_param_dump

                        # == NESTING OF NAMESPACES ==
                        for i in range(
                                len(curr_namespace_list_ph
                                    ) if len(curr_namespace_list_ph) <
                                len(namespace_list) else len(namespace_list)):
                            # if there are any differences...
                            if curr_namespace_list_ph[i] != namespace_list[i]:
                                # ...we've changed namespaces
                                CNS_ph = True
                                break

                            # if we don't trip the above if-statement, there
                            # was no namespace change
                            CNS_ph = False

                        # now we check if we need to add namespace code to
                        # parameters.h by checking CNS_ph and the lengths
                        # of the old and new namespace lists since either
                        # condition could mean we've changed namespaces;
                        # both need to be checked because the loop above
                        # only went as fast as indexing of the shortest
                        # of curr_namespace_list and
                        # namespace_list (so there may be missed differences)
                        if CNS_ph or len(curr_namespace_list_ph) < len(
                                namespace_list):
                            # if this is the first namespace, we will write to
                            # parameters.h
                            if first_ph:
                                # update current namespace list
                                curr_namespace_list_ph = namespace_list
                                # update the indents
                                indent_ph = t * (len(curr_namespace_list_ph) -
                                                 1)

                                # and now add the namespace opener to
                                # parameters_h
                                paramh_str += "\n"
                                paramh_str += indent_ph
                                paramh_str += "namespace "
                                paramh_str += f"{namespace_list[-1]}\n"
                                paramh_str += indent_ph + "{\n"

                                # now it's in the parameters
                                first_ph = False

                            else:
                                # if it isn't the first namespace we've
                                # written to, then if
                                # len(curr_namespace_list.ph) >=
                                # len(namespace_list),
                                # we need to close the previous namespace
                                # (we already know we've changed)
                                # namespaces, so if the new namespace
                                # is the same length or shorter,
                                # we should know there's some closing to do
                                if len(curr_namespace_list_ph) >= len(
                                        namespace_list):
                                    # closing the previous namespace simply
                                    # requires writing out the appropriate
                                    # series of semicolons, newlines, and tabs:
                                    for j in range(
                                            len(curr_namespace_list_ph) - i):
                                        paramh_str += indent_ph + "}\n"
                                        indent_ph = indent_ph[:-len(t)]

                                # we then update the current namespace list...
                                curr_namespace_list_ph = namespace_list
                                # ... and update the indents
                                indent_ph = t * (len(curr_namespace_list_ph) -
                                                 1)

                                # we then write the new namespace opener to
                                # parameters.h
                                paramh_str += "\n" + indent_ph
                                paramh_str += "namespace "
                                paramh_str += f"{namespace_list[-1]}\n"
                                paramh_str += f"{indent_ph}{{\n"

                        # if we changed namespaces in params.h then we can
                        # check if we've changed namespaces in parameters.cpp
                        # using mostly the same method as above:
                        # by looping through the current and new namespace
                        # lists and running comparisons
                        for i in range(
                                len(curr_namespace_list_pc
                                    ) if len(curr_namespace_list_pc) <
                                len(namespace_list) else len(namespace_list)):
                            # if there are any differences...
                            if curr_namespace_list_pc[i] != namespace_list[i]:
                                # ...we've changed namespaces
                                CNS_pc = True
                                break

                            # if we don't trip the above if-statement, there
                            # was no namespace change
                            CNS_pc = False

                        # here we check if we need to add namespace code to
                        # parameters.cpp; we might not need to because
                        # we only need to write the opener for a namespace if
                        # we haven't previously written a parameter
                        # out to that namespace in parameters.cpp; this is
                        # technically true for parameters.h as well but for
                        # parameters.cpp the matter is complicated by the fact
                        # that not every parameter is written out to
                        # parameters.cpp (only the non-invariant ones are);
                        # this requires different logic than that used
                        # for parameters.h
                        if v["class"][0] != 'i' and (
                                CNS_pc or len(curr_namespace_list_pc) <
                                len(namespace_list)):
                            # if this is the first namespace we're using
                            if first_pc:
                                # update the current list
                                curr_namespace_list_pc = namespace_list

                                # update the indents
                                indent_pc = t * (len(curr_namespace_list_pc) -
                                                 1)

                                # and now add the namespace opener to params_c
                                paramc_str += "\n"
                                paramc_str += indent_pc
                                paramc_str += "namespace "
                                paramc_str += f"{namespace_list[-1]}\n"
                                paramc_str += indent_pc + "{\n"

                                # with a namespace now included, we update
                                first_pc = False

                                # "previous non-invariant parameter" in this NS
                                PNIP = True

                            # if this is not the first namespace
                            else:
                                # ...we check if we need to close the previous
                                # namespace; we only need to do this if
                                # there was a Previous Non-Invariant Parameter
                                # (PNIP) and the namespace change we
                                # we already know there was didn't involve
                                # a level deeper in the nesting structure
                                if len(curr_namespace_list_pc) >= len(
                                        namespace_list):
                                    # here we follow the same procedure as .h
                                    for n in range(
                                            len(curr_namespace_list_pc) - i):
                                        paramc_str += indent_pc + "}\n"
                                        indent_pc = indent_pc[:-len(t)]

                                # then update the current namespace lise
                                curr_namespace_list_pc = namespace_list

                                # and then update the indents
                                indent_pc = t * (len(curr_namespace_list_pc) -
                                                 1)

                                # we then write the new namespace opener .c
                                paramc_str += "\n" + indent_pc
                                paramc_str += "namespace "
                                paramc_str += f"{namespace_list[-1]}\n"
                                paramc_str += f"{indent_ph}{{\n"

                                # and then we know that we've seen a "previous
                                # non-invariant parameter"
                                PNIP = True

                        # now we add an extra indent needed for the declaration
                        # /definition code
                        indent_pc += t
                        indent_ph += t

                        # in this next section, we write the declarations and
                        # definitions required for the current parameter we're
                        # looking at; the way this is done varies based on the
                        # "class" of the parameter, of which there are 3:
                        # variant, semivariant, and invariant

                        if v["class"][0] == 'v':  # If parameter is variant
                            # In parameters.h: declare extern
                            # In parameters.cpp: declare but don't define
                            # In grUtils.cpp: define from user input; throw
                            # error if value input not given

                            temp_ph, temp_pc, temp_read = get_variant_text(
                                namespace_list, table_names, k, v, indent_ph,
                                indent_pc, curr_namespace_list_ph, 3)

                            paramh_str += temp_ph
                            paramc_str += temp_pc
                            param_read += temp_read

                        # if the parameter is semivariant
                        elif v["class"][0] == 's':

                            temp_ph, temp_pc, temp_read = get_semivariant_text(
                                namespace_list, table_names, k, v, indent_ph,
                                indent_pc, curr_namespace_list_ph, 3)

                            paramh_str += temp_ph
                            paramc_str += temp_pc
                            param_read += temp_read

                        # if the parameter is invariant
                        elif v["class"][0] == "i":
                            temp_ph = get_invariant_text(
                                table_names, k, v, indent_ph)

                            paramh_str += temp_ph

                        # if a parameter is hyperinvariant
                        elif v["class"][0] == "h":
                            temp_ph = get_hyperinvariant_text(
                                table_names, k, v, indent_ph)

                            paramh_str += temp_ph

                        # if a parameter is dependent
                        elif v["class"][0] == "d":
                            temp_ph, temp_bcast = get_dependent_text(
                                namespace_list, table_names, k, v, indent_ph,
                                curr_namespace_list_ph, 3)

                            paramh_str += temp_ph
                            # add a "bcast" string because it's easy to slot here
                            bcasts.insert(0, temp_bcast)

                        # now, for all but invariant, we write broadcast code
                        if v["class"][0] != "i":
                            bcasts.append(
                                get_broadcast(table_names, k, v,
                                              curr_namespace_list_ph, 2))

                        # now update indent_ph and indent_pc
                        indent_ph = indent_ph[:-len(TAB)]
                        indent_pc = indent_pc[:-len(TAB)]

                        break

    # now we close up all namespaces
    for ii in range(len(curr_namespace_list_ph)):
        paramh_str += indent_ph
        paramh_str += "}\n"
        indent_ph = indent_ph[:-len(TAB)]

    # and close up all namespaces in the cpp file
    for ii in range(len(curr_namespace_list_pc)):
        paramc_str += indent_pc
        paramc_str += "}\n"
        indent_pc = indent_pc[:-len(TAB)]

    # blackhole 1
    param_read += f"{TAB * 3}{ps}::BH1 = BH({ps_u}_BH1_MASS,\n"
    for ii, bhvar in enumerate(BHOLE_VARS[1:]):
        param_read += f"  {TAB * 6}{ps_u}_BH1_{bhvar.upper()}"
        if ii < len(BHOLE_VARS) - 2:
            param_read += ",\n"
        else:
            param_read += ");\n\n"

    # blackhole 2
    param_read += f"{TAB * 3}{ps}::BH2 = BH({ps_u}_BH2_MASS,\n"
    for ii, bhvar in enumerate(BHOLE_VARS[1:]):
        param_read += f"  {TAB * 6}{ps_u}_BH2_{bhvar.upper()}"
        if ii < len(BHOLE_VARS) - 2:
            param_read += ",\n"
        else:
            param_read += ");\n\n"
    
    # then we close out of rank0-only code that reads things in
    param_read += f"{TAB*2}}}\n\n"

    # combine the brodcasts to the param_read
    param_read += f"{TAB*2}// Broadcast code to send parameters to other processes\n"
    param_read += "\n".join(bcasts)
    param_read += f"\n{TAB*1}}}\n"

    # then add param_read to the c file
    paramc_str += param_read
    paramc_str += "\n\n"
    paramc_str += param_dump

    # then close it out with braces
    paramc_str += f"{TAB*2}}}\n{TAB}}}\n}}\n"

    return paramh_str, paramc_str


def recurse_param_dict_for_toml(in_dict: toml.table,
                                namespaces: list = [],
                                table_name: str = ""):
    """Recurses an input dictionary for parameter information

    This function takes in the input dictionary and recursively
    iterates through it to find all of the parameters contained
    within the dictionary. Then, it can generate the lines of
    text for the sample parameter TOML file based on that information.

    Parameters
    ----------
    in_dict : dict or tomlkit.table
        The input dictionary based on the template file. The structure
        of this dictionary depends entirely on what stage the recursion
        is on and on the template file.
    namespaces : list
        The list of encountered namespaces.
    table_name : str
        The name of the current table (often nothing)

    Returns
    -------
    str
        The currently built-up string with all of the TOML code
        for the parameters.
    """
    out_str = ""

    for the_key, the_var in in_dict.items():

        if isinstance(the_var, dict):
            # if we've got "class" in the variable, it's an item
            if "class" in the_var.keys():
                if the_var["class"] not in [
                        "invariant", "dependent", "hypervariant"
                ]:
                    out_str += get_toml_lines(the_key, the_var, namespaces,
                                              table_name)

            else:
                namespace_name = the_key.split("-NMSPC")[0]
                if namespace_name != the_key:
                    namespaces.append(namespace_name)
                else:
                    table_name = namespace_name + "_"

                out_str += "\n# ========\n# === "
                out_str += "::".join(namespaces)
                if table_name:
                    out_str += "::" + table_name[:-1]
                out_str += " PARAMETERS\n"
                out_str += "# ========\n"
                out_str += recurse_param_dict_for_toml(the_var, namespaces,
                                                       table_name)

            # if there was no table name, we pop off the last namespace
            if the_key.endswith("-NMSPC"):
                try:
                    namespaces.pop()
                except IndexError:
                    pass

    return out_str


def get_toml_lines(var_name: str, in_dict: dict, namespaces: list,
                   table_name: str):
    """Generates the lines of TOML for a given parameter

    With the input information, it generates the lines for the
    sample TOML file. This includes a brief description
    of the parameter (based on a "desc" field in the template file)
    as well as setting the value equal to the default.

    Parameters
    ----------
    var_name : str
        The name of the variable (not the full name)
    in_dict : dict or tomlkit.table
        The set of information corresponding to the parameter
    namespaces : list
        The list of namespaces that correspond to this variable
    table_name : str
        The corresponding table name

    Returns
    -------
    str
        The generated TOML code for the parameter
    """
    # basic parameter information
    out_str = "# @brief: "
    out_str += in_dict.get("desc", "No description given.")
    out_str += f"\n# param type: {in_dict['class']} "
    out_str += f"| data type: {in_dict['dtype']} "

    if "default" in in_dict.keys():
        out_str += f"| default: {in_dict['default']} "
    if "min" in in_dict.keys():
        out_str += f"| min: {in_dict['min']} "
    if "max" in in_dict.keys():
        out_str += f"| max: {in_dict['max']}"

    full_v_name = get_full_vname(namespaces, table_name, var_name)

    out_str += "\n\"" + full_v_name + "\" = "
    if in_dict['dtype'] == "string":
        out_str += "\"" + str(in_dict["default"]) + "\""
    elif in_dict['dtype'] == "bool":
        out_str += str(in_dict["default"]).lower()
    else:
        if isinstance(in_dict["default"], str):
            print(f"\nWARNING: The parameter {full_v_name} was identified",
                  "as a number but the default value is not.",
                  "The C++ code generation will take that value,",
                  "as we accept expressions.",
                  "However, THIS WILL NEED TO BE UPDATED",
                  "BY HAND BEFORE RUNNING.",
                  "Either remove this parameter or fix the",
                  "values in the generated file.\n")
        out_str += str(in_dict["default"])

    out_str += "\n\n"

    return out_str


def generate_sample_config_file_text(project_short: str, filename: str):
    """Generate TOML file text from the template Parameters

    This function will take the template parameter file
    just as the `generate_all_parameter_text` function does,
    but it will generate the text needed to create a sample
    configuration file that specifies parameters. This file
    can then be used in conjunction with the compiled project
    for setting the various parameters.

    It will create entries for every single parameter in the
    template file, except for the invariant parameters. Invariant
    parameters can only be changed in the source code and the project
    needs to be recompiled for the changes to persist.

    The output text will also assign all of the default values
    to the different variables. This means that if C++
    expressions are used for the defaults, they will need to be
    changed by hand before running. If such a problem occurs,
    a warning will be printed in the console with the variable
    name that needs to be adjusted.

    Parameters
    ----------
    project_short : str
        A short string for prefixing various things. "bssn" or
        "ccz4" for example.
    filename : str
        The input filename that has the template for all parameters.

    Returns
    -------
    str
        The string that can be written or inserted into the
        sample parameters TOML file.
    """

    ps = project_short.lower()
    ps_u = project_short.upper()

    setup_dict = get_toml_data(filename)

    # start with the output string
    out_str = ""
    out_str += "# ===========================================\n"
    out_str += "# ==== PARAMETER FILE FOR PROJECT: " + ps_u + "\n"
    out_str += "# ===========================================\n"
    out_str += "#\n# This file contains all of the configurable parameters" + \
        " for running this Dendro project.\n"
    out_str += "#\n# This file was initially generated from the file: " + \
        os.path.basename(filename) + "\n"
    out_str += "#\n# Feel free to edit any of the values listed in this " + \
        " file or make copies for individual runs.\n"
    out_str += "#\n# Each of the parameters listed in this file should" + \
        " contain information about the parameter according to the template.\n"
    out_str += "#\n# Please note that all \"invariant\", \"dependent\", and " + \
        "\"hyperinvariant\" parameters were not included.\n"
    out_str += "#\n# NOTE: What follows is an explanation on parameter types:"
    out_str += "\n# A variant parameter is required for execution."
    out_str += "\n# A semivariant parameter is not required for execution" + \
        " and can be removed from this file safely."
    out_str += "\n\n"

    # time to iterate through the dictionary
    out_str += recurse_param_dict_for_toml(setup_dict)

    out_str += "# ===========================================\n"
    out_str += "# END PARAMTER FILE GENERATION\n"
    out_str += "# ===========================================\n"

    return out_str
