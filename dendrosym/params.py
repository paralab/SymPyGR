'''
@brief parameter file generation should go here.
'''

import tomlkit as toml

TAB = "    "

BHOLE_VARS = [
    "MASS", "X", "Y", "Z", "V_X", "V_Y", "V_Z", "SPIN", "SPIN_THETA",
    "SPIN_PHI"
]


def get_toml_data(filename):
    """Function that can read a TOML file and create a dictionary

    """
    with open(filename, "r") as in_file:
        toml_data = toml.parse(in_file.read())

    return toml_data


# NOTE: i'm ignoring all of the "includes, jumping straight to the params"


def get_rank_npes(n_tabs=2, tab_char="    "):
    t = tab_char * n_tabs
    outstr = f"{t}int rank, npes;\n"
    outstr += f"{t}MPI_Comm_rank(comm, &rank);\n"
    outstr += f"{t}MPI_Comm_size(comm, &npes);\n\n"
    return outstr


def params_file_data(project_short: str, filename: str):

    namespace_tag = "-NMSPC"

    ps = project_short.lower()
    ps_u = project_short.upper()

    setup_dict = get_toml_data(filename)
    # tab delimiter
    t = "    "

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
    param_dump = f"{t}void dumpParamFile(std::osstream& sout, " + \
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
                            tmp_param_dump += f"\n{t*3}for (unsigned " + \
                                " int i = 0; i < "
                            tmp_param_dump += str(v["size"])
                            tmp_param_dump += f"; ++i)\n{t*3}{{\n{t*4}sout << "
                            tmp_param_dump += get_full_vname(
                                namespace_list, table_names, k)
                            tmp_param_dump += "[i] << (i<"
                            tmp_param_dump += str(v["size"])
                            tmp_param_dump += "-1"
                            tmp_param_dump += "?',',:']');\n"
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

                            temp_ph, temp_pc, temp_read = get_variant(
                                namespace_list, table_names, k, v, indent_ph,
                                indent_pc, curr_namespace_list_ph, 3)

                            paramh_str += temp_ph
                            paramc_str += temp_pc
                            param_read += temp_read

                        # if the parameter is semivariant
                        elif v["class"][0] == 's':

                            temp_ph, temp_pc, temp_read = get_semiinvariant(
                                namespace_list, table_names, k, v, indent_ph,
                                indent_pc, curr_namespace_list_ph, 3)

                            paramh_str += temp_ph
                            paramc_str += temp_pc
                            param_read += temp_read

                        # if the parameter is invariant
                        elif v["class"][0] == "i":
                            temp_ph = get_invariant(table_names, k, v,
                                                    indent_ph)
                            
                            paramh_str += temp_ph

                        # now, for all but invariant, we write broadcast code
                        if v["class"][0] != "i":
                            bcasts.append(
                                get_broadcast(table_names, k, v,
                                              curr_namespace_list_ph, 3))

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
    param_read += f"{TAB * 3}{ps}:BH1 = BH({ps_u}_BH1_MASS,\n"
    for bhvar in BHOLE_VARS[1:]:
        param_read += f"{TAB * 6}{ps_u}_BH1_{bhvar.upper()}\n"

    # blackhole 2
    param_read += f"{TAB * 3}{ps}:BH2 = BH({ps_u}_BH2_MASS,\n"
    for bhvar in BHOLE_VARS[1:]:
        param_read += f"{TAB * 6}{ps_u}_BH2_{bhvar.upper()}\n"

    # combine the brodcasts to the param_read
    param_read += "\n".join(bcasts)
    param_read += f"\n{TAB*2}}}\n{TAB}}}"

    # then add param_read to the c file
    paramc_str += param_read
    paramc_str += "\n"
    paramc_str += param_dump

    # then close it out with braces
    paramc_str += f"{TAB*2}}}\n{TAB}}}\n}}\n"

    return paramh_str, paramc_str


def get_variant(namespaces,
                table_name,
                vname,
                vinfo,
                indent_ph,
                indent_pc,
                curr_namespace_list_ph,
                base_t=3):

    # In parameters.h: declare extern
    # In parameters.cpp: declare but don't define
    # In grUtils.cpp: define from user input; throw
    # error if value input not given

    paramh_str = indent_pc
    if vinfo.get("desc", "") != "":
        paramh_str += "/** @brief: " + vinfo.get("desc", "No description given") + " */\n"
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


def get_semiinvariant(namespaces,
                      table_name,
                      vname,
                      vinfo,
                      indent_ph,
                      indent_pc,
                      curr_namespace_list_ph,
                      base_t=3):

    paramh_str = ""
    paramc_str = ""
    param_read = ""

    # header needs extern declaration
    # source needs declaration and define
    # grUtils.cpp define from user input
    paramh_str += indent_ph
    if vinfo.get("desc", "") != "":
        paramh_str += "/** @brief: " + vinfo.get("desc", "No description given") + " */\n"
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


def get_invariant(table_name, vname, vinfo, indent_ph):

    paramh_str = indent_ph
    if vinfo.get("desc", "") != "":
        paramh_str += "/** @brief: " + vinfo.get("desc", "No description given") + " */\n"
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


def get_broadcast(table_name, vname, vinfo, curr_namespace_list_ph, base_t=3):

    bcasts = f"{TAB * base_t}"

    # if we have an array
    if vinfo["dtype"][-1] == ']':
        bcasts += "MPI_Bcast(&("
        bcasts += get_full_vname(curr_namespace_list_ph, table_name, vname)
        bcasts += "), " + str(vinfo["size"]) + ", "
        bcasts += "MPI_DOUBLE" if vinfo["dtype"][0] == "d" else "MPI_INT"
        bcasts += ", 0, comm);"

    elif vinfo["dtype"] == "s":
        bcasts += "MPI_Bcast(const_cast<char*>("
        bcasts += get_full_vname(curr_namespace_list_ph, table_name, vname)
        bcasts += ".c_str()), "
        bcasts += get_full_vname(curr_namespace_list_ph, table_name, vname)
        bcasts += ".size() + 1, MPI_CHAR, 0, comm);"

    else:
        bcasts += "par::Mpi_Bcast(&"
        bcasts += get_full_vname(curr_namespace_list_ph, table_name, vname)
        bcasts += ", 1, 0, comm);"

    return bcasts


def get_full_vname(namespaces, table_name, vname):
    out_str = "::".join(namespaces)
    out_str += "::" + table_name + vname
    return out_str


def get_array_assignment_start(vinfo, base_t=4):

    out_str = "for (int i = 0; i < "
    out_str += str(vinfo["size"])
    out_str += f"; ++i)\n{TAB*base_t}{{\n{TAB*(base_t+1)}"
    return out_str


def get_bound_text(namespaces, table_name, vname, vinfo, base_t=4):
    out_str = "if (" + str(vinfo["min"]) + " > file[\""
    out_str += get_full_vname(namespaces, table_name, vname)
    out_str += "\"].as_"

    out_str += "floating" if vinfo["dtype"][0] == 'd' else "integer"

    out_str += "() || " + str(vinfo["max"]) + " < file[\""

    out_str += get_full_vname(namespaces, table_name, vname)

    out_str += "\".as_"
    out_str += "floating" if vinfo["dtype"][0] == 'd' else "integer"

    # now for the text that declares it's an invalid value
    out_str += f"())\n{TAB * base_t}{{\n"
    out_str += f"{TAB * (base_t + 1)}std::cerr << R\"(Invalid value for \""

    out_str += get_full_vname(namespaces, table_name, vname)

    out_str += f"\")\" << std::endl;\n{TAB * (base_t + 1)}"
    out_str += f"exit(-1);\n{TAB * base_t}}}\n\n"

    return out_str + f"{TAB * base_t}"
