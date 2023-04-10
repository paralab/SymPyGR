# This script generates the c++ files needed for the BSSN code to access the parameters set up
# in setup.toml and parameters.toml

import tomlkit as toml
import sys

# GET INPUT FROM SETUP FILE ############################################################################################
with open(sys.argv[1]) as in_file:
    # this is where we will store the setup file while we read from it
    setup_dict = toml.parse(in_file.read())

# CREATE parameters.h and  parameters.cpp ##############################################################################
with open("parameters.h", "w") as params_h, open("parameters.cpp",
                                                 "w") as params_cpp:
    # this will store the definition of the readParamFile function
    param_read_def = ""

    # this will store the definition of the dumpParamFile function
    param_dump_def = ""

    # First we write the headers for each of the files
    params_h.write('''#ifndef SFCSORTBENCH_PARAMETERS_H
#define SFCSORTBENCH_PARAMETERS_H

#include "bh.h"
#include "grDef.h"
#include <string.h>
#include <iostream>
#include "dendro.h"
#include "memory_pool.h"
#include "toml.hpp"

namespace bssn
{
    void readParamFile(const char* inFile, MPI_Comm comm);
    void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm);

    extern mem::memory_pool<double> BSSN_MEM_POOL;
    //extern RefinementMode BSSN_REFINEMENT_MODE;
    extern BH BH1;
    extern BH BH2;
    extern Point BSSN_BH_LOC[2];
}

''')

    params_cpp.write('''#include "parameters.h"

namespace bssn
{
    mem::memory_pool<double> BSSN_MEM_POOL = mem::memory_pool<double>(0,16);
    //RefinementMode BSSN_REFINEMENT_MODE = RefinementMode::WAMR;
    BH BH1;
    BH BH2;
    Point BSSN_BH_LOC[2];
}
''')

    # then we write the opening lines for the definitions of readParamFile and dumpParamFile
    param_read_def = "".join([
        param_read_def, '''
namespace bssn
{
    void readParamFile(const char* inFile, MPI_Comm comm)
    {
        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        auto file = toml::parse(inFile);

        if(!rank)
        {
'''
    ])

    param_dump_def = "".join([
        param_dump_def, '''    
    void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm)
    {

        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        if(rank==root)
        {
            sout << "parameters read: " << std::endl;
'''
    ])

    # now we're going to use an iterative algorithm to traverse the contents of setup_dict, extracting the desired
    # information and writing it out to parameters.h, parameters.cpp, param_read_def, and param_dump_def (as the case
    # requires); below are helper variables that will be used in this process:

    stack = [
    ]  # this will be used much as the call stack is in a recursive algorithm
    table_names = ""  # this will store the names of the tables containing the parameter we are looking at
    namespace_list = [
    ]  # this will store the list of namespaces the parameter we are looking at is in
    curr_namespace_list_ph = [
        ""
    ]  # this will store the namespace from which we are writing to parameters.h
    curr_namespace_list_pc = [
        ""
    ]  # this will store the namespace from which we are writing to parameters.cpp
    indent_ph = ""  # this will store the characters to appropriately indent the text that we write to parameters.h
    indent_pc = ""  # this will store the characters to appropriately indent the text that we write to parameters.cpp
    bcasts = [
    ]  # this will store the info for the code for the mpi broadcasts for the parameter we're looking at
    first_ph = True  # indicates if we are dealing with the first namespace to be written to parameters.h
    first_pc = True  # indicates if we are dealing with the first namespace to be written to parameters.cpp
    PNIP = False  # tells if we've had a previous non-invariant parameter in current namespace (True = yes, False = no)
    CNS_ph = False  # indicates if we've changed namespaces in parameters.h
    CNS_pc = False  # indicates if we've changed namespaces in parameters.cpp

    # we initialize by putting setup_dict on the stack along with the empty table_names and namespace_list
    stack.append((setup_dict, table_names, namespace_list))

    # iterate so long as stack isn't empty
    while stack:
        # pop the top item off the stack
        (dikt, table_names, namespace_list) = stack.pop()

        # Note: both "tables" and "namespaces" correspond (in the toml file) to toml tables, it's just that some tables
        # correspond to namespaces in parameters.h and/or parameters.cpp and some are just used for organizational
        # purposes in the toml file; to indicate that a table designates a namespace, the tag "-NMSPC" is appended to
        # its name (in the toml file). The names of purely organizational tables are prepended to that of the parameter
        # (when it is declared/defined in parameters.h and parameters.cpp). So, if there is a parameter whose name (in
        # the toml file) is "param1" that is inside "table2", which is inside "table1 " (both organizational tables)
        # then in parameters.h and parameters.cpp it is given the name "table1_table2_param1"; if it is furthermore
        # inside the table "nam2-NMSPC", which is itself inside the table "nam1-NMSPC", it would be defined/declared
        # inside the namespace nam2, which would be defined within the namespace nam1.

        # check what kind of stuff is inside dikt (tables and namespaces get put on the stack; parameters have their
        # info written to parameters.h, parameters.cpp, param_read_def, and param_dump_def).
        for k, v in dikt.items():
            # v is a dict, then k is a table/namespace
            if isinstance(v, dict):
                for kk, vv in v.items():
                    # if vv is a dict, then k has nested tables/namespaces
                    if isinstance(vv, dict):
                        # if k has the "-NMSPC" tag, it is a namespace
                        if k[-6:] == "-NMSPC":
                            namespace_list.append(k[:-6])
                            stack.append((v, table_names, namespace_list[:]))
                            namespace_list.pop()

                        # If k doesn't have the "-NMPSC" tag, k is a normal table
                        else:
                            stack.append((v, "".join([table_names, k,
                                                      "_"]), namespace_list))

                        break

                    # if k doesn't have nested tables or namespaces, it's a parameter
                    else:
                        # For parameters, first write code to output parameter name and value to console
                        param_dump_def = "".join(
                            [param_dump_def, "\t\t\tsout << \"\\t"])
                        param_dump_def = "".join(
                            [param_dump_def, "::".join(namespace_list)])
                        param_dump_def = "".join([param_dump_def, "::"])
                        param_dump_def = "".join([param_dump_def, table_names])
                        param_dump_def = "".join([param_dump_def, k])
                        if v['type'][-2] == '[':  # for arrays
                            param_dump_def = "".join(
                                [param_dump_def, ": [\";"])
                            param_dump_def = "".join([
                                param_dump_def,
                                "\n\t\t\tfor (unsigned int i = 0; i < "
                            ])
                            param_dump_def = "".join(
                                [param_dump_def,
                                 str(v["size"])])
                            param_dump_def = "".join([
                                param_dump_def,
                                "; ++i)\n\t\t\t{\n\t\t\t\tsout << "
                            ])
                            param_dump_def = "".join(
                                [param_dump_def, "::".join(namespace_list)])
                            param_dump_def = "".join([param_dump_def, "::"])
                            param_dump_def = "".join(
                                [param_dump_def, table_names])
                            param_dump_def = "".join([param_dump_def, k])
                            param_dump_def = "".join(
                                [param_dump_def, "[i] << (i<"])
                            param_dump_def = "".join(
                                [param_dump_def,
                                 str(v["size"])])
                            param_dump_def = "".join([param_dump_def, "-1"])
                            param_dump_def = "".join([
                                param_dump_def,
                                "?',':']');\n\t\t\t}\n\t\t\tsout "
                            ])
                        else:  # for non-arrays
                            param_dump_def = "".join(
                                [param_dump_def, ": \" << "])
                            param_dump_def = "".join(
                                [param_dump_def, "::".join(namespace_list)])
                            param_dump_def = "".join([param_dump_def, "::"])
                            param_dump_def = "".join(
                                [param_dump_def, table_names])
                            param_dump_def = "".join([param_dump_def, k])
                        param_dump_def = "".join(
                            [param_dump_def, " << std::endl;\n"])

                        # note: this next section contains code that handles the nesting of namespaces
                        # in parameters.h and parameters.cpp

                        # First we check if we've changed namespaces in parameters.h; we begin by looping through the
                        # current and new namespace lists and comparing their elements:
                        for i in range(
                                len(curr_namespace_list_ph
                                    ) if len(curr_namespace_list_ph) <
                                len(namespace_list) else len(namespace_list)):
                            # if there are any differences...
                            if curr_namespace_list_ph[i] != namespace_list[i]:
                                # ...we've changed namespaces
                                CNS_ph = True
                                break

                            # if we don't trip the above if-statement, there was no namespace change
                            CNS_ph = False

                        # here we check if we need to add namespace code to parameters.h by checking CNS_ph and the lengths
                        # of the old and new namespace lists since either condition could mean we've changed namespaces;
                        # we need to check both because the loop we did above only went as far as the indexing of the
                        # shortest of curr_namespace_list and namespace_list (so we might have missed differences)
                        if CNS_ph or len(curr_namespace_list_ph) < len(
                                namespace_list):
                            # if this is the first namespace we will write to parameters.h...
                            if first_ph:
                                # ...update current namespace list...
                                curr_namespace_list_ph = namespace_list
                                # ...update indents...
                                indent_ph = '\t' * (
                                    len(curr_namespace_list_ph) - 1)

                                # ...and write the new namespace opener to parameters.h
                                params_h.write('\n')
                                params_h.write(indent_ph)
                                params_h.write("namespace ")
                                params_h.write(namespace_list[-1])
                                params_h.write('\n')
                                params_h.write(indent_ph)
                                params_h.write("{\n")

                                # we've now written a namespace to parameters.h, so...
                                first_ph = False

                            # if this is not the first namespace we've written to parameters.h...
                            else:
                                # ...then, if len(curr_namespace_list_ph) >= len(namespace_list), we need to close the
                                # previous namespace (note: we already know we've changed namespaces, so if the new
                                # namespace is the same length or shorter, we know there's some closing to do):
                                if len(curr_namespace_list_ph) >= len(
                                        namespace_list):
                                    # closing the previous namespace simply requires writing out the appropriate
                                    # series of semicolons, newlines, and tabs:
                                    for j in range(
                                            len(curr_namespace_list_ph) - i):
                                        params_h.write(indent_ph)
                                        params_h.write('}\n')
                                        indent_ph = indent_ph[:-1]

                                # we then update the current namespace list...
                                curr_namespace_list_ph = namespace_list
                                # ... and update the indents
                                indent_ph = '\t' * (
                                    len(curr_namespace_list_ph) - 1)

                                # we then write the new namespace opener to parameters.h
                                params_h.write('\n')
                                params_h.write(indent_ph)
                                params_h.write("namespace ")
                                params_h.write(namespace_list[-1])
                                params_h.write('\n')
                                params_h.write(indent_ph)
                                params_h.write("{\n")

                        # having checked if we changed namespaces in parameters.h, we now check if we've changed namespaces
                        # in parameters.cpp using mostly the same method as above: by looping through the current and new
                        # namespace lists and comparing their elements:
                        for i in range(
                                len(curr_namespace_list_pc
                                    ) if len(curr_namespace_list_pc) <
                                len(namespace_list) else len(namespace_list)):
                            # if there are any differences...
                            if curr_namespace_list_pc[i] != namespace_list[i]:
                                # ...we've changed namespaces
                                CNS_pc = True
                                break

                            # if we don't trip the above if-statement, there was no namespace change
                            CNS_pc = False

                        # here we check if we need to add namespace code to parameters.cpp; we might not need to because
                        # we only need to write the opener for a namespace if we haven't previously written a parameter
                        # out to that namespace in parameters.cpp; this is technically true for parameters.h as well but for
                        # parameters.cpp the matter is complicated by the fact that not every parameter is written out to
                        # parameters.cpp (only the non-invariant ones are); this requires different logic than that used
                        # for parameters.h
                        if v["class"][0] != 'i' and (
                                CNS_pc or len(curr_namespace_list_pc) <
                                len(namespace_list)):

                            # if this is the first namespace we're writing to parameters.cpp...
                            if first_pc:
                                # ... just update current namespace list...
                                curr_namespace_list_pc = namespace_list

                                # ...update indents...
                                indent_pc = '\t' * (
                                    len(curr_namespace_list_pc) - 1)

                                # ... and write namespace opener to parameters.cpp
                                params_cpp.write('\n')
                                params_cpp.write(indent_pc)
                                params_cpp.write("namespace ")
                                params_cpp.write(namespace_list[-1])
                                params_cpp.write('\n')
                                params_cpp.write(indent_pc)
                                params_cpp.write("{\n")

                                # since we've now written a namespace to parameters.cpp:
                                first_pc = False

                                # we now know that we have seen a "Previous Non-Invariant Parameter" in this namespace
                                PNIP = True

                            # if this is not the first namespace we've written to parameters.cpp...
                            else:
                                # ...we check if we need to close the previous namespace; we only need to do this if
                                # there was a "Previous Non-Invariant Parameter" (PNIP) and the namespace change we
                                # we already know there was didn't involve going a level deeper in the nesting structure
                                if len(curr_namespace_list_pc) >= len(
                                        namespace_list):
                                    # here we follow the same procedure as with parameters.h
                                    for n in range(
                                            len(curr_namespace_list_pc) - i):
                                        params_cpp.write(indent_pc)
                                        params_cpp.write('}\n')
                                        indent_pc = indent_pc[:-1]

                                # then we update the current namespace list...
                                curr_namespace_list_pc = namespace_list

                                # ... and then update the indents...
                                indent_pc = '\t' * (
                                    len(curr_namespace_list_pc) - 1)

                                # ... and then write the new namespace opener to parameters.h
                                params_cpp.write('\n')
                                params_cpp.write(indent_pc)
                                params_cpp.write("namespace ")
                                params_cpp.write(namespace_list[-1])
                                params_cpp.write('\n')
                                params_cpp.write(indent_pc)
                                params_cpp.write("{\n")

                                # we now know that we have seen a "Previous Non-Invariant Parameter" in this namespace
                                PNIP = True

                        # here we add an extra indent needed for the parameter declaration/definition code
                        indent_ph += '\t'
                        indent_pc += '\t'

                        # in this next section, we write the declarations and definitions required for
                        # the current parameter we're looking at; the way this is done varies based on the "class" of
                        # the parameter, of which there are 3: variant, semivariant, and invariant

                        if v["class"][0] == 'v':  # If parameter is variant
                            # In parameters.h: declare extern
                            # In parameters.cpp: declare but don't define
                            # In grUtils.cpp: define from user input; throw error if value input not given
                            params_h.write(indent_ph)
                            params_h.write(
                                "extern std::" if v["type"][0] == 's' else
                                "extern ")  # strings need "std::"
                            params_h.write(
                                v["type"][:-2] if v["type"][-2] == '[' else
                                v["type"])  # arrays: chop "[]"
                            params_h.write(' ')
                            params_h.write(table_names)
                            params_h.write(k)
                            if v["type"][
                                    -2] == '[':  # this writes the size declaration needed for arrays
                                params_h.write("[")
                                params_h.write(str(v["size"]))
                                params_h.write("]")
                            params_h.write(';\n')

                            params_cpp.write(indent_pc)
                            if v["type"][0] == 's':
                                params_cpp.write(
                                    "std::")  # strings need "std::"
                            params_cpp.write(
                                v["type"][:-2] if v["type"][-2] == '[' else
                                v["type"])  # arrays: chop "[]"
                            params_cpp.write(' ')
                            params_cpp.write(table_names)
                            params_cpp.write(k)
                            # this writes the size declaration needed for arrays
                            if v["type"][-2] == '[':
                                params_cpp.write("[")
                                params_cpp.write(str(v["size"]))
                                params_cpp.write("]")
                            params_cpp.write(';\n')

                            param_read_def = "".join(
                                [param_read_def, "\t\t\tif(file.contains(\""])
                            param_read_def = "".join(
                                [param_read_def, "::".join(namespace_list)])
                            param_read_def = "".join([param_read_def, "::"])
                            param_read_def = "".join(
                                [param_read_def, table_names])
                            param_read_def = "".join([param_read_def, k])
                            param_read_def = "".join(
                                [param_read_def, "\"))\n\t\t\t{\n\t\t\t\t"])
                            if (v["type"][0] == 'i' or v["type"][0]
                                    == 'd') and v["type"][-2] != '[':
                                param_read_def = "".join(
                                    [param_read_def, "if ("])
                                param_read_def = "".join(
                                    [param_read_def,
                                     str(v["min"])])
                                param_read_def = "".join(
                                    [param_read_def, " > file[\""])
                                param_read_def = "".join([
                                    param_read_def, "::".join(namespace_list)
                                ])
                                param_read_def = "".join(
                                    [param_read_def, "::"])
                                param_read_def = "".join(
                                    [param_read_def, table_names])
                                param_read_def = "".join([param_read_def, k])
                                param_read_def = "".join(
                                    [param_read_def, "\"].as_"])
                                param_read_def = "".join([
                                    param_read_def, "floating"
                                    if v["type"][0] == 'd' else "integer"
                                ])
                                param_read_def = "".join(
                                    [param_read_def, "() || "])
                                param_read_def = "".join(
                                    [param_read_def,
                                     str(v["max"])])
                                param_read_def = "".join(
                                    [param_read_def, " < file[\""])
                                param_read_def = "".join([
                                    param_read_def, "::".join(namespace_list)
                                ])
                                param_read_def = "".join(
                                    [param_read_def, "::"])
                                param_read_def = "".join(
                                    [param_read_def, table_names])
                                param_read_def = "".join([param_read_def, k])
                                param_read_def = "".join(
                                    [param_read_def, "\"].as_"])
                                param_read_def = "".join([
                                    param_read_def, "floating"
                                    if v["type"][0] == 'd' else "integer"
                                ])
                                param_read_def = "".join([
                                    param_read_def,
                                    "())\n\t\t\t\t{\n\t\t\t\t\tstd::cerr << R\"(Invalid value for \""
                                ])
                                param_read_def = "".join([
                                    param_read_def, "::".join(namespace_list)
                                ])
                                param_read_def = "".join(
                                    [param_read_def, "::"])
                                param_read_def = "".join(
                                    [param_read_def, table_names])
                                param_read_def = "".join([param_read_def, k])
                                param_read_def = "".join([
                                    param_read_def,
                                    "\")\" << std::endl;\n\t\t\t\t\texit(-1);\n\t\t\t\t}\n\n\t\t\t\t"
                                ])
                            if v["type"][
                                    -2] == '[':  # only for arrays: write a loop to assign values
                                param_read_def = "".join(
                                    [param_read_def, "for (int i = 0; i < "])
                                param_read_def = "".join(
                                    [param_read_def,
                                     str(v["size"])])
                                param_read_def = "".join([
                                    param_read_def,
                                    "; ++i)\n\t\t\t\t{\n\t\t\t\t\t"
                                ])
                            param_read_def = "".join(
                                [param_read_def, "::".join(namespace_list)])
                            param_read_def = "".join([param_read_def, "::"])
                            param_read_def = "".join(
                                [param_read_def, table_names])
                            param_read_def = "".join([param_read_def, k])
                            param_read_def = "".join([
                                param_read_def, "[i] = file[\""
                                if v["type"][-2] == '[' else " = file[\""
                            ])  # array loop stuff
                            param_read_def = "".join(
                                [param_read_def, "::".join(namespace_list)])
                            param_read_def = "".join([param_read_def, "::"])
                            param_read_def = "".join(
                                [param_read_def, table_names])
                            param_read_def = "".join([param_read_def, k])
                            param_read_def = "".join([
                                param_read_def, "\"][i].as_"
                                if v["type"][-2] == '[' else "\"].as_"
                            ])  # array loop stuff
                            param_read_def = "".join([
                                param_read_def,
                                "floating" if v["type"][0] == 'd'
                                else  # pick correct type for parameter
                                "integer" if (v["type"][0] == 'i'
                                              or v["type"][0] == 'u') else
                                "string" if v["type"][0] == 's' else "boolean"
                            ])
                            param_read_def = "".join([
                                param_read_def, "();\n\t\t\t\t}\n\t\t\t}\n\n"
                                if v["type"][-2] == '[' else "();\n\t\t\t}\n\n"
                            ])
                            param_read_def = "".join([
                                param_read_def,
                                "\t\t\telse\n\t\t\t{\n\t\t\t\tstd::cerr << R\"(No value for \""
                            ])
                            param_read_def = "".join(
                                [param_read_def, table_names])
                            param_read_def = "".join([param_read_def, k])
                            param_read_def = "".join(
                                [param_read_def, "\"; \""])
                            param_read_def = "".join(
                                [param_read_def, table_names])
                            param_read_def = "".join([param_read_def, k])
                            param_read_def = "".join([
                                param_read_def,
                                "\" must be given a value)\" << std::endl;\n\t\t\t\texit(-1);\n\t\t\t}\n\n"
                            ])

                        elif v["class"][
                                0] == 's':  # If parameter is semivariant
                            # In parameters.h: declare extern
                            # In parameters.cpp: declare and define (using default)
                            # In grUtils.cpp: define from user input if given
                            params_h.write(indent_ph)
                            params_h.write(
                                "extern std::" if v["type"][0] == 's' else
                                "extern ")  # strings need "std::"
                            params_h.write(
                                v["type"][:-2] if v["type"][-2] == '[' else
                                v["type"])  # arrays: chop "[]"
                            params_h.write(' ')
                            params_h.write(table_names)
                            params_h.write(k)
                            if v["type"][
                                    -2] == '[':  # this writes the size declaration for arrays
                                params_h.write("[")
                                params_h.write(str(v["size"]))
                                params_h.write("]")
                            params_h.write(';\n')

                            params_cpp.write(indent_pc)
                            if v["type"][0] == 's':
                                params_cpp.write(
                                    "std::")  # strings need "std::"
                            params_cpp.write(
                                v["type"][:-2] if v["type"][-2] == '[' else
                                v["type"])  # arrays: chop "[]"
                            params_cpp.write(' ')
                            params_cpp.write(table_names)
                            params_cpp.write(k)
                            if v["type"][
                                    -2] == '[':  # this writes the special declaration code for arrays
                                params_cpp.write("[")
                                params_cpp.write(str(v["size"]))
                                params_cpp.write("] = {")
                                params_cpp.write(str(v["default"])[1:-1])
                                params_cpp.write("}")
                            else:  # if parameter is not an array, we just do a normal declaration
                                params_cpp.write(
                                    " = \"" if v["type"][0] == 's' else
                                    " = ")  # strings need quotes
                                params_cpp.write(
                                    str(v["default"]).lower() if v["type"][0]
                                    == 'b' else str(v["default"]))
                            params_cpp.write("\";\n" if v["type"][0] == 's'
                                             else ";\n")  # strings need quotes

                            param_read_def = "".join(
                                [param_read_def, "\t\t\tif(file.contains(\""])
                            param_read_def = "".join(
                                [param_read_def, "::".join(namespace_list)])
                            param_read_def = "".join([param_read_def, "::"])
                            param_read_def = "".join(
                                [param_read_def, table_names])
                            param_read_def = "".join([param_read_def, k])
                            param_read_def = "".join(
                                [param_read_def, "\"))\n\t\t\t{\n\t\t\t\t"])
                            if (v["type"][0] == 'i' or v["type"][0]
                                    == 'd') and v["type"][-2] != '[':
                                param_read_def = "".join(
                                    [param_read_def, "if ("])
                                param_read_def = "".join(
                                    [param_read_def,
                                     str(v["min"])])
                                param_read_def = "".join(
                                    [param_read_def, " > file[\""])
                                param_read_def = "".join([
                                    param_read_def, "::".join(namespace_list)
                                ])
                                param_read_def = "".join(
                                    [param_read_def, "::"])
                                param_read_def = "".join(
                                    [param_read_def, table_names])
                                param_read_def = "".join([param_read_def, k])
                                param_read_def = "".join(
                                    [param_read_def, "\"].as_"])
                                param_read_def = "".join([
                                    param_read_def, "floating"
                                    if v["type"][0] == 'd' else "integer"
                                ])
                                param_read_def = "".join(
                                    [param_read_def, "() || "])
                                param_read_def = "".join(
                                    [param_read_def,
                                     str(v["max"])])
                                param_read_def = "".join(
                                    [param_read_def, " < file[\""])
                                param_read_def = "".join([
                                    param_read_def, "::".join(namespace_list)
                                ])
                                param_read_def = "".join(
                                    [param_read_def, "::"])
                                param_read_def = "".join(
                                    [param_read_def, table_names])
                                param_read_def = "".join([param_read_def, k])
                                param_read_def = "".join(
                                    [param_read_def, "\"].as_"])
                                param_read_def = "".join([
                                    param_read_def, "floating"
                                    if v["type"][0] == 'd' else "integer"
                                ])
                                param_read_def = "".join([
                                    param_read_def,
                                    "())\n\t\t\t\t{\n\t\t\t\t\tstd::cerr << R\"(Invalid value for \""
                                ])
                                param_read_def = "".join([
                                    param_read_def, "::".join(namespace_list)
                                ])
                                param_read_def = "".join(
                                    [param_read_def, "::"])
                                param_read_def = "".join(
                                    [param_read_def, table_names])
                                param_read_def = "".join([param_read_def, k])
                                param_read_def = "".join([
                                    param_read_def,
                                    "\")\" << std::endl;\n\t\t\t\t\texit(-1);\n\t\t\t\t}\n\n"
                                ])
                                param_read_def = "".join(
                                    [param_read_def, "\t\t\t\t"])
                            if v["type"][
                                    -2] == '[':  # only for arrays: write a loop to assign values
                                param_read_def = "".join(
                                    [param_read_def, "for (int i = 0; i < "])
                                param_read_def = "".join(
                                    [param_read_def,
                                     str(v["size"])])
                                param_read_def = "".join([
                                    param_read_def,
                                    "; ++i)\n\t\t\t\t{\n\t\t\t\t\t"
                                ])
                            param_read_def = "".join([
                                param_read_def,
                                "::".join(curr_namespace_list_ph)
                            ])
                            param_read_def = "".join([param_read_def, "::"])
                            param_read_def = "".join(
                                [param_read_def, table_names])
                            param_read_def = "".join([param_read_def, k])
                            param_read_def = "".join([
                                param_read_def, "[i] = file[\""
                                if v["type"][-2] == '[' else " = file[\""
                            ])  # array loop stuff
                            param_read_def = "".join(
                                [param_read_def, "::".join(namespace_list)])
                            param_read_def = "".join([param_read_def, "::"])
                            param_read_def = "".join(
                                [param_read_def, table_names])
                            param_read_def = "".join([param_read_def, k])
                            param_read_def = "".join([
                                param_read_def, "\"][i].as_"
                                if v["type"][-2] == '[' else "\"].as_"
                            ])  # array loop stuff
                            param_read_def = "".join([
                                param_read_def,
                                "floating" if v["type"][0] == 'd'
                                else  # pick correct type for parameter
                                "integer" if (v["type"][0] == 'i'
                                              or v["type"][0] == 'u') else
                                "string" if v["type"][0] == 's' else "boolean"
                            ])
                            param_read_def = "".join([
                                param_read_def, "();\n\t\t\t\t}\n\t\t\t}\n\n"
                                if v["type"][-2] == '[' else "();\n\t\t\t}\n\n"
                            ])

                        elif v["class"][0] == 'i':  # If parameter is invariant
                            # In parameters.h: declare static const and define
                            # In parameters.cpp: do nothing
                            # In grUtils.cpp: do nothing
                            params_h.write(indent_ph)
                            params_h.write(
                                "static const std::" if v["type"][0] == 's'
                                else "static const ")  # strings need "std::"
                            params_h.write(
                                v["type"][:-2] if v["type"][-2] == '[' else
                                v["type"])  # arrays chop "[]"
                            params_h.write(' ')
                            params_h.write(table_names)
                            params_h.write(k)
                            if v["type"][
                                    -2] == '[':  # this writes special declaration code for arrays
                                params_h.write("[")
                                params_h.write(str(v["size"]))
                                params_h.write("] = {")
                                params_h.write(str(v["default"])[1:-1])
                                params_h.write('}')

                            else:
                                params_h.write(
                                    " = \"" if v["type"][0] == 's' else
                                    " = ")  # strings need quotes
                                params_h.write(str(v["default"]))
                            params_h.write("\";\n" if v["type"][0] == 's' else
                                           ";\n")  # strings need quotes

                        # here we write the mpi broadcast code (variant and semivariant parameters only)
                        if v["class"][0] != 'i':
                            if v["type"][
                                    -1] == ']':  # if parameter is an array
                                bcasts.append("".join([
                                    "\t\t\t",
                                    "MPI_Bcast(&(",
                                    "::".join(curr_namespace_list_ph),
                                    "::",
                                    table_names,
                                    k,
                                    "),",
                                    str(v["size"]),
                                    ",",
                                    "MPI_DOUBLE" if v["type"][0] == "d" else
                                    "MPI_INT",  # pick correct type
                                    ",0,comm);"
                                ]))

                            elif v["type"][
                                    0] == 's':  # if parameter is a string
                                bcasts.append("".join([
                                    "\t\t\t", "MPI_Bcast(const_cast<char*>(",
                                    "::".join(curr_namespace_list_ph), "::",
                                    table_names, k, ".c_str()),",
                                    "::".join(curr_namespace_list_ph), "::",
                                    table_names, k,
                                    ".size() + 1,MPI_CHAR,0,comm);"
                                ]))

                            else:  # if parameter is bool, double, or int
                                bcasts.append("".join([
                                    "\t\t\t", "par::Mpi_Bcast(&(",
                                    "::".join(curr_namespace_list_ph), "::",
                                    table_names, k, "),1,0,comm);"
                                ]))

                        indent_ph = indent_ph[:-1]
                        indent_pc = indent_pc[:-1]

                        break

    # close up all namespaces in parameters.h
    for i in range(len(curr_namespace_list_ph)):
        params_h.write(indent_ph)
        params_h.write("}\n")
        indent_ph = indent_ph[:-1]

    # close up all namespaces in paramters.cpp
    for i in range(len(curr_namespace_list_pc)):
        params_cpp.write(indent_pc)
        params_cpp.write("}\n")
        indent_pc = indent_pc[:-1]

    # The code this adds define the specifications for the blackholes using the special BH type
    param_read_def = "".join([
        param_read_def, '''            bssn::BH1 = BH(BSSN_BH1_MASS,
                           BSSN_BH1_X,
                           BSSN_BH1_Y,
                           BSSN_BH1_Z,
                           BSSN_BH1_V_X,
                           BSSN_BH1_V_Y,
                           BSSN_BH1_V_Z,
                           BSSN_BH1_SPIN,
                           BSSN_BH1_SPIN_THETA,
                           BSSN_BH1_SPIN_PHI);

            bssn::BH2 = BH(BSSN_BH2_MASS,
                           BSSN_BH2_X,
                           BSSN_BH2_Y,
                           BSSN_BH2_Z,
                           BSSN_BH2_V_X,
                           BSSN_BH2_V_Y,
                           BSSN_BH2_V_Z,
                           BSSN_BH2_SPIN,
                           BSSN_BH2_SPIN_THETA,
                           BSSN_BH2_SPIN_PHI);

            bssn::BSSN_BH_LOC[0]=Point(BH1.getBHCoordX(),BH1.getBHCoordY(),BH1.getBHCoordZ());
            bssn::BSSN_BH_LOC[1]=Point(BH2.getBHCoordX(),BH2.getBHCoordY(),BH2.getBHCoordZ());

'''
    ])

    # add mpi bcasts
    param_read_def = "".join([param_read_def, "\n".join(bcasts)])
    param_read_def = "".join([param_read_def, "\n\t\t}\n\t}"])

    # add code for readParamFile and dumpParamFile
    params_cpp.write(param_read_def)
    params_cpp.write('\n')
    params_cpp.write(param_dump_def)

    # Close some braces in parameters.cpp, add footer for parameters.h
    params_cpp.write("\t\t}\n\t}\n}")
    params_h.write(
        "#include \"grUtils.tcc\"\n\n#endif  //SFCSORTBENCH_PARAMETERS_H")
