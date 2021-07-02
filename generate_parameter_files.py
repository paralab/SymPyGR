"""Simple file for running all of the parameter generation

This basically is an example for how the parameter code is generated.
"""

import argparse
import dendrosym

PARAM_H_BEGINNING = '''#ifndef SFCSORTBENCH_PARAMETERS_H
#define SFCSORTBENCH_PARAMETERS_H

#include "bh.h"
#include "grDef.h"
#include <string.h>
#include <iostream>
#include "dendro.h"
#include "memory_pool.h"
#include "toml.hpp"

'''

PARAM_C_BEGINNING = '''#include "parameters.h"

'''

PARAM_H_ENDING = '''
#endif //SFCSORTBENCH_PARAMETERS_H

'''


def main(template_file, project_name, output_param_h, output_param_c,
         output_param_sample):
    """Main function for this basic script

    Creates the parameter header, source, and sample files given
    the template file, a project name, and the names of the files
    that should be saved.
    """

    print("Now generating the text for:", template_file)

    paramh_str, paramc_str = dendrosym.params.generate_all_parameter_text(
        project_name, template_file)

    param_sample = dendrosym.params.generate_sample_config_file_text(
        project_name, template_file)

    print("Finished generating the text, now saving to output files...")
    # now we can start by opening the h file
    with open(output_param_h, "w") as f:
        # write the beginning of the file
        f.write(PARAM_H_BEGINNING)

        # write the generated string
        f.write(paramh_str)

        # write the end of the file
        f.write(PARAM_H_ENDING)

    # now we move on to the cpp source file
    with open(output_param_c, "w") as f:
        # write the beginning of the file
        f.write(PARAM_C_BEGINNING)

        # write the generated string
        f.write(paramc_str)

    # then the parameter sample file
    with open(output_param_sample, "w") as f:
        # no beginning here, the whole thing is generated
        f.write(param_sample)

    print("Finished generating the files!")
    print("The header file has been written to:", output_param_h)
    print("The cpp file has been written to:", output_param_c)
    print("The sample parameter TOML file has been written to:",
          output_param_sample)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Create the C++ code for parameters " +
        "based on a single template file")

    parser.add_argument("template_file", type=str)

    parser.add_argument("project_name", type=str)

    parser.add_argument(
        '-hpp',
        "--output_param_h",
        default="parameters.h",
        type=str,
        help="File name (with optional path) for the header file")

    parser.add_argument(
        "-cpp",
        "--output_param_c",
        default="parameters.cpp",
        type=str,
        help="File name (with optional path) for the source file")

    parser.add_argument(
        "-s",
        "--output_param_sample",
        default="params-sample.toml",
        type=str,
        help="File name (with optional path) for the sample parameter file")

    args = parser.parse_args()

    main(**vars(args))
