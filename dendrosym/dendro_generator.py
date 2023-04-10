"""dendro_generator.py

This file contains the class implementation for the Dendro
Generator. This class provides the functionality for 
"""

import datetime
import sympy
from sympy.printing import print_ccode
import git

from pathlib import Path

GIT_REPOSITORY_URL = "https://github.com/dfvankomen/dendrogr-main.git"


class DendroGenerator:
    def __init__(
        self, variable_list, evolution_expressions, configuration_name="solver"
    ) -> None:

        self.folder = None
        self.configuration_name = configuration_name

        self.variable_list = variable_list
        self.evolution_expressions = evolution_expressions

        pass

    def generate_code(self, folder=""):

        # convert incoming folder to a Path object
        self.folder = Path(folder)
        # create the directory so that we have a place it resides
        self.folder.mkdir(parents=True, exist_ok=True)

        # start by pulling the repository to the folder
        # self._pull_git_code()

        # then go ahead and start generating the code
        out_code = (
            "##### GENERATED PYTHON SCRIPT CODE FOR USE WITH DENDRO\n"
            f"# Code Generated on {datetime.date.today()}\n\n"
            "# ==== DENDRO INITIALIZATION"
            "import dendrosym\nimport sympy as sym\n"
            f'dendroConfigs = dendrosym.NRConfig("{self.configuration_name}")\n\n'
        )

        # TODO: parameter variables

        # TODO: constraint variables, if given

        # now for the physical variables
        out_code += "# ==== VARIABLE INITIALIZATION\n"
        temp_str = "dendroConfigs.add_evolution_variables(["
        for ii, var_name in enumerate(self.variable_list):
            out_code += f'{var_name} = dendrosym.dtypes.scalar("{var_name}")\n'
            temp_str += f"{var_name}" + (
                ", " if ii != len(self.variable_list) - 1 else ""
            )
        out_code += temp_str + "])\n\n"

        # set metric? or are we expecting people to do this whole thing themselves?
        # out_code += "dendroConfigs.set_metric()"

        # RHS equation calculation
        out_code += "# ==== EVOLUTION RIGHT HAND SIDE EQUATIONS\n"
        out_code += "def evolution_rhs_equations():\n"
        tb = "    "
        temp_str = f"{tb}rhs_list = ["
        temp_str2 = f"{tb}var_list = ["
        # TODO: modify for if the derivatives are second order!
        for ii, (var, eqn) in enumerate(
            zip(self.variable_list, self.evolution_expressions)
        ):
            # get the pure sympy code
            sym_code = sympy.sstr(eqn)

            # TODO: any transformations? Like converting derivatives
            out_code += f"{tb}{var}_rhs = {sym_code}\n"

            temp_str += f"{var}_rhs" + (
                ", " if ii != len(self.variable_list) - 1 else ""
            )
            temp_str2 += f"{var}" + (
                ", " if ii != len(self.variable_list) - 1 else ""
            )

        out_code += "\n"
        out_code += temp_str + "]\n"
        out_code += temp_str2 + "]\n"
        out_code += f"\n{tb}return rhs_list, var_list\n\n"

        # TODO: physical constraint equations function

        out_code += 'dendroConfigs.set_rhs_equation_function("evolution", evolution_rhs_equations)\n'

        # TODO: boundaries and evolution constraints like trace_zero, pos_floor

        print(out_code)

        pass

    def _pull_git_code(self):
        # From GitPython's documentation:
        empty_repo = git.Repo.init(self.folder, "empty")

        origin = empty_repo.create_remote("origin", GIT_REPOSITORY_URL)

        assert origin.exists()
        assert (
            origin == empty_repo.remotes.origin == empty_repo.remotes["origin"]
        )

        origin.fetch()

        # create the head for main, get the tracking branch, and checkout
        empty_repo.create_head("main", origin.refs.main).set_tracking_branch(
            origin.refs.main
        ).checkout()
