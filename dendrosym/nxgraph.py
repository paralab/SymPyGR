"""
@author : Milinda Fernando, David Van Komen
@brief  : Compute a directed graph from a sympy expression.
    1. Rename custom functions to symbols
    2. Graphs are merged from multiple expressions
NetworkX is required, and used to store the graph.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use,copy, modify, merge, publish, distribute, sublicense,and/or sell copies
of the Software,and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
"""

import sympy as sympy
import networkx as nx
import matplotlib.pyplot as plt


class ExpressionGraph:
    def __init__(self):
        self._nx_graphs = dict()
        self._sympy_expr = dict()

        self._G_ = None

    def __pre_traversal_1(self, expr, node_list, edge_list):
        """Preorder traversal of the expression converting undefined functions
        to sympy symbols
        """
        if isinstance(expr.func, sympy.core.function.UndefinedFunction):
            sym_name = str(expr.func)
            for a in expr.args:
                sym_name = sym_name + "_" + str(a)

            node_list.append(sympy.Symbol(sym_name))
        else:
            node_list.append(expr)

        for arg in expr.args:
            if isinstance(arg.func, sympy.core.function.UndefinedFunction):
                f = arg.func
                sym_name = str(f)
                for a in arg.args:
                    sym_name = sym_name + "_" + str(a)

                node_list.append(sympy.Symbol(sym_name))
                edge_list.append((expr, sympy.Symbol(sym_name)))
            else:
                edge_list.append((expr, arg))
                self.__pre_traversal_1(arg, node_list, edge_list)

    def __pre_traversal_2(self, expr, node_list, edge_list):
        """Perform the pretraversal, second depth

        Keep undefined function references as it is pruning
        but not renaming.
        """
        if isinstance(expr.func, sympy.core.function.UndefinedFunction):
            # sym_name=str(expr.func)
            # for a in expr.args:
            #     sym_name = sym_name + '_' + str(a)
            node_list.append(expr)
        else:
            node_list.append(expr)

        for arg in expr.args:
            if isinstance(arg.func, sympy.core.function.UndefinedFunction):
                f = arg.func
                sym_name = str(f)
                for a in arg.args:
                    sym_name = sym_name + "_" + str(a)

                node_list.append(sympy.Symbol(sym_name))
                edge_list.append((expr, sympy.Symbol(sym_name)))
            else:
                edge_list.append((expr, arg))
                self.__pre_traversal_2(arg, node_list, edge_list)

    def __preorder_traversal__(self, expr):
        """Pre order traversal and returns the node list and edge list"""

        expr_list = list()
        edge_list = list()
        self.__pre_traversal_2(expr, expr_list, edge_list)
        return [expr_list, edge_list]

    def add_expression(self, expr, expr_name):
        """Generate a networkx graph for a given expression"""

        G = nx.DiGraph(vname=str(expr_name))
        G.name = str(expr_name)

        # total nodes on the graph.
        [node_list, edge_list] = self.__preorder_traversal__(expr)

        G.add_nodes_from(node_list)
        G.add_edges_from(edge_list)

        self._nx_graphs[str(expr_name)] = G
        self._sympy_expr[str(expr_name)] = expr

    def add_expressions(self, outs, vnames, suffix_idx="", verbose=False):
        """Adds list of sympy expressions"""
        mi = [0, 1, 2, 4, 5, 8]
        midx = ["00", "01", "02", "11", "12", "22"]

        num_e = 0
        for i, e in enumerate(outs):
            if type(e) == list:
                num_e = num_e + len(e)
                for j, ev in enumerate(e):
                    expr_name = vnames[i] + "" + str(j) + str(suffix_idx)
                    if verbose:
                        print(
                            "processing expr : %d var name %s[%s]"
                            % (i, vnames[i], str(j))
                        )
                    self.add_expression(ev, expr_name)
            elif type(e) == sympy.Matrix:
                num_e = num_e + len(e)
                for j, k in enumerate(mi):
                    expr_name = vnames[i] + "" + str(midx[j]) + str(suffix_idx)
                    if verbose:
                        print(
                            "processing expr : %d var name %s[%s]"
                            % (i, vnames[i], midx[j])
                        )
                    self.add_expression(e[k], expr_name)
            else:
                num_e = num_e + 1
                if verbose:
                    print("processing expr : %d var name %s" % (i, vnames[i]))
                expr_name = vnames[i] + str(suffix_idx)
                self.add_expression(e, expr_name)

    def composed_graph(self, verbose=False):
        """Calculates and computes the composed graph

        This is no longer supported, it is preferred to use the
        get_composed_graph function.
        """

        print(
            "WARNING: composed_graph() is a depreciated function, "
            + "please update to get_composed_graph"
        )

        if self._G_ is None:
            self.compute_composed_graph(verbose=verbose)

        return self._G_

    def compute_composed_graph(self, verbose=False):
        """Computes the composed graph, can be called any time
        there is an update
        """
        g_list = list()
        for (v, g) in self._nx_graphs.items():

            if verbose:
                print("graph for var: %s" % str(v))
                # print(nx.info(g))

            g_list.append(g)

        self._G_ = nx.compose_all(g_list)
        self._G_.name = "full graph added expressions"
        if verbose:
            print("full graph")
            # print(nx.info(self._G_))

    def get_composed_graph(self, verbose=False):
        """Get's the currently stored/generated composed graph"""

        if self._G_ is None:
            self.compute_composed_graph(verbose=verbose)

        return self._G_

    def plot_adjmatrix(self):
        """plots the adjacency matrix for the composed big G"""

        A = nx.adjacency_matrix(self._G_)
        plt.figure()
        plt.spy(A, markersize=3)
        plt.savefig("testfig.png")

    def draw_graph(self, expr_name, save_name=None):
        """plots the graph for a given sympy expression"""

        g = self._nx_graphs[expr_name]

        is_planar, embedding = nx.check_planarity(g)

        plt.figure()

        if not is_planar:
            print(
                f"WARNING: {expr_name} is not planar! "
                + "Will draw with the Kamada Kawai algorithm, "
                + "this may take a while!"
            )
            nx.draw_kamada_kawai(g)
        else:
            nx.draw_networkx(g, pos=nx.planar_layout(g), font_size=6)

        plt.title(f"Graph for {expr_name}")

        if save_name is not None:
            plt.savefig(save_name)

        plt.show()

    def draw_all_graphs(
        self, save_prefix="graph", save_format="png", verbose=False
    ):
        """Plots all individual graphs automatically"""

        for key in self._nx_graphs.keys():
            if verbose:
                print("NOW DRAWING GRAPH FOR", key)
            self.draw_graph(
                key, save_name=f"{save_prefix}_{key}.{save_format}"
            )

    def get_graph(self, expr_name):
        """
        Get the computational graph for a given expression
        """
        g = self._nx_graphs.get(expr_name, None)
        if g is None:
            raise NameError(f"Graph by name {expr_name} not found!")
        return g
