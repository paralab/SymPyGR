import sympy as sym
import dendrosym
from emda_configs import dendroConfigs
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy as sp

import pathlib

import graph_coarsening

import pygsp
import warnings
import matplotlib.cbook

warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

method = "variation_neighborhood"

r = 0.6
k = 5
kmax = int(3 * k)
max_depth = 3

coordSaveName = "emda_coords.npy"

eqns, vnames = dendroConfigs.get_rhs_eqns_flat("evolution")

# now build up the actual expression graph
g = dendrosym.nxgraph.ExpressionGraph()

g.add_expressions(eqns, vnames, verbose=True)

g1 = g.get_composed_graph()

print(f"The composed graph has {len(g1.nodes)} nodes, {len(g1.edges)} edges")

g1_sparse_directed = nx.to_scipy_sparse_array(g1)

g1_undirected = g1.to_undirected()
g1_sparse = nx.to_scipy_sparse_array(g1_undirected)

# convert to pygsp graph
g1_pygsp = pygsp.graphs.Graph(g1_sparse)

# PLOT
plt.figure()
plt.spy(g1_pygsp.W, markersize=1.0)
plt.savefig("EMDA_composed_undirected.png")

plt.figure()
plt.spy(g1_sparse_directed, markersize=1.0)
plt.savefig("EMDA_composed_directed.png")

print(f"pygsp graph: {g1_pygsp.N} nodes, {g1_pygsp.Ne} edges")

# ==========
# CALCULATE THE COORDINATES
coordfile = pathlib.Path(coordSaveName)

if coordfile.exists():
    print("Found coordinate file, now loading...")
    coords = np.load(coordfile)
else:
    print("Calculating coordinates...")
    pos = nx.nx_agraph.graphviz_layout(g1, prog="sfdp")
    coords = np.array(list(pos.values()))
    np.save(coordfile, coords)

# # set the coordinates to the g1_pygsp graph
g1_pygsp.set_coordinates(coords)


# =========
# graph coarsening
C, Gc, Call, Gall = graph_coarsening.coarsen(g1_pygsp, K=k, r=r, method=method)

plt.figure()
plt.spy(Gc.W, markersize=1.0)
plt.savefig("EMDA_coarsened.png")

fig = graph_coarsening.plot_coarsening(Gall, Call, size=5, title=method)
fig.savefig("EMDA_coarsened_nodeview.png")


g1_pygsp.compute_fourier_basis()
N = g1_pygsp.N
L = g1_pygsp.L.toarray()
S = graph_coarsening.graph_utils.get_S(g1_pygsp).T

metrics = graph_coarsening.coarsening_quality(g1_pygsp, C, kmax=kmax)
n = Gc.N

print("Now making error figure")
plt.figure()
size = 2.04; fig, axes = plt.subplots(1, 3, figsize=(4*size*3, 3*size)); lineWidth = 1

axes[0].plot(np.arange(kmax), np.abs(metrics['error_subspace']), 'or-') 
axes[0].set_xlabel('$k$'); axes[0].set_ylabel('$\epsilon$')
axes[0].plot( [k, k], [0, max(metrics['error_subspace'])], ':k') 

axes[1].boxplot(np.abs(metrics['error_eigenvalue'])) 
axes[1].set_ylabel('relative eigenvalue error')

axes[2].imshow(abs(metrics['angle_matrix'][:,0:kmax]) )
axes[2].plot( [k, k], [0, kmax], ':w') 
axes[2].plot( [0, kmax], [k, k], ':w') 
axes[2].plot( [0, N], [n-1, n-1], ':w') 
axes[2].set_xlim([0, kmax-1])
axes[2].set_ylim([0, kmax-1])
axes[2].set_xlabel('Lc eigenvectors lifted'); axes[2].set_ylabel('L eigenvectors');
plt.savefig("EMDA_graph_stats.png")


# get the number of levels
n_levels = len(Gall) - 1

print(n_levels)

colors = ["k", "g", "b", "r", "y"]

for level in range(n_levels):
    print(level)
    G_temp = Gall[level]
    edges = np.array(G_temp.get_edge_list()[0:2])

    Gc_temp = Gall[level + 1]

    edges_c = np.array(Gc_temp.get_edge_list()[0:2])

    Cs = Call[level]
    Cs = Cs.toarray()

    print(edges.shape)

    print(Gc_temp.N)

    # gather the nodes that we want to try and look at

    node_depth = {}

    for i in range(Gc_temp.N):
        subgraph = np.arange(G_temp.N)[Cs[i, :] > 0]

        list_node_depth = node_depth.get(len(subgraph) - 1, [])

        list_node_depth.append(i)

        node_depth[len(subgraph) - 1] = list_node_depth

        # clipped_subgraph = colors[np.clip(len(subgraph) - 1, 0, 4)]

    # print(node_depth)

    all_depths = list(node_depth.keys())
    all_depths.sort()

    print(all_depths)

    # so then we want to start with the nodes that are deepest

    # then match up the node IDs from before with coordinates (since the original graph loses it all for some reason?)
    # TODO: there might be a smarter way with labels?

    for depth in reversed(all_depths):
        # print(depth)
        # get the one we want
        print(
            "For depth", depth, " There are ", len(node_depth[depth]), "values"
        )
        for idx in node_depth[depth]:
            # print("  ", idx)

            # compare G_temp coords with the first graph
            # g1_pygsp.coords vs G_temp.coords
            if not np.all(g1_pygsp.coords[idx] == G_temp.coords[idx]):
                print("FAILURE FOUND")

            # okay, so now we have it on the first round of coarsening, and we have a nice list!

# now with the list, we can pull out the various nodes that are identified as most important

# convert the graph to a dictionary so we can get the list of nodes
dictofnodes = nx.to_dict_of_dicts(g1)

nodes_sorted_by_idx = list(dictofnodes.keys())

# then we can just go through our node depth reversed

current_expr_idx = 0
exprs_to_compute_first = []
graph_simpl_names = []
updated_exprs = eqns.copy()

for depth in reversed(all_depths):
    if depth > max_depth:
        for idx in node_depth[depth]:
            # print(nodes_sorted_by_idx[idx])

            expr_to_replace = nodes_sorted_by_idx[idx]

            expr_temp_var = sym.Symbol(
                f"DENDRO_GRAPH_SIMPL_{current_expr_idx:0>4}"
            )

            exprs_to_compute_first.append(expr_to_replace)
            graph_simpl_names.append(expr_temp_var)

            print(expr_temp_var)

            # then we go through our equations and try to replace it
            for eq_idx, eqn in enumerate(updated_exprs):
                eqn = eqn.xreplace({expr_to_replace: expr_temp_var})
                # eqn = eqn.replace(expr_to_replace, expr_temp_var)

                updated_exprs[eq_idx] = eqn

            current_expr_idx += 1

            # done
    # done again

# now with updated_exprs and exprs_to_compute_first we can generate the C++ code

# start with the exprs_to_compute_first code

# count all operations
orig_num_ops = sym.count_ops(eqns)

# then start by constructing the cse from expressions
coarsened_ops = sym.count_ops(exprs_to_compute_first)
cse_exp = dendrosym.codegen.construct_cse_from_list(exprs_to_compute_first)
# then we can generate CPU preextracted
from_coarsened = dendrosym.codegen.generate_cpu_preextracted(cse_exp, graph_simpl_names, "", coarsened_ops)

with open("graph_simplified_portion.cpp", "w") as f:
    f.write(from_coarsened)

# then we can generate CSE for the rest of it
orig_ops_after_coarsen = sym.count_ops(updated_exprs)


# cse time
cse_exp = dendrosym.codegen.construct_cse_from_list(updated_exprs)

after_coarsened = dendrosym.codegen.generate_cpu_preextracted(cse_exp, vnames, "", orig_ops_after_coarsen)

with open("post_graph_simplified_portion.cpp", "w") as f:
    f.write(after_coarsened)

print("")
print(f"NOTE: original number of operations was ", orig_num_ops)
