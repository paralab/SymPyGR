import sympy as sym
import dendrosym
from bssn_eqns import dendroConfigs
import matplotlib.pyplot as plt
import networkx as nx

idx = ""
iters = 1000

eqns, vnames = dendroConfigs.get_rhs_eqns_flat("evolution")

# now build up the actual expression graph
g = dendrosym.nxgraph.ExpressionGraph()

g.add_expressions(eqns, vnames, verbose=True)

G = g.get_composed_graph(verbose=False)

g.plot_adjmatrix()

print(g._nx_graphs.keys())

g.draw_all_graphs(verbose=True)

plt.figure()
nx.draw_kamada_kawai(G)
plt.title("COMPOSED GRAPH")
plt.savefig("composed_graph.png")
plt.show()

print("|V|= %d |E| =%d" % (G.number_of_nodes(), G.number_of_edges()))
