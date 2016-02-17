
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


er=nx.erdos_renyi_graph(20,0.05) # 1st arg: nb of node, 2nd arg: proba to create an edge
nx.draw(er)
#plt.show()
print er.get_edge_attributes(er)