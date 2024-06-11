# -*- coding: utf-8 -*-
"""
Created on Wed May  8 00:59:47 2024

@author: rodri
"""

from gurobipy import *
import math
import networkx as nx
import matplotlib.pyplot as plt
import random
import pandas as pd
import re
import io
from random import randrange
import numpy as np
import itertools
import time
graph=nx.Graph()

V=nx.read_gml('./atlanta.gml',destringizer=int)
N=3
print('This are the node locations for the selected topology:')
print(V.nodes(data='name'))
print('\\')
print('This are the existing links for the selected topology:')
print(V.edges(data='name'))
node_positions = {}
print(V.number_of_edges())
for node, attributes in V.nodes(data=True):
    
    graphics = attributes.get('graphics', {})
    x = graphics.get('x', 0)
    y = graphics.get('y', 0)
    
    node_positions[node] = (x, y)
#nx.draw_networkx(V, pos=node_positions, with_labels=True)
min_vals_arb=[]
combinations = list(itertools.combinations(V.nodes(), N))
print("First few combinations:", combinations[:5])
print("Total number of combinations: ", len(combinations))


start = time.time()
#print("hello")

for k in range(len(combinations)):
    print("combination " +str(k)+" out of "+str(len(combinations)))
    #We arbitrarly select DCNs
    #DCN_sel=[]
    #i=0
    #while i<N:
     #   tmp=randrange(V.number_of_nodes())
      #  if tmp not in DCN_sel:
       #     DCN_sel.append(tmp)
        #    i=i+1
    #print(DCN_sel)
    DCN_sel=combinations[k]
    for i in V.nodes():
        if i in DCN_sel:
            V.nodes[i]['S']=1
            
        else:
            V.nodes[i]['S']=0


    nw_mod=Model(name="nw")


    C=100
    c=10


    delta_ij = nw_mod.addVars(V.number_of_nodes(),V.number_of_nodes(),name='delta_{i}_{j}',vtype=GRB.BINARY)
    x=nw_mod.addVars(V.number_of_nodes(),name='x',vtype=GRB.BINARY)
    y=nw_mod.addVars(V.number_of_nodes(),name='y',vtype=GRB.BINARY)
    z=nw_mod.addVars(V.number_of_nodes(),name='z',vtype=GRB.BINARY)

    obj_fun=sum(c*delta_ij[i,j] for i,j in V.edges())
    nw_mod.setObjective(obj_fun, GRB.MINIMIZE)

    for i in V.nodes():
        nw_mod.addConstr(x[i] + y[i] + z[i] == 1, 'constraint_one_subd{i}')
        

    nw_mod.addConstr(sum(V.nodes[i]['S'] * x[i] for i in V.nodes()) <= (math.ceil(N/2)) - 1, 'constraint_sum_S_x')
    nw_mod.addConstr(sum(V.nodes[i]['S'] * y[i] for i in V.nodes()) <= (math.ceil(N/2)) - 1, 'constraint_sum_S_y')
    nw_mod.addConstr(sum(V.nodes[i]['S'] * z[i] for i in V.nodes()) <= (math.ceil(N/2)) - 1, 'constraint_sum_S_z')


    for i, j in V.edges():
        nw_mod.addConstr(delta_ij[i, j] <= 2 - x[i] - x[j], f'constraint_delta_ij_x_{i}_{j}_same_x')
        nw_mod.addConstr(delta_ij[j, i] <= 2 - x[i] - x[j], f'constraint_delta_ij_x_{i}_{j}_same_x_2')
        nw_mod.addConstr(delta_ij[i, j] <= 2 - y[i] - y[j], f'constraint_delta_ij_x_{i}_{j}_same_y')
        nw_mod.addConstr(delta_ij[j, i] <= 2 - y[i] - y[j], f'constraint_delta_ij_x_{i}_{j}_same_y_2')
        nw_mod.addConstr(delta_ij[i, j] <= 2 - z[i] - z[j], f'constraint_delta_ij_x_{i}_{j}_same_z')
        nw_mod.addConstr(delta_ij[j, i] <= 2 - z[i] - z[j], f'constraint_delta_ij_x_{i}_{j}_same_z_2')
        nw_mod.addConstr(delta_ij[i, j] >= x[i] - x[j], f'constraint_delta_ij_x_{i}_{j}_not_same_x')
        nw_mod.addConstr(delta_ij[j, i] >= x[i] - x[j], f'constraint_delta_ij_x_{i}_{j}_not_same_x_2')
        nw_mod.addConstr(delta_ij[i, j] >= y[i] - y[j], f'constraint_delta_ij_x_{i}_{j}_not_same_y')
        nw_mod.addConstr(delta_ij[j, i] >= y[i] - y[j], f'constraint_delta_ij_x_{i}_{j}_not_same_y_2')
        nw_mod.addConstr(delta_ij[i, j] >= z[i] - z[j], f'constraint_delta_ij_x_{i}_{j}_not_same_z')
        nw_mod.addConstr(delta_ij[j, i] >= z[i] - z[j], f'constraint_delta_ij_x_{i}_{j}_not_same_z_2')

    nw_mod.setParam('OutputFlag',False)
    nw_mod.optimize()
    #nw_mod.computeIIS()
    #print('Optimization is done. Objective function value: %.2f' % nw_mod.objVal)
    min_vals_arb.append(nw_mod.objVal)

end = time.time()
print("Execution time in seconds")
print(end - start)
(uniq, freq) = (np.unique(min_vals_arb, return_counts=True))
print(np.column_stack((uniq,freq)))

plt.xlabel('Objective Function Value', fontsize=12)
plt.ylabel('Number of occurrences [1000 iterations]', fontsize=12)
plt.title('US-CA Arbitrary DCN positioning [N=5, c_ij=10]')
plt.bar(uniq,freq)
plt.savefig('./foo.png')

