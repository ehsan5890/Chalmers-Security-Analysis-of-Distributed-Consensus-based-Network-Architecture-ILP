# -*- coding: utf-8 -*-
"""
Created on Wed May  8 00:59:47 2024

@author: rodri
"""

from gurobipy import *
import math
import random
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


def optimization_result(V, N: int, c: int, C:int):

    nw_mod = Model(name="nw")


    delta_ij = nw_mod.addVars(V.number_of_nodes(), V.number_of_nodes(), name='delta_{i}_{j}', vtype=GRB.BINARY)
    x = nw_mod.addVars(V.number_of_nodes(), name='x', vtype=GRB.BINARY)
    y = nw_mod.addVars(V.number_of_nodes(), name='y', vtype=GRB.BINARY)
    z = nw_mod.addVars(V.number_of_nodes(), name='z', vtype=GRB.BINARY)

    obj_fun = sum(c * delta_ij[i, j] for i, j in V.edges())
    nw_mod.setObjective(obj_fun, GRB.MINIMIZE)

    for i in V.nodes():
        nw_mod.addConstr(x[i] + y[i] + z[i] == 1, 'constraint_one_subd{i}')

    # nw_mod.addConstr(sum(V.nodes[i]['S'] * x[i] for i in V.nodes()) <= (math.ceil(N / 2)) - 1, 'constraint_sum_S_x')
    # nw_mod.addConstr(sum(V.nodes[i]['S'] * y[i] for i in V.nodes()) <= (math.ceil(N / 2)) - 1, 'constraint_sum_S_y')
    # nw_mod.addConstr(sum(V.nodes[i]['S'] * z[i] for i in V.nodes()) <= (math.ceil(N / 2)) - 1, 'constraint_sum_S_z')
    nw_mod.addConstr(sum(V.nodes[i]['S'] * x[i] for i in V.nodes()) <= N//2, 'constraint_sum_S_x')
    nw_mod.addConstr(sum(V.nodes[i]['S'] * y[i] for i in V.nodes()) <= N//2, 'constraint_sum_S_y')
    nw_mod.addConstr(sum(V.nodes[i]['S'] * z[i] for i in V.nodes()) <= N//2, 'constraint_sum_S_z')

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

    nw_mod.setParam('OutputFlag', False)
    nw_mod.optimize()
    deltas = {}
    for v in nw_mod.getVars():
        if re.match('delta', v.varName) is not None and v.x != 0:
            deltas[v.varName] = v.x

    edges_to_remove = []
    for i, j in V.edges():
        if delta_ij[i, j].X == 1:
            edges_to_remove.append((i, j))
    # nw_mod.computeIIS()
    # print('Optimization is done. Objective function value: %.2f' % nw_mod.objVal)

    return (nw_mod.objVal, deltas, edges_to_remove)


graph = nx.Graph()

V = nx.read_gml('./atlanta.gml', destringizer=int)
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
# nx.draw_networkx(V, pos=node_positions, with_labels=True)
min_vals_arb = []

C = 100
c = 10

# end=0
# if V.number_of_nodes()%2==0:
#     end=V.number_of_nodes()-1
# else:
#     end=V.number_of_nodes()
# N=5
# costs = []
# while N<=end:
#     meth = 0
#     number_dcn_sharden1 = N // 2
#     number_dcn_sharden2 = N//2 + 1
#     V_list = list(V.nodes())
#     set1 = random.sample(V_list, number_dcn_sharden1)
#     remaining_nodes = [node for node in V_list if node not in set1]
#     # Select 3 unique nodes for the second set
#     set2 = random.sample(remaining_nodes, number_dcn_sharden2)
#     while meth < 2:
#         if meth == 1:
#             for i in V.nodes():
#                 if i in set1:
#                     V.nodes[i]['S'] = 1
#                 else:
#                     V.nodes[i]['S'] = 0
#             objVal, deltas, edge_to_remove1 = optimization_result(V, number_dcn_sharden1, c, C)
#         else:
#             for i in V.nodes():
#                 if i in set2:
#                     V.nodes[i]['S'] = 1
#                 else:
#                     V.nodes[i]['S'] = 0
#             objVal1, deltas1, edge_to_remove2 = optimization_result(V, number_dcn_sharden2, c, C)
#         meth +=1
#
#
#     common_elements = list(set(edge_to_remove1) & set(edge_to_remove2))
#     non_common_elements = list(set(edge_to_remove1) ^ set(edge_to_remove2))
#     cost = (len(non_common_elements) + len(common_elements)) *c
#     costs.append(cost)
#     N+=2
#
#
# print(costs)


# Define the number of iterations
iterations = 300
all_costs = []

# Run the code 300 times
for _ in range(iterations):
    costs = []
    end = 0
    if V.number_of_nodes() % 2 == 0:
        end = V.number_of_nodes() - 1
    else:
        end = V.number_of_nodes()
    N = 5
    while N <= end:
        meth = 0
        number_dcn_sharden1 = N // 2
        number_dcn_sharden2 = N // 2 + 1
        V_list = list(V.nodes())
        set1 = random.sample(V_list, number_dcn_sharden1)
        remaining_nodes = [node for node in V_list if node not in set1]
        # Select unique nodes for the second set
        set2 = random.sample(remaining_nodes, number_dcn_sharden2)
        while meth < 2:
            if meth == 1:
                for i in V.nodes():
                    if i in set1:
                        V.nodes[i]['S'] = 1
                    else:
                        V.nodes[i]['S'] = 0
                objVal, deltas, edge_to_remove1 = optimization_result(V, number_dcn_sharden1, c, C)
            else:
                for i in V.nodes():
                    if i in set2:
                        V.nodes[i]['S'] = 1
                    else:
                        V.nodes[i]['S'] = 0
                objVal1, deltas1, edge_to_remove2 = optimization_result(V, number_dcn_sharden2, c, C)
            meth += 1

        common_elements = list(set(edge_to_remove1) & set(edge_to_remove2))
        non_common_elements = list(set(edge_to_remove1) ^ set(edge_to_remove2))
        cost = (len(non_common_elements) + len(common_elements)) * c
        costs.append(cost)
        N += 2

    # Store this run's costs
    all_costs.append(costs)

# Calculate the element-wise average
average_costs = np.mean(all_costs, axis=0)
print("Average Costs for each element over 300 runs:", average_costs)



