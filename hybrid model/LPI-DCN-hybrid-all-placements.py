# -*- coding: utf-8 -*-
"""
Created on Thu May  9 22:39:49 2024

@author: rodri
"""
from collections import Counter

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
N=5
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

u_avg=[]
delta_avg=[]
c=10
t=30
u_zero_count = 0
delta_zero_count = 0
comb=[]
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
    '''print("Nodes selected arbitrarily to work as DCN")
    for i in V.nodes():
        if V.nodes[i]['S']==1:
            print(V.nodes[i]["name"])
    print(V.nodes(data='S'))'''

    nw_mod=Model(name="nw")
    #nw_mod=Model(name="nw")


    
    delta_ij = nw_mod.addVars(V.number_of_nodes(),V.number_of_nodes(),name='delta_{i}_{j}',vtype=GRB.BINARY)
    x=nw_mod.addVars(V.number_of_nodes(),name='x',vtype=GRB.BINARY)
    y=nw_mod.addVars(V.number_of_nodes(),name='y',vtype=GRB.BINARY)
    z=nw_mod.addVars(V.number_of_nodes(),name='z',vtype=GRB.BINARY)
    u=nw_mod.addVars(V.number_of_nodes(),name='u',vtype=GRB.BINARY)
    a=nw_mod.addVars(V.number_of_nodes(),name='a',vtype=GRB.BINARY)
    b=nw_mod.addVars(V.number_of_nodes(),name='b',vtype=GRB.BINARY)
    w=nw_mod.addVars(V.number_of_nodes(),name='w',vtype=GRB.BINARY)
    
    obj_fun=sum(c*delta_ij[i,j] for i,j in V.edges())+sum(t*u[k] for k in V.nodes())
    nw_mod.setObjective(obj_fun, GRB.MINIMIZE)
    
    for i in V.nodes():
        nw_mod.addConstr(x[i] + y[i] + z[i] == 1, 'constraint_one_subd{i}')
        
    for i in V.nodes():
        nw_mod.addConstr(2*a[i]<=u[i]+x[i],'constraint_ui_xi')
        nw_mod.addConstr(2*b[i]<=u[i]+y[i],'constraint_ui_yi')
        nw_mod.addConstr(2*w[i]<=u[i]+z[i],'constraint_ui_zi')
    
        nw_mod.addConstr(1+2*a[i]>=u[i]+x[i],'constraint_ui_xi_2')
        nw_mod.addConstr(1+2*b[i]>=u[i]+y[i],'constraint_ui_yi_2')
        nw_mod.addConstr(1+2*w[i]>=u[i]+z[i],'constraint_ui_zi_2')
    
   
    nw_mod.addConstr(sum(V.nodes[i]['S'] * (x[i]-a[i]) for i in V.nodes()) <= (math.ceil(N/2)) - 1, 'constraint_sum_S_x')
    nw_mod.addConstr(sum(V.nodes[i]['S'] * (y[i]-b[i]) for i in V.nodes()) <= (math.ceil(N/2)) - 1, 'constraint_sum_S_y')
    nw_mod.addConstr(sum(V.nodes[i]['S'] * (z[i]-w[i]) for i in V.nodes()) <= (math.ceil(N/2)) - 1, 'constraint_sum_S_z')
    for i, j in V.edges():
        nw_mod.addConstr(delta_ij[i, j] <= 2 - x[i] - x[j], f'constraint_delta_ij_x_{i}_{j}_same_x')
        nw_mod.addConstr(delta_ij[j, i] <= 2 - x[i] - x[j], f'constraint_delta_ij_x_{i}_{j}_same_x_2')
        nw_mod.addConstr(delta_ij[i, j] <= 2 - y[i] - y[j], f'constraint_delta_ij_x_{i}_{j}_same_y')
        nw_mod.addConstr(delta_ij[j, i] <= 2 - y[i] - y[j], f'constraint_delta_ij_x_{i}_{j}_same_y_2')
        nw_mod.addConstr(delta_ij[i, j] <= 2 - z[i] - z[j], f'constraint_delta_ij_x_{i}_{j}_same_z')
        nw_mod.addConstr(delta_ij[j, i] <= 2 - z[i] - z[j], f'constraint_delta_ij_x_{i}_{j}_same_z_2')
        nw_mod.addConstr(delta_ij[i, j] >= x[i] - x[j], f'constraint_delta_ij_x_{i}_{j}_not_same_x')
        nw_mod.addConstr(delta_ij[i, j] >= x[j] - x[i], f'constraint_delta_ij_x_{i}_{j}_not_same_x_2')
        nw_mod.addConstr(delta_ij[i, j] >= y[i] - y[j], f'constraint_delta_ij_x_{i}_{j}_not_same_y')
        nw_mod.addConstr(delta_ij[i, j] >= y[j] - y[i], f'constraint_delta_ij_x_{i}_{j}_not_same_y_2')
        nw_mod.addConstr(delta_ij[i, j] >= z[i] - z[j], f'constraint_delta_ij_x_{i}_{j}_not_same_z')
        nw_mod.addConstr(delta_ij[i, j] >= z[j] - z[i], f'constraint_delta_ij_x_{i}_{j}_not_same_z_2')
    nw_mod.setParam('OutputFlag',False)
    nw_mod.optimize()
    print('Optimization is done. Objective function value: %.2f' % nw_mod.objVal)
    min_vals_arb.append(nw_mod.objVal)
    deltas={}
    x_arr={}
    y_arr={}
    z_arr={}
    u_arr={}
    for v in nw_mod.getVars():
        if re.match('delta',v.varName) is not None and v.x!=0:
            deltas[v.varName]=v.x
        if re.match('u',v.varName) is not None and v.x!=0:
            u_arr[v.varName]=v.x
    tmp_d=len(deltas)
    tmp_u=len(u_arr)
    u_avg.append(tmp_u)
    delta_avg.append(tmp_d)
    comb.append((tmp_u,tmp_d))

end = time.time()
print("Execution time in seconds")
print(end - start)
(uniq, freq) = (np.unique(min_vals_arb, return_counts=True))
print(np.column_stack((uniq,freq)))


plt.xlabel('Objective function value', fontsize=12)
plt.ylabel('Number of occurrences', fontsize=12)
plt.title('Atlanta Arbitrary DCN positioning with default value N='+str(N)+', c='+str(c)+',and t='+str(t))
plt.bar(uniq,freq)
plt.xlim(10, 100)
plt.savefig(f'./value_occur_eu_N{N}c{c}t{t}.png')

plt.show()
f = plt.figure()
f.set_figwidth(13)
f.set_figheight(8)
plt.clf()
width =0.5
#plt.plot(runs,min_vals_arb,'-.', label='Optimized minimum cost per run')
plt.xlabel('Iteration', fontsize=24)
plt.ylabel('Number of occurrences', fontsize=24)
plt.plot(np.arange(len(u_avg)), u_avg,label='Number of DCNs disabled for c='+str(c)+' and t='+str(t))
plt.plot(np.arange(len(delta_avg)), delta_avg,label='Number of links cut for c='+str(c)+' and t='+str(t))
plt.legend(loc='best',fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.savefig(f'./value_choice_eu_N{N}c{c}t{t}.png')

plt.show()
plt.clf()
u_zero_count = len([val for val in u_avg if val == 0])
delta_zero_count = len([val for val in delta_avg if val == 0])

# Print zero counts
print("Total number of iterations with zero disabled DCNs:", u_zero_count)
print("Total number of iterations with zero cut links:", delta_zero_count)
f = plt.figure()
f.set_figwidth(13)
f.set_figheight(8)
plt.clf()
plt.xlabel('Iteration', fontsize=12)
plt.ylabel('Objective function value', fontsize=12)
plt.plot(np.arange(len(min_vals_arb)),min_vals_arb,'-.', label='Objective value variation per iteration')
plt.legend(loc='best',fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.savefig(f'./value_fluc_eu_N{N}c{c}t{t}.png')
plt.clf()
plt.show()
#print(comb)
count = Counter(comb)
print("Combination of attack occurrences reported as (link cuts, nodes disrupted): "+str(count))