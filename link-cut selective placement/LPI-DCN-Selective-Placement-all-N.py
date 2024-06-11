# -*- coding: utf-8 -*-
"""
Created on Fri May 10 15:41:12 2024

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
import numpy as np

graph=nx.Graph()
V=nx.read_gml('./us-ca.gml',destringizer=int)
N=3
print('This are the node locations for the selected topology:')
print(V.nodes(data='name'))
print('\\')
print('This are the existing links for the selected topology:')
print(V.edges(data='name'))
node_positions = {}

for node, attributes in V.nodes(data=True):
    
    graphics = attributes.get('graphics', {})
    x = graphics.get('x', 0)
    y = graphics.get('y', 0)
    
    node_positions[node] = (x, y)


end=0
if V.number_of_nodes()%2==0:
    end=V.number_of_nodes()-1
else:
    end=V.number_of_nodes()

reads={}
while (N<=end):
    meth=0
    nw_mod=Model(name="nw")
    tmp1={}
    while(meth<4):
        nw_mod.reset(0) 
        


        C=100
        c=10
                
                
        delta_ij = nw_mod.addVars(V.number_of_nodes(),V.number_of_nodes(),name='delta_{i}_{j}',vtype=GRB.BINARY)
        x=nw_mod.addVars(V.number_of_nodes(),name='x',vtype=GRB.BINARY)
        y=nw_mod.addVars(V.number_of_nodes(),name='y',vtype=GRB.BINARY)
        z=nw_mod.addVars(V.number_of_nodes(),name='z',vtype=GRB.BINARY)
                
        obj_fun=sum(c*delta_ij[i,j] for i,j in V.edges())
        nw_mod.setObjective(obj_fun, GRB.MINIMIZE)
        d_central=nx.degree_centrality(V)

        top_dc={}
        min_dc={}
        for key in sorted(d_central, key=d_central.get):
            min_dc[key]=d_central[key]
        for key in sorted(d_central, key=d_central.get, reverse=True):
            top_dc[key]=d_central[key]

        c_central=nx.closeness_centrality(V)

        top_cc={}
        min_cc={}
        for key in sorted(c_central, key=c_central.get):
            min_cc[key]=c_central[key]
        for key in sorted(c_central, key=c_central.get, reverse=True):
            top_cc[key]=c_central[key]

        b_central=nx.degree_centrality(V)

            

        top_degree=[]
        top_close=[]
        #top_between=[]
        min_degree=[]
        min_close=[]
        #min_between=[]

        S_likely_deg=[]
        S_likely_deg_low=[]
        for key in  top_dc.keys():
            top_degree.append(key)
        for key in  min_dc.keys():
            min_degree.append(key)



        for i in range(0,N):
            S_likely_deg.append(top_degree[i])
        for i in range(0,N):
            S_likely_deg_low.append(min_degree[i])
            

        print("Nodes with higher degree centrality")
        print(S_likely_deg)
        print("Nodes with lowest degree centrality")
        print(S_likely_deg_low)

        S_likely_cls=[]
        S_likely_cls_low=[]
        for key in  top_cc.keys():
            top_close.append(key)
        for key in  min_cc.keys():
            min_close.append(key)
            

        for i in range(0,N):
            S_likely_cls.append(top_close[i])
        for i in range(0,N):
            S_likely_cls_low.append(min_close[i])
        print("Nodes with higher closeness centrality")
        print(S_likely_cls)
        print("Nodes with lowest closeness centrality")
        print(S_likely_cls_low)
        
        if meth==0:
            print("Allocation method chosen: High degree centrality")
            print("Nodes chosen to be DCN:")
            meth_name='HDC'
            for i in V.nodes():
                
                if i in S_likely_deg:
                    V.nodes[i]['S']=1
                    print(V.nodes[i]["name"])
                else:
                    V.nodes[i]['S']=0
            print(V.nodes(data='S'))
            
        elif meth==1:
            meth_name='HCC'
            print("Allocation method chosen: High closeness centrality")
            print("Nodes chosen to be DCN:")
            for i in V.nodes():
                if i in S_likely_cls:
                    V.nodes[i]['S']=1
                    print(V.nodes[i]["name"])
                else:
                    V.nodes[i]['S']=0
            print(V.nodes(data='S'))
            
        elif meth==2:
            meth_name='LDC'
            print("Allocation method chosen: Low degree centrality")
            print("Nodes chosen to be DCN:")
            for i in V.nodes():
                if i in S_likely_deg_low:
                    V.nodes[i]['S']=1
                    print(V.nodes[i]["name"])
                else:
                    V.nodes[i]['S']=0
            print(V.nodes(data='S'))
        
        elif meth==3:
            meth_name='LCC'
            print("Allocation method chosen: Low closeness centrality")
            print("Nodes chosen to be DCN:")
            for i in V.nodes():
                if i in S_likely_cls_low:
                    V.nodes[i]['S']=1
                    print(V.nodes[i]["name"])
                else:
                    V.nodes[i]['S']=0
            print(V.nodes(data='S'))

        nw_mod.update()
        for i in V.nodes():
            nw_mod.addConstr(x[i] + y[i] + z[i] == 1, 'constraint_one_subd{i}')
            
        nw_mod.addConstr(sum(V.nodes[i]['S'] * x[i] for i in V.nodes()) <= (math.ceil(N/2)) - 1, 'constraint_sum_S_x')
        nw_mod.addConstr(sum(V.nodes[i]['S'] * y[i] for i in V.nodes()) <= (math.ceil(N/2)) - 1, 'constraint_sum_S_y')
        nw_mod.addConstr(sum(V.nodes[i]['S'] * z[i] for i in V.nodes()) <= (math.ceil(N/2)) - 1, 'constraint_sum_S_z')
        
        nw_mod.update()
        
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
        nw_mod.update()
        nw_mod.optimize()
        
        print('Optimization is done. Objective function value: %.2f' % nw_mod.objVal)
        #min_vals_arb.append(nw_mod.objVal)
        deltas={}

        for v in nw_mod.getVars():
            if re.match('delta',v.varName) is not None and v.x!=0:
                deltas[v.varName]=v.x
        tmp_d=len(deltas)

        tmp1[meth_name]=nw_mod.objVal
        
        meth+=1
    reads[N]=tmp1
    N+=2
print(reads)
N_values = list(reads.keys())
methods = list(reads[N_values[0]].keys())
f = plt.figure()
f.set_figwidth(15)
f.set_figheight(10)
plt.clf()
# Prepare the data for plotting
width = 0.2  # Width of each bar
num_methods = len(methods)
x_indexes = np.arange(len(N_values))

fig, ax = plt.subplots()

# Plotting each method as a group of bars
for i, method in enumerate(methods):
    # Calculate the bar positions
    bar_positions = x_indexes + i * width
    # Get the values for each method across different N values
    method_values = [reads[N][method] for N in N_values]
    # Plot the bars
    ax.bar(bar_positions, method_values, width=width, label=method)

#avg_cost_rand=[44.06,44.52,44.56,44.3,43.3,42.66,40]
#avg_cost_rand=[48.08,50.32,50.35,52.01,50.47,51.01,50.25,50.42,50.48,50.24,50,50,50]
avg_cost_rand=[51.41,50.76,53.46,53.022,49.688,54.152,60.432,64.004,61.596,62.714,62.62,61.26,60.666,60.2,60.37,60.224,60,60,60]
ax.plot(x_indexes, avg_cost_rand, color='k', linestyle='--', marker='o', markersize=5, label='Avg cost for arbitrary placement for N')
# Customize the graph
ax.set_xlabel('DCN density')
ax.set_ylabel('Optimized minimum cost')
ax.set_title('US-CA graph behaviour for N=3...|V| with c=10 and selective DCN placement')
ax.set_xticks(x_indexes + width * (num_methods - 1) / 2)
ax.set_xticklabels(N_values)
ax.legend()

# Display the graph
plt.show()