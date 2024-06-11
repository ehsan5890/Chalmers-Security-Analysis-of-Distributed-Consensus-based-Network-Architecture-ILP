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
V=nx.read_gml('./atlanta.gml',destringizer=int)
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


c=10
t=30
u_avg=[]
delta_avg=[]

reads={}
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

        #nw_mod.update()
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
        
        #nw_mod.update()
        
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
        #nw_mod.update()
        nw_mod.optimize()
        
        print('Optimization is done. Objective function value: %.2f' % nw_mod.objVal)
        #min_vals_arb.append(nw_mod.objVal)


        deltas={}
        u_arr={}
        for v in nw_mod.getVars():
            if re.match('delta',v.varName) is not None and v.x!=0:
                #print ('%s: %g' % (v.varName,v.x))
                deltas[v.varName]=v.x
                
            if re.match('u',v.varName) is not None and v.x!=0:
                u_arr[v.varName]=v.x
                #print ('%s: %g' % (v.varName,v.x))
        #print(u_arr)
        #print(deltas)
        tmp_d=len(deltas)
        tmp_u=len(u_arr)
        tmp1[meth_name]=nw_mod.objVal
        u_avg.append(tmp_u)
        delta_avg.append(tmp_d)
        meth+=1
    reads[N]=tmp1
    N+=2
print(reads)
print(u_avg)
print(delta_avg)
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

# Customize the graph
#avg_cost_rand=[20,30,39.414,40.0,40.0,40.0,40.0]
#avg_cost_rand=[30,41.01398601398601,42.056,42.095,41.69230769230769,41.333333333333336,40.0]
avg_cost_rand=[44.08791208791209,44.52547452547453,44.702,44.19,43.38461538461539,42.666666666666664,40.0]
ax.plot(x_indexes, avg_cost_rand, color='k', linestyle='--', marker='o', markersize=5, label='Avg cost for arbitrary placement for N')
ax.set_xlabel('DCN density')
ax.set_ylabel('Optimized minimum cost')
ax.set_title('Atlanta graph behaviour for N=5...13 with c='+str(c)+', t='+str(t)+' and selective DCN placement')
ax.set_xticks(x_indexes + width * (num_methods - 1) / 2)
ax.set_xticklabels(N_values)
#for c=t=10
#u_vals_sel=[2,3,4,0.75,0,0,0]
#d_vals_sel=[0,0,0,3.25,4,4,4]
#d_vals_rand=[0,0.011988011988011988,0.1248,3.0102,3.38021978021978,3.5047619047619047,4]
#u_vals_rand=[2,2.988011988011988,3.8166,0.9898,0.6197802197802198,0.49523809523809526,0]
#for c=10 t=15
#u_vals_sel=[2,2.25,0.5,0.5,0,0,0]
#d_vals_sel=[0,1,3.5,3.5,4,4,4]
#d_vals_rand=[0,1.7512487512487513,3.413,3.581,3.6615384615384614,3.7333333333333334,4]
#u_vals_rand=[2,1.5667665667665667,0.5284,0.419,0.3384615384615385,0.26666666666666666,0]
#for c=10 t=30
u_vals_sel=[0,0,0,0,0,0,0]
d_vals_sel=[4.5,4.75,4.5,4.5,4,4,4]
d_vals_rand=[3.940659340659341,4.442557442557442,4.469,4.419,4.338461538461538,4.266666666666667,4]
u_vals_rand=[0.15604395604395604,0.00333000333000333,0.0004,0,0,0,0]
# Display the graph
ax.legend()
plt.show()
plt.xlabel('DCN Density')
plt.ylabel('Quantity')
plt.plot(N_values, u_vals_sel,marker='o',label='Average number of DCNs disabled for c='+str(c)+' and t='+str(t)+' with selective method')
plt.plot(N_values, u_vals_rand,marker='s',label='Average number of DCNs disabled for c='+str(c)+' and t='+str(t)+' with arbitrary method')
#plt.legend(loc='best')
#plt.show()


#ax.legend()
#plt.show()
#plt.xlabel('DCN Density')
#plt.ylabel('Quantity')
plt.plot(N_values, d_vals_sel,marker='o',label='Average number of links cut for c='+str(c)+' and t='+str(t)+' with selective method')
plt.plot(N_values, d_vals_rand,marker='s',label='Average number of links_cut for c='+str(c)+' and t='+str(t)+' with arbitrary method')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()