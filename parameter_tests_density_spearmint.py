
# coding: utf-8

# In[117]:

import os
from shutil import copy, copyfile
import subprocess
from save_embedded_graph27 import main_binary as embed_main
from spearmint_ghsom import main as ghsom_main
import numpy as np
import pickle
from time import time
from __future__ import division

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
def fitness(x, gml, labels):
    
    params = {'w': 0.0001,
         'eta': 0.0001,
         'sigma': 1,
          'e_sg': x,
         'e_en': 0.8}
    
    return 1 - ghsom_main(params, gml, labels, 1000)[0][0]

def neighbour(x, p=0.1):
    return x + 2 * (np.random.rand() - 0.5) * p

def simulated_annealing(gml_filename, labels, initial_x=0.8, 
                        num_iter=50, start_temp=0.01, epsilon=0.05, no_improv_thres=10): 

    #load_progress
    progress_file = 'sa_progress'
    if os.path.isfile('{}.pkl'.format(progress_file)):
        start, xb, f_xb, temp, no_improvements = load_obj(progress_file)
        print 'loading simulated annealing progress'
    else:
        #new progress
        start = 0
        xb = initial_x
        f_xb = fitness(xb, gml_filename, labels)
        temp = start_temp
        no_improvements = 0
    print 'starting from iteration {}'.format(start)

    for i in range(start, num_iter):
        
        if no_improvements >= no_improv_thres:
            print 'no improvement for {} iterations, stopping'.format(no_improv_thres)
            break
        
        if f_xb < epsilon:
            print 'function minimised to tolerance of {}, stopping'.format(epsilon)
            break

        x = neighbour(xb)

        print 'trying x={}'.format(x)

        temp = 0.99 * temp

        f_x = fitness(x, gml_filename, labels)
        
        p = np.exp((f_xb - f_x) / temp)
        
        print 'f_xb={}, f_x={}, p={}'.format(f_xb, f_x, p)

        if f_x < f_xb or np.random.rand() < p:
            xb = x
            f_xb = f_x
            
            no_improvements = 0
        else:
            no_improvements += 1

        #save progress
        save_obj((i + 1, xb, f_xb, temp, no_improvements), progress_file)
                
        print 'epoch={}, xb={}, f_xb={}'.format(i, xb, f_xb)    
        print

    return xb

def graph_measures(gml_filename):
    
    G = nx.read_gml(gml_filename)
    
    density = nx.density(G)
    
    arr = np.array([val for key, val in nx.betweenness_centrality(G).iteritems()])
    arr = np.sort(arr)[::-1]
    
    k = len(arr)

    k1 = k

    var = 1

    while var > 0.95:

        k -= 1

        var = np.sum(arr[:k]) / sum 
    
    centrality = k / k1
    
    return density, centrality


# In[119]:

#root dir
os.chdir("C:\Miniconda3\Jupyter\GHSOM_simplex_dsd")

#save directory
dir = os.path.abspath("parameter_tests_density_spearmint")

#number of networks to generate
num_networks = 200

#number of times to repeat
num_jobs = 10

#number of nodes in the graph
N = 64

#make save directory
if not os.path.isdir(dir):
    os.mkdir(dir)

#change to dir
os.chdir(dir)    

#network file names -- output of network generator
network = "network.dat"
first_level = "community.dat"

#community labels
labels = 'firstlevelcommunity'

#mixing factors
mu = 0.3

densities = np.zeros(num_networks)
centralities = np.zeros(num_networks)
best_e_sgs = np.zeros(num_networks)

spearmint_dir = 'C:\Miniconda3\Jupyter\GHSOM_simplex_dsd\spearmint\spearmint'

for i in range(num_networks):
    
    #create directory
    dir_string = os.path.join(dir, str(i))
    if not os.path.isdir(dir_string):
        os.mkdir(dir_string)
    
    #change working directory    
    os.chdir(dir_string)
    
    #make benchmark parameter file
    filename = "benchmark_flags_{}.dat".format(i)

    if os.path.isfile('density.txt'):
        
        density = np.genfromtxt('density.txt')
        
    else:
    
        #number of edges
        num_edges = np.random.randint(256, 2017 * 0.8)

        #number of communities
        num_communities = np.random.randint(1, 5)

        #number of nodes in micro community
        minc = 1
        maxc = np.random.randint(16, 30)

        #average number of edges
        k = float(num_edges) / N

        #max number of edges
        maxk = 2 * k

        ##calculate density
        density = 2 * float(num_edges) / (N * (N-1))
        with open('density.txt','w') as f:
            f.write('{}\n'.format(density))

        if not os.path.isfile(filename):
            print 'density: {}'.format(density)
            print '-N {} -k {} -maxk {} -minc {} -maxc {} -mu {}'.format(N, k, maxk, minc, maxc, mu)
            with open(filename,"w") as f:
                f.write("-N {} -k {} -maxk {} -minc {} -maxc {} -mu {}".format(N, k, maxk, minc, maxc, mu))
            print 'written flag file: {}'.format(filename)
     
    print 'density of random network {}: {}'.format(i, density)
#     densities[i] = density
    
    #copy executable
    ex = "benchmark.exe"   
    if not os.path.isfile(ex):

        source = "C:\\Users\\davem\\Documents\\PhD\\Benchmark Graph Generators\\binary_networks\\benchmark.exe"
        copyfile(source, ex)
            
    #cmd strings
    change_dir_cmd = "cd {}".format(dir_string)
    generate_network_cmd = "benchmark -f {}".format(filename)

    #output of cmd
    output_file = open("cmd_output.out", 'w')

    #generate network and rename
    if not os.path.isfile(network):

        process = subprocess.Popen(change_dir_cmd + " && " + generate_network_cmd, 
                                stdout=output_file, 
                                stderr=output_file, 
                                shell=True)
        process.wait()

        print 'generated network {}'.format(i)

        
    gml_filename = 'embedded_network_{}.gml'.format(i)    
    
    #embed into gml file
    if not os.path.isfile(gml_filename):

        ##embed graph
        embed_main(network, first_level, gml_filename)

        print 'embedded network {} as {} in {}'.format(i, gml_filename, os.getcwd())
        
    densities[i], centralities[i] = graph_measures(gml_filename)
    print 'calculated graph measures'
    
    e_sg_file = 'best_e_sg.txt'
    
    if not os.path.isfile(e_sg_file):
        
        print 'using simulated annealing to optimise e_sg'
        
        best_e_sg = simulated_annealing(gml_filename, labels)
        
        print 'best setting: {}'.format(best_e_sg)
        
        with open(e_sg_file,'w') as f:
            f.write('{}\n'.format(best_e_sg))
    
    else:

        best_e_sg = np.genfromtxt(e_sg_file)
        
        print 'loading best setting: {}'.format(best_e_sg)
    
    best_e_sgs[i] = best_e_sg
    
    
print 'DONE'


# In[126]:

for i in range(num_networks):
    
    print i
    print densities[i]
    print best_e_sgs[i]
    print


# In[123]:

def trendline(xd, yd, order=1, c='r', alpha=1, Rval=False):
    """Make a line of best fit"""

    #Calculate trendline
    coeffs = np.polyfit(xd, yd, order)

    intercept = coeffs[-1]
    slope = coeffs[-2]
    power = coeffs[0] if order == 2 else 0

    minxd = np.min(xd)
    maxxd = np.max(xd)

    xl = np.array([minxd, maxxd])
    yl = power * xl ** 2 + slope * xl + intercept

    #Plot trendline
    plt.plot(xl, yl, c, alpha=alpha)

    #Calculate R Squared
    p = np.poly1d(coeffs)

    ybar = np.sum(yd) / len(yd)
    ssreg = np.sum((p(xd) - ybar) ** 2)
    sstot = np.sum((yd - ybar) ** 2)
    Rsqr = ssreg / sstot

    if not Rval:
        #Plot R^2 value
        plt.text(0.8 * maxxd + 0.2 * minxd, 0.8 * np.max(yd) + 0.2 * np.min(yd),
                 '$R^2 = %0.2f$' % Rsqr)
    else:
        #Return the R^2 value:
        return Rsqr

import matplotlib.pyplot as plt

plt.plot(densities, best_e_sgs, 'x')
trendline(densities, best_e_sgs)
plt.xlabel('density')
plt.ylabel('best $\epsilon_{sg}$')
plt.show()


# In[67]:

import scipy
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(densities, best_e_sgs)
print 'slope={}'.format(slope)
print 'intercept={}'.format(intercept)
print 'r_value={}'.format(r_value)
print 'p_value={}'.format(p_value)
print 'std_err={}'.format(std_err)


# In[124]:

##y = mx + c
#e_sg = 0.3333 density + 0.59
import networkx as nx
import numpy as np

os.chdir("C:\Miniconda3\Jupyter\GHSOM_simplex_dsd")

G = nx.read_gml('embedded_polbooks.gml')



# In[125]:

from  __future__ import division

# print nx.density(G)
# print nx.degree_centrality(G)
arr = np.array([val for key, val in nx.betweenness_centrality(G).iteritems()])
arr = np.sort(arr)[::-1]
print arr
sum = np.sum(arr)
# print np.sort(arr)[::-1]
# print np.mean(arr)

k = len(arr)

k1 = k
print 'k1={}'.format(k1) 

var = 1

while var > 0.95:
    
    k -= 1
    
    var = np.sum(arr[:k]) / sum

print 'k={}'.format(k)    

k2 = k / k1    

print 'k2={}'.format(k2)


# In[ ]:

def graph_measures(G):
    
    density = nx.density(G)
    
    arr = np.array([val for key, val in nx.betweenness_centrality(G).iteritems()])
    arr = np.sort(arr)[::-1]
    
    k = len(arr)

    k1 = k

    var = 1

    while var > 0.95:

        k -= 1

        var = np.sum(arr[:k]) / sum 
    
    centrality = k / k1
    
    return density, centrality

