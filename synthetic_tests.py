
# coding: utf-8

# In[1]:

import os
from shutil import copyfile
import subprocess
from save_embedded_graph27 import main_hierarchical as embed_main
from spearmint_ghsom import main as ghsom_main
import numpy as np
import pickle
from time import time

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

#root dir
os.chdir("C:\Miniconda3\Jupyter\GHSOM_simplex_dsd")

#save directory
dir = os.path.abspath("synthetic_benchmarks")

#number of times to repeat
num_repeats = 100

#number of micro communities
k1 = 16
#number of macro communities
k2 = 4
#number of nodes in micro communitiy
s1 = 32
#number of nodes in macro community = s1 * s2
s2 = 4
#number of nodes in the network
N = s1 * s2 * k2
#number of links to same micro community
z1 = 16
#number of links to same macro community
z2 = 16
#nuber of nodes in micro community
minc = s1
maxc = s1
#number of nodes in macro community
minC = s1 * s2
maxC = s1 * s2

#make save directory
if not os.path.isdir(dir):
    os.mkdir(dir)

#change to dir
os.chdir(dir)    

#network file names -- output of network generator
network = "network.dat"
first_level = "community_first_level.dat"
second_level = "community_second_level.dat"

#community labels
labels = 'firstlevelcommunity,secondlevelcommunity'

#ghsom parameters
params = {'w': 0.0001,
         'eta': 0.0001,
         'sigma': 1,
         'e_sg': 0.8,
         'e_en': 0.8}

mixing_factors = [16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36]

overall_nmi_scores = np.zeros((len(mixing_factors), num_repeats, 2))

for i in range(len(mixing_factors)):
    
    z3 = mixing_factors[i]
    
    #node degree
    k = z1 + z2 + z3
    maxk = k
    
    #mixing factors
    mu1 = float(z3) / k
    mu2 = float(z2) / k 
    
    #create directory
    dir_string = os.path.join(dir, str(z3))
    if not os.path.isdir(dir_string):
        os.mkdir(dir_string)
    
    #change working directory    
    os.chdir(dir_string)
    
    if os.path.isfile('nmi_scores.csv'):
        print 'already completed {}, loading nmi scores and continuing'.format(z3)
        nmi_scores = np.genfromtxt('nmi_scores.csv', delimiter=',')
        overall_nmi_scores[i] = nmi_scores
        continue
    
    #copy executable
    ex = "hbenchmark.exe"   
    if not os.path.isfile(ex):
        
        source = "C:\Users\davem\Documents\PhD\Benchmark Graph Generators\hierarchical_bench2_2\hbenchmark.exe"
        copyfile(source, ex)
        
    #make benchmark parameter file
    filename = "benchmark_flags_{}.dat".format(z3)
    if not os.path.isfile(filename):
        with open(filename,"w") as f:
            f.write("-N {} -k {} -maxk {} -minc {} -maxc {} -minC {} -maxC {} -mu1 {} -mu2 {}".format(N, k, maxk, minc, maxc, minC, maxC, mu1, mu2))
    
    #cmd strings
    change_dir_cmd = "cd {}".format(dir_string)
    generate_network_cmd = "hbenchmark -f {}".format(filename)
    
    #output of cmd
    output_file = open("cmd_output.out", 'w')
    
    #record NMI scores
    if not os.path.isfile('nmi_scores.pkl'):
        print 'creating new nmi scores array'
        nmi_scores = np.zeros((num_repeats, len(labels.split(','))))
    else:
        print 'loading nmi score progress'
        nmi_scores = load_obj('nmi_scores')
        
    #record running times
    if not os.path.isfile('running_times.pkl'):
        print 'creating new running time array'
        running_times = np.zeros(num_repeats)
    else:
        print 'loading running time progress'
        running_times = load_obj('running_times')
    
    #generate networks
    for r in range(1, num_repeats+1):
        
        network_rename = "{}_{}".format(r,network)
        first_level_rename = "{}_{}".format(r,first_level)
        second_level_rename = "{}_{}".format(r,second_level)
        gml_filename = 'embedded_network_{}.gml'.format(r)
        
        if not os.path.isfile(network_rename):
        
            process = subprocess.Popen(change_dir_cmd + " && " + generate_network_cmd, 
                                    stdout=output_file, 
                                    stderr=output_file, 
                                    shell=True)
            process.wait()

            os.rename(network, network_rename)
            os.rename(first_level, first_level_rename)
            os.rename(second_level, second_level_rename)
            
        if not os.path.isfile(gml_filename):
            
            ##embed graph
            embed_main(network_rename, first_level_rename, second_level_rename)
            
            print 'created {} in {}'.format(gml_filename, os.getcwd())
            
        ##score for this network
        if not np.all(nmi_scores[r-1]):
            
            start_time = time()
            
            print 'starting ghsom for: {}/{}'.format(z3, gml_filename)
            nmi_score, communities_detected = ghsom_main(params, gml_filename, labels, 10000)
            nmi_scores[r-1] = nmi_score
            
            running_time = time() - start_time
            print 'running time of algorithm: {}'.format(running_time)
            running_times[r-1] = running_time
            
            #save
            save_obj(nmi_scores, 'nmi_scores')
            save_obj(running_times, 'running_times')
            
            print 'saved nmi score for network {}: {}'.format(gml_filename, nmi_score)
            print
            
    ##output nmi scores to csv file
    print 'writing nmi scores and running times to file'
    np.savetxt('nmi_scores.csv',nmi_scores,delimiter=',')
    np.savetxt('running_times.csv',running_times,delimiter=',')
    
print 'DONE'

print 'OVERALL NMI SCORES'
print overall_nmi_scores


# In[13]:

import matplotlib.pyplot as plt
import numpy as np


os.chdir("C:\Miniconda3\Jupyter\GHSOM_simplex_dsd")

first_level = np.genfromtxt('first_level.csv',delimiter=',')
first_level = first_level[1:]
first_level[:,0] = np.rint(first_level[:,0])

first_level_m = np.zeros((len(first_level) / 3, 4))
for i in range(len(first_level) / 3):
    first_level_m[i] = np.mean(first_level[3*i:3*i+2], axis=0)
first_level_m[first_level_m > 1] = 1

second_level = np.genfromtxt('second_level.csv',delimiter=',')
second_level = second_level[1:]
second_level[:,0] = np.rint(second_level[:,0])

second_level_m = np.zeros((len(second_level) / 3, 4))
for i in range(len(second_level) / 3):
    second_level_m[i] = np.mean(second_level[3*i:3*i+2], axis=0)
second_level_m[second_level_m > 1] = 1

means = np.zeros((len(mixing_factors), 2))
ses = np.zeros((len(mixing_factors), 2))

for i in range(len(mixing_factors)):
# for score in overall_nmi_scores:
    score = overall_nmi_scores[i]
    m = np.mean(score, axis=0)
    means[i] = m
    print m
    sd = np.std(score, axis=0)
    se = sd / np.sqrt(num_repeats)
    ses[i] = se
    print se
    print

plt.errorbar(mixing_factors, means[:, 0], yerr=ses[:, 0], fmt='-o')
plt.errorbar(mixing_factors, means[:, 1], yerr=ses[:, 1], fmt='-o')
plt.axis([16, 36, 0.95, 1])
plt.legend(['First level','Second level'], loc=3)
plt.xlabel('Mixing parameter $z_3$')
plt.ylabel('Normalized mutual information')

plt.show()

##first level
plt.errorbar(mixing_factors, means[:, 0], yerr=ses[:, 0], fmt='-o')
plt.plot(mixing_factors, first_level_m[:, 1], '-rx')
plt.plot(mixing_factors, first_level_m[:, 2], '-gx')
plt.plot(mixing_factors, first_level_m[:, 3], '-bx')
plt.axis([16, 36, 0.4, 1])
plt.legend(['FM', 'FUC', 'PMC', 'GHSOM'], loc=3)
plt.xlabel('Mixing parameter $z_3$')
plt.ylabel('Normalized mutual information')
plt.title('First level')
# plt.title('mixing factor vs. NMI score for both levels of community')

plt.show()

#second level
plt.errorbar(mixing_factors, means[:, 1], yerr=ses[:, 0], fmt='-o')
plt.plot(mixing_factors, second_level_m[:, 1], '-rx')
plt.plot(mixing_factors, second_level_m[:, 2], '-gx')
plt.plot(mixing_factors, second_level_m[:, 3], '-bx')
plt.axis([16, 36, 0.4, 1])
plt.legend(['FM', 'FUC', 'PMC', 'GHSOM'], loc=3)
plt.xlabel('Mixing parameter $z_3$')
plt.ylabel('Normalized mutual information')
plt.title('Second level')



plt.show()


# In[17]:

import numpy as np

first_level = np.genfromtxt('first_level.csv',delimiter=',')
first_level = first_level[1:]
first_level[:,0] = np.rint(first_level[:,0])

print first_level

first_level_m = np.zeros((len(first_level) / 3, 4))
for i in range(len(first_level) / 3):
    first_level_m[i] = np.mean(first_level[3*i:3*i+2], axis=0)
first_level_m[first_level_m > 1] = 1
print first_level_m


# In[22]:

for score in overall_nmi_scores:
    m = np.mean(score, axis=0)
    print m
    sd = np.std(score, axis=0)
    se = sd / np.sqrt(100)
    print se
    print


# In[2]:

os.chdir("C:\Miniconda3\Jupyter\GHSOM_simplex_dsd")
gml_filename = 'embedded_network_69.gml'
nmi_scores, num_communities = ghsom_main(params, gml_filename, labels, 10000)
print nmi_scores


# In[ ]:



