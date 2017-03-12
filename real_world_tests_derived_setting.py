
# coding: utf-8

# In[3]:

import os
from shutil import copyfile
import subprocess
from spearmint_ghsom import main as ghsom_main
import numpy as np
import pickle
from time import time
import networkx as nx

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

#root dir
os.chdir("C:\Miniconda3\Jupyter\GHSOM_simplex_dsd")

#save directory
dir = os.path.abspath("real_world_benchmarks_derived")

#number of times to repeat
num_repeats = 100

#make save directory
if not os.path.isdir(dir):
    os.mkdir(dir)

#change to dir
os.chdir(dir)    



#network names
network_names = ['karate','dolphin','polbooks','football']

#community labels
labels = ['club','group','value','value']

overall_nmi_scores = np.zeros((len(network_names), num_repeats))
overall_communities_detected = np.zeros((len(network_names), num_repeats))

for i in range(len(network_names)):
    
    #name of current network
    network_name = network_names[i]

    #label of current network
    label = labels[i]
    
    #create directory
    dir_string = os.path.join(dir, network_name)
    if not os.path.isdir(dir_string):
        os.mkdir(dir_string)
    
    #change working directory    
    os.chdir(dir_string)
    
    gml_filename = 'embedded_{}.gml'.format(network_name)  
    if not os.path.isfile(gml_filename):
        
        source = "C:\Miniconda3\Jupyter\GHSOM_simplex_dsd\{}".format(gml_filename)
        copyfile(source, gml_filename)
    ##calculate density and derive parameter setting
    
    #load graph and calculate density
    G = nx.read_gml(gml_filename)
    density = nx.density(G)
    
    #derive parameter setting -- from scipy
    e_sg = 0.377746404462 * density + 0.590217653032
    
    print 'density of network={}'.format(density)
    print 'e_sg={}'.format(e_sg)
    
    if os.path.isfile('nmi_scores.csv'):
        print 'already completed {} network, loading nmi scores and continuing'.format(network_name)
        nmi_scores = np.genfromtxt('nmi_scores.csv', delimiter=',')
        overall_nmi_scores[i] = nmi_scores
        communities_detected = np.genfromtxt('communties_detected.csv', delimiter=',')
        overall_communities_detected[i] = communities_detected
        continue
    
    #record NMI scores
    if not os.path.isfile('nmi_scores.pkl'):
        print 'creating new nmi scores array'
        nmi_scores = np.zeros(num_repeats)
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
        
    #record communities detected
    if not os.path.isfile('communities_detected.pkl'):
        print 'creating new communites detected array'
        communities_detected = np.zeros(num_repeats)
    else:
        print 'loading communites detected progress'
        communities_detected = load_obj('communities_detected')
        
    #copy embedded gml
    gml_filename = 'embedded_{}.gml'.format(network_name)  
    if not os.path.isfile(gml_filename):
        
        source = "C:\Miniconda3\Jupyter\GHSOM_simplex_dsd\{}".format(gml_filename)
        copyfile(source, gml_filename)
    ##calculate density and derive parameter setting
    
    #load graph and calculate density
    G = nx.read_gml(gml_filename)
    density = nx.density(G)
    
    #derive parameter setting -- from scipy
    e_sg = 0.377746404462 * density + 0.590217653032
    
    print 'density of network={}'.format(density)
    print 'e_sg={}'.format(e_sg)
    
    #ghsom parameters
    params = {'w': 0.0001,
         'eta': 0.0001,
         'sigma': 1,
         'e_sg': e_sg,
         'e_en': 0.8}
    
    #generate networks
    for r in range(1,num_repeats+1):
            
        ##score for this network
        if not np.all(nmi_scores[r-1]):
            
            start_time = time()
            
            print 'starting ghsom for: {}, repeat: {}'.format(gml_filename, r)
            nmi_score, comm_det = ghsom_main(params, gml_filename, label, 10000)
            nmi_scores[r-1] = nmi_score
            communities_detected[r-1] = comm_det
            
            running_time = time() - start_time
            print 'running time of algorithm: {}'.format(running_time)
            running_times[r-1] = running_time
            
            #save
            save_obj(nmi_scores, 'nmi_scores')
            save_obj(running_times, 'running_times')
            save_obj(communities_detected, 'communities_detected')
            
            print 'saved nmi score for network {}: {}'.format(gml_filename, nmi_score)
            print 'saved communities detected for network {}: {}'.format(gml_filename, comm_det)
            print
            
    ##output nmi scores to csv file
    print 'writing nmi scores and running times to file'
    np.savetxt('nmi_scores.csv',nmi_scores,delimiter=',')
    np.savetxt('running_times.csv',running_times,delimiter=',')
    np.savetxt('communties_detected.csv',communities_detected,delimiter=',')
    
    overall_nmi_scores[i] = nmi_scores
    overall_communities_detected[i] = communities_detected
    
print 'DONE'

print 'OVERALL NMI SCORES'
print overall_nmi_scores
print overall_communities_detected


# In[2]:

for score in overall_nmi_scores:
    
    mean = np.mean(score)
    print mean
    se = np.std(score) / np.sqrt(num_repeats)
    print se
    print


# In[3]:

import networkx as nx

G = nx.read_gml('embedded_karate.gml')

print nx.density(G)
print 2.0 * nx.number_of_edges(G) / (nx.number_of_nodes(G) * (nx.number_of_nodes(G) - 1))


# In[ ]:



