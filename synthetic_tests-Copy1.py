
# coding: utf-8

# In[ ]:

import os
from shutil import copyfile
import subprocess
from save_embedded_graph27 import main as embed_main
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

for z3 in [16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36]:
    
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
    
    #copy executable
    ex = "hbenchmark.exe"   
    if not os.path.isfile(ex):
        
        source = "C:\Users\davem\Documents\PhD\Benchmark Graph Generators\hierarchical_bench2_2\hbenchmark.exe"
        copyfile(source, ex)
        
    #make benchmark parameter file
    filename = "benchmark_flags_{}.dat".format(z3)
    if not os.path.isfile(filename):
        with open(filename,"w") as f:
            f.write("-N {} -k {} -maxk {} -minc {} -maxc {} -minC {} -maxC {} -mu1 {} -mu2 {}".format(N, iik, maxk, minc, maxc, minC, maxC, mu1, mu2))
    
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
            nmi_score = ghsom_main(params, gml_filename, labels)
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


# In[ ]:



