
# coding: utf-8

# In[ ]:

import os
from shutil import copyfile
import subprocess
from save_embedded_graph27 import main_binary as embed_main
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
dir = os.path.abspath("parameter_tests")

#number of times to repeat
num_repeats = 10

#number of nodes in communitiy
s1 = 32

#number of links to same community
z1 = 16

#number of nodes in micro community
minc = s1
maxc = s1

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

#mixing parameter
z2 = 16

#node degree
k = z1 + z2
maxk = k

#mixing factors
mu = float(z2) / k 

num_communities = [3, 4, 5, 6]
parameter_settings = [0.5, 0.6, 0.7, 0.8, 0.9, 1][::-1]

overall_nmi_scores = np.zeros((len(num_communities), len(parameter_settings)))

for i in range(len(num_communities)):
# for k1 in num_communities:
    k1 = num_communities[i]
    
    #number of nodes in the network
    N = k1 * s1
    
    #create directory
    dir_string = os.path.join(dir, str(k1))
    if not os.path.isdir(dir_string):
        os.mkdir(dir_string)
    
    #change working directory    
    os.chdir(dir_string)
    
    for j in range(len(parameter_settings)):
#     for p in parameter_settings:
        
        p = parameter_settings[j]
        
        #ghsom parameters
        params = {'w': 0.0001,
                 'eta': 0.0001,
                 'sigma': 1,
                  'e_sg': p,
                 'e_en': 0.8}
    
        #create directory
        dir_string_p = os.path.join(dir_string, str(p))
        if not os.path.isdir(dir_string_p):
            os.mkdir(dir_string_p)
    
        #change working directory    
        os.chdir(dir_string_p)
        
        if os.path.isfile('nmi_scores.csv'):
            print 'already completed {}/{}, loading scores and continuing'.format(k1, p)
            nmi_scores = np.genfromtxt('nmi_scores.csv', delimiter=',')
            overall_nmi_scores[i,j] = np.mean(nmi_scores, axis=0)
            continue
        
        #copy executable
        ex = "benchmark.exe"   
        if not os.path.isfile(ex):

            source = "C:\\Users\\davem\\Documents\\PhD\\Benchmark Graph Generators\\binary_networks\\benchmark.exe"
            copyfile(source, ex)

        #make benchmark parameter file
        filename = "benchmark_flags_{}_{}.dat".format(k1,p)
        if not os.path.isfile(filename):
            with open(filename,"w") as f:
                f.write("-N {} -k {} -maxk {} -minc {} -maxc {} -mu {}".format(N, k, maxk, minc, maxc, mu))
            print 'written flag file: {}'.format(filename)
            
        #cmd strings
        change_dir_cmd = "cd {}".format(dir_string_p)
        generate_network_cmd = "benchmark -f {}".format(filename)

        #output of cmd
        output_file = open("cmd_output.out", 'w')

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
            
        print

        #generate networks
        for r in range(1, num_repeats+1):

            network_rename = "{}_{}".format(r,network)
            first_level_rename = "{}_{}".format(r,first_level)
            gml_filename = 'embedded_network_{}.gml'.format(r)

            if not os.path.isfile(network_rename):

                process = subprocess.Popen(change_dir_cmd + " && " + generate_network_cmd, 
                                        stdout=output_file, 
                                        stderr=output_file, 
                                        shell=True)
                process.wait()
                
                print 'generated graph {}'.format(r)

                os.rename(network, network_rename)
                os.rename(first_level, first_level_rename)
                
                print 'renamed graph {}'.format(r)

            if not os.path.isfile(gml_filename):

                ##embed graph
                embed_main(network_rename, first_level_rename)

                print 'embedded graph {} as {} in {}'.format(r, gml_filename, os.getcwd())

            ##score for this network
            if not np.all(nmi_scores[r-1]):

                start_time = time()

                print 'starting ghsom for: {}/{}/{}'.format(k1, p, gml_filename)
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
        print
    
print 'DONE'

print 'OVERALL NMI SCORES'
print overall_nmi_scores


# In[3]:

for scores in overall_nmi_scores:
    
    print scores
    idx = np.argsort(scores)[::-1]
    
    print parameter_settings[idx[0]]


# In[ ]:



