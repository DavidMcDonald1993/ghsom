
# coding: utf-8

# In[5]:

import numpy as np

for i in np.arange(0.5, 1, 0.1)[::-1]:
    print i


# In[4]:

import os
import sys
import time
import optparse
import subprocess
import numpy as np

def main():
    parser = optparse.OptionParser(usage="usage: %prog [options]")
    parser.add_option("--experiment-name", "-e", dest="experiment_name",
                      help="The name of the experiment in the database.",
                      type="string")
    
    
    kargs, args = parser.parse_args()
    
    print kargs
    
if __name__ == '__main__':
    main()


# In[10]:

import json

def write_config_file(exp_name):

    data = {"language": "PYTHON",
            "main-file": "speamrin_ghsom.py",
            "experiment-name": "{}_exp".format(exp_name),
            "likelihood": "GAUSSIAN",
            "variables": {
            "w": {"type": "FLOAT","size": 1,"min": 0.0001,"max": 1},
            "eta": {"type": "FLOAT","size": 1,"min": 0.0001,"max": 1},
            "sigma": {"type": "FLOAT","size": 1,"min": 0.001,"max": 1},
            "e_sg": {"type": "FLOAT","size": 1,"min": 0.3,"max": 1},
            "e_en": {"type": "FLOAT","size": 1,"min": 0.3,"max": 1}}}

    with open('test_config.json', 'w') as outfile:
        json.dump(data, outfile)


# In[19]:

def double(x = 2, y = 3):
    return 2 * x + y

print double(y=9)


# In[20]:

print str(len('value'.split(',')))


# In[9]:

import networkx as nx
import sklearn.metrics as met

G = nx.Graph()
G.add_edges_from([(1,2),(2,3)])

G.node[1]['act'] = 'a'
G.node[2]['act'] = 'a'
G.node[3]['act'] = 'b'

G.node[1]['pred'] = 'a'
G.node[3]['pred'] = 'b'

num_nodes = nx.number_of_nodes(G)

actual_community = nx.get_node_attributes(G, 'act')
predicted_community = nx.get_node_attributes(G, 'pred')

# assigned_nodes = 

labels_true = [v for k,v in actual_community.items() if k in predicted_community]
labels_pred = [v for k,v in predicted_community.items()]


score = met.normalized_mutual_info_score(labels_true, labels_pred) * len(labels_pred) / num_nodes
print score

