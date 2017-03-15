
# coding: utf-8

# In[1]:

import goenrich
import numpy as np
import pandas as pd
import pickle
import networkx as nx
from spearmint_ghsom import main_no_labels as ghsom_main


# In[2]:

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


# In[29]:

import sys
from time import sleep
from __future__ import division

len = 100

for i in range(100):
    line = "[" + "=" * int(i * len / 100) + ">" + "-" * int((100 - i - 1) * len / 100) + "]"
    sys.stdout.write("\r{}".format(line))
    sys.stdout.flush()
#     print line
    sleep(0.1)
sys.stdout.write(" DONE")


# In[4]:

import os

os.chdir("/home/david/Documents/ghsom")

G, map = load_obj('yeast_communities')
print 'num communities: {}'.format(len(map))


# In[5]:

import os

os.chdir("/home/david/Documents/ghsom")

dir_name = "uetz_communities"

if not os.path.isdir(dir_name):
    os.mkdir(dir_name)
    print 'made directory {}'.format(dir_name)
    
os.chdir(dir_name)

shortest_path = nx.floyd_warshall_numpy(map).astype(np.int)
np.savetxt("shortest_path.csv", shortest_path, fmt='%i', delimiter=",")
print 'written shortest path matrix'

c = 0
for n, d in map.nodes(data=True):
    ls = d['ls']
    with open('community_{}.txt'.format(c),'w') as f:
        for l in ls:
            f.write('{}\n'.format(G.node[l]['label']))
    print 'written community_{}.txt'.format(c)
    c += 1


# In[27]:

nx.adjacency_matrix(map).toarray().astype(np.int)


# In[28]:

n1, d1 = map.nodes(data=True)[1]

for n in map.neighbors(n1):
    print n, len(map.node[n]['ls'])


# In[29]:

n1, d1 = map.nodes(data=True)[1]
n1, d2 = map.nodes(data=True)[8]

H1 = G.subgraph(d1['ls'])
H2 = G.subgraph(d2['ls'])

l1 = [v for k,v in nx.get_node_attributes(H1, 'label').items()]
l2 = [v for k,v in nx.get_node_attributes(H2, 'label').items()]
print l2 


# In[30]:

O = goenrich.obo.ontology('db/go-basic.obo')

annot = goenrich.read.sgd('db/gene_association.sgd.gz')
gene2go = goenrich.read.gene2go('db/gene2go.gz')
values = {k: set(v) for k,v in annot.groupby('go_id')['db_object_symbol']}

# propagate the background through the ontology
background_attribute = 'gene2go'
goenrich.enrich.propagate(O, values, background_attribute)


# In[6]:

annot.loc[annot['db_object_synonym'].str.contains('YOL020W')==True]


# In[31]:

q1 = np.array(l1)
q2 = np.array(l2)


# In[32]:

for q in q1:
    print q
print
for q in q2:
    print q


# In[33]:

# for additional export to graphviz just specify the gvfile argument
# the show argument keeps the graph reasonably small
df1 = goenrich.enrich.analyze(O, q1, background_attribute)
df1 = df1.sort_values('q').dropna()
df1 = df1.loc[df1['rejected'] == 0.0]
df1 = df1.loc[df1['x'] > 0]

df2 = goenrich.enrich.analyze(O, q2, background_attribute)
df2 = df2.sort_values('q').dropna()
df2 = df2.loc[df2['rejected'] == 0.0]
df2 = df2.loc[df2['x'] > 0]


# In[34]:

df1


# In[16]:

pd.merge(df1, df2, on=['M', 'N', 'n', 'name', 'namespace', 'term'], how='inner')


# In[71]:

df


# In[72]:

# generate html
df.dropna().head(n = 20).to_html('q2.html')
# df.head(n = 20).to_html('example.html')


# In[ ]:



