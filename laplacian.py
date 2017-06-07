
# coding: utf-8

# In[9]:

import numpy as np
import networkx as nx


# In[2]:

G = nx.karate_club_graph()


# In[3]:

L = nx.laplacian_matrix(G)


# In[4]:

L


# In[11]:

d = np.array([v for k,v in nx.degree(G).items()])


# In[12]:

d


# In[14]:

L.dot(d)


# In[16]:

e = np.array([v for k,v in nx.betweenness_centrality(G).items()])


# In[17]:

L.dot(e)


# In[ ]:



