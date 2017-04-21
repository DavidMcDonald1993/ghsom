
# coding: utf-8

# In[1]:

from multiprocessing import Pool
from itertools import repeat


# In[19]:

def f((x, y, z)):
    if x + y < 100:
        return x + y * z


# In[20]:

inputs = zip(range(10000), repeat(6), reversed(range(10000)))


# In[21]:

inputs[:10]


# In[22]:

time a=map(f, inputs)


# In[23]:

p = Pool(processes=4)


# In[24]:

get_ipython().run_cell_magic(u'time', u'', u'b = p.map(f, inputs)\np.close()')


# In[17]:

zip(range(10), repeat(2), range(10), [5,6])


# In[ ]:



