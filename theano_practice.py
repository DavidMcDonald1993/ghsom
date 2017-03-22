
# coding: utf-8

# In[ ]:

from theano import function, config, shared, sandbox
import theano.tensor as T
import numpy as np
import time

vlen = 10 * 30 * 768  # 10 x #cores x # threads per core
iters = 1000

rng = np.random.RandomState(22)
x = shared(np.asarray(rng.rand(vlen), config.floatX))
f = function([], T.exp(x))
print(f.maker.fgraph.toposort())
t0 = time.time()
for i in range(iters):
    r = f()
t1 = time.time()
print("Looping %d times took %f seconds" % (iters, t1 - t0))
print("Result is %s" % (r,))
if np.any([isinstance(x.op, T.Elemwise) for x in f.maker.fgraph.toposort()]):
    print('Used the cpu')
else:
    print('Used the gpu')


# In[1]:

import os
print os.environ['USERPROFILE']


# In[9]:

import networkx as nx
import numpy as np

G = nx.read_gml('embedded_karate.gml')
print np.sum(list(G.degree().values()))


# In[ ]:

import theano.tensor as T


# In[ ]:

import theano
print 'imported theano'


# In[1]:

import numpy as np
import time
import theano
A = np.random.rand(1000,10000).astype(theano.config.floatX)
B = np.random.rand(10000,1000).astype(theano.config.floatX)
np_start = time.time()
AB = A.dot(B)
np_end = time.time()
X,Y = theano.tensor.matrices('XY')
mf = theano.function([X,Y],X.dot(Y))
t_start = time.time()
tAB = mf(A,B)
t_end = time.time()
print("NP time: %f[s], theano time: %f[s] (times should be close when run on CPU!)" %(
                                           np_end-np_start, t_end-t_start))
print("Result difference: %f" % (np.abs(AB-tAB).max(), ))


# In[2]:

import numpy as np

density = np.genfromtxt('density.txt')

print density


# In[2]:

get_ipython().run_cell_magic(u'R', u'', u'\nx <- c(1,2,3)\nx')


# In[ ]:



