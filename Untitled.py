
# coding: utf-8

# In[6]:

from urllib.request import urlopen
html = urlopen("http://www.google.com/")
# print(html.read())

for i in html.read():
    
    print(i)

