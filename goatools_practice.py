
# coding: utf-8

# In[1]:

from goatools import obo_parser

go = obo_parser.GODag("./db/go-basic.obo")


# In[7]:

go_id3 = 'GO:0048584'
go_id4 = 'GO:0040008'
print(go[go_id3])
print(go[go_id4])


# In[5]:

from goatools.associations import read_gaf

associations = read_gaf("http://geneontology.org/gene-associations/gene_association.sgd.gz")


# In[8]:

from goatools.semantic import semantic_similarity

sim = semantic_similarity(go_id3, go_id4, go)
print('The semantic similarity between terms {} and {} is {}.'.format(go_id3, go_id4, sim))


# In[10]:

from goatools.semantic import TermCounts, ic

# First get the counts of each GO term.
termcounts = TermCounts(go, associations)

# Calculate the information content
go_id = "GO:0048584"
infocontent = ic(go_id, termcounts)
print('Information content ({}) = {}'.format(go_id, infocontent))


# In[11]:

from goatools.semantic import resnik_sim

sim_r = resnik_sim(go_id3, go_id4, go, termcounts)
print('Resnik similarity score ({}, {}) = {}'.format(go_id3, go_id4, sim_r))


# In[12]:

from goatools.semantic import lin_sim

sim_l = lin_sim(go_id3, go_id4, go, termcounts)
print('Lin similarity score ({}, {}) = {}'.format(go_id3, go_id4, sim_l))


# In[17]:

l = []
for line in open('gene_list.txt','r'):
    l.append(line.rstrip())
print l


# In[19]:

from goatools.go_enrichment import GOEnrichmentStudy

goeaobj = GOEnrichmentStudy(
        l, # List of mouse protein-coding genes
        associations, # geneid/GO associations
        go, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # default multipletest correction method


# In[ ]:



