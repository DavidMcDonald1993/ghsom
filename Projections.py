
# coding: utf-8

# In[ ]:

##function to compute G = (I - P + W)^-1
#N: normalised laplacian
#d: list of degrees
#v: eigenvector corresponding to smallest eigenvalue
def compute_gr(N, B, d, v):
    
    #number of nodes in the graph
    n = len(N)
    
    #initialise Gr
    Gr = np.zeros((n, n))
    
    for i in range(n):
        
        #b2
        b2 = (np.transpose(v) @ B[:,i]) * v
        
        #b1 
        b1 = B[:,i] - b2
        
        #x1
        x1 = np.linalg.lstsq(N, b1)[0]
        
        #add to Gr
        Gr[:,i] = np.diag(d ** 0.5) @ (x1 + b2)
    
    #return Gr
    return Gr


# In[ ]:

##function to compute DSD matrix using AMG method
#N: laplacian matrix
#d: list of degrees
#v: smallest eigenvector
def dsd(N, d, v):
    
    #number of nodes in G
    n = len(N)
    
    #initialize dsd matrix
    dsd = np.zeros((n, n))
    
    #B
    B = np.diag(d ** -0.5) @ np.identity(n)
    
    #compute G
    G = compute_gr(N, B, d, v)
    
    print('computed greens matrix')
    
    #compute dsd for each pair
    for i in range(n):
        for j in range(i, n):
            
            #compute distance
            dis = np.linalg.norm(np.transpose(B[:,i] - B[:,j]) @ G, ord=1)
            
            #add to dsd matrix
            dsd[i,j] = dis
            dsd[j,i] = dis
    
    #return
    return dsd


# In[ ]:

##dsd embedding
def dsd_embedding(G):
    
    #number of nodes in the graph
    n = nx.number_of_nodes(G)
    
    #adjacency matrix
    A = nx.adjacency_matrix(G).toarray()
    
    #list of degrees
    deg = np.sum(A, axis=1)
    
    #normalised graph laplacian
    N = np.identity(n) - np.diag(deg ** -0.5) @ A @ np.diag(deg ** -0.5)
    
    print('constructed normalised laplacian')
    
    #eigen decompose normalised laplacian
    l, U = np.linalg.eigh(N)
    
    ##sort eigenvalues (and eigenvectors) into ascending order
    idx = l.argsort()
    l = l[idx]
    U = U[:,idx]
    
    #compute dsd matrix as metric
    D = dsd(N, deg, U[:,0]) 
    
    print('computed dsd matrix')
    
    #centreing matrix
    C = np.identity(n) - np.ones((n, n)) / n
    
    #similarity matrix
    K = - 1/2 * C @ D ** 2 @ C
    
    #eigen decompose K
    l, U = np.linalg.eigh(K)
    
    ##sort eigenvalues (and eigenvectors) into descending order
    idx = l.argsort()[::-1]
    l = l[idx]
    U = U[:,idx]
    
    #sum of all eigen values
    s = np.sum(l)
    
    #estimate the number of dimensions to keep
    k = len(l)
    var = 1
    
    while var > 0.95:
        k -= 1
        var = np.sum(l[:k]) / s
    
    k += 1
   
    print('reduced dimension of data',k)

    #run mds
    mds = man.MDS(n_components=k, max_iter=30000, dissimilarity="precomputed", n_jobs=1)
    emb = mds.fit(K)
    
    return emb.embedding_


# In[ ]:

##GHSOM2 algorithm

#H: graph
#lam: lambda -- the number of epochs to train before assessing error
#eta: learning rate
#sigma: initial neighbourhood range
#e_0: error of previous layer
#e_sg: error must reduce by this much for growth to stop
#e_en: error must be greater than this to expand neuron
#layer: layer of som
def ghsom2(G, lam, w, eta, sigma, e_0, e_sg, e_en, layer):
    
    print('num nodes',len(G))
    print('connected',nx.is_connected(G))
    
    #embedding
#     X = dsd_embedding(G)
    X = floyd_embedding(G)
    
    #save embedding to graph
    set_embedding(G, X, layer)
    
    print('embedded graph')
    
    #number of nodes in G
    num_nodes = nx.number_of_nodes(G)
    
    ##number of training patterns to visit
    N = min(num_nodes, 100)
    
    #create som for this neuron
    network = som.initialise_network(X, 1, w)
    
    #initialise MQE
    MQE = math.inf
    
    #train for l epochs
    som.train_network(X, network, lam, eta, sigma, N)
    
    #classify nodes
    som.assign_nodes(G, X, network, layer)
    
    while MQE >= e_sg * e_0:
        
        #save current error
        prev_MQE = MQE
    
        #find neuron with greatest error
        error_unit = som.identify_error_unit(network)
        
        if layer > 0:
            
            #expand network
            som.expand_network2(G, network, error_unit, layer)
                    
            print('ghsom has expanded som',layer,'error',MQE)
        
        #train for l epochs
        som.train_network(X, network, lam, eta, sigma, N)
        
        #delete superfluoous neurons
        som.delete_neurons(network)

        #classify nodes
        som.assign_nodes(G, X, network, layer)

        #calculate mean network error
        MQE = som.update_errors(network)
        
        if np.linalg.norm(MQE - prev_MQE) < 1e-3:
            
            print('no improvement, stopping growth')
            
            break
        
    print('ghsom has terminated expansion',layer)
    print('error',MQE)
    
    #recalculate error after neuron expansion
    MQE = 0
    
    ##neuron expansion phase
    #iterate thorugh all neruons and find neurons with error great enough to expand
    for i in network.nodes():
        
        #unpack
        ls = network.node[i]['ls']
        e = network.node[i]['e']
        
        #check error
        if e > e_en * e_0 and layer < len(labels) or e_0 == math.inf:
#         if layer < len(labels) and len(ls) > 0:

            if e_0 == math.inf:
                e_0 = e
        
            #subgraph
            H = G.subgraph(ls)
            
            #recursively run algorithm to create new network for subgraph of this neurons nodes
            net, e = ghsom2(H, lam, w, eta, sigma, e_0, e_sg, e_en, layer + 1)
            
            #repack
            network.node[i]['e'] = e
            network.node[i]['n'] = net
            
            print('ghsom has built new layer',layer+1)
            
        #increase overall network error
        MQE += e
    
    #mean MQE
    MQE /= nx.number_of_nodes(network)
    
    #return network
    return network, MQE


# In[ ]:

##LLE embed graph to euclidean space using graph laplacian
#G: graph
#k: number of nearest neighbours to look for
def LLE(G, k):
    
    #number of nodes in the graph
    n = nx.number_of_nodes(G)
    
    #adjacency graph matrix
    A = nx.floyd_warshall(G)
    
    #initialise weight matrix
    W = np.zeros((n, n))
    
    #iterate over all nodes to find neighbours
    for i in range(n):
        
        #node of graph
        node = G.nodes()[i]
        
        #distances from node
        dist = np.array([v for k,v in A[node].items()])
        
        #get id of k nearest nodes
        sorted_ids = np.argsort(dist)
        
        #k nearest neighbors
        knn = sorted_ids[1:k+1]
        
        #sum
        s = sum(dist[knn])
        
        #normalise weights for k nearest nodes
        for j in knn:
            W[i,j] = dist[j] / s
    
    #similarity matrix
    K = (k - 1) * np.identity(n) - k / n * np.ones((n,n)) + W + np.transpose(W) - np.transpose(W) @ W
    
    #eigen decompose K
    l, U = np.linalg.eig(W)
    
    #embedding into euclidean space
    X = U
    
    #return
    return X


# In[ ]:

##function for isomap embedding
#G: graph
#k: number of nearest neighbours
def isomap_embedding(G, k):
    
    #number of nodes in the graph
    n = nx.number_of_nodes(G)
    
    #distance matrix
    fl = nx.floyd_warshall(G)
    
    ##use MDS to compute kernel and embed graph into euclidean space

    #distance matrix
    D = np.zeros((n, n))
    
    for i in range(n):
        n1 = G.nodes()[i]
        for j in range(n):
            n2 = G.nodes()[j]
            D[i,j] = fl[n1][n2]
    
    #centering matrix
    C = np.identity(n) - np.ones((n, n)) / n
    
    #similarity matrix
    K = -0.5 * C * D * C
    
    ##eigen decompose K (real values)
    l, U = np.linalg.eigh(K)
    
    ##sort eigenvalues (and eigenvectors) into ascending order
    idx = l.argsort()
    l = l[idx]
    U = U[:,idx]
    
    ##dimensions to keep
    k = len(l[l < 1e-12])
    
    l = l[k:]
    U = U[:,k:]
    
    #diagonal array of eigen values
    lam = np.diag(l)
    lam_inv_sqrt = np.diag(l ** 0.5)
    
    #position in euclidean space
    X = U @ lam_inv_sqrt
    
    #return
    return X


# In[ ]:

#perform principle component analysis
# def pca(X, preserved_var):
def my_pca(X, k):

    #number of data points
    n = len(X)
    d = len(X[0])
    
    #centre
    X = X - np.ones((n, n)) @ X / n
    
    #estimate co-variance matrix
    C = 1 / n * np.transpose(X) @ X
    
    #eigen decompose C
    l, U = np.linalg.eigh(C)
    
    ##sort eigenvalues (and eigenvectors) into descending order
    idx = l.argsort()[::-1]
    l = l[idx]
    U = U[:,idx]
    
    #number of eigenvalues
    num_eig = len(l)
    
    #normalise U
    for j in range(num_eig):
        U[:,j] = U[:,j] / np.linalg.norm(U[:,j])
    
#     #sum of eigenvalues
#     sl = sum(l)
    
#     #determine number of dimensions to keep
#     k = num_eig - 1
#     var = 1
    
#     while var > preserved_var:
        
#         k -= 1
#         var = sum(l[:k]) / sl
        
#     k += 1
    
#     print("pca dim lost",d-k)
    
    #project X
    Y = X @ U[:,:k]
    
    #return
    return Y


# In[ ]:

##laplacian eigenmap embedding
def laplacian_eigenmap(G, k, sigma):
    
    #number of nodes in the graph
    n = nx.number_of_nodes(G)
    
    #distance matrix
    fl = nx.floyd_warshall(G)
    
    #intitialise distance matrix
    d = np.zeros((n, n))
    
    #initialise adjacency matric
    A = np.zeros((n, n))
    
    #find closest k neighbours
    for i in range(n):
        n1 = G.nodes()[i]
        for j in range(n):
            n2 = G.nodes()[j]
            d[i,j] = fl[n1][n2]
            
        #get id of k nearest nodes
        sorted_ids = np.argsort(d[i])
        
        #k nearest neighbors
        knn = sorted_ids[:k]
        
        #gaussian similarity -- heat kernel
        for j in knn:
            A[i,j] = np.exp( -d[i,j] ** 2 / (2 * sigma ** 2) )   
        
    #degree matrix
    D = np.sum(A, axis=1)
    
    #normalised graph laplacian
    N = np.identity(n) - np.diag(D ** -0.5) @ A @ np.diag(D ** -0.5)
    
    #eigen decompose
    l, U = np.linalg.eigh(N)
    
    ##sort eigenvalues (and eigenvectors) into ascending order
    idx = l.argsort()
    l = l[idx]
    U = U[:,idx]
    
    ##dimensions to keep
    k = len(l[l < 1e-12])
    
    print('embedding dim lost',k)
    
    #position matrix
    X = np.diag(D ** -0.5) @ U
    
    #return 
    return X


# In[ ]:

from sklearn import manifold

##function to read in dsd file and embed graph
def embed_from_dsd_file(file, n):

    #initialise D -- distance matrix
    D = np.zeros((n + 1, n + 1))
    
    #initialise line counter i
    i = 0
    
    #read dsd file
    with open(file) as dsd_file:
        
        #read each line
        for line in dsd_file:
            
            #append -1 to first line
            l = np.array([])
            
            if i == 0:
                l = np.append(l,-1)
            
            #save to D
            D[i,:] = np.append(l,np.fromstring(line,sep='\t'))
            
    
    #remove frst row and column
    D = D[1:,1:]
    
     #centreing matrix
#     C = np.identity(n) - np.ones((n, n)) / n
    
#     #similarity matrix
#     K = - 1/2 * C @ D @ C
#     K = np.exp(- D ** 2 / (2 * nx.radius(G) ** 2))
    #eigen decompose K
#     l, U = np.linalg.eigh(K)
    
#     ##sort eigenvalues (and eigenvectors) into descending order
#     idx = l.argsort()[::-1]
#     l = l[idx]
#     U = U[:,idx]
    
#     k = 3
    
# #     print('dsd embedding dim lost',k)
    
#     #position matrix
# #     X = U[:,k:] @ np.diag(l[k:] ** 0.5)
#     X = U[:,:k] @ np.diag(l[:k] ** 0.5)
    
# #     print(D[1,2])
# #     print(np.linalg.norm(X[1] - X[2]))

    mds = manifold.MDS(n_components=3, max_iter=300000, eps=1e-9, dissimilarity="precomputed", n_jobs=-1)

    emb = mds.fit(D)
    
    print('stress',emb.stress_)
    
    return emb.embedding_


# In[ ]:

from sklearn import manifold

##floyds embedding
def floyd_embedding(G):
    
    n = len(G)
    
    fl = nx.floyd_warshall(G)
    
    #intitialise distance matrix
    D = np.zeros((n, n))
    
    #find closest k neighbours
    for i in range(n):
        n1 = G.nodes()[i]
        for j in range(n):
            n2 = G.nodes()[j]
            D[i,j] = fl[n1][n2]
    
    C = np.identity(n) - np.ones((n, n)) / n
    
    #similarity matrix
#     K = - 1/2 * C @ D ** 2 @ C
    K = np.exp(- D ** 2 / (2 * nx.radius(G) ** 2))
    
    #eigen decompose K
    l, U = np.linalg.eigh(K)
    
    ##sort eigenvalues (and eigenvectors) into ascending order
    idx = l.argsort()[::-1]
    l = l[idx]
    U = U[:,idx]
    
    s = sum(l)
    
    k = len(l)
    var = 1
    
    while var > 0.95:
        k -= 1
        var = sum(l[:k]) / s
    
    k += 1
    print('number of dimensions',k)
    
#     #position matrix
#     X = U[:,:k] @ np.diag(l[:k] ** 0.5)
    
#     return X
    
#     mds = manifold.MDS(n_components=k, max_iter=300, dissimilarity="precomputed", n_jobs=1)
    mds = manifold.SpectralEmbedding(n_components=k, affinity="precomputed")
#     mds = manifold.TSNE(n_components=k, n_iter=10000000, metric="precomputed")
   
    emb = mds.fit(K)
    
#     print('stress',emb.stress_)
#     print(emb.affinity_matrix_)
#     print(emb.kl_divergence_)
    
    return emb.embedding_

