from sklearn import *
import numpy as np
import matplotlib.pyplot as plt
from copy import *

if __name__=="__main__":
    #n_individual  = 1066
    #n_gene = 585
    n_factor = 40

    X = np.load('data/real_data_brain_c22_quantile.npy')
    print X.shape

    '''
    # shuffle test
    arr = np.arange(1066)
    np.random.shuffle(arr)

    result = []
    for i in range(1066):
        ind = arr[i]
        a = X[ind,:].tolist()
        result.append(a)

    result = np.matrix(result)
    '''

    pca = decomposition.PCA(n_components = n_factor)
    pca_ind = pca.fit_transform(X)
    pca_gene = pca.components_.T

    '''
    pca_ind1 = deepcopy(pca_ind)

    for i in range(1066):
        pca_ind1[arr[i], :] = pca_ind[i, :]

    pca_ind = pca_ind1
    '''
    print pca_ind.shape
    print pca_gene.shape

    #var_list = [np.var(pca_ind[:,i]) for i in range(n_factor)]
    #print var_list[0]
    #print var_list[1]
    #plt.plot(var_list)
    #plt.ylabel('Var of components')
    #plt.show()

    arr = pca.explained_variance_ratio_
    print arr
    print sum(arr)


    #plt.plot(arr)
    #plt.xlabel("Principle Components")
    #plt.ylabel("Explained ratio")
    #plt.title("PCA Explained Variance Ratio")
    #plt.show()
    np.save('data/real_individual_brain_c22_quantile', pca_ind)
    np.save('data/real_gene_c22_quantile', pca_gene)


    mu = np.mean(pca_ind, axis = 0)
    T = np.cov(pca_ind, rowvar=0)
    v = n_factor + 1
    kappa = 2


    np.save('data/mu', mu)
    np.save('data/kappa', kappa)
    np.save('data/precision', T)
    np.save('data/v', v)

    print mu
    print T
    print v
    print kappa
