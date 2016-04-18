from sklearn import *
import numpy as np

if __name__=="__main__":
    n_individual  = 200
    n_gene = 2000
    n_factor = 40

    X = np.load('data/mat_simulation.npy')
    pca = decomposition.PCA(n_components = n_factor)
    pca_ind = pca.fit_transform(X)
    pca_gene = pca.components_.T
    print pca_ind.shape
    print pca_gene.shape


    np.save('data/pca_individual', pca_ind)
    np.save('data/pca_gene', pca_gene)



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
