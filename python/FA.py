from sklearn import *
import numpy as np

if __name__=="__main__":
    n_individual  = 200
    n_gene = 2000
    n_factor = 40
    fmlist = []

    X = np.load('data/mat_simulation.npy')
    fa = decomposition.FactorAnalysis(n_components = n_factor)
    test_ind = fa.fit_transform(X)
    test_gene = fa.components_
    print test_ind.shape
    print test_gene.shape


    np.save('data/fa_individual', test_ind)
    np.save('data/fa_gene', test_gene)


    fmlist.append(test_ind)
    fmlist.append(test_gene)


    mu = np.mean(test_ind, axis = 0)
    T = np.cov(test_ind, rowvar=0)
    v = n_factor + 1
    kappa = 2


    np.save('data/mu', mu)
    np.save('data/kappa', kappa)
    np.save('data/precision', T)
    np.save('data/v', v)
