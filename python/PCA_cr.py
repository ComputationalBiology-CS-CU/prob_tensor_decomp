from sklearn import *
import numpy as np
import matplotlib.pyplot as plt
from copy import *
from scipy import stats

if __name__=="__main__":
    #n_individual  = 1066
    #n_gene = 585
    n_factor = 20

    X = np.load('data/real_data/real_data_brain_c22_quantile.npy')
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
    np.save('data/real_data/real_individual_brain_c22_quantile', pca_ind)
    np.save('data/real_data/real_gene_brain_c22_quantile', pca_gene)


    mu = np.mean(pca_ind, axis = 0)
    T = np.cov(pca_ind, rowvar=0)
    v = n_factor + 1
    kappa = 2


    np.save('data/real_data/mu', mu)
    np.save('data/real_data/kappa', kappa)
    np.save('data/real_data/precision', T)
    np.save('data/real_data/v', v)

    print mu
    print T
    print v
    print kappa

    #====Sparsity ======
    gene_t = pca_gene.T
    z_t = []
    mu = []
    lamb = []
    pi = []
    for row in xrange(len(gene_t)):
        non_zero = []
        z_row = []
        for col in xrange(len(gene_t[row])):
            if abs(gene_t[row][col] - 0.0) > 0.02:
                non_zero.append(gene_t[row][col])
                z_row.append(1)
            else:
                z_row.append(0)
        mu.append(np.mean(non_zero))
        lamb.append(1/np.var(non_zero))
        z_t.append(z_row)
        pi.append(len(non_zero)/float(len(z_row)))

    np.save('data/real_data/pi_array', pi)
    print "pi is ", pi

    mu_NG = np.mean(mu)
    kappa_NG = 2
    # alpha_NG, loc, beta_NG = stats.gamma.fit(lamb, loc=0)
    alpha_NG, beta_NG = 1, 1

    np.save('data/real_data/mu_NG', mu_NG)
    np.save('data/real_data/kappa_NG', kappa_NG)
    np.save('data/real_data/alpha_NG', alpha_NG)
    np.save('data/real_data/beta_NG', beta_NG)

    print 'mu ', mu_NG
    print 'kappa ',kappa_NG
    print 'alpha_NG ',alpha_NG
    print 'beta_NG ',beta_NG
    # print 'loc ', loc

    # beta_alpha, beta_loc, beta_beta = stats.gamma.fit(pi, loc=0)
    beta_alpha, beta_beta = 2,2

    np.save('data/real_data/beta_alpha', beta_alpha)
    np.save('data/real_data/beta_beta', beta_beta)

    print beta_alpha
    print beta_beta
    # print beta_loc

    z = np.array(z_t).T
    np.save('data/real_data/z', z)
    print z







