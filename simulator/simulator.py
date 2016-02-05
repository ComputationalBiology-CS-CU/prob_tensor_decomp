##===============
##==== libraries
##===============
import numpy as np
from numpy.linalg import inv
from scipy.stats import wishart
from scipy.stats import bernoulli
import math
from numpy import linalg as LA
import matplotlib.pyplot as plt

##=====================
##==== global variables
##=====================
n_factor = 0
n_individual = 0
n_gene = 0
#The following parameters need to be determined by test-and-trials
beta = [1,2]                #parameters of beta[alpha, beta]
normalGamma = [1,2,1,2]     #parameters of NG[mu, kappa, alpha, beta]
normalWishart = [[2,2],2,[[1,0],[0,1]],1]   #parameters of NW[mu, kappa, Lambda, v]

##=================================================
##===The Following code is adapted from Shuo's code
##=================================================

##=====================
##==== subroutines
##=====================
# get the "xxx-yyy" from "xxx-yyy-zzz-aaa-qqq", which is defined as the individual ID of the GTEx samples
def get_individual_id(s):
    ## naively find the second '-'
    id = ''
    count = 0
    for i in range(len(s)):
        if s[i] == '-':
            count += 1

        if count == 2:
            break

        id += s[i]

    return id


##==== load and format the target dataset
## TODO: the current dir are all in my local Mac (Shuo)
def data_prepare():
    global n_individual
    global n_gene

    # == get all the individuals
    individual_rep = {}  # map individual ID into its index in the tensor
    count = 0
    for i in range(n_tissue):
        index_tissue = i + 1
        file = open(
            "/Users/shuoyang/Desktop/Genetics_GeneExpression/GTEx/workbench_v.6/data_processed/expression_by_etissue/tissue_" + str(
                index_tissue) + ".txt", 'r')
        line = (file.readline()).strip()
        line = line.split('\t')[1:]
        for sample in line:
            id = get_individual_id(sample)
            if id not in individual_rep:
                individual_rep[id] = count
                count += 1
        file.close()
    n_individual = len(individual_rep)

    # == get the number of genes
    file = open(
        "/Users/shuoyang/Desktop/Genetics_GeneExpression/GTEx/workbench_v.6/data_processed/expression_by_etissue/tissue_1.txt",
        'r')
    file.readline()
    count = 0
    while 1:
        line = file.readline()
        if not line:
            break

        count += 1
    file.close()
    n_gene = count


##==== sampling from Gaussian (with mean and std)
def sampler_Normal(mu, sigma):
    sample = np.random.normal(mu, sigma)
    return sample


##==== sampling from Multivariate Gaussian
def sampler_MVN(mean, cov):
    array = [0] * len(mean)
    array = np.array(array)
    #
    x = np.random.multivariate_normal(mean, cov, 1).T
    for i in range(len(x)):
        y = x[i][0]
        array[i] = y
    return array


##==== sampling from Wishart
def sampler_W(df, scale):
    #
    sample = wishart.rvs(df, scale, size=1, random_state=None)
    matrix = sample[0]
    return matrix


##==== sampling from Gamma
def sampler_Gamma(alpha, beta):
    precision = 0
    x = np.random.gamma(alpha, beta, 1)
    precision = x[0]
    return precision


## ==== End of adaptation

## ==== sampling Beta
def sampler_beta(a, b):
    return np.random.beta(a, b)



## ==== Start to simulate
def simulator_spike_slab():
    #First, simulate Beta-Bernoulli conjugate prior
    z = []     #n_factor * n_gene array
    v = []     #n_factor * n_gene array

    #initiate v to all zeros
    for i in range(n_factor):
        temp = []
        for j in range(n_gene):
            temp = temp + [0]
        v.append(temp)

    #simulate pi and z
    for i in range(n_factor):
        pi = sampler_beta(beta[0], beta[1])
        z.append(bernoulli.rvs(pi, size = n_gene))

    #simulate v
    for i in range(n_factor):
        for j in range(n_gene):
            if (z[i][j] != 0):
                sigma = sampler_Gamma(normalGamma[2], normalGamma[3])
                v[i][j] = sampler_Normal(normalGamma[0], sigma/normalGamma[1])

    return v

def simulator_MVN(n_sample, isTrans):
    u = []

    for sample in range(n_sample):
        precisionMatrix = sampler_W(normalWishart[2], normalWishart[3])
        precisionMatrix_scaled = []
        for i in range(precisionMatrix):
            temp = []
            for j in range(precisionMatrix[0]):
                temp.append(precisionMatrix[i][j]/normalWishart[1])
            precisionMatrix_scaled.append(temp)
        mu = np.random.multivariate_normal(normalWishart[0], precisionMatrix_scaled)
        u.append(np.random.multivariate_normal(mu, precisionMatrix))

    if (isTrans):
        return np.array(u).T
    else:
        return u

if __name__ == '__main__':



    # DEBUG
    print "enter program..."


    # DEBUG
    print "now start preparing the data..."



    ##==================================
    ##==== loading and preparing dataset
    ##==================================
    data_prepare()			# prepare the "dataset" and "markerset"



    # DEBUG
    print "finish data preparation..."


    # DEBUG
    print "now initializing all the variables..."

    ##================================
    ##==== initialize global variables
    ##================================
    n_factor = 400			# TODO: this is tunable, and the number 400 comes from results of other more complex methods
    useSpike = True

    #initialize normal wishart parameter
    mu = []
    precision = []
    for i in range(n_factor):
        mu.append(2)

    for i in range(n_factor):
        temp = []
        for j in range(n_factor):
            if (i == j):
                temp.append(1)
            else:
                temp.append(0)
        precision.append(temp)

    normalWishart[0] = mu
    normalWishart[2] = precision

    if (useSpike):
        v = simulator_spike_slab()
        u = simulator_MVN(n_individual, False)

    else:
        v = simulator_MVN(n_gene, True)
        u = simulator_MVN(n_individual, False)















