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
import seaborn as sns


##=====================
##==== global variables
##=====================
n_factor = 40
n_individual = 200
n_gene = 2000
n_tissue = 30
n_genotype = 10000
#The following parameters need to be determined by test-and-trials
#According to Barbara, they used alpha=beta=1 for the uniform on sparsity
#alpha = 1 beta = 2 is a line of y = -2x + 2
beta = [1,2]                #parameters of beta[alpha, beta]
normalGamma = [1,2,1,2]     #parameters of NG[mu, kappa, alpha, beta]
normalWishart = [[2,2],2,[[10,5],[5,10]],3]   #parameters of NW[mu, kappa, Lambda, v]

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
## TODO: the current dir are all in my local Mac
'''
def data_prepare():
    global n_individual
    global n_gene
    global n_tissue
    global n_genotype

    #== get the number of tissue types
    file = open("/Users/rachelren/Documents/Columbia_CS_Research/expression_by_etissue/tissue_list.txt", 'r')
    count = 0
    while 1:
        line = file.readline()
        if not line:
            break
        count += 1
    file.close()
    n_tissue = count

    # == get all the individuals
    individual_rep = {}  # map individual ID into its index in the tensor
    count = 0
    for i in range(n_tissue):
        index_tissue = i + 1
        file = open(
            "/Users/rachelren/Documents/Columbia_CS_Research/expression_by_etissue/tissue_" + str(
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
        "/Users/rachelren/Documents/Columbia_CS_Research/expression_by_etissue/tissue_1.txt",
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


    # == get the number of genotypes

'''

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
#    matrix = sample[0]
    return sample


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
        precisionMatrix = sampler_W(normalWishart[3], normalWishart[2])
        precisionMatrix_scaled = []
        for i in range(len(precisionMatrix)):
            temp = []
            for j in range(len(precisionMatrix[0])):
                temp.append(precisionMatrix[i][j]/normalWishart[1])
            precisionMatrix_scaled.append(temp)
        mu = np.random.multivariate_normal(normalWishart[0], precisionMatrix_scaled)
        u.append(np.random.multivariate_normal(mu, precisionMatrix))

    if (isTrans):
        return np.array(u).T
    else:
        return u

def compare_sparsity (v1, v2):
    k1 = len(n_gene)
    k2 = len(n_gene)

    sigma = []   #correlation matrix K1-by-K2
    for i in range(n_factor):
        sigma_row = []
        for j in range(n_factor):
            sigma_row.append(calc_corr(v1[i], v2[j]))
        sigma.append(sigma_row)

    sum_col = 0
    for row in range(len(sigma)):
        for col in range(len(sigma[row])):
            #TODO!
            continue

    sum_col = 0

def calc_corr (col1, col2):
    #TODO!
    return


if __name__ == '__main__':



    # DEBUG
    print "enter program..."


    # DEBUG
    print "now start preparing the data..."



    ##==================================
    ##==== loading and preparing dataset
    ##==================================
    '''
    data_prepare()          # prepare the "dataset" and "markerset"
    '''



    # DEBUG
    print "finish data preparation..."


    # DEBUG
    print "now initializing all the variables..."

    ##================================
    ##==== initialize global variables
    ##================================
    #n_factor = 400          # TODO: this is tunable, and the number 400 comes from results of other more complex methods
    useSpike = True

    #initialize normal wishart parameter
    mu = []
    precision = []
    for i in range(n_factor):
        mu.append(0)

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
    normalWishart[3] = n_factor + 1

    np.save("mu", normalWishart[0])
    np.save("kappa", normalWishart[1])
    np.save("precision", normalWishart[2])
    np.save("v", normalWishart[3])

    #print normalWishart[0]
    #print normalWishart[2]
    #print normalWishart[3]

    #DEBUG
    print "finish initializing all the variables..."

    #DEBUG
    print "now start simulation..."

    v = simulator_MVN(n_gene, False)
    u = simulator_MVN(n_individual, False)
    
    np.save("Individual", u)
    np.save("Gene", v)

    print len(u), len(u[0])
    print len(v), len(v[0])

    #DEBUG
    print "matrix multiplication..."

    #product = np.dot(u, v).T
    #np.save("IxG", product)
    '''
    #DEBUG
    print "prepare to draw matrix..."
    
    sns.set(context="paper", font="monospace")
    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(12, 9))

    
    #DEBUG
    print "draw product using seaborn..."
    
    # Draw the heatmap using seaborn
#    sns_plot = sns.heatmap(product)
    sns_plot = sns.heatmap(v)

    
#    f.tight_layout()
    
    #DEBUG
    print "finish drawing product..."
    
    #DEBUG
    print "saving figures..."
    
    fig = sns_plot.get_figure()
    fig.savefig("test.png")
    '''
    
    
    
    #DEBUG
    print "End of program..."
    
    


