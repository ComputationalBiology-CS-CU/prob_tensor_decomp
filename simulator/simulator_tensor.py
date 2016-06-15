## this simulation is for gene expression tensor, without genotype
## we will probably try the Sparsity prior for tensor decomposition later on (L_p, ARD, Spike and Slab); NOTE: here I commented all the Spike and Slab code


##===============
##==== libraries
##===============
import numpy as np
from numpy.linalg import inv
from scipy.stats import wishart
from scipy.stats import bernoulli
import math
from numpy import linalg as LA
from numpy.linalg import inv
#import matplotlib.pyplot as plt				# not supported in C2B2 cluster
#import seaborn as sns						# not supported in C2B2 cluster
import re

# for individual ID: get the "xxx-yyy" from "xxx-yyy-zzz-aaa-qqq", which is defined as the individual ID of the GTEx samples
pattern_indiv = re.compile(r'^(\w)+([\-])(\w)+')





##=====================
##==== global variables
##=====================
n_factor = 40			# TODO: this is tunable, and the number 400 comes from results of other more complex methods; 40 is for one chr
n_individual = 100		# the same magnitude
n_gene = 2000			# 10% of the original
n_tissue = 30			# the same magnitude
'''
#The following parameters need to be determined by test-and-trials
#According to Barbara, they used alpha=beta=1 for the uniform on sparsity
#alpha = 1 beta = 2 is a line of y = -2x + 2
beta = [1,2]                #parameters of beta[alpha, beta]
normalGamma = [1,2,1,2]     #parameters of NG[mu, kappa, alpha, beta]
'''
normalWishart = [[2,2],2,[[10,5],[5,10]],3]   #priors of NW[mu, kappa, Lambda, v]





##=====================
##==== subroutines
##=====================
def get_individual_id(s):
	match = p.match(s)
	if match:
		return match.group()
	else:
		print "!!! no individual ID is found..."
		return ""


##==== sampling from Gaussian (with mean and std)
def sampler_Normal(mu, sigma):
	sample = np.random.normal(mu, sigma)
	return sample


##==== sampling from Multivariate Gaussian
def sampler_MVN(mean, cov):
	array = np.zeros(len(mean))
	x = np.random.multivariate_normal(mean, cov, 1).T
	for i in range(len(x)):
		y = x[i][0]
		array[i] = y
	return array


##==== sampling from Wishart
def sampler_W(df, scale):
	sample = wishart.rvs(df, scale, size=1, random_state=None)
	return sample


##==== sampling from Gamma (for precision of the uni-Gaussian)
def sampler_Gamma(alpha, beta):
	return np.random.gamma(alpha, beta, 1)[0]


'''
## ==== sampling from Beta
def sampler_beta(a, b):
	return np.random.beta(a, b)
'''





## ==== Start to simulate
'''
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
	z = np.array(z)

	#simulate v
	for i in range(n_factor):
		for j in range(n_gene):
			if (z[i][j] != 0):
				lamb = sampler_Gamma(normalGamma[2], normalGamma[3])
				v[i][j] = sampler_Normal( normalGamma[0], 1.0/(lamb*normalGamma[1]) )
	v = np.array(v)
	return v
'''


def simulator_MVN(mean, cov, n_sample):
	M = []
	for i in range(n_sample):
		sample = sampler_MVN(mean, cov)
		M.append(sample)
	M = np.array(M)
	return M





if __name__ == '__main__':


	##================================
	##==== initialize global variables
	##================================
	print "now initializing all the variables..."

	#initialize normal wishart parameter
	#normalWishart = [[2,2],2,[[10,5],[5,10]],3]   #parameters of NW[mu, kappa, Lambda, v]
	mu = np.zeros(n_factor)
	precision = np.identity(n_factor)

	normalWishart[0] = mu
	normalWishart[2] = precision
	normalWishart[3] = n_factor		# TODO: greater than or equal to n_factor

	np.save("./para_data/mu", normalWishart[0])
	np.save("./para_data/kappa", normalWishart[1])
	np.save("./para_data/precision", normalWishart[2])
	np.save("./para_data/v", normalWishart[3])



	##================================
	##==== simulation
	##================================
	print "now start simulating..."

	#==== individual factor matrix
	precisionMatrix = sampler_W(normalWishart[3], normalWishart[2])
	precisionMatrix_scaled = precisionMatrix * normalWishart[1]
	cov = inv(precisionMatrix_scaled)
	mu = sampler_MVN(normalWishart[0], cov)
	#
	np.save("./para_data/Mu_indiv", mu)
	np.save("./para_data/Lambda_indiv", precisionMatrix)
	#
	cov = inv(precisionMatrix)
	U = simulator_MVN(mu, cov, n_individual)
	np.save("./para_data/Individual", U)

	print "individual factor matrix:",
	print len(U), len(U[0])


	#==== gene factor matrix
	precisionMatrix = sampler_W(normalWishart[3], normalWishart[2])
	precisionMatrix_scaled = precisionMatrix * normalWishart[1]
	cov = inv(precisionMatrix_scaled)
	mu = sampler_MVN(normalWishart[0], cov)
	#
	np.save("./para_data/Mu_gene", mu)
	np.save("./para_data/Lambda_gene", precisionMatrix)
	#
	cov = inv(precisionMatrix)
	V = simulator_MVN(mu, cov, n_gene)
	np.save("./para_data/Gene", V)

	print "gene factor matrix:",
	print len(V), len(V[0])


	#==== tissue factor matrix
	precisionMatrix = sampler_W(normalWishart[3], normalWishart[2])
	precisionMatrix_scaled = precisionMatrix * normalWishart[1]
	cov = inv(precisionMatrix_scaled)
	mu = sampler_MVN(normalWishart[0], cov)
	#
	np.save("./para_data/Mu_tissue", mu)
	np.save("./para_data/Lambda_tissue", precisionMatrix)
	#
	cov = inv(precisionMatrix)
	W = simulator_MVN(mu, cov, n_tissue)
	np.save("./para_data/Tissue", W)

	print "tissue factor matrix:",
	print len(W), len(W[0])




	'''
	#==== the tensor
	tensor = np.zeros((n_individual, n_gene))
	U = U.T
	V = V.T
	for i in range(n_factor):
		array_indiv = np.array(U[i])
		array_gene = np.array(V[i])
		tensor += np.outer(array_indiv, array_gene)
	np.save("./para_data/Tensor", tensor)

	print "tensor dimension:",
	print len(tensor), len(tensor[0])
	'''




	#==== the tensor
	U = U.T
	V = V.T
	for i in range(n_tissue):
		tensor = np.zeros((n_individual, n_gene))
		array_tissue = W[i]
		for j in range(n_factor):
			array_indiv = np.array(U[j])
			array_gene = np.array(V[j])
			factor_tissue = array_tissue[j]
			tensor += factor_tissue * np.outer(array_indiv, array_gene)
		np.save("./para_data/Tensor_tissue_" + str(i), tensor)

		print "tissue#",
		print i,
		print "tensor dimension:",
		print len(tensor), len(tensor[0])




