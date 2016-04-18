## this is the course project of FOGM15Fall, Bayesian probabilistic tensor decomposition
## this was designed with:
##	1. the Spike and Slab prior, on the gene factor matrix;
##	2. the tissue hierarchical regulation (on the tissue factor matrix) iteration scheme;
## however, here we only achieve the very baseline without the above, to:
##	1. stress more on a flexible framework, and probably with some OOP and Parallel schemes;
##	2. achieve some modules, and check the convergence of the basic model on our processed genomic dataset;
##	3. check the speed of the Gibbs sampler;
##	4. to be compared with later on for the move developed models with above attributes;
##	5. ...



##===============
##==== libraries
##===============
import numpy as np
from numpy.linalg import inv
from scipy.stats import wishart
import math
from numpy import linalg as LA
from copy import *
import cProfile
import timeit

##=====================
##==== global variables
##=====================
n_factor = 40
n_individual = 200
n_gene = 2000
n_tissue = 30
dimension = (n_individual, n_gene)
factor_name = {}
dataset = np.zeros(shape=dimension) # individual x gene x tissue
markerset = np.ones(shape=dimension) # mark the position where there are data	
individual_rep = {}
##==== mapping: #1: individual; #2: gene; #3: tissue
fmlist = []
prior1 = []		# 0: mean array; 1: precision matrix
prior2 = []		# as above
prior3 = []		# as above
prior = []

hyper_prior1 = []	# 0: xxx; 1: xxx; 2: xxx; 3: xxx
hyper_prior2 = []	# as above
hyper_prior3 = []	# as above
hyper_prior = []

alpha = 0		# the precision of the final observation
alpha_prior = []	# 0: xxx; 1: xxx


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

def cal_product(i, j):
	global n_factor
	global fmlist

	product = 0
	for count in range(n_factor):
		product += fmlist[0][i][count] * fmlist[1][j][count]
	return product

def factor_combination(dataset, fmlist, dim_depth, n_factor, fm_depth, prod, path):
	if fm_depth==len(fmlist):
		ind_tup = tuple(path)
		dataset[ind_tup] = np.sum(prod)
		return 
	
	for i in range(dim_depth[fm_depth]):
		new_prod = np.multiply(prod, fmlist[fm_depth][i])
		path.append(i)
		factor_combination(dataset, fmlist, dim_depth, n_factor, fm_depth+1, new_prod, path)
		path.pop()



##==== this function will load 
def load_dataset():
	global dataset
	global fmlist
	global dimension
	global n_factor
	# load fmlist from simulated data
	fmlist.append(np.load('data/true_individual.npy'))
	fmlist.append(np.load('data/true_gene.npy'))
	prod = np.ones(n_factor)
	factor_combination(dataset, fmlist, dimension, n_factor, 0, prod, [])
	print dataset

##==== sampling from Gaussian (with mean and std)
def sampler_Normal(mu, sigma):
	sample = np.random.normal(mu, sigma)
	return sample

##==== sampling from Multivariate Gaussian
def sampler_MVN(mean, cov):
	'''
	array = [0] * len(mean)
	array = np.array(array)
	
	x = np.random.multivariate_normal(mean, cov, 1)
	for i in range(len(x[0])):
		y = x[0][i]
		array[i] = y
	return array
	'''
	x = np.random.multivariate_normal(mean, cov, 1)
	return x[0]

##==== sampling from Wishart
def sampler_W(df, scale):
	#
	sample = wishart.rvs(df, scale, size=1, random_state=None)
	return sample


##==== sampling from Gamma, for the variance
def sampler_Gamma(para1, para2):
	para2 = 1.0/para2
	x = np.random.gamma(para1, para2, 1)
	return x[0]


def sampler_factor_helper(dataset, markerset, fmlist, dim_depth, n_factor, prod, path, factor_id, precision_matrix):
	if fm_depth==len(fmlist):
		if markerset[tuple(path)] == 0:
			return
		precision_matrix = np.add(precision_matrix, alpha*np.dot(np.array([array]).T,np.array([array]))) 
		return 

	if fm_depth!=factor_id:
		for i in range(dim_depth[fm_depth]):
			new_prod = np.multiply(prod, fmlist[fm_depth][i])
			path.append(i)
			factor_combination(dataset, fmlist, dim_depth, n_factor, fm_depth+1, new_prod, path, factor_id, precision_matrix)
			path.pop()
	else:
		path.append(dim_depth[factor_id])
		factor_combination(dataset, fmlist, dim_depth, n_factor, fm_depth+1, prod, path, factor_id, precision_matrix)
		path.pop()


##==== sample factor relevant components (the factor matrix, and its Gaussian-Wishart prior)
def sampler_factor(factor_id):
	global dataset
	global fmlist
	global prior
	global hyper_prior
	global n_factor
	global alpha

	# We define that
	# 	dimension 1: tissue
	# 	dimension 2: individual
	# 	dimension 3: gene


	cur_dimension = dimension[factor_id]


	# DEBUG
	#print "we'll sample factor matrix first..."


	#==== sample factor matrix
	for i in range(cur_dimension):	# there are n factor array, that are independent with each other --> parallel


		# DEBUG
		#print "now sampling factor#",
		#print i+1,
		#print "out of",
		#print cur_dimension

		precision_matrix = prior[factor_id][1]
		ids = [0,1]
		ids.remove(factor_id)
		dimension1 = dimension[ids[0]]
		#dimension2 = dimension[ids[1]]

		#prod = np.ones(n_factor)
		#sampler_factor_helper(dataset, markerset, fmlist, dimension, n_factor, prod, [], factor_id, precision_matrix)

		Q = []
		for j in range(dimension1):
			hash_temp = {factor_id: i, ids[0]: j}
			index1 = hash_temp[0]
			index2 = hash_temp[1]
			if markerset[(index1, index2)] == 0:
				continue
			#array = np.multiply(fmlist[ids[0]][j], fmlist[ids[1]][k])
			#precision_matrix = np.add(precision_matrix, alpha*np.dot(np.array([fmlist[ids[0]][j]]).T, np.array([fmlist[ids[0]][j]])))

			#vt_vector = np.multiply(v[i], t[j])
			Q.append(fmlist[ids[0]][j])

		Q = np.array(Q)
		precision_matrix = np.add(precision_matrix, np.multiply(alpha, np.dot(Q.T, Q)))
		cov = inv(precision_matrix)
		mean = np.dot(prior[factor_id][0], prior[factor_id][1])


		for j in range(dimension1):
			hash_temp = {factor_id: i, ids[0]: j}
			index1 = hash_temp[0]
			index2 = hash_temp[1]
			if markerset[index1][index2] == 0:
				continue

			#array = np.multiply(fmlist[ids[0]][j], fmlist[ids[1]][k])
			mean = np.add(alpha * dataset[(index1, index2)] * fmlist[ids[0]][j], mean)
		
		mean = np.dot(mean, cov)

		#== sampling

		fmlist[factor_id][i] = sampler_MVN(mean, cov)
	
	# DEBUG
	#print "now we are sampling the prior..."

	#==== sample Gaussian-Wishart (Wishart first, then Gaussian) prior
	#== sample Wishart first

	# DEBUG
	#print "preparing and sampling Wishart first..."


	# factor_mean
	# factor_var
	# not sure
	'''
	factor_mean = np.array([0] * n_factor)
	for i in range(dimension[factor_id]):
		for j in range(n_factor):
			factor_mean[j] += fmlist[factor_id][i][j]
	
	for i in range(n_factor):
		factor_mean[i] = factor_mean[i] * 1.0 / dimension[factor_id]

	factor_var = np.zeros((n_factor, n_factor))
	for i in range(dimension[factor_id]):
		for count1 in range(n_factor):
			for count2 in range(n_factor):
				factor_var[count1][count2] += (fmlist[factor_id][i][count1] - factor_mean[count1]) * (fmlist[factor_id][i][count2] - factor_mean[count2])
	for count1 in range(n_factor):
		for count2 in range(n_factor):
			factor_var[count1][count2] = factor_var[count1][count2] * 1.0 / dimension[factor_id]
	
	print factor_var
	'''
	# mw
	factor_mean = [0] * n_factor
	for i in range(dimension[factor_id]):
		for j in range(n_factor):
			factor_mean[j] += fmlist[factor_id][i][j]
	
	for i in range(n_factor):
		factor_mean[i] = factor_mean[i] * 1.0 / dimension[factor_id]

	factor_mean_t = np.array([factor_mean]).T
	mean_matrix = np.repeat(factor_mean_t, dimension[factor_id], axis=1)

	factor_var = np.dot(fmlist[factor_id].T - mean_matrix, (fmlist[factor_id].T - mean_matrix).T)

	#print factor_var

	# cov_matrix
	'''
	cov_matrix = inv(hyper_prior[factor_id][0])
	for count1 in range(n_factor):
		for count2 in range(n_factor):
			cov_matrix[count1][count2] += factor_var[count1][count2]
	temp = hyper_prior[factor_id][3] * dimension[factor_id] / ( hyper_prior[factor_id][3] + dimension[factor_id] )
	for count1 in range(n_factor):
		for count2 in range(n_factor):
			cov_matrix[count1][count2] += temp * (hyper_prior[factor_id][2][count1] - factor_mean[count1]) * (hyper_prior[factor_id][2][count2] - factor_mean[count2])
	precision_matrix = inv(cov_matrix)
	'''

	# mw
	cov_matrix = inv(hyper_prior[factor_id][0])
	cov_matrix = np.add(cov_matrix, factor_var)
	temp = hyper_prior[factor_id][3] * dimension[factor_id] / ( hyper_prior[factor_id][3] + dimension[factor_id] )
	mean_diff = hyper_prior[factor_id][2]-factor_mean
	mean_diff = np.array([mean_diff])
	cov_matrix = np.add(cov_matrix, temp * np.dot(mean_diff.T, mean_diff))
	precision_matrix = inv(cov_matrix)

	# df new
	df = hyper_prior[factor_id][1] + dimension[factor_id]

	## sampling Wishart
	prior[factor_id][1] = sampler_W(df, precision_matrix)


	# DEBUG
	#print "now sampling MVN..."

	#== sample Gaussian
	# beta new
	beta = hyper_prior[factor_id][3] + dimension[factor_id]
	precision_matrix = beta * prior[factor_id][1]
	cov = inv(precision_matrix)

	# mean
	mean = (hyper_prior[factor_id][3] * hyper_prior[factor_id][2] + dimension[factor_id] * np.array(factor_mean)) / (hyper_prior[factor_id][3] + dimension[factor_id])
	
	# sampling MVN
	prior[factor_id][0] = sampler_MVN(mean, cov)

	#return


##==== sampling precision
def sampler_precision():
	global alpha_prior
	global n_individual
	global n_gene
	global n_tissue
	global dataset
	global markerset
	global alpha

	print "old alpha" + str(alpha)

	para1_old = alpha_prior[0]
	para2_old = alpha_prior[1]
	n = 0
	s = 0

	for i in range(n_individual):
		for j in range(n_gene):
			if markerset[i][j] == 0:
				continue
			#
			R_real = dataset[i][j]
			R_exp = cal_product(i, j)
			s += math.pow(R_real - R_exp, 2)
			n += 1

	para1_new = para1_old + 0.5 * n/2
	para2_new = para2_old + 0.5 * s/2

	print "para2_new: " + str(para2_new)
	alpha = sampler_Gamma(para1_new, para2_new)

	print "alpha" + str(alpha)
	return



##==== log likelihood for univariate Gaussian
def loglike_Gaussian(obs, mean, var):
	#print obs, mean
	#print "obs" + str(obs)
	#print "mean" + str(mean)
	#print "var" + str(var)
	like_log = 0
	like_log += (- math.log( math.sqrt(var) ))	# removed the \pi relevant constant item
	like_log += (- (obs - mean) * (obs - mean) / (2 * var))
	#print "ll" + str(like_log)
	return like_log


##==== log likelihood for multivariate Gaussian, TODO: removed constant terms
def loglike_MVN(obs, mean, precision):
	#
	like_log = 0

	# precision term
	cov = inv(precision)
	cov_norm = LA.norm(cov)
	like_log += (- 0.5 * math.log(cov_norm))

	# distance term
	array1 = map(lambda x1,x2: x1-x2, obs, mean)
	array2 = []
	for i in range(n_factor):
		value = 0
		for j in range(n_factor):
			value += array1[j] * precision[j][i]
		array2.append(value)
	array3 = []
	temp = 0
	for i in range(n_factor):
		temp += array1[i] * array2[i]
	like_log += (- 0.5 * temp)
	return like_log


##==== log likelihood for Gaussian-Wishart
def loglike_GW(obs1, obs2, scale, df, mean, scaler):
	like_log = 0
	# Wishart loglike
	like_log += wishart.logpdf(obs2, df, scale)
	# MVN loglike
	# scaler the precision first
	precision = []
	for i in range(n_factor):
		precision.append([])
		for j in range(n_factor):
			temp = obs2[i][j] * scaler
			precision[i].append(temp)
	like_log += loglike_MVN(obs1, mean, precision)
	return like_log


##==== log likelihood for Gamma, TODO: removed the constant terms
def loglike_Gamma(obs, para1, para2):
	like_log = 0
	like_log += (para1 - 1) * math.log(obs)
	like_log += (- para2 * obs)
	return like_log

##==== calculate the joint log likelihood
def loglike_joint():
	global n_individual, n_gene, n_tissue
	global dataset, markerset
	global fmlist
	global prior		# prior[0, 1, 2]: 0: mean array; 1: precision matrix
	global hyper_prior	# hyper_prior[0, 1, 2]: 0: scale; 1: df; 2: mean; 3: scaler
	global alpha
	global alpha_prior	# 0: shape parameter; 1: rate parameter

	like_log = 0


	#==== observation likelihood
	for i in range(n_individual):
		for j in range(n_gene):
			if markerset[i][j] == 0:
				continue

			obs = dataset[i][j]
			mean = cal_product(i, j)
			var = 1.0 / alpha
			like_log += loglike_Gaussian(obs, mean, var)
	print "gaussian: " + str(like_log)

	#print "LL is ",like_log

	#==== factor matrix likelihood
	a = 0
	for i in range(n_individual):
		a += loglike_MVN(fmlist[0][i], prior[0][0], prior[0][1])
	like_log += a
	print "mvn_individual: " + str(a)
	#for j in range(n_gene):
	#	like_log += loglike_MVN(fm[1][j], prior[1][0], prior[1][1])
	a = 0
	for k in range(n_gene):
		a += loglike_MVN(fmlist[1][k], prior[1][0], prior[1][1])
	like_log += a
	print "mvn_gene: " + str(a)

	#==== factor prior likelihood
	a = 0
	for i in range(2):
		a += loglike_GW(prior[i][0], prior[i][1], hyper_prior[i][0], hyper_prior[i][1], hyper_prior[i][2], hyper_prior[i][3])
	like_log += a
	print "gw: " + str(a)

	#==== precision/variance likelihood
	a = loglike_Gamma(alpha, alpha_prior[0], alpha_prior[1])
	like_log += a
	print "gamma: " + str(a)


	return like_log


if __name__ == '__main__':

	# DEBUG
	print "enter program..."


	# DEBUG
	print "now start preparing the data..."



	##==================================
	##==== loading and preparing dataset
	##==================================

	# prepare the "dataset" and "markerset"
	# data_prepare()
	load_dataset()		

	# DEBUG
	print "finish data preparation..."


	# DEBUG
	print "now initializing all the variables..."


	##================================
	##==== initialize global variables
	##================================
	#n_factor = 40			# TODO: this is tunable, and the number 400 comes from results of other more complex methods

	#== set the hyper-prior
	hyper_prior1 = []
	hyper_prior2 = []
	hyper_prior3 = []
	hyper_prior = []

	hyper_prior.append(hyper_prior1)
	hyper_prior.append(hyper_prior2)
	hyper_prior.append(hyper_prior3)


	
	for n in range(3):
		# 4 parts here: scale matrix, df, mean, scaler of the precision matrix
		
		scale = np.identity(n_factor)
		hyper_prior[n].append(np.load("data/precision.npy"))    # lambda
		hyper_prior[n].append(np.int(np.load("data/v.npy")))		# TODO: tunable   v_0
		hyper_prior[n].append(np.load("data/mu.npy"))		# TODO: tunable	 mu_0
		hyper_prior[n].append(np.int(np.load("data/kappa.npy")))		# TODO: tunable  kappa_0
	

	'''
	mu_0 = np.zeros(n_factor)
	precision_0 = np.identity(n_factor)

	# test some random parameters
	for n in range(3):
		# 4 parts here: scale matrix, df, mean, scaler of the precision matrix
		
		scale = np.identity(n_factor)
		hyper_prior[n].append(precision_0)    # lambda precision
		hyper_prior[n].append(60)		# TODO: tunable   v_0
		hyper_prior[n].append(mu_0)		# TODO: tunable	 mu_0
		hyper_prior[n].append(3)		# TODO: tunable  kappa_0
	'''

	#== the prior of MVN (mean and precision matrix)
	prior1 = []
	prior2 = []
	prior3 = []
	prior = []
	prior.append(prior1)
	prior.append(prior2)
	prior.append(prior3)
	for n in range(3):
		mean = np.array([0] * n_factor)
		precision = np.identity(n_factor)
		prior[n].append(deepcopy(mean))
		prior[n].append(deepcopy(precision))


	#== the MVN drawing (mean 0, cov 1) for factorized matrices
	mean = np.array([0] * n_factor)

	cov = np.identity(n_factor)

	#== set the prior for precision (Gamma)
	alpha_prior = [1, 0.5]		# shape parameter and rate parameter, TODO: tunable


	#== drawing precision from Gaussian
	alpha = sampler_Normal(0, 1)



	# DEBUG
	print "finish variables initialization..."



	##==============================
	##==== sampler calling iteration
	##==============================
	ITER = 100
	ll_result = []
	for i in range(ITER):
		print "current iteration#",
		print i+1
		for j in range(2):
			print "sampling factor#",
			print j+1
			#print "..."
			sampler_factor(j)
		#print "sample precision..."
		sampler_precision()
		like_log = loglike_joint()	# monitor the log joint likelihood
		#print "sampling done. the log joint likelihood is",
		print like_log
		ll_result.append(like_log)


	fo = open("result/true_param.txt","w+")
	for ele in ll_result:
		fo.write(str(ele)+"\n")

	fo.close()



