## this is the course project of FOGM15Fall, Bayesian probabilistic tensor decomposition
## this was designed with:
##  1. the Spike and Slab prior, on the gene factor matrix;
##  2. the tissue hierarchical regulation (on the tissue factor matrix) iteration scheme;
## however, here we only achieve the very baseline without the above, to:
##  1. stress more on a flexible framework, and probably with some OOP and Parallel schemes;
##  2. achieve some modules, and check the convergence of the basic model on our processed genomic dataset;
##  3. check the speed of the Gibbs sampler;
##  4. to be compared with later on for the move developed models with above attributes;
##  5. ...


## Jul.19: MVN doesn't work on real data, need to re-model



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
from scipy.stats import multivariate_normal


##=====================
##==== global variables
##=====================
#n_factor = 40
#n_individual = 100						#NOTE mw: 1066
#n_gene = 2000							#NOTE mw: 22473
#n_tissue = 30
n_factor = 40
n_tissue = 15
n_individual = 100
n_gene = 600

#dimension = (n_individual, n_gene, n_tissue)			#NOTE mw
dimension = (n_tissue, n_individual, n_gene)
factor_name = {}
dataset = np.zeros(shape=dimension) # individual x gene x tissue
markerset = np.ones(shape=dimension) # mark the position where there are data
individual_rep = {}
##==== mapping: #1: individual; #2: gene; #3: tissue
fmlist = []

prior1 = []     # 0: mean array; 1: precision matrix
prior2 = []     # as above
prior3 = []     # as above
prior = []

hyper_prior1 = []   # 0: xxx; 1: xxx; 2: xxx; 3: xxx
hyper_prior2 = []   # as above
hyper_prior3 = []   # as above
hyper_prior = []

alpha = 0       # the precision of the final observation
alpha_prior = []    # 0: xxx; 1: xxx


##==== NOTE: for re-using variables and saving computation
tensor_mean = []
tensor_null_n = 0



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

## re-write cal_product with numpy
'''
def cal_product(i, j, k):
	global n_factor
	global fmlist

	product = 0
	for count in range(n_factor):
		product += fmlist[0][i][count] * fmlist[1][j][count] * fmlist[2][k][count]
	return product
'''
def cal_product(i, j, k):
	global n_factor
	global fmlist

	array = np.multiply(fmlist[0][i], fmlist[1][j])
	product = np.inner(array, fmlist[2][k])
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





##==== this function will load simu data
def load_dataset_simu():
	global fmlist
	global prior1
	global prior2
	global prior3
	global prior

	global n_factor
	global n_tissue
	global n_individual
	global n_gene
	global dimension

	global dataset
	global markerset


	#== load fmlist from real data init	NOTE: TODO
	#factor_tissue = np.load('data/Tissue.npy')
	factor_tissue = np.load('data_init/Tissue.npy')
	fmlist.append(factor_tissue)
	#factor_indiv = np.load('data/Individual.npy')
	factor_indiv = np.load('data_init/Individual.npy')
	fmlist.append(factor_indiv)
	#factor_gene = np.load('data/Gene.npy')
	factor_gene = np.load('data_init/Gene.npy')
	fmlist.append(factor_gene)
	#== global variable init
	n_factor = len(factor_tissue[0])
	n_tissue = len(factor_tissue)
	n_individual = len(factor_indiv)
	n_gene = len(factor_gene)
	dimension = (n_tissue, n_individual, n_gene)

	#== calculate mean and cov (sample) for MVN prior; we use sample mean and cov as the init for this prior
	prior1 = []
	prior2 = []
	prior3 = []
	prior = []
	prior.append(prior1)
	prior.append(prior2)
	prior.append(prior3)
	for n in range(3):
		mean = np.array(np.mean(fmlist[n], axis = 0))
		cov = np.cov(fmlist[n], rowvar=0)
		precision = np.array(inv(cov))
		prior[n].append(mean)
		prior[n].append(precision)


	##========================
	##==== tensor data without indexing (e.g., full tensor data)
	##========================
	"""
	#== load tensor data, NOTE: for simu data, we don't need to index them
	dataset = []								# TODO: to load incomplete tensor, and also the markerset
	for i in range(n_tissue):
		dataset.append(np.load("data/Tensor_tissue_" + str(i) + ".npy"))
	dataset = np.array(dataset)
	print dataset.shape

	markerset = np.ones(shape=dimension)

	## NOTE: make it incomplete for the simu data
	'''
	markerset = np.ones(shape=dimension)
	for i in range(n_tissue):
		list_random = np.random.permutation(n_individual)
		for j in range(60):	# TODO: maybe tune this number later on: 10, 25, 40, 60
			index = list_random[j]
			markerset[i][index] = np.zeros(n_gene)


	## DEBUG
	print n_tissue * n_individual * n_gene
	print np.sum(markerset)
	print "# of factors:",
	print n_factor
	'''
	"""

	##========================
	##==== tensor data with indexing
	##========================
	#== load tensor data (and markerset)
	dataset = np.zeros(shape=dimension)
	markerset = np.zeros(shape=dimension)
	list_tissue = np.load("data_init/Tissue_list.npy")
	list_indiv = np.load("data_init/Individual_list.npy")
	rep_indiv = {}
	for i in range(len(list_indiv)):
		individual = list_indiv[i]
		rep_indiv[individual] = i

	## DEBUG
	print len(list_indiv)
	print len(list_tissue)
	count = 0

	for i in range(len(list_tissue)):
		tissue = list_tissue[i]
		list_sample = np.load("data_init/Tensor_tissue_" + str(i) + "_list.npy")

		## DEBUG
		count += len(list_sample)

		data = np.load("data_init/Tensor_tissue_" + str(i) + ".npy")
		for j in range(len(list_sample)):
			sample = list_sample[j]
			individual = (sample.split('-'))[1]		# NOTE: IndividualID of simulated samples
			index = rep_indiv[individual]
			dataset[i][index] = data[j]
			markerset[i][index] = np.ones(n_gene)

	## DEBUG
	print count

	dataset = np.array(dataset)
	markerset = np.array(markerset)
	print "dataset and markerset shape:",
	print dataset.shape,
	print markerset.shape

	## DEBUG
	np.save("temp/dataset", dataset)
	np.save("temp/markerset", markerset)








##==== this function will load real data
def load_dataset_real():
	global fmlist
	global prior1
	global prior2
	global prior3
	global prior

	global n_factor
	global n_tissue
	global n_individual
	global n_gene
	global dimension

	global dataset
	global markerset

	#== load fmlist from real data init
	factor_tissue = np.load('data/Tissue.npy')
	fmlist.append(factor_tissue)
	factor_indiv = np.load('data/Individual.npy')
	fmlist.append(factor_indiv)
	factor_gene = np.load('data/Gene.npy')
	fmlist.append(factor_gene)
	#== global variable init
	n_factor = len(factor_tissue[0])
	n_tissue = len(factor_tissue)
	n_individual = len(factor_indiv)
	n_gene = len(factor_gene)
	dimension = (n_tissue, n_individual, n_gene)

	#== calculate mean and cov (sample) for MVN prior; we use sample mean and cov as the init for this prior
	prior1 = []
	prior2 = []
	prior3 = []
	prior = []
	prior.append(prior1)
	prior.append(prior2)
	prior.append(prior3)
	for n in range(3):
		mean = np.array(np.mean(fmlist[n], axis = 0))
		cov = np.cov(fmlist[n], rowvar=0)
		precision = np.array(inv(cov))
		prior[n].append(mean)
		prior[n].append(precision)

	#== load tensor data (and markerset)
	dataset = np.zeros(shape=dimension)
	markerset = np.zeros(shape=dimension)
	list_tissue = np.load("data/Tissue_list.npy")
	list_indiv = np.load("data/Individual_list.npy")
	rep_indiv = {}
	for i in range(len(list_indiv)):
		individual = list_indiv[i]
		rep_indiv[individual] = i


	## DEBUG
	print len(list_indiv)
	print len(list_tissue)
	count = 0

	for i in range(len(list_tissue)):
		tissue = list_tissue[i]
		list_sample = np.load("data/Tensor_tissue_" + str(i) + "_list.npy")

		## DEBUG
		count += len(list_sample)

		data = np.load("data/Tensor_tissue_" + str(i) + ".npy")
		for j in range(len(list_sample)):
			sample = list_sample[j]
			individual = get_individual_id(sample)
			index = rep_indiv[individual]
			dataset[i][index] = data[j]
			markerset[i][index] = np.ones(n_gene)

	## DEBUG
	print count


	dataset = np.array(dataset)
	markerset = np.array(markerset)
	print "dataset and markerset shape:",
	print dataset.shape,
	print markerset.shape

	## DEBUG
	np.save("temp/dataset", dataset)
	np.save("temp/markerset", markerset)









##==== sampling from Gaussian (with mean and std)
def sampler_Normal(mu, sigma):
	sample = np.random.normal(mu, sigma)
	return sample

##==== sampling from Multivariate Gaussian
def sampler_MVN(mean, cov):
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
	global markerset
	global fmlist
	global prior
	global hyper_prior
	global alpha
	global alpha_prior
	global n_factor
	global dimension

	# We define that
	#   dimension 1: tissue
	#   dimension 2: individual
	#   dimension 3: gene

	cur_dimension = dimension[factor_id]
	print "cur_dimension:"
	print cur_dimension

	# DEBUG
	#print "we'll sample factor matrix first..."


	#==== sample factor matrix
	start = timeit.default_timer()
	for i in range(cur_dimension):  # there are n factor array, that are independent with each other --> parallel


		# DEBUG
		#print "now sampling factor#",
		#print i+1,
		#print "out of",
		#print cur_dimension

		precision_matrix = prior[factor_id][1]
		ids = [0,1,2]
		ids.remove(factor_id)
		dimension1 = dimension[ids[0]]
		dimension2 = dimension[ids[1]]

		#prod = np.ones(n_factor)
		#sampler_factor_helper(dataset, markerset, fmlist, dimension, n_factor, prod, [], factor_id, precision_matrix)

		Q = []
		for j in range(dimension1):
			for k in range(dimension2):
				hash_temp = {factor_id: i, ids[0]: j, ids[1]: k}
				index1 = hash_temp[0]
				index2 = hash_temp[1]
				index3 = hash_temp[2]
				if markerset[(index1, index2, index3)] == 0:
					continue

				#array = np.multiply(fmlist[ids[0]][j], fmlist[ids[1]][k])
				#precision_matrix = np.add(precision_matrix, alpha*np.dot(np.array([fmlist[ids[0]][j]]).T, np.array([fmlist[ids[0]][j]])))

				#vt_vector = np.multiply(v[i], t[j])
				Q.append(np.multiply(fmlist[ids[0]][j], fmlist[ids[1]][k]))

		Q = np.array(Q)
		precision_matrix = np.add(precision_matrix, np.multiply(alpha, np.dot(Q.T, Q)))
		cov = inv(precision_matrix)
		mean = np.dot(prior[factor_id][0], prior[factor_id][1])

		for j in range(dimension1):
			for k in range(dimension2):
				hash_temp = {factor_id: i, ids[0]: j, ids[1]: k}
				index1 = hash_temp[0]
				index2 = hash_temp[1]
				index3 = hash_temp[2]
				if markerset[(index1, index2, index3)] == 0:
					continue

				#array = np.multiply(fmlist[ids[0]][j], fmlist[ids[1]][k])
				mean = np.add(alpha * dataset[(index1, index2, index3)] * np.multiply(fmlist[ids[0]][j], fmlist[ids[1]][k]), mean)

		mean = np.dot(mean, cov)

		#== sampling
		fmlist[factor_id][i] = sampler_MVN(mean, cov)

	print "Time elapsed for sampling is ", timeit.default_timer() - start



	#==== sample Gaussian-Wishart (Wishart first, then Gaussian) prior
	start = timeit.default_timer()
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
	# TODO: to opt the following
	##====
	'''
	factor_mean = [0] * n_factor
	for i in range(dimension[factor_id]):
		for j in range(n_factor):
			factor_mean[j] += fmlist[factor_id][i][j]

	for i in range(n_factor):
		factor_mean[i] = factor_mean[i] * 1.0 / dimension[factor_id]

	factor_mean_t = np.array([factor_mean]).T
	mean_matrix = np.repeat(factor_mean_t, dimension[factor_id], axis=1)

	factor_var = np.dot(fmlist[factor_id].T - mean_matrix, (fmlist[factor_id].T - mean_matrix).T)
	'''
	##====
	## NOTE: the new routine is as followed:
	factor_mean = np.mean(fmlist[factor_id], axis=0)
	#factor_var = np.cov(fmlist[factor_id], rowvar=0)	# NOTE: seems not work (as we need to time (n-1) )
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

	
	print "Time elapsed for sampling is ", timeit.default_timer() - start




##==== sampling precision
def sampler_precision():
	global alpha_prior
	global dimension
	global n_individual
	global n_gene
	global n_tissue
	global dataset
	global markerset
	global alpha
	global tensor_mean		# NOTE: for re-using
	global tensor_null_n		# NOTE: for re-using

	print "old alpha" + str(alpha)

	para1_old = alpha_prior[0]
	para2_old = alpha_prior[1]


	'''
	n = 0
	s = 0

	start = timeit.default_timer()
	## new testing routine for tensor multi:
	tensor_mean = []
	for i in range(n_tissue):
		tensor_mean.append(np.zeros((n_individual, n_gene)))
		factor_indiv = fmlist[1].T
		factor_gene = fmlist[2].T
		for k in range(n_factor):
			factor = fmlist[0][i][k]
			tensor_mean[-1] += factor * np.outer(factor_indiv[k], factor_gene[k])
	tensor_mean = np.array(tensor_mean)
	for i in range(n_tissue):
		for j in range(n_individual):
			for k in range(n_gene):
				if markerset[(i,j,k)] == 0:
					continue
				#
				R_real = dataset[(i,j,k)]
				#R_exp = cal_product(i, j, k)
				R_exp = tensor_mean[(i,j,k)]
				s += math.pow(R_real - R_exp, 2)
				n += 1
	print "n, s:"
	print s
	print n
	print "Time elapsed:", timeit.default_timer() - start
	'''
	## re-write the following with numpy, but with limited improvement (say, from 4.86396002769 secs down to 4.06863713264 secs)
	start = timeit.default_timer()
	##====
	n = 0
	s = 0

	tensor_mean = []
	for i in range(n_tissue):
		tensor_mean.append(np.zeros((n_individual, n_gene)))
		factor_indiv = fmlist[1].T
		factor_gene = fmlist[2].T
		for k in range(n_factor):
			factor = fmlist[0][i][k]
			tensor_mean[-1] += factor * np.outer(factor_indiv[k], factor_gene[k])
	tensor_mean = np.array(tensor_mean)
	for i in range(n_tissue):
		for j in range(n_individual):
			for k in range(n_gene):
				if markerset[(i,j,k)] == 0:
					tensor_mean[(i,j,k)] = dataset[(i,j,k)]		# in order to make the log-likelihood as 0
					continue
				n += 1
	s = np.sum(np.power(np.add(dataset, -tensor_mean), 2))
	tensor_null_n = n
	##====
	print "s, n:"
	print s
	print n
	print "Time elapsed:", timeit.default_timer() - start



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
	like_log += (- math.log( math.sqrt(var) ))  # removed the \pi relevant constant item
	like_log += (- (obs - mean) * (obs - mean) / (2 * var))
	#print "ll" + str(like_log)
	return like_log


##==== log likelihood for multivariate Gaussian, NOTE: removed constant terms
def loglike_MVN(obs, mean, precision):
	#
	like_log = 0

	# precision term
	cov = inv(precision)
	'''
	cov_norm = LA.norm(cov)			## NOTE: what's this?
	like_log += (- 0.5 * math.log(cov_norm))
	'''
	#cov_det = abs(LA.det(cov))
	#like_log += (- 0.5 * math.log(cov_det))
	cov_log_det = LA.slogdet(cov)[1]
	print "cov_log_det:",
	print cov_log_det
	like_log += (- 0.5 * cov_log_det)
	'''
	cov_det = LA.det(cov)
	print "cov_det:",
	print cov_det
	like_log += (- 0.5 * math.log(cov_det))
	'''

	# distance term
	'''
	array1 = map(lambda x1,x2: x1-x2, obs, mean)
	array2 = []
	for i in range(n_factor):			# NOTE: this nested loop is problematic when n_factor is very large
		value = 0
		for j in range(n_factor):
			value += array1[j] * precision[j][i]
		array2.append(value)
	array3 = []
	temp = 0
	for i in range(n_factor):
		temp += array1[i] * array2[i]
	'''
	diff = obs - mean
	temp = np.dot(diff, precision)
	temp = np.dot(temp, diff)

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
	## NOTE:
	'''
	cov = inv(precision)
	like_log += multivariate_normal.logpdf(obs1, mean, cov)
	'''
	return like_log


##==== log likelihood for Gamma, NOTE: removed the constant terms
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
	global prior        # prior[0, 1, 2]: 0: mean array; 1: precision matrix
	global hyper_prior  # hyper_prior[0, 1, 2]: 0: scale; 1: df; 2: mean; 3: scaler
	global alpha
	global alpha_prior  # 0: shape parameter; 1: rate parameter
	global dimension
	## NOTE: I will re-use the "tensor_mean" array from global for the computing here, this depends on the current sampling order!!!
	global tensor_mean
	global tensor_null_n


	like_log = 0

	#==== observation likelihood
	start = timeit.default_timer()
	'''
	var = 1.0 / alpha
	for i in range(n_tissue):
		for j in range(n_individual):
			for k in range(n_gene):
				if markerset[(i,j,k)] == 0:
					continue

				obs = dataset[(i,j,k)]
				mean = cal_product(i, j, k)
				like_log += loglike_Gaussian(obs, mean, var)
	'''
	## new testing routine for tensor multi:
	'''
	tensor_mean = []
	for i in range(n_tissue):
		tensor_mean.append(np.zeros((n_individual, n_gene)))
		factor_indiv = fmlist[1].T
		factor_gene = fmlist[2].T
		for k in range(n_factor):
			factor = fmlist[0][i][k]
			tensor_mean[-1] += factor * np.outer(factor_indiv[k], factor_gene[k])
	tensor_mean = np.array(tensor_mean)
	'''
	'''
	##========
	var = 1.0 / alpha
	for i in range(n_tissue):
		for j in range(n_individual):
			for k in range(n_gene):
				if markerset[(i,j,k)] == 0:
					continue
				obs = dataset[(i,j,k)]
				mean = tensor_mean[(i,j,k)]
				like_log += loglike_Gaussian(obs, mean, var)
	print "gaussian:",
	print like_log
	print "Time elapsed for cal obs loglike is ", timeit.default_timer() - start
	##========
	'''
	## new routine with numpy
	var = 1.0 / alpha
	like_log += (- alpha / 2) * np.sum(np.power(np.add(dataset, -tensor_mean), 2)) + tensor_null_n * (- math.log( math.sqrt(var) ))
	print "gaussian:",
	print like_log
	print "Time elapsed for cal obs loglike is ", timeit.default_timer() - start



	#==== factor matrix likelihood
	##============
	a = 0
	cur_dimension = dimension[0]
	for i in range(cur_dimension):
		## NOTE: optimize the following one:
		#a += loglike_MVN(fmlist[0][i], prior[0][0], prior[0][1])

		# distance term			# NOTE: we can also optimize this do to matrix operation only once
		obs = fmlist[0][i]
		mean = prior[0][0]
		precision = prior[0][1]

		diff = obs - mean
		temp = np.dot(diff, precision)
		temp = np.dot(temp, diff)

		a += (- 0.5 * temp)
	# precision term
	precision = prior[0][1]
	cov = inv(precision)
	#cov_norm = LA.norm(cov)		## NOTE: what's this?
	#cov_det = abs(LA.det(cov))
	cov_log_det = LA.slogdet(cov)[1]
	print "cov_log_det:",
	print cov_log_det
	#a += cur_dimension * (- 0.5 * math.log(cov_det))
	a += cur_dimension * (- 0.5 * cov_log_det)
	'''
	cov_det = LA.det(cov)
	print "cov_det:",
	print cov_det
	a += cur_dimension * (- 0.5 * math.log(cov_det))
	'''

	like_log += a
	print "mvn_tissue:",
	print a

	##============
	a = 0
	cur_dimension = dimension[1]
	for j in range(cur_dimension):
		#a += loglike_MVN(fmlist[1][j], prior[1][0], prior[1][1])

		# distance term
		obs = fmlist[1][j]
		mean = prior[1][0]
		precision = prior[1][1]

		diff = obs - mean
		temp = np.dot(diff, precision)
		temp = np.dot(temp, diff)

		a += (- 0.5 * temp)
	# precision term
	precision = prior[1][1]
	cov = inv(precision)
	#cov_norm = LA.norm(cov)		## NOTE: what's this?
	#cov_det = abs(LA.det(cov))
	cov_log_det = LA.slogdet(cov)[1]
	print "cov_log_det:",
	print cov_log_det
	#a += cur_dimension * (- 0.5 * math.log(cov_det))
	a += cur_dimension * (- 0.5 * cov_log_det)
	'''
	cov_det = LA.det(cov)
	print "cov_det:",
	print cov_det
	a += cur_dimension * (- 0.5 * math.log(cov_det))
	'''

	like_log += a
	print "mvn_individual:",
	print a

	##============
	a = 0
	cur_dimension = dimension[2]
	for k in range(cur_dimension):
		#a += loglike_MVN(fmlist[2][k], prior[2][0], prior[2][1])

		# distance term
		obs = fmlist[2][k]
		mean = prior[2][0]
		precision = prior[2][1]

		diff = obs - mean
		temp = np.dot(diff, precision)
		temp = np.dot(temp, diff)

		a += (- 0.5 * temp)
	# precision term
	precision = prior[2][1]
	cov = inv(precision)
	#cov_norm = LA.norm(cov)		## NOTE: what's this?
	#cov_det = abs(LA.det(cov))
	cov_log_det = LA.slogdet(cov)[1]
	print "cov_log_det:",
	print cov_log_det
	#a += cur_dimension * (- 0.5 * math.log(cov_det))
	a += cur_dimension * (- 0.5 * cov_log_det)
	'''
	cov_det = LA.det(cov)
	print "cov_det:",
	print cov_det
	a += cur_dimension * (- 0.5 * math.log(cov_det))
	'''

	like_log += a
	print "mvn_gene:",
	print a



	"""
	## NOTE: test Numpy logpdf cal
	##============
	a = 0
	mean = prior[0][0]
	precision = prior[0][1]
	cov = inv(precision)
	cur_dimension = dimension[0]
	for i in range(cur_dimension):
		#a += loglike_MVN(fmlist[0][i], prior[0][0], prior[0][1])
		obs = fmlist[0][i]
		a += multivariate_normal.logpdf(obs, mean, cov)
	like_log += a
	print "mvn_tissue:",
	print a

	##============
	a = 0
	mean = prior[1][0]
	precision = prior[1][1]
	cov = inv(precision)
	cur_dimension = dimension[1]
	for j in range(cur_dimension):
		#a += loglike_MVN(fmlist[1][j], prior[1][0], prior[1][1])
		obs = fmlist[1][j]
		a += multivariate_normal.logpdf(obs, mean, cov)
	like_log += a
	print "mvn_individual:",
	print a

	##============
	a = 0
	mean = prior[2][0]
	precision = prior[2][1]
	cov = inv(precision)
	cur_dimension = dimension[2]
	for k in range(cur_dimension):
		#a += loglike_MVN(fmlist[2][k], prior[2][0], prior[2][1])
		obs = fmlist[2][k]
		a += multivariate_normal.logpdf(obs, mean, cov)
	like_log += a
	print "mvn_gene:",
	print a
	"""




	#==== factor prior likelihood
	a = 0
	for i in range(3):
		a += loglike_GW(prior[i][0], prior[i][1], hyper_prior[i][0], hyper_prior[i][1], hyper_prior[i][2], hyper_prior[i][3])
	like_log += a
	print "gw:",
	print a


	#==== precision/variance likelihood
	a = loglike_Gamma(alpha, alpha_prior[0], alpha_prior[1])
	like_log += a
	print "gamma:",
	print a


	return like_log





if __name__ == '__main__':


	print "enter program..."


	##==================================
	##==== loading and preparing dataset
	##==================================
	print "now start preparing the data..."
	#-> prepare the "dataset" and "markerset"		TODO: markerset is not properly loaded
	#-> prepare the initialized factor matrices (PCA), and their mean/cov as the MVN prior
	#load_dataset_simu()		# NOTE: load the simulated dataset
	load_dataset_real()		# NOTE: and also fill in the dimension
	print "finish data preparation..."



	##================================
	##==== initialize global variables
	##================================
	print "now initializing all the hyper-prior..."
	#== set the hyper-prior
	hyper_prior1 = []
	hyper_prior2 = []
	hyper_prior3 = []
	hyper_prior = []
	hyper_prior.append(hyper_prior1)
	hyper_prior.append(hyper_prior2)
	hyper_prior.append(hyper_prior3)
	'''
	##========
	for n in range(3):
		# 4 parts here: scale matrix, df, mean, scaler of the precision matrix
		hyper_prior[n].append(np.load("data/precision.npy"))    # lambda
		hyper_prior[n].append(np.int(np.load("data/v.npy")))        # TODO: tunable   v_0
		hyper_prior[n].append(np.load("data/mu.npy"))       # TODO: tunable  mu_0
		hyper_prior[n].append(np.int(np.load("data/kappa.npy")))        # TODO: tunable  kappa_0
	##========
	'''
	## NOTE: for those 400 factors, we need new dimension for these hyper priors
	for n in range(3):
		mu = np.zeros(n_factor)
		kappa = 2
		precision = np.identity(n_factor)
		v = n_factor			# TODO: greater than or equal to n_factor

		# 4 parts here: scale matrix, df, mean, scaler of the precision matrix
		hyper_prior[n].append(precision)    				# lambda
		hyper_prior[n].append(v)        				# TODO: tunable   v_0
		hyper_prior[n].append(mu)					# TODO: tunable  mu_0
		hyper_prior[n].append(kappa)					# TODO: tunable  kappa_0




	#== set the prior for precision (Gamma)
	alpha_prior = [1, 0.5]      # shape parameter and rate parameter, TODO: tunable


	#== drawing precision from Gaussian
	alpha = sampler_Normal(0, 1)


	print "finish variables initialization..."





	##==============================
	##==== sampler calling iteration
	##==============================
	start_time = timeit.default_timer()

	ITER = 500
	ll_result = []
	for i in range(ITER):
		print "@@@@current iteration#",
		print i+1

		#==== factor matrix (and prior)
		print "@@sample factors..."
		for j in range(3):
			start = timeit.default_timer()
			print "sampling factor#",
			print j+1
			sampler_factor(j)
			print "Time elapsed for sampling is ", timeit.default_timer() - start

		#==== precision
		print "@@sample precision..."
		start = timeit.default_timer()
		sampler_precision()
		print "Time elapsed for sampling is ", timeit.default_timer() - start

		#==== loglike
		print "@@calculate logolike..."
		start = timeit.default_timer()
		like_log = loglike_joint()  # monitor the log joint likelihood
		print "sampling done. the log joint likelihood is",
		print like_log
		print "Time elapsed for cal loglike is ", timeit.default_timer() - start


		# log the factor/loading matrix
		#np.save("result/ind_loading_brain_full_quantile", fmlist[0])			NOTE mw
		#np.save("result/gene_loading_brain_full_quantile", fmlist[1])			NOTE mw
		np.save("result/fm_tissue", fmlist[0])
		np.save("result/fm_indiv", fmlist[1])
		np.save("result/fm_gene", fmlist[2])

		ll_result.append(like_log)


	elapsed = timeit.default_timer() - start_time
	print "Time elapsed for sampling is ",
	print elapsed

	fo = open("result/loglike.txt","w+")
	for ele in ll_result:
		fo.write(str(ele)+"\n")
	fo.close()





