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
from scipy.stats import t
from scipy.stats import norm
from scipy.stats import bernoulli
import math
from numpy import linalg as LA
import simulator



##=====================
##==== global variables
##=====================
n_factor = 0
n_individual = 0
n_gene = 0
n_tissue = 0
dataset = []		# individual x gene x tissue
markerset = []		# mark the position where there are data
individual_rep = {}
##==== mapping: #1(index = 0): individual; #2(index = 1): gene; #3(index = 2): tissue
fm1 = []
fm2 = []
fm3 = []
fm = []
prior1 = []		# 0: mean array; 1: precision matrix
prior2 = []		# as above
prior3 = []		# as above
prior = [] 
hyper_prior1 = []	# 0: xxx; 1: xxx; 2: xxx; 3: xxx
hyper_prior2 = []	# as above
hyper_prior3 = []	# as above
hyper_prior = []

sparsity_prior = [] # 0: mean; 1: precision; 2: z; 3: pi
sparsity_hyper_prior = [] # 0: alpha_0; 1: beta_0; 2: mu_0; 3: kappa_0; 4: alpha_pi; 5: gamma_pi


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


##==== load and format the target dataset
## TODO: the current dir are all in my local Mac (Shuo)
def data_prepare():
	global dataset
	global markerset
	global n_individual
	global n_gene
	global n_tissue
	global individual_rep


	#== get the number of tissue types
	file = open("/Users/shuoyang/Desktop/Genetics_GeneExpression/GTEx/workbench_v.6/data_processed/expression_by_etissue/tissue_list.txt", 'r')
	count = 0
	while 1:
		line = file.readline()
		if not line:
			break

		count += 1
	file.close()
	n_tissue = count


	#== get all the individuals
	individual_rep = {}	# map individual ID into its index in the tensor
	count = 0
	for i in range(n_tissue):
		index_tissue = i + 1
		file = open("/Users/shuoyang/Desktop/Genetics_GeneExpression/GTEx/workbench_v.6/data_processed/expression_by_etissue/tissue_" + str(index_tissue) + ".txt", 'r')
		line = (file.readline()).strip()
		line = line.split('\t')[1:]
		for sample in line:
			id = get_individual_id(sample)
			if id not in individual_rep:
				individual_rep[id] = count
				count += 1
		file.close()
	n_individual = len(individual_rep)


	#== get the number of genes
	file = open("/Users/shuoyang/Desktop/Genetics_GeneExpression/GTEx/workbench_v.6/data_processed/expression_by_etissue/tissue_1.txt", 'r')
	file.readline()
	count = 0
	while 1:
		line = file.readline()
		if not line:
			break

		count += 1
	file.close()
	n_gene = count


	#== initialize the empty tensor first of all
	dataset = []
	markerset = []
	for i in range(n_individual):
		dataset.append([])
		markerset.append([])
		for j in range(n_gene):
			dataset[i].append([])
			markerset[i].append([])
			for k in range(n_tissue):
				dataset[i][j].append(0)
				markerset[i][j].append(0)


	#== get all samples and fill them into the empty tensor (and the marker tensor)
	for k in range(n_tissue):
		index_tissue = k + 1
		file = open("/Users/shuoyang/Desktop/Genetics_GeneExpression/GTEx/workbench_v.6/data_processed/expression_by_etissue/tissue_" + str(index_tissue) + ".txt", 'r')
		# get the rep for gene array of individuals first
		rep_temp = {}		# all samples of this tissue
		index_individual_map = {}
		line = (file.readline()).strip()
		line = line.split('\t')
		for count in range(1, len(line)):
			sample = line[count]
			individual = get_individual_id(sample)
			if individual in rep_temp:			# there seems to be more than one samples in one tissue for some individuals
				continue

			index_individual_map[count] = individual
			rep_temp[individual] = []
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			for count in range(1, len(line)):
				if count not in index_individual_map:	# there seems to be more than one samples in one tissue for some individuals
					continue

				rpkm = float(line[count])
				individual = index_individual_map[count]
				rep_temp[individual].append(rpkm)
		file.close()

		# then map them to the tensor
		for individual in rep_temp:
			index = individual_rep[individual]
			for count in range(len(rep_temp[individual])):
				rpkm = rep_temp[individual][count]
				dataset[index][count][k] = rpkm
				markerset[index][count][k] = 1

	return

def data_prepare_simulation(useSpike, smallMatrix):
	#This function will prepare data from simulation data

	dataset = simulator.simulation(useSpike, smallMatrix, True)



##==== calculate < U_i, V_j, T_k >
def cal_product(i, j, k):
	global n_factor
	global fm

	product = 0
	for count in range(n_factor):
		product += fm[0][i][count] * fm[1][j][count] * fm[2][k][count]
	return product

##==== Deep copy of gene matrix (used when calculate V~nk)
def deep_copy(v):
	v_new = []
	for n in range(len(v)):
		v_row = []
		for k in range(len(v[0])):
			v_row.append(v[n][k])
		v_new.append(v_row)

	return v_new

def cal_pR(n,k,num_factor):
	global markerset
	global alpha

	result = 1

	v = deep_copy(fm[num_factor])

	v[n][k] = 0

	for i in range(n_individual):
		for j in range(n_gene):
			for t in range(n_tissue):
				if (markerset[i][j][t] != 0):
					#Assuming each R_{ij}^{k} is independent to others
					result = result * norm(cal_product(i,j,t), alpha**(-1)).pdf(dataset[i][j][t])

	return result

def cal_pR_1(n,k,num_factor):
	global alpha
	global dataset
	global fm
	global markerset
	global sparsity_prior
	global n_factor
	global n_gene

	v = deep_copy(fm[num_factor])

	N = n_gene

	term1 = 0
	term2 = 0
	term3 = 0
	for i in range(n_individual):
		for t in range(n_tissue):
			if (markerset[i][n][t] != 0):
				term1 += dataset[i][n][t] * fm[0][i][k] * fm[2][t][k]
				term3 += (fm[0][i][k] * fm[2][t][k])**2
				term2_rest = 0
				for m in range(n_factor):
					if (m != k):
						term2_rest += fm[0][i][m] * fm[2][t][m] * fm[1][n][m]
				term2 += fm[0][i][k] * fm[2][t][k] * term2_rest
	term1 = sparsity_prior[1][k][k]**2 * term1
	term3 = sparsity_prior[1][k][k]**2 * term3

	v[n][k] = (alpha**(-2) * sparsity_prior[0][k] + term1 - term2)/(term3 + alpha**(-2))

	likelihood = 1

	for i in range(n_individual):
		for j in range(n_gene):
			for t in range(n_tissue):
				if (markerset[i][j][t] != 0):
					#Assuming each R_{ij}^{k} is independent to others
					likelihood = likelihood * norm(cal_product(i,j,t), alpha**(-1)).pdf(dataset[i][j][t])

	sigma2 = N/(alpha**(2) * term3 + sparsity_hyper_prior[1][k][k]**(-2))
	result = likelihood * norm(sparsity_prior[0][k], sparsity_hyper_prior[1][k][k]).pdf(v[n][k]) * math.sqrt(2 * math.pi) * math.sqrt(sigma2) * (N**(-0.5))

	return result

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


##==== sampling from Gamma, for the variance
def sampler_Gamma(para1, para2):
	precision = 0
	#
	x = np.random.gamma(para1, para2, 1)
	precision = x[0]
	return precision

## ==== sampling Beta
def sampler_beta(a, b):
	return np.random.beta(a, b)

## ==== sampling Bernoulli
def sampler_bernoulli(p, size):
	return bernoulli.rvs(p, loc=0, size=size)

## ==== sampling student-t
def sampler_student_t(df, loc, scale):
	return t.rvs(df, loc = loc, scale = scale)


##==== sample factor relevant components (the factor matrix, and its Gaussian-Wishart prior)
def sampler_factor(num_factor):
	global dataset
	global fm
	global prior
	global hyper_prior
	global n_factor
	global alpha

	# We define that
	# 	dimension 1: tissue
	# 	dimension 2: individual
	# 	dimension 3: gene

	# num_factor
	dimension = 0
	dimension1 = 0		# for the other dimension
	dimension2 = 0		# as above
	num_factor_1 = 0	# the index of other factor
	num_factor_2 = 0	# as above
	if num_factor == 0:
		dimension = n_individual
		dimension1 = n_gene
		dimension2 = n_tissue
		num_factor_1 = 1
		num_factor_2 = 2
	elif num_factor == 1:
		dimension = n_gene
		dimension1 = n_individual
		dimension2 = n_tissue
		num_factor_1 = 0
		num_factor_2 = 2
	else:
		dimension = n_tissue
		dimension1 = n_individual
		dimension2 = n_gene
		num_factor_1 = 0
		num_factor_2 = 1



	# DEBUG
	print "we'll sample factor matrix first..."


	#==== sample factor matrix
	for i in range(dimension):	# there are n factor array, that are independent with each other --> parallel


		# DEBUG
		print "now sampling factor#",
		print i+1,
		print "out of",
		print dimension



		#== precision_matrix
		precision_matrix = []
		for count1 in range(n_factor):
			precision_matrix.append([])
			for count2 in range(n_factor):
				value = prior[num_factor][1][count1][count2]	# initialize with the prior precision matrix, then later update
				precision_matrix[count1].append(value)
		precision_matrix = np.array(precision_matrix)

		for j in range(dimension1):
			for k in range(dimension2):
				# re-arrange the three dimension, for querying original dataset
				hash_temp = {num_factor: i, num_factor_1: j, num_factor_2: k}
				index1 = hash_temp[0]
				index2 = hash_temp[1]
				index3 = hash_temp[2]
				if markerset[index1][index2][index3] == 0:
					continue
				array = np.multiply(fm[num_factor_1][j], fm[num_factor_2][k])
				# array = [0] * n_factor
				# for count in range(n_factor):
				# 	array[count] = fm[num_factor_1][j][count] * fm[num_factor_2][k][count]
				# for count1 in range(n_factor):
				# 	for count2 in range(n_factor):
				# 		precision_matrix[count1][count2] += alpha * array[count1] * array[count2]
		precision_matrix = np.add(precision_matrix, alpha*np.dot(np.array([array]).T, np.array([array])))
		cov = inv(precision_matrix)


		#== mean array
		# mean = [0] * n_factor

		# for count1 in range(n_factor):
		# 	for count2 in range(n_factor):
		# 		mean[count1] += prior[num_factor][0][count2] * prior[num_factor][1][count1][count2]

		#mean = np.dot(prior[num_factor][0], prior[num_factor][1].T)
		mean = np.dot(prior[num_factor][0], prior[num_factor][1])

		for j in range(dimension1):
			for k in range(dimension2):
				# re-arrange the three dimension, for querying original dataset
				hash_temp = {num_factor: i, num_factor_1: j, num_factor_2: k}
				index1 = hash_temp[0]
				index2 = hash_temp[1]
				index3 = hash_temp[2]
				if markerset[index1][index2][index3] == 0:
					continue
				# array = [0] * n_factor
				# for count in range(n_factor):
				# 	array[count] = fm[num_factor_1][j][count] * fm[num_factor_2][k][count]
				# R = dataset[index1][index2][index3]
				# for count in range(n_factor):
				# 	mean[count] += alpha * R * array[count]
				array = np.multiply(fm[num_factor_1][j], fm[num_factor_2][k])
				mean = np.add(alpha*dataset[index1][index2][index3] * array, mean)
		
		# array = [0] * n_factor
		# array = np.array(array)
		# for count1 in range(n_factor):
		# 	for count2 in range(n_factor):
		# 		array[count1] += mean[count2] * cov[count1][count2]
		# mean = array
		mean = np.dot(mean, cov.T)


		#== sampling
		fm[num_factor][i] = sampler_MVN(mean, cov)



	# DEBUG
	print "now we are sampling the prior..."



	#==== sample Gaussian-Wishart (Wishart first, then Gaussian) prior
	#== sample Wishart first



	# DEBUG
	print "preparing and sampling Wishart first..."



	# factor_mean

	# factor_mean = np.array([0] * n_factor)
	# for i in range(dimension):
	# 	for j in range(n_factor):
	# 		factor_mean[j] += fm[num_factor][i][j]
	# for i in range(n_factor):
	# 	factor_mean[i] = factor_mean[i] * 1.0 / dimension
	factor_mean = np.average(fm[num_factor], axis = 0)

	# factor_var
	factor_var = np.zeros((n_factor, n_factor))
	for i in range(dimension):
		# for count1 in range(n_factor):
		# 	for count2 in range(n_factor):
		# 		factor_var[count1][count2] += (fm[num_factor][i][count1] - factor_mean[count1]) * (fm[num_factor][i][count2] - factor_mean[count2])
		factor_var += np.dot((fm[num_factor][i] - factor_mean).T, (fm[num_factor][i] - factor_mean))

	#I think the following code is redundant
	'''
	for count1 in range(n_factor):
		for count2 in range(n_factor):
				factor_var[count1][count2] = factor_var[count1][count2] * 1.0 / dimension
	'''
	# cov_matrix
	cov_matrix = inv(hyper_prior[num_factor][0])
	# for count1 in range(n_factor):
	# 	for count2 in range(n_factor):
	# 			cov_matrix[count1][count2] += dimension * factor_var[count1][count2]
	cov_matrix = cov_matrix + factor_var
	temp = hyper_prior[num_factor][3] * dimension / ( hyper_prior[num_factor][3] + dimension )
	# for count1 in range(n_factor):
	# 	for count2 in range(n_factor):
	# 			cov_matrix[count1][count2] += temp * (hyper_prior[num_factor][2][count1] - factor_mean[count1]) * (hyper_prior[num_factor][2][count2] - factor_mean[count2])
	cov_matrix = cov_matrix + temp * np.dot((hyper_prior[num_factor][2] - factor_mean).T, (hyper_prior[num_factor][2] - factor_mean))
	precision_matrix = inv(cov_matrix)

	# df new
	df = hyper_prior[num_factor][1] + dimension

	## sampling Wishart
	prior[num_factor][1] = sampler_W(df, precision_matrix)



	# DEBUG
	print "now sampling MVN then..."



	#== sample Gaussian then
	# beta new
	beta = hyper_prior[num_factor][3] + dimension
	precision_matrix = beta * prior[num_factor][1]
	cov = inv(precision_matrix)

	# mean
	mean = (hyper_prior[num_factor][3] * hyper_prior[num_factor][2] + dimension * factor_mean) / (hyper_prior[num_factor][3] + dimension)
	
	# sampling MVN
	prior[num_factor][0] = sampler_MVN(mean, cov)

	return




# the following is implemented by Chuqiao Ren

##==== sample factor relevant components (when special prior involved, such like Spike and Slab prior)
def sampler_factor_sparsity():

	global dataset
	global fm
	global prior
	global hyper_prior
	global n_factor
	global alpha
	global sparsity_prior
	global sparsity_hyper_prior
	global maxIter

	#Only gene dimension will have sparsity prior

	dimension = n_gene
	dimension1 = n_individual
	dimension2 = n_tissue
	num_factor = 1
	num_factor_1 = 0
	num_factor_2 = 2

	# DEBUG
	print "we'll sample factor matrix first..."

	#DEBUG
	print "preparing to do LaPlace approximation..."
	p_0 = []     #A n-by-k matrix contains all the probability of z_nk = 0
	p_1 = []     #A n-by-k matrix contains all the probability of z_nk = 1

	for n in range(dimension):
		p_k0 = []
		p_k1 = []
		for k in range(n_factor):
			p_k0.append((1 - sparsity_prior[3][k])*cal_pR(n,k, num_factor))
			p_k1.append(sparsity_prior[3][k] * cal_pR_1(n,k,num_factor))
		p_0.append(p_k0)
		p_1.append(p_k1)

	print "now sampling z"
	#Note: bernoulli.rvs(p) takes p_1 (prob of getting a 1)
	z = []
	for i in range(dimension):
		#sampler_Bernoulli will return np.array object
		z.append(sampler_bernoulli(p_1/(p_0 + p_1), loc=0, size=n_factor))

	sparsity_prior[2] = np.array(z)

	#DEBUG
	print "now sampling V"


	#==== sample factor matrix
	for dim in range(dimension):
		# DEBUG
		print "now sampling factor#",
		print dim+1,
		print "out of",
		print dimension
        
		factor_ind = []
		z = sparsity_prior[2]
		for i in range(len(z[dim])):
			if (z[dim][i] == 1):
				factor_ind.append(i)

		#== precision_matrix
		precision_matrix = []
		for count1 in range(len(factor_ind)):
			count1_new = factor_ind[count1]
			precision_matrix.append([])
			for count2 in range(len(factor_ind)):
				count2_new = factor_ind[count2]
				value = sparsity_prior[1][count1_new][count2_new]	# initialize with the prior precision matrix, then later update
				precision_matrix[count1_new].append(value)
		precision_matrix = np.array(precision_matrix)

		for j in range(dimension1):
			for k in range(dimension2):
				# re-arrange the three dimension, for querying original dataset
				hash_temp = {num_factor: dim, num_factor_1: j, num_factor_2: k}
				index1 = hash_temp[0]
				index2 = hash_temp[1]
				index3 = hash_temp[2]
				if markerset[index1][index2][index3] == 0:
					continue
				array = [0] * len(factor_ind)
				for count in range(len(factor_ind)):
					count_new = factor_ind[count]
					array[count_new] = fm[num_factor_1][j][count_new] * fm[num_factor_2][k][count_new]
				for count1 in range(len(factor_ind)):
					count1_new = factor_ind[count1]
					for count2 in range(len(factor_ind)):
						count2_new = factor_ind[count2]
						precision_matrix[count1_new][count2_new] += alpha * array[count1_new] * array[count2_new]
		cov = inv(precision_matrix)


		#== mean array
		mean = [0] * len(factor_ind)
		mean = np.array(mean)
		for count1 in range(len(factor_ind)):
			count1_new = factor_ind[count1]
			for count2 in range(len(factor_ind)):
				count2_new = factor_ind[count2]
				mean[count1_new] += sparsity_prior[0][count2_new] * sparsity_prior[1][count1_new][count2_new]
		for j in range(dimension1):
			for k in range(dimension2):
				# re-arrange the three dimension, for querying original dataset
				hash_temp = {num_factor: i, num_factor_1: j, num_factor_2: k}
				index1 = hash_temp[0]
				index2 = hash_temp[1]
				index3 = hash_temp[2]
				if markerset[index1][index2][index3] == 0:
					continue
				array = [0] * len(factor_ind)
				for count in range(len(factor_ind)):
					count_new = factor_ind[count]
					array[count_new] = fm[num_factor_1][j][count_new] * fm[num_factor_2][k][count_new]
				R = dataset[index1][index2][index3]
				for count in range(len(factor_ind)):
					count_new = factor_ind[count]
					mean[count_new] += alpha * R * array[count_new]

		array = [0] * len(factor_ind)
		array = np.array(array)
		for count1 in range(len(factor_ind)):
			count1_new = factor_ind[count1]
			for count2 in range(len(factor_ind)):
				count2_new = factor_ind[count2]
				array[count1_new] += mean[count2_new] * cov[count1_new][count2_new]
		mean = array


		#== sampling
		fm[num_factor][dim] = sampler_MVN(mean, cov)


	# for dim in range(dimension):
	# 	# DEBUG
	# 	print "now sampling factor#",
	# 	print dim+1,
	# 	print "out of",
	# 	print dimension
    	
	# 	result_array = [0] * n_factor

	# 	for m in range(n_factor):
	# 		if (sparsity_prior[2][dim][m] != 0):
	# 			result_array[m] = sampler_Normal(sparsity_prior[0][m], sparsity_prior[1][m])

	# 	#== sampling
	# 	fm[num_factor][dim] = result_array

	# # DEBUG
	# print "now we are sampling the prior..."


	#==== sample Gaussian-Gamma (Gamma first, then Gaussian) prior
	#== sample Gamma first


	# DEBUG
	print "preparing and sampling Normal-Gamma..."
	precision = []
	mean = []
	for k in range(n_factor):
		v_n = []  #This is 1-dimensional vector V
		for n in range(dimension):
			v_n.append(fm[num_factor][n][k])

		#Now calcualte \bar{V}
		v_bar = sum(v_n)/dimension

		#Calcualte alpha_n
		alpha_n = sparsity_hyper_prior[0] + dimension/2

		#Calculate beta_n

		diff_v2 = []  #This is (V_n - \bar{V})^2
		for n in range(dimension):
			diff_v2.append((v_n[n] - v_bar)**2)

		beta_n = sparsity_hyper_prior[1] + sum(diff_v2)/2 + (sparsity_hyper_prior[3]*dimension*(v_bar - sparsity_hyper_prior[2])**2)/(2*(sparsity_hyper_prior[3] + dimension))

		#Calculate kappa_n

		kappa_n = sparsity_hyper_prior[3] + dimension

		#Calculate mu_n

		mu_n = (sparsity_hyper_prior[3] * sparsity_hyper_prior[2] + dimension * v_bar)/(kappa_n)

		#Update mean

		mean.append(sampler_student_t(2*alpha_n, mu_n, beta_n/(alpha_n * kappa_n)))

		#Update precision

		precision.append(sampler_Gamma(alpha_n, beta_n))

	#DEBUG
	print "Update sparsity prior (mean)..."

	sparsity_prior[0] = mean

	print "Make precision matrix"
	precision_m = []
	for i in range(len(precision)):
		precision_each = i * [0] + [precision[i]] + (len(precision) - 1 - i)*[0]
		precision_m.append(precision_each)

	sparsity_prior[1] = precision_m


	#DEBUG
	print "preparing and sampling pi..."

	pi = []

	for k in range(n_factor):
		z_k = []    #Each column of z
		for n in range(dimension):
			z_k.append(sparsity_prior[2][n][k])

		alpha_bar = sparsity_hyper_prior[4] + sum(z_k)
		gamma_bar = sparsity_hyper_prior[5] + dimension - sum(z_k)

		pi.append(sampler_beta(alpha_bar,gamma_bar))


	sparsity_prior[3] = pi

	return




##==== sampling precision
def sampler_precision():
	global alpha_prior
	global fm
	global n_individual
	global n_gene
	global n_tissue
	global dataset
	global markerset
	global alpha

	para1_old = alpha_prior[0]
	para2_old = alpha_prior[1]
	n = 0
	sum = 0

	for i in range(n_individual):
		for j in range(n_gene):
			for k in range(n_tissue):
				if markerset[i][j][k] == 0:
					continue
				#
				R_real = dataset[i][j][k]
				R_exp = cal_product(i, j, k)
				sum += math.pow(R_real - R_exp, 2)
				n += 1

	para1_new = para1_old + 0.5 * n
	para2_new = para2_old + 0.5 * sum

	alpha = sampler_Gamma(para1_new, para2_new)
	return



##==== log likelihood for univariate Gaussian
def loglike_Gaussian(obs, mean, var):
	like_log = 0
	like_log += (- math.log( math.sqrt(var) ))	# removed the \pi relevant constant item
	like_log += (- (obs - mean) * (obs - mean) / (2 * var))
	return like_log


##==== log likelihood for multivariate Gaussian, TODO: removed constant terms
def loglike_MVN(obs, mean, precision):
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
	global fm
	global prior		# prior[0, 1, 2]: 0: mean array; 1: precision matrix
	global hyper_prior	# hyper_prior[0, 1, 2]: 0: scale; 1: df; 2: mean; 3: scaler
	global alpha
	global alpha_prior	# 0: shape parameter; 1: rate parameter

	like_log = 0

	#==== edml
	for i in range(n_individual):
		for j in range(n_gene):
			for k in range(n_tissue):
				if markerset[i][j][k] == 0:
					continue

				obs = dataset[i][j][k]
				mean = cal_product(i, j, k)
				var = 1.0 / alpha
				like_log += loglike_Gaussian(obs, mean, var)


	#==== factor matrix likelihood
	for i in range(n_individual):
		like_log += loglike_MVN(fm[0][i], prior[0][0], prior[0][1])
	for j in range(n_gene):
		like_log += loglike_MVN(fm[1][j], prior[1][0], prior[1][1])
	for k in range(n_tissue):
		like_log += loglike_MVN(fm[2][k], prior[2][0], prior[2][1])


	#==== factor prior likelihood
	for i in range(3):
		like_log += loglike_GW(prior[i][0], prior[i][1], hyper_prior[i][0], hyper_prior[i][1], hyper_prior[i][2], hyper_prior[i][3])


	#==== precision/variance likelihood
	like_log += loglike_Gamma(alpha, alpha_prior[0], alpha_prior[1])


	return like_log
		




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
		scale = []				# TODO: not sure this is a good initialization, say, an appropriate scale matrix
		for i in range(n_factor):
			scale.append([])
			for j in range(n_factor):
				if j == i:
					scale[i].append(1)
				else:
					scale[i].append(0)
		scale = np.array(scale)
		hyper_prior[n].append(scale)
		hyper_prior[n].append(n_factor)		# TODO: tunable
		hyper_prior[n].append(0)		# TODO: tunable
		hyper_prior[n].append(1)		# TODO: tunable



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
		precision = []		# will be initialized with a diagonal 1 matrix
		for i in range(n_factor):
			precision.append([])
			for j in range(n_factor):
				if j == i:
					precision[i].append(1)
				else:
					precision[i].append(0)
		precision = np.array(precision)
		prior[n].append(mean)
		prior[n].append(precision)



	#== the MVN drawing (mean 0, cov 1) for factorized matrices
	fm1 = []
	fm2 = []
	fm3 = []
	fm = []
	mean = np.array([0] * n_factor)
	cov = []		# will be initialized with a diagonal 1 matrix
	for i in range(n_factor):
		cov.append([])
		for j in range(n_factor):
			if j == i:
				cov[i].append(1)
			else:
				cov[i].append(0)
	cov = np.array(cov)
	for i in range(n_individual):
		array = sampler_MVN(mean, cov)
		fm1.append(array)
	for j in range(n_gene):
		array = sampler_MVN(mean, cov)
		fm2.append(array)
	for k in range(n_tissue):
		array = sampler_MVN(mean, cov)
		fm3.append(array)
	fm.append(fm1)
	fm.append(fm2)
	fm.append(fm3)


	#== set the prior for precision (Gamma)
	alpha_prior = [1, 0.5]		# shape parameter and rate parameter, TODO: tunable


	#== drawing precision from Gaussian
	alpha = sampler_Normal(0, 1)



	# DEBUG
	print "finish variables initialization..."



	##==============================
	##==== sampler calling iteration
	##==============================
	ITER = 1000
	for i in range(ITER):
		print "current iteration#",
		print i+1
		for j in range(3):
			print "sampling factor#",
			print j+1,
			print "..."
			sampler_factor(j)
		print "sample precision..."
		sampler_precision()
		like_log = loglike_joint()	# monitor the log joint likelihood
		print "sampling done. the log joint likelihood is",
		print like_log





	##=====================
	##==== parameter saving
	##=====================
	# probably need this





