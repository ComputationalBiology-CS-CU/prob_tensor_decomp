## init the simulated tensor for sampling test (other than using the real factor matrices for init)

## we can keep the complete tensor, or make it incomplete

## the input is the output of the tensor simulator

## the final purpose of this script: we know the data is generated this way, we try to test whether init complete/incomplete tensor with the way we designed will give the sampling enough meaning




## NOTE: this script will work on incomplete tensor (first make the tensor incomplete, then init the model)
##	we will also need the R script





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
from sklearn.decomposition import PCA



n_factor = 0
n_tissue = 0
n_individual = 0
n_gene = 0
dimension = (n_tissue, n_individual, n_gene)



##===============
##==== sub-routines
##===============
# get the "yyy" from "xxx-yyy"
def get_individual_id(s):
	return (s.split('-'))[1]


def get_tissue_id(s):
	return (s.split('-'))[0]





if __name__ == "__main__":



	##==== factor matrix
	fm_tissue = np.load("./data/Tissue.npy")
	fm_individual = np.load("./data/Individual.npy")
	fm_gene = np.load("./data/Gene.npy")
	n_factor = len(fm_tissue[0])
	n_tissue = len(fm_tissue)
	n_individual = len(fm_individual)
	n_gene = len(fm_gene)
	dimension = (n_tissue, n_individual, n_gene)

	print dimension
	print n_factor

	Tensor = []
	for i in range(n_tissue):
		print i,
		matrix = np.load("./data/Tensor_tissue_" + str(i) + ".npy")
		print matrix.shape
		Tensor.append(matrix)
	Tensor = np.array(Tensor)
	print Tensor.shape







	######## SESSION I ########
	"""
	##=============
	##==== make the tensor incomplete
	##====	to save:
	##====		Tissue_list, Individual_list, Gene_list, Tensor_tissue_X.npy, Tensor_tissue_X_list.npy
	##====		list_sample.npy
	##====	where:
	##====		./data_init/
	##====	special NOTE:
	##====		sample ID is defined by "TissueID-IndividualID"
	##=============
	##
	list_tissue = []
	for i in range(n_tissue):
		list_tissue.append(str(i))
	list_tissue = np.array(list_tissue)
	np.save("./data_init/Tissue_list", list_tissue)
	##
	list_individual = []
	for i in range(n_individual):
		list_individual.append(str(i))
	list_individual = np.array(list_individual)
	np.save("./data_init/Individual_list", list_individual)
	##
	list_gene = []
	for i in range(n_gene):
		list_gene.append(str(i))
	list_gene = np.array(list_gene)
	np.save("./data_init/Gene_list", list_gene)


	##
	Data = []
	list_sample = []
	for i in range(n_tissue):

		list_random = np.random.permutation(n_individual)[:50]		# NOTE: we take half the individuals for each tissue

		list_sample_cur = []
		matrix_sample_cur = []
		for index in list_random:
			sample = Tensor[i][index]
			list_sample.append(str(i) + '-' + str(index))
			Data.append(sample)

			list_sample_cur.append(str(i) + '-' + str(index))
			matrix_sample_cur.append(sample)
		list_sample_cur = np.array(list_sample_cur)
		matrix_sample_cur = np.array(matrix_sample_cur)
		np.save("./data_init/Tensor_tissue_" + str(i) + "_list", list_sample_cur)
		np.save("./data_init/Tensor_tissue_" + str(i), matrix_sample_cur)

	Data = np.array(Data)
	list_sample = np.array(list_sample)
	np.save("./data_init/list_sample", list_sample)

	print "Data and list_sample shape:",
	print Data.shape,
	print list_sample.shape








	##=============
	##==== do PCA for Sample x Gene matrix
	##====	to save:
	##====		Gene.npy
	##=============
	print "performing PCA..."
	pca = PCA(n_components=n_factor)
	pca.fit(Data)
	Y2 = (pca.components_).T
	Y1 = pca.transform(Data)
	variance = pca.explained_variance_ratio_

	print variance
	print "and the cumulative variance are:"
	for i in range(len(variance)):
		print i,
		print np.sum(variance[:i+1]),
	print ""

	# DEBUG
	print "sample factor matrix:",
	print len(Y1),
	print len(Y1[0])
	print "gene factor matrix:",
	print len(Y2),
	print len(Y2[0])

	##==== save PCA results (two matrices, for coefficient matrix and factor matrix; and also the sample_list)
	#np.save("./data_processed/pca_sample", Y1)
	#np.save("./data_processed/pca_gene", Y2)
	#np.save("./data_processed/pca_variance", variance)

	np.save("./data_init/Gene", Y2)










	##=============
	##==== re-shape the Sample x Factor matrix back to tensor (incomplete), and prepare the data for incomplete PCA
	##==== save the Individual x Tissue matrix (with Nan in) under "./data_inter/"
	##=============
	rep_tissue_index = {}
	for i in range(len(list_tissue)):
		tissue = list_tissue[i]
		rep_tissue_index[tissue] = i
	rep_individual_index = {}
	for i in range(len(list_individual)):
		individual = list_individual[i]
		rep_individual_index[individual] = i

	Data = np.zeros((n_tissue, n_individual, n_factor))
	for i in range(n_tissue):
		for j in range(n_individual):
			for k in range(n_factor):
				Data[(i,j,k)] = float("Nan")

	print "shape of list_sample:",
	print list_sample.shape
	print "shape of Sample x Factor matrix:",
	print Y1.shape
	for i in range(len(list_sample)):
		sample = list_sample[i]
		tissue = get_tissue_id(sample)
		index_tissue = rep_tissue_index[tissue]
		individual = get_individual_id(sample)
		index_individual = rep_individual_index[individual]

		Data[index_tissue][index_individual] = Y1[i]
	print "the Tissue x Individual x Factor tensor has the dimension:",
	print Data.shape



	for k in range(n_factor):
		m_factor = np.zeros((n_tissue, n_individual))
		for i in range(n_tissue):
			for j in range(n_individual):
				m_factor[i][j] = Data[i][j][k]
		np.save("./data_inter/f" + str(k) + "_tissue_indiv", m_factor)
	print "save done..."


	##== need to save the results in tsv file (including Nan), in order to load in R
	for k in range(n_factor):
		m = np.load("./data_inter/f" + str(k) + "_tissue_indiv.npy")
		file = open("./data_inter/f" + str(k) + "_tissue_indiv.txt", 'w')
		for i in range(len(m)):
			for j in range(len(m[i])):
				value = m[i][j]
				file.write(str(value))
				if j != len(m[i])-1:
					file.write('\t')
			file.write('\n')
		file.close()
	"""











	######## SESSION II ########
	##=============
	##==== to save:
	##====		Tissue.npy, Individual.npy
	##=============
	factor_tissue = []
	factor_indiv = []
	for k in range(n_factor):
		#
		factor_tissue.append([])
		file = open("./data_inter/f" + str(k) + "_tissue.txt", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break
			factor_tissue[-1].append(float(line))
		file.close()

		#
		factor_indiv.append([])
		file = open("./data_inter/f" + str(k) + "_indiv.txt", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break
			factor_indiv[-1].append(float(line))
		file.close()

	factor_tissue = (np.array(factor_tissue)).T
	factor_indiv = (np.array(factor_indiv)).T

	print "factor tissue:",
	print factor_tissue.shape

	print "factor indiv:",
	print factor_indiv.shape


	np.save("./data_init/Tissue", factor_tissue)
	np.save("./data_init/Individual", factor_indiv)







