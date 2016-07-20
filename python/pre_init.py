## initialize the tensor factor matrices
## targets:
##	Tissue.npy		, in-directly
##	Individual.npy		, in-directly
##	Gene.npy		, directly
##


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


n_factor = 40		# TODO: to decide
n_tissue = 0
n_individual = 0
n_gene = 0
dimension = (n_tissue, n_individual, n_gene)



##===============
##==== sub-routines
##===============
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



if __name__ == "__main__":



	######## SESSION I ########
	"""
	##=============
	##==== loading
	##=============
	list_tissue = np.load("./data_processed/Tissue_list.npy")
	list_individual = np.load("./data_processed/Individual_list.npy")
	list_gene = np.load("./data_processed/Gene_list.npy")
	list_sample = np.load("./data_raw/list_sample.npy")
	Data = np.load("./data_raw/Data.npy")
	## sample_tissue_rep
	file = open("./data_raw/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_sample_tissue_type", 'r')
	sample_tissue_rep = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		if len(line) < 3:
			print line
			continue

		sample = line[0]
		tissue = line[2]
		sample_tissue_rep[sample] = tissue
	file.close()



	## TODO: we probably need to subtract the Chr22 genes for convergence testing purpose
	''' will do the following in "pre_process.py"
	##================
	print "now we are picking up a subset of all the genes..."
	## to update the following;
	##list_gene = np.load("./data_processed/Gene_list.npy")
	##Data = np.load("./data_raw/Data.npy")
	file = open("./data_raw/gene_tss.txt", 'r')
	rep_gene_chr = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		chr = line[1]
		rep_gene_chr[gene] = chr
	file.close()

	list_gene_new = []
	rep_index_gene = {}
	for i in range(len(list_gene)):
		gene = list_gene[i]
		if rep_gene_chr[gene] == '22':
			list_gene_new.append(gene)
			rep_index_gene[i] = 1
	list_gene = np.array(list_gene_new)
	print "length of new gene list:",
	print len(list_gene)
	np.save("./data_processed/Gene_list.npy", list_gene)
	Data_new = []
	for i in range(len(Data)):
		Data_new.append([])
		for j in range(len(Data[i])):
			if j in rep_index_gene:
				value = Data[i][j]
				Data_new[-1].append(value)
	Data = np.array(Data_new)
	print "shape of new Data matrix (Sample x Gene):",
	print Data.shape
	##================
	'''


	n_tissue = len(list_tissue)
	n_individual = len(list_individual)
	n_gene = len(list_gene)
	dimension = (n_tissue, n_individual, n_gene)

	print dimension





	##=============
	##==== do PCA for Sample x Gene matrix
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
	print "saving the .npy data done..."

	np.save("./data_processed/Gene", Y2)



	##=============
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
	for i in range(len(list_sample)):
		sample = list_sample[i]
		tissue = sample_tissue_rep[sample]
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


	np.save("./data_processed/Tissue", factor_tissue)
	np.save("./data_processed/Individual", factor_indiv)









