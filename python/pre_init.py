## initialize the tensor factor matrices


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


n_factor = 40
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


	##=============
	##==== loading
	##=============
	list_tissue = np.load("./data_processed/Tissue_list.npy")
	list_individual = np.load("./data_processed/Individual_list.npy")
	list_gene = np.load("./data_processed/Gene_list.npy")
	list_sample = np.load("./data_raw/list_sample.npy")
	Data = np.load("./data_processed/Data.npy")
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
	n_tissue = len(list_tissue)
	n_individual = len(list_individual)
	n_gene = len(list_gene)

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

	# DEBUG
	print "sample factor matrix:",
	print len(Y1),
	print len(Y1[0])
	print "gene factor matrix:",
	print len(Y2),
	print len(Y2[0])

	##==== save PCA results (two matrices, for coefficient matrix and factor matrix; and also the sample_list)
	np.save("./data_processed/pca_sample", Y1)
	np.save("./data_processed/pca_gene", Y2)
	np.save("./data_processed/pca_variance", variance)
	print "saving the .npy data done..."





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






