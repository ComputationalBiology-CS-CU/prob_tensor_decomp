## init the simulated tensor for sampling test (other than using the real factor matrices for init)

## we can keep the complete tensor, or make it incomplete

## the input is the output of the tensor simulator

## the final purpose of this script: we know the data is generated this way, we try to test whether init complete/incomplete tensor with the way we designed will give the sampling enough meaning




## NOTE: this script will work on full tensor





import numpy as np
from sklearn.decomposition import PCA




n_factor = 0
n_tissue = 0
n_individual = 0
n_gene = 0
dimension = (n_tissue, n_individual, n_gene)



if __name__ == "__main__":


	print "hello world..."
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

	Data = []
	for i in range(n_tissue):
		print i,
		matrix = np.load("./data/Tensor_tissue_" + str(i) + ".npy")
		print matrix.shape
		for i in range(len(matrix)):
			Data.append(matrix[i])
	Data = np.array(Data)	# NTOE: this is a Sample x Gene matrix
	print Data.shape



	##=============
	##==== do PCA for Sample x Gene matrix
	##=============
	## TODO: tunable
	#n_factor = 40

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
	##==== save the whole tensor
	##=============
	Tensor = []
	for i in range(n_tissue):
		Tensor.append([])
		for j in range(n_individual):
			index = i * n_individual + j
			Tensor[-1].append(Y1[index])
	Tensor = np.array(Tensor)
	print Tensor.shape



	##=============
	##==== factorize Tissue x Individual matrix for each factor
	##=============
	fm_tissue = []
	fm_individual = []
	for k in range(n_factor):
		print "factor#:",
		print k

		factor_m = np.zeros((n_tissue, n_individual))
		for i in range(n_tissue):
			for j in range(n_individual):
				factor_m[i][j] = Tensor[i][j][k]
		## pca
		pca = PCA(n_components=1)
		pca.fit(factor_m)
		Y2 = (pca.components_).T
		Y1 = pca.transform(factor_m)
		variance = pca.explained_variance_ratio_
		print variance
		loading_tissue = Y1.T[0]
		loading_individual = Y2.T[0]

		fm_tissue.append(loading_tissue)
		fm_individual.append(loading_individual)

	fm_tissue = np.array(fm_tissue).T
	fm_individual = np.array(fm_individual).T

	print "tissue and individual factor matrices shapes:"
	print fm_tissue.shape
	print fm_individual.shape

		


	np.save("./data_init/Tissue", fm_tissue)
	np.save("./data_init/Individual", fm_individual)





