## initialize the matrix factor matrices
## targets:
##	Individual.npy		, directly
##	Gene.npy		, directly
##

## NOTE: this script will depend on the results from "pre_process.py", which picks up the samples, processes the genes, and normalizes the genes



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




n_factor = 200		# TODO: to decide




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
	list_sample = np.load("./data_raw/list_sample.npy")
	Data = np.load("./data_raw/Data.npy")



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
	np.save("./data_mf/Individual", Y1)
	np.save("./data_mf/Gene", Y2)
	np.save("./data_mf/pca_variance", variance)
	print "saving the .npy data done..."





