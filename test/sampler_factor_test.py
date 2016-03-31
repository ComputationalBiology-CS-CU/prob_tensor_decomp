#This is the test code for sampler_factor
import numpy as np
from scipy.stats import wishart

#====================
#==Global variables==
#====================
n_factor = 30
n_v = 20000
n_t = 200
normalWishart = [[2,2],2,[[10,5],[5,10]],3]

def test_cr(v, t, alpha, lambda_u, precision_matrix):
	print "hello world"


def test_mw(v, t, alpha, lambda_u, precision_matrix):
	d1 = v.shape
	d2 = t.shape

	for i in range(d1[0]):
		for j in range(d2[0]):
			vt_vector = np.multiply(v[i], t[j])
			precision_matrix = np.add(precision_matrix, alpha * np.dot(np.array([vt_vector]).T, np.array([vt_vector])))
	return precision_matrix



def test_main():
	global n_factor
	global n_v
	global n_t


	#DEBUG
	print "Initiate precision_matrix"
	precision = []		# will be initialized with a diagonal 1 matrix
	for i in range(n_factor):
		precision.append([])
		for j in range(n_factor):
			if j == i:
				precision[i].append(1)
			else:
				precision[i].append(0)
	precision = np.array(precision)

	simulator_MVN()



#========================
#==The helper functions==
#========================

def simulator_MVN(n_sample, isTrans):
	global normalWishart

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

def sampler_W(df, scale):
	#
    sample = wishart.rvs(df, scale, size=1, random_state=None)
#    matrix = sample[0]
    return sample