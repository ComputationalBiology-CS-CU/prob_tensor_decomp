#This is the test code for sampler_factor
import numpy as np
from scipy.stats import wishart
import timeit

#====================
#==Global variables==
#====================
n_factor = 10
n_v = 200
n_t = 20
alpha = 0.5
normalWishart = [[2,2],2,[[10,5],[5,10]],3]

def test_cr(v, t, lambda_u, precision_matrix):
	global alpha
	dimension1 = n_v
	dimension2 = n_t

	Q = []

	for i in range(dimension1):
		for j in range(dimension2):
			vt_vector = np.multiply(v[i], t[j])
			Q.append(vt_vector)

	Q = np.array(Q)
	precision_matrix = np.add(precision_matrix,alpha * np.dot(Q.T, Q))
	return precision_matrix


def test_mw(v, t,  lambda_u, precision_matrix):
	global alpha
	d1 = n_v
	d2 = n_t

	for i in range(d1):
		for j in range(d2):
			vt_vector = np.multiply(v[i], t[j])
			precision_matrix = np.add(precision_matrix, alpha * np.dot(np.array([vt_vector]).T, np.array([vt_vector])))

	return precision_matrix



def test_main():
	global n_factor
	global n_v
	global n_t
	global normalWishart

	#DEBUG
	print "Initiate normal Wishart prior"
	mu = []
	precision = []
	for i in range(n_factor):
		mu.append(2)
	for i in range(n_factor):
		temp = []
		for j in range(n_factor):
			if (i == j):
				temp.append(1)
			else:
				temp.append(0)
		precision.append(temp)

	normalWishart[0] = mu
	normalWishart[2] = precision
	normalWishart[3] = n_factor + 1


	#DEBUG
	print "Initiate precision_matrix and lambda_U"
	precision_matrix = []		# will be initialized with a diagonal 1 matrix
	lambda_U = []
	for i in range(n_factor):
		precision_matrix.append([])
		lambda_U.append([])
		for j in range(n_factor):
			if j == i:
				precision_matrix[i].append(1)
				lambda_U[i].append(1)
			else:
				precision_matrix[i].append(0)
				lambda_U[i].append(0)
	precision_matrix = np.array(precision_matrix)
	lambda_U = np.array(lambda_U)


	#DEBUG
	print "simulate V and T"
	v = simulator_MVN(n_v, False)
	t = simulator_MVN(n_t, False)

	#DEBUG
	print "test Mengqing's method"
	start_time_mw = timeit.default_timer()
	mw = test_mw(v, t, lambda_U, precision_matrix)
	elapsed_mw = timeit.default_timer() - start_time_mw

	print "test Chuqiao's method"
	start_time_cr = timeit.default_timer()
	cr = test_cr(v, t, lambda_U,precision_matrix)
	elapsed_cr = timeit.default_timer() - start_time_cr

	isSame = matrix_equal(mw, cr)

	print
	print "---------RESULTS SECTION--------------"

	print "Have we got the same answer? ", isSame
	print "--------------------------------------"
	print "Mengqing's result is ", elapsed_mw
	print "--------------------------------------"
	print "Chuqiao's result is ", elapsed_cr


#========================
#==The helper functions==
#========================

def matrix_equal(m1, m2):
	if (m1.shape != m2.shape):
		print "Shape is not the same!"
		return False

	for i in range(m1.shape[0]):
		for j in range(m1.shape[1]):
			if abs(m1[i][j] - m2[i][j]) > 0.01:
				print m1[i][j]
				print m2[i][j]
				print "Entry ERROR"
				return False

	return True

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
		return np.array(u)

def sampler_W(df, scale):
	#
    sample = wishart.rvs(df, scale, size=1, random_state=None)
#    matrix = sample[0]
    return sample

if __name__ == '__main__':
	test_main()