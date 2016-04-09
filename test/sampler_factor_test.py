#This is the test code for sampler_factor
#Written by Chuqiao Ren and Mengqing Wang
import numpy as np
from scipy.stats import wishart
from numpy.linalg import inv
import timeit

#====================
#==Global variables==
#====================
n_factor = 10
n_v = 2000
n_t = 20
n_u = 100
alpha = 0.5
normalWishart = [[2,2],2,[[10,5],[5,10]],3]
testTensor = True

def test_cr(v, t, lambda_u, precision_matrix):
	global alpha
	global n_factor
	dimension1 = n_v
	dimension2 = n_t
	Q = []

	for i in range(dimension1):
		for j in range(dimension2):
			vt_vector = np.multiply(v[i], t[j])
			Q.append(vt_vector)

	Q = np.array(Q)
	precision_matrix = np.add(precision_matrix, np.multiply(alpha, np.dot(Q.T, Q)))
	return precision_matrix

def test_cr_tensor(v, t, lambda_u, precision_matrix, mean, product, u_i):
	global alpha
	global n_factor
	dimension1 = n_v
	dimension2 = n_t
	Q = []

	for i in range(dimension1):
		for j in range(dimension2):
			vt_vector = np.multiply(v[i], t[j])
			Q.append(vt_vector)
			mean = np.add(mean,np.multiply(alpha,np.multiply(vt_vector , product[i][u_i][j])))

	Q = np.array(Q)

	precision_matrix = np.add(precision_matrix, np.multiply(alpha, np.dot(Q.T, Q)))

	cov = inv(precision_matrix)

	mean = np.dot(mean, cov.T)

	return mean, precision_matrix

def test_mw_tensor(v, t, lambda_u, precision_matrix, mean, product, u_i):
	global alpha
	d1 = n_v
	d2 = n_t

	for i in range(d1):
		for j in range(d2):
			vt_vector = np.multiply(v[i], t[j])
			precision_matrix = np.add(precision_matrix, np.multiply(alpha, np.dot(np.array([vt_vector]).T, np.array([vt_vector]))))
			mean = np.add(alpha * np.multiply(product[(i, u_i, j)], vt_vector), mean)

	cov = inv(precision_matrix)

	mean = np.dot(mean, cov.T)

	return mean, precision_matrix


def test_mw(v, t,  lambda_u, precision_matrix):
	global alpha
	d1 = n_v
	d2 = n_t

	for i in range(d1):
		for j in range(d2):
			vt_vector = np.multiply(v[i], t[j])
			precision_matrix = np.add(precision_matrix, np.multiply(alpha, np.dot(np.array([vt_vector]).T, np.array([vt_vector]))))

	return precision_matrix



def test_main():
	global n_factor
	global n_v
	global n_t
	global n_u
	global normalWishart
	global testTensor

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
	v = simulator_MVN(n_v, True)
	t = simulator_MVN(n_t, True)

	if not testTensor:
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

	if testTensor:
		#DEBUG
		print "simulate U"
		u = simulator_MVN(n_u, True)

		#DEBUG
		print "tensor multiplication..."

		vs = [v[0], u[0], t[0]]
		product = reduce(np.multiply, np.ix_(*vs))

		for i in xrange(1,n_factor):
			vs = [v[i], u[i], t[i]]
			product += reduce(np.multiply, np.ix_(*vs))

		print product.shape
		v = v.T
		u = u.T
		t = t.T

		#DEBUG
		print "initialize mean array"
		mean = [1] * n_factor
		mean = np.array(mean)

		#DEBUG
		print "test Mengqing's method"
		start_time_mw = timeit.default_timer()
		mw = test_mw_tensor(v, t, lambda_U, precision_matrix, mean, product, 1)
		elapsed_mw = timeit.default_timer() - start_time_mw

		print "test Chuqiao's method"
		start_time_cr = timeit.default_timer()
		cr = test_cr_tensor(v, t, lambda_U, precision_matrix, mean, product, 1)
		elapsed_cr = timeit.default_timer() - start_time_cr

		isSame = matrix_equal(mw[1], cr[1]) and array_equal(mw[0], cr[0])

		print
		print "---------RESULTS SECTION FOR TENSOR--------------"

		print "Have we got the same answer? ", isSame
		print "-------------------------------------------------"
		print "Mengqing's result is ", elapsed_mw
		print "-------------------------------------------------"
		print "Chuqiao's result is ", elapsed_cr


#========================
#==The helper functions==
#========================

def matrix_equal(m1, m2):
	if (m1.shape != m2.shape):
		print "Matrix shape is not the same!"
		return False

	for i in range(m1.shape[0]):
		for j in range(m1.shape[1]):
			if abs(m1[i][j] - m2[i][j]) > 0.01:
				print m1[i][j]
				print m2[i][j]
				print "Matrix Entry ERROR"
				return False

	return True

def array_equal(a1, a2):
	if (len(a1) != len(a2)):
		print "Array shape error!"
		print "Mengqing has ", len(a1)
		print "Chuqiao has ", len(a2)
		return False

	for i in range(len(a1)):
		if abs(a1[i] - a2[i]) > 0.01:
			print a1[i]
			print a2[i]
			print "Array Entry ERROR"
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