#This is the test code for cal_product function in sampler
#Written by Chuqiao Ren
import numpy as np
from scipy.stats import wishart
import timeit

#======The global variables
n_factor = 10
fm = []
n_individual = 100
n_tissue = 20
n_gene = 200000
normalWishart = [[2,2],2,[[10,5],[5,10]],3]


#======The following function is written by Shuo Yang
def cal_product(i, j, k):
	global n_factor
	global fm

	product = 0
	for count in range(n_factor):
		product += fm[0][i][count] * fm[1][j][count] * fm[2][k][count]
	return product

#=======The following function is written by Chuqiao Ren
def cal_product_cr(i, j, k):
	global fm
	# print len(fm[0][i])
	return np.sum(np.multiply(np.multiply(fm[0][i], fm[1][j]),fm[2][k]))

def test_main():
	global n_individual
	global n_tissue
	global n_gene
	global fm

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

	print "Ready to simulate three matricies"
	#First simulate three matricies
	u = simulator_MVN(n_individual, False)
	print "Finished u"
	v = simulator_MVN(n_gene, False)
	print "Finished v"
	t = simulator_MVN(n_tissue, False)
	print "Finished t"

	print "Ready to construct the fm matrix"
	#Then calculate factor matricies
	fm = []
	fm.append(u)
	fm.append(v)
	fm.append(t)
	fm = np.array(fm)

	#DEBUG
	print "test Shuo's method"
	start_time_mw = timeit.default_timer()
	array_shuo = []
	for j in range(n_gene):
		for k in range(n_tissue):
				array_shuo.append(cal_product(1,j,k))
	elapsed_mw = timeit.default_timer() - start_time_mw

	print "test Chuqiao's method"
	start_time_cr = timeit.default_timer()
	array_cr = []
	for j in range(n_gene):
		for k in range(n_tissue):
			array_cr.append(cal_product_cr(1,j,k))
	elapsed_cr = timeit.default_timer() - start_time_cr

	isSame = array_equal(array_shuo, array_cr)

	print
	print "---------RESULTS SECTION FOR TENSOR--------------"

	print "Have we got the same answer? ", isSame
	print "-------------------------------------------------"
	print "Shuo's result is ", elapsed_mw
	print "-------------------------------------------------"
	print "Chuqiao's result is ", elapsed_cr





#=======The following methods are helper methods

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
			print "at ",i
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