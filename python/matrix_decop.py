# This is a script to do matrix decomposition of simulated gene data
# Assume gene data has two dimension: gene and individual
# Gene has spike and slab sparsity prior while individual has MVN

# Written by Chuqiao Ren 

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
import sys
from copy import *
import timeit


##=====================
##==== global variables
##=====================
n_factor = 10
n_individual = 40
n_gene = 50
dimension = (n_individual, n_gene)
factor_name = {}
dataset = np.zeros(shape=dimension) # individual x gene x tissue
markerset = np.ones(shape=dimension) # mark the position where there are data
individual_rep = {}
##==== mapping: #1: individual; #2: gene; #3: tissue
fmlist = []
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

time = []

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

def cal_product(i, j):
    global n_factor
    global fmlist

    product = 0
    for count in range(n_factor):
        product += fmlist[0][i][count] * fmlist[1][j][count]

    return product

def cal_product_k(i, j, k):
    global n_factor
    global fmlist

    return fmlist[0][i][k] * fmlist[1][j][k]

def factor_combination(dataset, fmlist, dim_depth, n_factor, fm_depth, prod, path):
    if fm_depth==len(fmlist):
        ind_tup = tuple(path)
        dataset[ind_tup] = np.sum(prod)
        return

    for i in range(dim_depth[fm_depth]):
        new_prod = np.multiply(prod, fmlist[fm_depth][i])
        path.append(i)
        factor_combination(dataset, fmlist, dim_depth, n_factor, fm_depth+1, new_prod, path)
        path.pop()

##==== this function will load
def load_dataset():
    global dataset
    global fmlist
    global dimension
    global n_factor
    # load fmlist from simulated data
    fmlist.append(np.load('./data/Individual.npy'))
    fmlist.append(np.load('./data/Gene.npy'))
    prod = np.ones(n_factor)
    factor_combination(dataset, fmlist, dimension, n_factor, 0, prod, [])

##==== Deep copy of gene matrix (used when calculate V~nk)
def deep_copy(v):
    v_new = []
    for n in range(len(v)):
        v_row = []
        for k in range(len(v[0])):
            v_row.append(v[n][k])
        v_new.append(v_row)

    return v_new


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
    return sample


##==== sampling from Gamma, for the variance
def sampler_Gamma(para1, para2):
    para2 = 1.0/para2
    x = np.random.gamma(para1, para2, 1)
    return x[0]

## ==== sampling Beta
def sampler_beta(a, b):
    return np.random.beta(a, b)


## ==== sampling Bernoulli
def sampler_bernoulli(p, size):
    return bernoulli.rvs(p, loc=0, size=size)


## ==== sampling student-t
def sampler_student_t(df, loc, scale):
    return t.rvs(df, loc = loc, scale = scale)

## ==== related function to sampler z

def norm_pdf(mean, variance, p):
    #This function is to calculate the pdf of a normal distribution given mean and variance at point p
    out = 1/(math.sqrt(variance*2*math.pi))
    return out*math.exp(-(p - mean)**2/(2 * variance))


def cal_pR_0(n, k):
    #TODO: Need to optimize here!
    #This function is to calculate the probability of z_0
    global markerset
    global alpha
    global n_factor

    # print "Deep copy v"
    # result = 1
    # v = deep_copy(fmlist[num_factor])
    # v[n][k] = 0

    result = 0

    #print "Ready to calculate z_0"
    for i in range(n_individual):
        if (markerset[i][n] != 0):
            #Assuming each R_{ij}^{k} is independent to others
            product = 0
            for m in range(n_factor):
                if m != k:
                    product += cal_product_k(i,n,m)


            im = norm_pdf(product, alpha**(-1), dataset[i][n])
            if im == 0:
                # print "i, j: ", i, j
                # print "im: ", im
                # print "n, k and num_factor :", n, k, num_factor
                # print "product: ",product
                # print "sigma^2: ",alpha**(-1)
                # print "dataset: ",dataset[i][j]
                continue
            else:
                result += math.log(im)
                # print "im is: ", im
                # print "result is: ", result
        # print "result: ", result


    # print "result: ", result
    # print "exp(result): ", math.exp(result)
    return result

def cal_pR_1(n,k):
    #This function is to calculate the probability of z_1
    global alpha
    global dataset
    global fmlist
    global markerset
    global sparsity_prior
    global n_factor
    global n_gene
    global n_individual

    #DEBUG
    #print 'Ready to calculate z_1'

    sigma2 = math.pow(alpha, -1)
    const = math.sqrt(2 * math.pi)
    C = 0
    frac = 1
    sum_frac_mu = 0
    sum_frac_mu2 = 0
    sum_frac_sigma2 = 0
    prod_sigma2 = 1
    sum_sigma2 = 0
    count = 0
    for m in range(n_individual):
        ut = fmlist[0][m][k]
        s = cal_product(m,n) - ut*fmlist[1][n][k]   #Note: cal_product takes individual * gene * tissue
        if (ut == 0):
            #Then update the C term
            int_temp1 = np.power((s - dataset[m][n]),2)/(2 * sigma2)
            temp_C = np.exp(0 - int_temp1)/(const * sigma2)
            C += math.log(temp_C)
        else:
            #First update the fraction 1/({u_{mk}t_{pk})^2
            frac += math.log(1/(ut * ut))

            #And calculate mu and sigma for this term
            mu_i = (dataset[m][n] - s)/ut
            sigma_i2 = sigma2/(ut * ut)

            #Then update overall term
            sum_frac_mu += mu_i / sigma_i2
            sum_frac_mu2 += mu_i * mu_i / sigma_i2
            sum_frac_sigma2 += 1 / sigma_i2

            #Then update the product of sigma^2
            # prod_sigma2 *= sigma_i2
            sum_sigma2 += math.log(sigma_i2)

            #update count
            count += 1

    sigma2_s = 1/(sum_frac_sigma2 + sparsity_prior[1][k][k])
    mu_s = (sum_frac_mu + sparsity_prior[0][k]*sparsity_prior[1][k][k]) * sigma2_s
    #Note: sparsity_prior[1] contains the precision matrix (sigma^2 = alpha^{-1}). The precision matrix is a diagonal matrix

    #Ready to calculate the final term
    # outter = math.pow(2 * math.pi, count/2)
    # outter2 = math.sqrt(sigma2_s*sparsity_prior[1][k][k]/prod_sigma2)
    # inner = sum_frac_mu2 + math.pow(sparsity_prior[0][k],2)*sparsity_prior[1][k][k] - mu_s * mu_s / sigma2_s
    # S = outter2 * math.exp(0 - inner/2)/outter
    outter = 0.5 * (math.log(sigma2_s) + math.log(sparsity_prior[1][k][k]) - sum_sigma2)
    S = 0.5*(sum_frac_mu2 + math.pow(sparsity_prior[0][k],2)*sparsity_prior[1][k][k] - mu_s * mu_s / sigma2_s)
    result = frac + C + math.log(sparsity_prior[3][k]) - count/2*math.log(2 * math.pi) + outter - S

    # print "sigma2_s", sigma2_s
    # print "mu_s", mu_s
    #
    # print "frac", frac
    # print "C", C
    # print "S", S
    # print "result", math.exp(result)
    # print "result: ", result
    return result

def cal_rest(n):
    #TODO: Need to optimize here!
    #This function is to calculate the probability of z_0
    global markerset
    global alpha
    global n_factor

    # print "Deep copy v"
    # result = 1
    # v = deep_copy(fmlist[num_factor])
    # v[n][k] = 0

    result = 0

    #print "Ready to calculate z_1 the rest"
    for i in range(n_individual):
        for j in range(n_gene):
            if (markerset[i][j] != 0 and j != n):
                #Assuming each R_{ij}^{k} is independent to others
                product = 0
                for m in range(n_factor):
                    product += cal_product_k(i,j,m)

                im = norm_pdf(product, alpha**(-1), dataset[i][j])
                if im == 0:
                    # print "i, j: ", i, j
                    # print "im: ", im
                    # print "n, k and num_factor :", n, k, num_factor
                    # print "product: ",product
                    # print "sigma^2: ",alpha**(-1)
                    # print "dataset: ",dataset[i][j]
                    continue
                else:
                    result += math.log(im)
                    # print "im is: ", im
                    # print "result is: ", result

        # print "result: ", result

    # print "result: ", result
    # print "exp(result): ", math.exp(result)
    return result


##==== sample factor relevant components (the factor matrix, and its Gaussian-Wishart prior)

def sampler_factor_sparsity():

    global dataset
    global fmlist
    global prior
    global hyper_prior
    global n_factor
    global alpha
    global sparsity_prior
    global sparsity_hyper_prior
    global time

    #Only gene dimension will have sparsity prior

    dimension = n_gene
    dimension1 = n_individual
    num_factor = 1
    num_factor_1 = 0

    # DEBUG
    print "we'll sample factor matrix first..."

    start_z = timeit.default_timer()
    #DEBUG
    print "preparing to calculate z_0 and z_1..."
    p_0 = []     #A n-by-k matrix contains all the probability of z_nk = 0
    p_1 = []     #A n-by-k matrix contains all the probability of z_nk = 1

    for n in range(dimension):
        p_k0 = []
        p_k1 = []
        rest = cal_rest(n)
        for k in range(n_factor):
            k0 = math.log(1 - sparsity_prior[3][k]) + cal_pR_0(n, k) + rest
            p_k0.append(k0)
            k1 = math.log(sparsity_prior[3][k]) + cal_pR_1(n, k) + rest
            p_k1.append(k1)

        print "p_k0 is ", p_k0
        print "p_k1 is ", p_k1

        p_0.append(p_k0)
        p_1.append(p_k1)


    end_z = timeit.default_timer() - start_z
    start_z2 = timeit.default_timer()

    print "now sampling z"
    #Note: bernoulli.rvs(p) takes p_1 (prob of getting a 1)
    z = []
    for i in range(dimension):
        z_row = []
        for k in range(n_factor):
            #sampler_Bernoulli will return np.array object
            logp1 = p_1[i][k]
            logp0 = p_0[i][k]
            diff = logp1 - logp0
            if diff < math.log(sys.float_info.max):
                ratio = math.exp(diff)
                z_row.append(sampler_bernoulli(1/(1 + ratio), size=1)[0]) #loc = 0 by default in the sampler_bernoulli function
            else:
                z_row.append(sampler_bernoulli(0, size = 1)[0])
        z.append(z_row)

    sparsity_prior[2] = np.array(z)
    print "z"
    print sparsity_prior[2]

    end_z2 = timeit.default_timer() - start_z2

    #DEBUG
    print "now sampling V"

    start_v = timeit.default_timer()

    #==== sample factor matrix
    for dim in range(dimension):
        # DEBUG
        print "now sampling factor#",
        print dim+1,
        print "out of",
        print dimension

        #First pick out those with z == 1
        '''
        factor_ind = []
        z = sparsity_prior[2]
        for i in range(len(z[dim])):
            if (z[dim][i] == 1):
                factor_ind.append(i)
        '''
        factor_ind = []
        factor_u = []
        sparsity_mean = []
        u = np.array(fmlist[num_factor_1])
        for i in range(len(z[dim])):
            if (z[dim][i] == 1):
                factor_u.append(u.T[i])
                factor_ind.append(i)
                sparsity_mean.append(sparsity_prior[0][i])
        factor_u = np.array(factor_u)
        factor_u = factor_u.T



        #== precision_matrix
        precision_matrix = []
        for count1 in range(len(factor_ind)):
            count1_new = factor_ind[count1]
            temp = []
            for count2 in range(len(factor_ind)):
                count2_new = factor_ind[count2]
                value = sparsity_prior[1][count1_new][count2_new]	# initialize with the prior precision matrix, then later update
                temp.append(value)
            precision_matrix.append(temp)

        precision_matrix = np.array(precision_matrix)
        mean = np.dot(sparsity_mean, precision_matrix)
        Q = []

        for j in range(dimension1):
            # re-arrange the three dimension, for querying original dataset
            # hash_temp = {num_factor: dim, num_factor_1: j}
            # index1 = hash_temp[0]
            # index2 = hash_temp[1]
            # index3 = hash_temp[2]
            if markerset[j][dim] == 0:
                continue
            array = factor_u[j]
            Q.append(array)
            mean = np.add(mean, np.multiply(alpha, np.multiply(array, dataset[j][dim])))

        Q = np.array(Q)
        precision_matrix = np.add(precision_matrix, np.multiply(alpha, np.dot(Q.T, Q)))
        cov = inv(precision_matrix)
        # mean = np.dot(mean, cov.T)
        mean = np.dot(mean, cov)

        #== sampling
        v_1 = sampler_MVN(mean, cov)
        count = 0
        re = []
        for i in range(n_factor):
            if i in factor_ind:
                re.append(v_1[count])
                count += 1
            else:
                re.append(0)

        fmlist[num_factor][dim] = re

    end_v = timeit.default_timer() - start_v

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

    start_ng = timeit.default_timer()

    #DEBUG
    print "factor matrix - gene"
    print fmlist[1]

    # DEBUG
    print "preparing and sampling Normal-Gamma..."
    precision = []
    mean = []
    for k in range(n_factor):
        v_n = []  #This is 1-dimensional vector V
        for n in range(dimension):
            if sparsity_prior[2][n][k] != 0:
                v_n.append(fmlist[num_factor][n][k])
        d = len(v_n)

        #Now calcualte \bar{V}
        v_bar = sum(v_n)/d

        #Calcualte alpha_n
        alpha_n = sparsity_hyper_prior[0] + d/2

        #Calculate beta_n

        diff_v2 = []  #This is (V_n - \bar{V})^2
        for n in range(d):
            diff_v2.append((v_n[n] - v_bar)**2)

        beta_n = sparsity_hyper_prior[1] + sum(diff_v2)/2 + (sparsity_hyper_prior[3]*d*(v_bar - sparsity_hyper_prior[2])**2)/(2*(sparsity_hyper_prior[3] + d))

        #Calculate kappa_n

        kappa_n = sparsity_hyper_prior[3] + d

        #Calculate mu_n

        mu_n = (sparsity_hyper_prior[3] * sparsity_hyper_prior[2] + d * v_bar)/(kappa_n)

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

    end_ng = timeit.default_timer() - start_ng

    start_pi = timeit.default_timer()
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

    end_pi = timeit.default_timer() - start_pi

    time.append([end_z, end_z2,end_v, end_ng, end_pi])
    return

def sampler_factor(factor_id):
    global dataset
    global fmlist
    global prior
    global hyper_prior
    global n_factor
    global alpha

    # We define that
    # 	dimension 1: tissue
    # 	dimension 2: individual
    # 	dimension 3: gene


    cur_dimension = dimension[factor_id]


    # DEBUG
    print "we'll sample factor matrix first..."


    #==== sample factor matrix
    for i in range(cur_dimension):	# there are n factor array, that are independent with each other --> parallel


        # DEBUG
        #print "now sampling factor#",
        #print i+1,
        #print "out of",
        #print cur_dimension

        precision_matrix = prior[factor_id][1]
        ids = [0,1]
        ids.remove(factor_id)
        dimension1 = dimension[ids[0]]
        #dimension2 = dimension[ids[1]]

        #prod = np.ones(n_factor)
        #sampler_factor_helper(dataset, markerset, fmlist, dimension, n_factor, prod, [], factor_id, precision_matrix)

        Q = []
        for j in range(dimension1):
            hash_temp = {factor_id: i, ids[0]: j}
            index1 = hash_temp[0]
            index2 = hash_temp[1]
            if markerset[(index1, index2)] == 0:
                continue
            #array = np.multiply(fmlist[ids[0]][j], fmlist[ids[1]][k])
            #precision_matrix = np.add(precision_matrix, alpha*np.dot(np.array([fmlist[ids[0]][j]]).T, np.array([fmlist[ids[0]][j]])))

            #vt_vector = np.multiply(v[i], t[j])
            Q.append(fmlist[ids[0]][j])

        Q = np.array(Q)
        precision_matrix = np.add(precision_matrix, np.multiply(alpha, np.dot(Q.T, Q)))
        cov = inv(precision_matrix)
        mean = np.dot(prior[factor_id][0], prior[factor_id][1])


        for j in range(dimension1):
            hash_temp = {factor_id: i, ids[0]: j}
            index1 = hash_temp[0]
            index2 = hash_temp[1]
            if markerset[index1][index2] == 0:
                continue

            #array = np.multiply(fmlist[ids[0]][j], fmlist[ids[1]][k])
            mean = np.add(alpha * dataset[(index1, index2)] * fmlist[ids[0]][j], mean)

        mean = np.dot(mean, cov)

        #== sampling

        fmlist[factor_id][i] = sampler_MVN(mean, cov)


    print fmlist[factor_id]

    # DEBUG
    print "now we are sampling the prior..."

    #==== sample Gaussian-Wishart (Wishart first, then Gaussian) prior
    #== sample Wishart first

    # DEBUG
    print "preparing and sampling Wishart first..."

    # factor_mean
    factor_mean = np.average(fmlist[factor_id], axis=0)

    # factor_var
    factor_var = np.zeros((n_factor, n_factor))
    for i in range(n_individual):
        factor_var += np.dot((fmlist[factor_id][i] - factor_mean).T, (fmlist[factor_id][i] - factor_mean))

    # cov_matrix
    cov_matrix = inv(hyper_prior[factor_id][0])
    cov_matrix = cov_matrix + factor_var
    temp = hyper_prior[factor_id][3] * n_individual / (hyper_prior[factor_id][3] + n_individual)
    cov_matrix = cov_matrix + temp * np.dot((hyper_prior[factor_id][2] - factor_mean).T,
                                            (hyper_prior[factor_id][2] - factor_mean))

    precision_matrix = inv(cov_matrix)

    # df new
    df = hyper_prior[factor_id][1] + n_individual

    ## sampling Wishart
    prior[factor_id][1] = np.array(sampler_W(df, precision_matrix))

    # DEBUG
    print "now sampling MVN then..."

    # == sample Gaussian then
    # beta new
    beta = hyper_prior[factor_id][3] + n_individual
    precision_matrix = np.multiply(beta,prior[factor_id][1])
    # print precision_matrix
    cov = inv(precision_matrix)

    # mean
    mean = (hyper_prior[factor_id][3] * hyper_prior[factor_id][2] + n_individual * factor_mean) / (
    hyper_prior[factor_id][3] + n_individual)

    # sampling MVN
    prior[factor_id][0] = sampler_MVN(mean, cov)

    #return

##==== sampling precision
def sampler_precision():
    global alpha_prior
    global n_individual
    global n_gene
    global dataset
    global markerset
    global alpha

    para1_old = alpha_prior[0]
    para2_old = alpha_prior[1]
    n = 0
    sum = 0

    for i in range(n_individual):
        for j in range(n_gene):
            if markerset[i][j] == 0:
                continue
            #
            R_real = dataset[i][j]
            R_exp = cal_product(i, j)
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

##==== log likelihood for Normal Gamma without constant terms
def loglike_Normal_Gamma(obs1, obs2, mu, kappa, alpha, beta):
    #Note: obs1 == mean and obs2 == precision (univariate)
    like_log = 0
    like_log += (alpha - 0.5) * math.log(obs2)
    like_log += 0 - obs2 * 0.5 * (kappa * (obs1 - mu) * (obs1 - mu) + 2 * beta)

    return like_log

##==== log likelihood for Beta without constant terms
def loglike_Beta(obs, alpha, beta):
    like_log = 0
    like_log += (alpha - 1) * math.log(obs) + (beta - 1) * math.log(1 - obs)
    return like_log

def loglike_joint_sparsity():
        global n_individual, n_gene, n_factor
        global dataset, markerset
        global prior  # prior[0, 1, 2]: 0: mean array; 1: precision matrix
        global hyper_prior  # hyper_prior[0, 1, 2]: 0: scale; 1: df; 2: mean; 3: scaler
        global alpha
        global alpha_prior  # 0: shape parameter; 1: rate parameter
        global sparsity_prior  # 0: mean; 1: precision; 2: z; 3: pi
        global sparsity_hyper_prior  # 0: alpha_0; 1: beta_0; 2: mu_0; 3: kappa_0; 4: alpha_pi; 5: gamma_pi

        like_log = 0

        # ==== edml
        for i in range(n_individual):
            for j in range(n_gene):
                if markerset[i][j] == 0:
                    continue

                obs = dataset[i][j]
                mean = cal_product(i, j)
                var = 1.0 / alpha
                like_log += loglike_Gaussian(obs, mean, var)

        # ==== factor matrix likelihood
        for i in range(n_individual):
            like_log += loglike_MVN(fmlist[0][i], prior[0][0], prior[0][1])
        for j in range(n_gene):
            like_log += loglike_MVN(fmlist[1][j], sparsity_prior[0], sparsity_prior[1])

        # ==== factor prior likelihood for individual
        like_log += loglike_GW(prior[0][0], prior[0][1], hyper_prior[0][0], hyper_prior[0][1], hyper_prior[0][2], hyper_prior[0][3])


        # ==== factor prior likelihood for gene
        for f in range(n_factor):
            # ==== Normal Gamma
            like_log += loglike_Normal_Gamma(sparsity_prior[0][f], sparsity_prior[1][f][f], sparsity_hyper_prior[2], sparsity_hyper_prior[3], sparsity_hyper_prior[0], sparsity_hyper_prior[1])
            # ==== Beta
            like_log += loglike_Beta(sparsity_prior[3][f], sparsity_hyper_prior[4], sparsity_hyper_prior[5])

        #TODO: Ask how if we need to include log-likelihood of z

        # ==== precision/variance likelihood
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

    # prepare the "dataset" and "markerset"
    # data_prepare()
    load_dataset()

    # DEBUG
    print "finish data preparation..."
    print "gene is ",
    print fmlist[1]


    # DEBUG
    print "now initializing all the variables..."


    ##================================
    ##==== initialize global variables
    ##================================
    #n_factor = 40			# TODO: this is tunable, and the number 400 comes from results of other more complex methods

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

        scale = np.identity(n_factor)
        hyper_prior[n].append(np.load("./data/precision.npy"))    # lambda
        hyper_prior[n].append(np.int(np.load("./data/v.npy")))		# TODO: tunable   v_0
        hyper_prior[n].append(np.load("./data/mu.npy"))		# TODO: tunable	 mu_0
        hyper_prior[n].append(np.int(np.load("./data/kappa.npy")))		# TODO: tunable  kappa_0



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
        precision = np.identity(n_factor)
        prior[n].append(deepcopy(mean))
        prior[n].append(deepcopy(precision))

    #== the prior of sparsity
    sparsity_prior = []  # 0: mean; 1: precision; 2: z; 3: pi
    mean = np.array([0] * n_factor)
    precision = np.identity(n_factor)
    sparsity_prior.append(deepcopy(mean))
    sparsity_prior.append(deepcopy(precision))
    sparsity_prior.append(np.load("./data/z.npy"))
    sparsity_prior.append(np.load("./data/pi_array.npy"))

    sparsity_hyper_prior = []  # 0: alpha_0; 1: beta_0; 2: mu_0; 3: kappa_0; 4: alpha_pi; 5: gamma_pi
    sparsity_hyper_prior.append(np.int(np.load("./data/alpha_NG.npy")))
    sparsity_hyper_prior.append(np.int(np.load("./data/beta_NG.npy")))
    sparsity_hyper_prior.append(np.int(np.load("./data/mu_NG.npy")))
    sparsity_hyper_prior.append(np.int(np.load("./data/kappa_NG.npy")))
    sparsity_hyper_prior.append(np.int(np.load("./data/beta_alpha.npy")))
    sparsity_hyper_prior.append(np.int(np.load("./data/beta_beta.npy")))

    #== the MVN drawing (mean 0, cov 1) for factorized matrices
    mean = np.array([0] * n_factor)

    cov = np.identity(n_factor)

    #== set the prior for precision (Gamma)
    alpha_prior = [1, 0.5]		# shape parameter and rate parameter, TODO: tunable


    #== drawing precision from Gaussian
    # alpha = sampler_Normal(0, 1)
    # while (alpha <= 0):
    #     alpha = sampler_Normal(0,1)
    alpha = 0.5


    # DEBUG
    print "finish variables initialization..."



    ##==============================
    ##==== sampler calling iteration
    ##==============================
    ITER = 100
    # ll_result = []
    fo = open("./result/test.txt", "w+")
    f_time = open("./result/time.txt", "w+")
    for i in range(ITER):
        time = []
        print "current iteration#",
        print i+1
        print "start to sample individual..."
        sampler_factor(0)
        print "start to sample gene..."
        sampler_factor_sparsity()
        print "sample precision..."
        sampler_precision()
        like_log = loglike_joint_sparsity()	# monitor the log joint likelihood
        #print "sampling done. the log joint likelihood is",
        print like_log

        turn = i + 1
        fo.write(str(turn) + ": " + str(like_log) + "\n")
        f_time.write(str(turn) + ": " + str(time) + "\n")
        # ll_result.append(like_log)

    # for t in time:
    #     line = str(t) + ": "
    #     for i in t:
    #         line+=str(i) + ","
    #     line = line[:-1]
    #     line += "\n"
    #     f_time.write(line)

    f_time.close()
    fo.close()










