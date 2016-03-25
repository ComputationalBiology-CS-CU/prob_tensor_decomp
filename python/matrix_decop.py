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
import sampler
import sys
sys.path.insert(0, '../simulator')
import simulator

##=====================
##==== global variables
##=====================
USE_SPARSITY = False
SMALL_DATA = True   # If true, we will simulate a smaller version of matrix
USE_SIMULATION = True

if __name__ == '__main__':

    # DEBUG
    print "enter program..."

    # DEBUG
    print "Do simulation"
    product = simulator.simulation(USE_SPARSITY, SMALL_DATA)

    print "Do Sampling"










