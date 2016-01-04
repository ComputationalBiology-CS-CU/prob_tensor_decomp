//
//  tensor.hpp
//  BayesianProbTensorDecomp
//
//  Created by Rachel Ren on 12/27/15.
//  Copyright © 2015 Rachel Ren. All rights reserved.
//

#ifndef tensor_hpp
#define tensor_hpp

#include <stdio.h>

#include "tensor.hpp"
#include "main.hpp"
#include "global_var.hpp"
#include "basic.hpp"

//Include Standard Libraries
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <algorithm>
#include <string>
#include <array>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <time.h>

//Include matrix headers
#include <gsl/gsl_matrix.h>

//Include vector headers
#include <gsl/gsl_vector.h>

//Include math headers
#include <gsl/gsl_math.h>

//Include random distributions
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

//Include random number generator
#include <gsl/gsl_rng.h>

using namespace std;

string get_individual_id(string);

int data_prepare();



#endif /* tensor_hpp */
