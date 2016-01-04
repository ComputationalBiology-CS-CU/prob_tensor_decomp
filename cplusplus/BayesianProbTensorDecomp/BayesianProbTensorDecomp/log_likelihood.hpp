//
//  log_likelihood.hpp
//  BayesianProbTensorDecomp
//
//  Created by Rachel Ren on 12/31/15.
//  Copyright © 2015 Rachel Ren. All rights reserved.
//

#ifndef log_likelihood_hpp
#define log_likelihood_hpp

#include <stdio.h>

#include "global_var.hpp"
//Include matrix headers
#include <gsl/gsl_matrix.h>

//Include linear algebra headers
#include <gsl/gsl_linalg.h>

//Include vector headers
#include <gsl/gsl_vector.h>

//Include math headers
#include <gsl/gsl_math.h>

//Include random distributions
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

//Include random number generator
#include <gsl/gsl_rng.h>

//Include cblas
#include <gsl/gsl_cblas.h>

//Include gamma distribution
#include <gsl/gsl_sf_gamma.h>

//Include blas
#include <gsl/gsl_blas.h>

using namespace std;

double loglike_Gaussian(double, double, double);

double loglike_MVN(individual_d, gsl_vector *, gsl_matrix *);

double loglike_GW(const gsl_vector *, const gsl_matrix *, const gsl_matrix *, const int, const gsl_vector *, const double);

double loglike_Gamma(double, double, double);

double loglike_joint();


#endif /* log_likelihood_hpp */
