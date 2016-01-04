//
//  sampler.hpp
//  BayesianProbTensorDecomp
//
//  Created by Rachel Ren on 12/29/15.
//  Copyright © 2015 Rachel Ren. All rights reserved.
//

#ifndef sampler_hpp
#define sampler_hpp

#include <stdio.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

using namespace std;

double cal_product(int, int , int);

double sampler_Normal(const gsl_rng *, const double, const double);

int sampler_MVN(const gsl_rng *, const int, const gsl_vector *, const gsl_matrix *, gsl_vector *);

int sampler_W(const gsl_rng *, const int, const int, const gsl_matrix *, gsl_matrix *);

double sampler_Gamma(const gsl_rng *, double, double);

int sampler_factor(int);

int sampler_precision();

#endif /* sampler_hpp */
