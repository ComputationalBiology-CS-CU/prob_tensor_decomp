//
//  main.cpp
//  BayesianProbTensorDecomp
//
//  Created by Rachel Ren on 12/25/15.
//  Copyright © 2015 Rachel Ren. All rights reserved.
//

#include "main.hpp"
#include "tensor.hpp"
#include "global_var.hpp"
#include "sampler.hpp"
#include "log_likelihood.hpp"

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

//Declare Global Variables
int n_factor = 0;
int n_individual = 0;
int n_gene = 0;
int n_tissue = 0;
//Not sure if dataset, markerset and fm should be initialized here
tensor_d dataset(100, gene_d(100000, individual_d(100)));		// individual x gene x tissue
tensor_d markerset(100, gene_d(100000, individual_d(100)));	// mark the position where there are data
unordered_map<string, int> individual_rep;
tensor_d fm(100, gene_d(100000, individual_d(100)));
vector<myPrior> prior (3); //0: mean array 1: precision matrix
vector<myArray> hyper_prior (3); //0: scale_matrix 1: df 2: mean 3: scaler of the precision matrix
double alpha = 0.0;		// the precision of the final observation
vector<double> alpha_prior(2);


int main() {
    cout<<"Start of the program";
    // prepare the "dataset" and "markerset"
    data_prepare();
    //Initializing the global varaiables
    n_factor = 400;
    
    for (int n = 0; n < 3; n++){
        gsl_matrix *scale = gsl_matrix_alloc(n_factor, n_factor);
        for (int i = 0; i < n_factor; i++) {
            for (int j = 0; j< n_factor; j++){
                if (i == j) {
                    gsl_matrix_set(scale, i, j, 1.0);
                }else{
                    gsl_matrix_set(scale, i, j, 0.0);
                }
            }
        }
        hyper_prior[n].scale_matrix = scale;
        hyper_prior[n].df = n_factor;
        gsl_vector *mean = gsl_vector_alloc(n_factor);
        gsl_vector_set_zero(mean);
        hyper_prior[n].mean = mean;
        hyper_prior[n].scaler = 1;
    }
    
    //the prior of MVN (mean and precision matrix)
    for (int n = 0; n < 3; n++){
        gsl_vector *mean = gsl_vector_alloc(n_factor);
        gsl_vector_set_zero(mean);
        gsl_matrix *precision = gsl_matrix_alloc(n_factor, n_factor);
        for (int i = 0; i < n_factor; i++){
            for (int j = 0; j < n_factor; j++){
                if (j == i){
                    gsl_matrix_set(precision, i, j, 1.0);
                }else{
                    gsl_matrix_set(precision, i, j, 0.0);
                }
            }
        }
        prior[n].mean = mean;
        prior[n].precision_matrix = precision;
    }
    
    //the MVN drawing (mean 0, cov 1) for factorized matrices
    
    gsl_vector *mean = gsl_vector_alloc(n_factor);
    gsl_vector_set_zero(mean);
    gsl_matrix *cov = gsl_matrix_alloc(n_factor, n_factor);
    
    for (int i = 0; i < n_factor; i++){
        for (int j = 0; j < n_factor; j++){
            if (j == i){
                gsl_matrix_set(cov, i, j, 1.0);
            }else{
                gsl_matrix_set(cov, i, j, 0.0);
            }
        }
    }

    //To be completed!
    
    //set the prior for precision (Gamma)
    
    alpha_prior[0] = 1;
    alpha_prior[1] = 0.5;
    
    //drawing precision from Gaussian.
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    alpha = sampler_Normal(r, 0, 1);
    
    gsl_rng_free(r);
    
    int ITER = 1000;
    for (int i = 0; i < ITER; i++){
        for (int j = 0; j < 3; j++){
            sampler_factor(j);
        }
        sampler_precision();
        double like_log = loglike_joint();
        cout << "Loglikihood is " << like_log;
    }
    
    gsl_vector_free(mean);
    gsl_matrix_free(cov);

    return 0;
}
