//
//  log_likelihood.cpp
//  BayesianProbTensorDecomp
//  This file contains functions to calculate the log-likelihood
//
//  Created by Rachel Ren on 12/31/15.
//  Copyright © 2015 Rachel Ren. All rights reserved.
//

#include "log_likelihood.hpp"
#include "sampler.hpp"
//#include "tensor.hpp"
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
#include <math.h>

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

//log likelihood for univariate Gaussian
double loglike_Gaussian(double obs, double mean, double var){
    double like_log = 0.0;
    like_log += 0 - (log(sqrt(var)));
    like_log += 0 - (pow((obs - mean), 2.0)/(2.0 * var));
    return like_log;
}

double loglike_MVN(individual_d obs, gsl_vector * mean, gsl_matrix * precision){
    double like_log = 0.0;
    
    //Precision term
    gsl_matrix * cov = gsl_matrix_alloc(n_factor, n_factor);
    gsl_permutation * perm = gsl_permutation_alloc(n_factor);
    gsl_linalg_LU_invert(precision, perm, cov);
    
    //Calculate the Frobenius norm of a gsl_matrix
    double sum = 0.0;
    double count = 0.0;
    for (int i = 0; i < n_factor; i++){
        for (int j = 0; j < n_factor; j++){
            sum += pow(gsl_matrix_get(cov, i, j),2.0);
            count += 1.0;
        }
    }
    double cov_norm = sqrt(sum/count);
    like_log += 0 - (0.5 * log(cov_norm));
    
    //Distance term
    vector<double> array1(n_factor);
    for (int i = 0; i < n_factor; i++){
        array1[i] = obs[i] - gsl_vector_get(mean, i);
    }
    vector<double> array2(n_factor);
    for (int i = 0; i < n_factor; i++){
        double value = 0.0;
        for (int j = 0; j < n_factor; j++){
            value += array1[j] * gsl_matrix_get(precision, j, i);
        }
        array2[i] = value;
    }
    double temp = 0.0;
    for (int i = 0; i < n_factor; i++){
        temp += array1[i] * array2[i];
    }
    like_log += 0 - (0.5* temp);
    return like_log;
}

//Log-likelihood for Gaussian Wishart
double loglike_GW(const gsl_vector *obs1, const gsl_matrix *obs2, const gsl_matrix *scale, const int df, const gsl_vector *mean, const double scaler){
    double like_log = 0.0;
    //To be completed!!
    return like_log;
}

//Log-likelihood for Gamma
double loglike_Gamma(double obs, double para1, double para2){
    double like_log = 0.0;
    like_log += (para1 - 1) * log(obs);
    like_log += 0 - (para2 * obs);
    return like_log;
}

//Calculate the joint log likelihood
double loglike_joint(){
    double like_log = 0.0;
    
    //observation likelihood
    for (int i = 0; i < n_individual; i++){
        for (int j = 0; j < n_gene; j++){
            for (int k = 0; k < n_tissue; k++){
                if (markerset[i][j][k] == 0){
                    continue;
                }else{
                    double obs = dataset[i][j][k];
                    double mean = cal_product(i, j, k);
                    double var = 1.0 / alpha;
                    like_log += loglike_Gamma(obs, mean, var);
                }
            }
        }
    }
    
    //Factor matrix likelihood
    for (int i = 0; i < n_individual; i++){
        like_log += loglike_MVN(fm[0][i], prior[0].mean, prior[0].precision_matrix);
    }
    for (int j = 0; j < n_gene; j++){
        like_log += loglike_MVN(fm[1][j], prior[1].mean, prior[1].precision_matrix);
    }
    for (int k = 0; k < n_tissue; k++) {
        like_log += loglike_MVN(fm[2][k], prior[2].mean, prior[2].precision_matrix);
    }
    
    //Factor prior likelihood
    for (int i = 0; i < 3; i++){
        like_log += loglike_GW(prior[i].mean, prior[i].precision_matrix, hyper_prior[i].scale_matrix, hyper_prior[i].df, hyper_prior[i].mean, hyper_prior[i].scaler);
    }
    
    //precision/variance likelihood
    like_log += loglike_Gamma(alpha, alpha_prior[0], alpha_prior[1]);
    
    return like_log;
}

/* log_likelihood_cpp */
