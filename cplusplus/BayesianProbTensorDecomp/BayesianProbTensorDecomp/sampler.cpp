//
//  sampler.cpp
//  BayesianProbTensorDecomp
//
//  Created by Rachel Ren on 12/29/15.
//  Copyright © 2015 Rachel Ren. All rights reserved.
//

#include "sampler.hpp"
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

double cal_product(int i, int j, int k){
    double product = 0.0;
    
    for (int m = 0; m < n_factor; m++){
        product += fm[0][i][m] * fm[1][j][m] * fm[2][k][m];
    }
    return product;
}

double sampler_Normal(const gsl_rng *r, const double mu, const double sigma){
//    const gsl_rng_type * T;
//    gsl_rng * r;
//    
//    gsl_rng_env_setup();
//    
//    T = gsl_rng_default;
//    r = gsl_rng_alloc (T);
//    
    return gsl_ran_gaussian(r, sigma) + mu;
}

/* The following function is adapted from https://lists.gnu.org/archive/html/help-gsl/2006-04/txtdb8Hdlx9uA.txt */
int sampler_MVN(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result){
    /* multivariate normal distribution random number generator */
    /*
     *	n	dimension of the random vetor
     *	mean	vector of means of size n
     *	var	variance matrix of dimension n x n
     *	result	output variable with a sigle random vector normal distribution generation
     */
    int k;
    gsl_matrix *work = gsl_matrix_alloc(n,n);
    
    gsl_matrix_memcpy(work,var);
    gsl_linalg_cholesky_decomp(work);
    
    for(k=0; k<n; k++)
        gsl_vector_set( result, k, gsl_ran_ugaussian(r) );
    
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
    gsl_vector_add(result,mean);
    
    gsl_matrix_free(work);
    
    return 0;
}

/* The following function is adapted from https://lists.gnu.org/archive/html/help-gsl/2006-04/txtdb8Hdlx9uA.txt */
int sampler_W(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, gsl_matrix *result){
    /* Wishart distribution random number generator */
    /*
     *	n	 gives the dimension of the random matrix
     *	dof	 degrees of freedom
     *	scale	 scale matrix of dimension n x n
     *	result	 output variable with a single random matrix Wishart distribution generation
     */
    int k,l;
    gsl_matrix *work = gsl_matrix_calloc(n,n);
    
    for(k=0; k<n; k++){
        gsl_matrix_set( work, k, k, sqrt( gsl_ran_chisq( r, (dof-k) ) ) );
        for(l=0; l<k; l++){
            gsl_matrix_set( work, k, l, gsl_ran_ugaussian(r) );
        }
    }
    gsl_matrix_memcpy(result,scale);
    gsl_linalg_cholesky_decomp(result);
    gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,result,work);
    gsl_blas_dsyrk(CblasUpper,CblasNoTrans,1.0,work,0.0,result);
    
    return 0;
}

double sampler_Gamma(const gsl_rng * r, double a, double b){
    return gsl_ran_gamma(r, a, b);
}

/* Sample factor relevant components*/
int sampler_factor(int num_factor){
    int dimension = 0;
    int dimension1 = 0;
    int dimension2 = 0;
    int num_factor_1 = 0;
    int num_factor_2 = 0;
    
    if (num_factor == 0){
        dimension = n_individual;
        dimension1 = n_gene;
        dimension2 = n_tissue;
        num_factor_1 = 1;
        num_factor_2 = 2;
    }else if (num_factor == 1){
        dimension = n_gene;
        dimension1 = n_individual;
        dimension2 = n_tissue;
        num_factor_1 = 0;
        num_factor_2 = 2;
    }else{
        dimension = n_tissue;
        dimension1 = n_individual;
        dimension2 = n_gene;
        num_factor_1 = 0;
        num_factor_2 = 1;
    }
    
    /* Sample factor matrix */
    for (int i = 0; i < dimension; i++){
        /* TODO: can use parallel here */
        
        // Precision Matrix
        gsl_matrix * precision_matrix = gsl_matrix_alloc(n_factor, n_factor);
        for (int count1 = 0; count1 < n_factor; count1++){
            for (int count2 = 0; count2 < n_factor; count2++){
                double value = gsl_matrix_get(prior[num_factor].precision_matrix,count1,count2);
                gsl_matrix_set(precision_matrix, count1, count2, value);
            }
        }
        
        for (int j = 0; j < dimension1; j++){
            for (int k = 0; k < dimension2; k++){
                unordered_map<int, int> hash_temp;
                hash_temp.clear(); //In case the initialization is not equivalent to an empty map
                hash_temp.emplace(num_factor, i);
                hash_temp.emplace(num_factor_1, j);
                hash_temp.emplace(num_factor_2, k);
                unordered_map<int,int>::const_iterator got1 = hash_temp.find (num_factor);
                unordered_map<int,int>::const_iterator got2 = hash_temp.find (num_factor_1);
                unordered_map<int, int>::const_iterator got3 = hash_temp.find (num_factor_2);
                
                int index1 = got1->second;
                int index2 = got2->second;
                int index3 = got3->second;
                
                if (markerset[index1][index2][index3] == 0){
                    continue;
                }
                double array[n_factor];
                // Initialize array
                for (int count = 0; count < n_factor; count++){
                    array[count] = fm[num_factor_1][j][count] * fm[num_factor_2][k][count];
                }
                for (int count1 = 0; count1 < n_factor; count1++){
                    for (int count2 = 0; count2 < n_factor; count2++){
                        double previous = gsl_matrix_get (precision_matrix, count1, count2);
                        gsl_matrix_set(precision_matrix, count1, count2, previous + alpha * array[count1] * array[count2]);
                    }
                }
            }
        }
        gsl_matrix * cov = gsl_matrix_alloc(n_factor, n_factor);
        gsl_permutation * perm = gsl_permutation_alloc(n_factor);
        gsl_linalg_LU_invert (precision_matrix, perm, cov);
        
        //Release perm and precision_matrix
        gsl_permutation_free(perm);
        gsl_matrix_free(precision_matrix);
        
        /* Mean array */
        gsl_vector * mean = gsl_vector_alloc(n_factor);
        for (int count1; count1 < n_factor; count1++){
            for (int count2; count2 < n_factor; count2++){
                double prev = gsl_vector_get(mean, count1);
                double new_val = prev + gsl_vector_get(prior[num_factor].mean,count2) * gsl_matrix_get(prior[num_factor].precision_matrix,count1,count2);
                gsl_vector_set(mean, count1, new_val);
            }
        }
        for (int j = 0; j < dimension1; j++){
            for (int k = 0; k < dimension2; k++){
                unordered_map<int, int> hash_temp;
                hash_temp.clear();  //In case the initialization above will not give us an empty map
                hash_temp.emplace(num_factor, i);
                hash_temp.emplace(num_factor_1, j);
                hash_temp.emplace(num_factor_2, k);
                unordered_map<int,int>::const_iterator got1 = hash_temp.find (num_factor);
                unordered_map<int,int>::const_iterator got2 = hash_temp.find (num_factor_1);
                unordered_map<int, int>::const_iterator got3 = hash_temp.find (num_factor_2);
                
                int index1 = got1->second;
                int index2 = got2->second;
                int index3 = got3->second;
                
                if (markerset[index1][index2][index3] == 0){
                    continue;
                }
                double array[n_factor];
                /* Since we assign each element in this array with new values, we don't need to worry about the initialization*/
                for (int count = 0; count < n_factor; count++){
                    array[count] = fm[num_factor_1][j][count] * fm[num_factor_2][k][count];
                }
                double R = dataset[index1][index2][index3];
                for (int count = 0; count < n_factor; count++){
                    double prev = gsl_vector_get(mean, count);
                    gsl_vector_set(mean, count, prev+alpha * R * array[count]);
                }
            }
        }
        // Clear array to full of 0's
        gsl_vector * new_mean = gsl_vector_alloc(n_factor);
        gsl_vector_set_zero(new_mean);
        
        for (int count1 = 0; count1 < n_factor; count1++){
            for (int count2 = 0; count2 < n_factor; count2++){
                double prev = gsl_vector_get(new_mean, count1);
                gsl_vector_set(new_mean, count1, prev+gsl_vector_get(mean,count2) * gsl_matrix_get(cov,count1,count2));
            }
        }
        
        /* Sampling */
        //Initialize the random number
        const gsl_rng_type * T;
        gsl_rng * r;
        
        gsl_rng_env_setup();
        
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        
        //Ready to sample
        gsl_vector * sample_result = gsl_vector_alloc(n_factor);
        sampler_MVN(r, n_factor, new_mean, cov, sample_result);
        
        // Update factorized matrix
        for (int iter = 0; iter < n_factor; iter++){
            fm[num_factor][i][iter] = gsl_vector_get(sample_result, iter);
        }
        
        // Release all the gsl data
        gsl_rng_free(r);
        gsl_vector_free(sample_result);
        gsl_vector_free(mean);
        gsl_vector_free(new_mean);
        gsl_matrix_free(cov);
        
    }
    
    /* Now we are ready to sample the priors */
    /* Wishart first */
    
    //Factor mean
    gsl_vector * factor_mean = gsl_vector_alloc(n_factor);
    gsl_vector_set_zero(factor_mean); //Set factor_mean to be a vector of 0's
    for (int i = 0; i < dimension; i++){
        for (int j = 0; j < n_factor; j++){
            double prev = gsl_vector_get(factor_mean, j);
            gsl_vector_set(factor_mean, j, prev+fm[num_factor][i][j]);
        }
    }
    for (int i = 0; i < n_factor; i++){
        gsl_vector_set(factor_mean, i, gsl_vector_get(factor_mean, i)/dimension);
    }
    
    //Factor varialbe
    gsl_matrix * factor_var = gsl_matrix_alloc(n_factor, n_factor);
    gsl_matrix_set_zero(factor_var); //Set factor_var to be a matrix of 0's
    for (int i = 0; i < dimension; i++){
        for (int count1 = 0; count1 < n_factor; count1++){
            for (int count2 = 0; count2 < n_factor; count2++){
                double prev = gsl_matrix_get(factor_var, count1, count2);
                gsl_matrix_set(factor_var, count1, count2, prev+(fm[num_factor][i][count1] - gsl_vector_get(factor_mean, count1)) * (fm[num_factor][i][count2] - gsl_vector_get(factor_mean, count2)));
            }
        }
    }
    for (int count1 = 0; count1 < n_factor; count1++){
        for (int count2 = 0; count2 < n_factor; count2++){
            gsl_matrix_set(factor_var, count1, count2, gsl_matrix_get(factor_var, count1, count2)/dimension);
        }
    }
    
    //Covariance matrix
    //First do the inverse
    gsl_matrix * cov_matrix = gsl_matrix_alloc(n_factor, n_factor);
    gsl_permutation * perm = gsl_permutation_alloc(n_factor);
    gsl_linalg_LU_invert(hyper_prior[num_factor].scale_matrix, perm, cov_matrix);
    
    //Release perm
    gsl_permutation_free(perm);
    
    //Operate on cov_matrix
    for (int count1 = 0; count1 < n_factor; count1++){
        for (int count2 = 0; count2 < n_factor; count2++){
            double prev = gsl_matrix_get(cov_matrix, count1, count2);
            gsl_matrix_set(cov_matrix, count1, count2, prev + dimension * gsl_matrix_get(factor_var, count1, count2));
        }
    }
    double temp = hyper_prior[num_factor].scaler * dimension / (hyper_prior[num_factor].scaler + dimension);
    for (int count1 = 0; count1 < n_factor; count1++){
        for (int count2 = 0; count2 < n_factor; count2++){
            double toAdd = temp * (gsl_vector_get(hyper_prior[num_factor].mean, count1) - gsl_vector_get(factor_mean, count1)) * (gsl_vector_get(hyper_prior[num_factor].mean,count2) - gsl_vector_get(factor_mean,count2));
            gsl_matrix_set(cov_matrix, count1, count2, gsl_matrix_get(cov_matrix, count1, count2) + toAdd);
        }
    }
    perm  = gsl_permutation_alloc(n_factor);
    gsl_matrix * precision_matrix = gsl_matrix_alloc(n_factor, n_factor);
    gsl_linalg_LU_invert(cov_matrix, perm, precision_matrix);
    
    //Release perm
    gsl_permutation_free(perm);
    
    //Df new
    int df = hyper_prior[num_factor].df + dimension;
    
    //Initialize the random number
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    //Ready to sample Wishart
    gsl_matrix * sample_wishart = gsl_matrix_alloc(n_factor, n_factor);
    
    sampler_W(r, n_factor, df, precision_matrix, sample_wishart);
    prior[num_factor].precision_matrix = sample_wishart;
    
    //Release memory
    gsl_matrix_free(sample_wishart);
    gsl_rng_free(r);
    
    /* Then MVN */
    //Beta new
    double beta = hyper_prior[num_factor].scaler + dimension;
    //gsl_matrix * precision_matrix = gsl_matrix_alloc(n_factor, n_factor);
    //Reuse the initialization in above code
    //Update the values of precision_matrix below
    gsl_matrix_memcpy(precision_matrix, prior[num_factor].precision_matrix);
    gsl_matrix_scale(precision_matrix, beta);
    
    perm  = gsl_permutation_alloc(n_factor);
    gsl_matrix * cov = gsl_matrix_alloc(n_factor, n_factor);
    gsl_linalg_LU_invert(precision_matrix, perm, cov);
    
    gsl_vector * mean = gsl_vector_alloc(n_factor);
    gsl_vector_memcpy(mean, hyper_prior[num_factor].mean);
    gsl_vector_scale(mean, hyper_prior[num_factor].scaler);
    
    gsl_vector * factor_mean_cp = gsl_vector_alloc(n_factor);
    gsl_vector_memcpy(factor_mean_cp, factor_mean);
    gsl_vector_scale(factor_mean_cp, dimension);
    
    gsl_vector_add(mean, factor_mean_cp);
    
    double divisor = hyper_prior[num_factor].scaler + dimension;
    gsl_vector_scale(mean, 1/divisor);
    
    r = gsl_rng_alloc (T);
    gsl_vector * sample_mvn = gsl_vector_alloc(n_factor);
    sampler_MVN(r, n_factor, mean, cov, sample_mvn);
    prior[num_factor].mean = sample_mvn;
    
    //Release memories
    gsl_vector_free(factor_mean);
    gsl_matrix_free(factor_var);
    gsl_matrix_free(cov_matrix);
    gsl_matrix_free(precision_matrix);
    gsl_permutation_free(perm);
    gsl_matrix_free(cov);
    gsl_vector_free(mean);
    gsl_vector_free(factor_mean_cp);
    gsl_rng_free(r);
    gsl_vector_free(sample_mvn);
    
    return 0;
}

/* Sampling precision */
int sampler_precision(){
    //Get parameter
    double para1_old = alpha_prior[0];
    double para2_old = alpha_prior[1];
    
    int n = 0; //Number of non-zero entries
    double sum = 0.0;
    
    for (int i = 0; i < n_individual; i++){
        for (int j = 0; j < n_gene; j++){
            for (int k = 0; k < n_tissue; k++){
                if (markerset[i][j][k] == 0){
                    continue;
                }else{
                    double R_real = dataset[i][j][k];
                    double R_exp = cal_product(i, j, k);
                    sum += pow(R_real - R_exp, 2.0);
                    n++;
                }
            }
        }
    }
    
    //Update the parameter
    double para1_new = para1_old + 0.5 * n;
    double para2_new = para2_old + 0.5 * sum;
    
    //Prepare rng
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    //Ready to sample gamma
    alpha = sampler_Gamma(r, para1_new, para2_new);
    
    //Release memory
    gsl_rng_free(r);
    return 0;
}

/* End of sampler.cpp */
