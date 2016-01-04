//
//  global_var.hpp
//  BayesianProbTensorDecomp
//
//  Created by Rachel Ren on 12/27/15.
//  Copyright © 2015 Rachel Ren. All rights reserved.
//

#ifndef global_var_hpp
#define global_var_hpp

using namespace std;

#include <stdio.h>
#include <unordered_map>
#include <vector>

#include <gsl/gsl_matrix.h>

extern int n_factor;
extern int n_individual;
extern int n_gene;
extern int n_tissue;
extern unordered_map<string, int> individual_rep;

//Tensor related arrays (type = INT or DOUBLE)
typedef vector<int> individual;
typedef vector<individual> gene;
typedef vector<gene> tensor;

typedef vector<double> individual_d;
typedef vector<individual_d> gene_d;
typedef vector<gene_d> tensor_d;

extern tensor_d fm; //individual x gene x tissue
//extern tensor dataset(100, gene(100000, individual(100))); //individual x gene x tissue
extern tensor_d dataset; //individual x gene x tissue
extern tensor_d markerset; //individual x gene x tissue
//extern tensor markerset (100, gene(100000, individual(100))); //individual x gene x tissue

// Priors
typedef struct myPrior{
    gsl_vector * mean; //0: mean array
    gsl_matrix * precision_matrix; //1: precision matrix
}myPrior;

extern vector<myPrior> prior;

// Hyper-Priors
typedef struct myArray{
    gsl_matrix * scale_matrix; //0
    int df; //1: degree of freedom
    gsl_vector * mean; //2: mean
    double scaler; // 3: Scaler of the precision matrix
}myArray;

extern vector<myArray> hyper_prior;

//Alpha Prior
extern vector<double> alpha_prior;

//The precision of the final observations
extern double alpha;



#endif /* global_var_hpp */
