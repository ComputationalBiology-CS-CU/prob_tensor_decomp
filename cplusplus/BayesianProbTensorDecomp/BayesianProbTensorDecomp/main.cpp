//
//  main.cpp
//  BayesianProbTensorDecomp
//
//  Created by Rachel Ren on 12/25/15.
//  Copyright © 2015 Rachel Ren. All rights reserved.
//

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


//Declare Global Variables
int n_factor = 0;
int n_individual = 0;
int n_gene = 0;
int n_tissue = 0;
//dataset = []		// individual x gene x tissue
//markerset = []		// mark the position where there are data
//individual_rep = {}
//fm1 = []
//fm2 = []
//fm3 = []
//fm = []
//prior1 = []		// 0: mean array; 1: precision matrix
//prior2 = []		// as above
//prior3 = []		// as above
//prior = []
//hyper_prior1 = []	// 0: xxx; 1: xxx; 2: xxx; 3: xxx
//hyper_prior2 = []	// as above
//hyper_prior3 = []	// as above
//hyper_prior = []
//
//alpha = 0		// the precision of the final observation
//alpha_prior = []	// 0: xxx; 1: xxx

using namespace std;

string get_individual_id(string s){
    string id;
    size_t len = s.length();
    int count = 0;
    
    for (int i = 0; i < len; i++) {
        if (s.at(i) == '-'){
            count++;
        }
        if (count == 2){
            break;
        }
        id+=s.at(i);
    }
    return id;
}

int data_prepare(){
    int count = 0;
    return 0;
}



int main(int argc, const char * argv[]) {
    // insert code here...
    cout << "Start the program\n";
    string id = get_individual_id("xxx-yyy-zzz-aaa-qqq");
    cout << id << "\n";
    return 0;
}
