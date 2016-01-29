//
//  tensor.cpp
//  BayesianProbTensorDecomp
//
//  Created by Rachel Ren on 12/27/15.
//  Copyright © 2015 Rachel Ren. All rights reserved.
//
//  This file contians all the relavant functions to read data and construct the tensor

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
#include <fstream>
#include <sstream>

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
    /* This function will load and format tensor*/

    //Get the number of tissue types (LOCAL!!)
    char filename_type[100] = "/Users/rachelren/Documents/Columbia_CS_Research/expression_by_etissue/tissue_list.txt";
    int count = 0;
    
    FILE * file_in = fopen(filename_type, "r");
    if (file_in == NULL){
        fputs("File error\n", stderr);
        exit(1);
    }
    int input_length = 100;
    char input[input_length];
    while (fgets(input, input_length, file_in) != NULL) {
        count++;
    }
    fclose(file_in);
    n_tissue = count;
    
    //Get all the individuals
    count = 0;
    char filename_indiv[100];
    string header = "/Users/rachelren/Documents/Columbia_CS_Research/expression_by_etissue/tissue_";
    string end = ".txt";
    FILE * file_in_indiv;
    for (int i = 0; i < n_tissue; i++){
        int index_tissue = i + 1;
        string name_tmp = header + to_string(index_tissue) + end;
        //convert name_tmp to char array
        strncpy(filename_indiv, name_tmp.c_str(), sizeof(filename_indiv));
        filename_indiv[sizeof(filename_indiv) - 1] = 0;
        file_in_indiv = fopen(filename_indiv, "r");
        if (file_in_indiv == NULL){
            fputs("File error \n", stderr);
            exit(1);
        }
        int input_length_indiv = 1000;
        char input_indiv[input_length_indiv];
        fgets(input_indiv, input_length_indiv, file_in_indiv);
        trim(input_indiv); //NOTE: Need to adapt from Shuo's code!!!!
        const char * sep = "\t\n";
        char * p;
        p = strtok(input_indiv, sep);
        while (p != NULL) {
            string s(p);
            if (s.compare("Name") != 0){
                //Skip the first one
                string id = get_individual_id(p);
                unordered_map<string, int>::const_iterator got = individual_rep.find(id);
                if (got == individual_rep.end()){
                    individual_rep.emplace(id, count);
                    count++;
                }
            }
            p = strtok(NULL, sep);
        }
        delete[] p;
        fclose(file_in_indiv);
    }
    n_individual = count + 1;
    
    //Get the number of genes
    char filename_gene[100] = "/Users/rachelren/Documents/Columbia_CS_Research/expression_by_etissue/tissue_1.txt";
    count = 0;
    FILE * file_in_gene = fopen(filename_gene, "r");
    if (file_in_gene == NULL){
        fputs("File error\n", stderr);
        exit(1);
    }
    int input_length_gene = 100000;
    char input_gene[input_length_gene];
    while (fgets(input_gene, input_length_gene, file_in_gene) != NULL) {
        count++;
    }
    fclose(file_in_gene);
    n_gene = count;
    
    //Spring edited started
    
    //Initializa the empty tensor first of all
    tensor_d dataset(n_individual, gene_d(n_gene, tissue_d(n_tissue)));
    tensor_d markerset(n_individual, gene_d(n_gene, tissue_d(n_tissue)));
    
    //Get all samples and fill them into the empty tensor and the marker tensor
    
    //Now we try another way to handle file reading
    //string header = "/Users/rachelren/Documents/Columbia_CS_Research/expression_by_etissue/tissue_";
    //string end = ".txt";
    
    for (int k = 0; k < n_tissue; k++){
        int index_tissue = k + 1;
        string name_tmp = header + to_string(index_tissue) + end;
        ifstream inFile(name_tmp);
        if (!inFile) {
            cerr << "File not found." << endl;
            return -1;
        }
        unordered_map<string, vector<double> > rep_temp;
        unordered_map<int, string> index_individual_map;
        
        //Read the first line
        string line;
        getline(inFile, line);
        char *cline = new char[line.length() + 1];
        strcpy(cline, line.c_str());
        trim(cline);
        
        const char * sep = "\t\n";
        char * p;
        p = strtok(cline, sep);
        int iter = 1;
        while (p!= NULL) {
            string s(p);
            if (s.compare("Name") != 0){
                //Skip the first one
                string individual = get_individual_id(p);
                auto search = rep_temp.find(individual);
                if (search != rep_temp.end()) {
                    continue;
                }
                vector<double> list;
                rep_temp[individual] = list;
                index_individual_map[iter] = individual;
            }
            iter++;
        }
        
        delete [] cline;
        
        //Read the rest of the file line by line
        
    }
    
    
    
    return 0;
}
