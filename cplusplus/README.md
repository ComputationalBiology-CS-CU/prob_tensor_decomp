###C++ version of the __Bayesian Probabilistic Tensor Decompositions__

_Based on the code wrote by Shuo in the python folder_  
Note: This code is the Gaussian distribution portion of the overall graphical model. The Spike-and-Slab Sparsity prior on gene dimension will be added (hopefully) in the future.

####Hierarchy of the folder
- __BayesianProbTensorDecomp__: This folder contains all the cpp files for this model
- __gsl_2.1__: This folder is the cpp version of the Gnu Scientific Computation Library.
- __Bayesian_Probabilistic_Tensor_Decomposition_with_Spike_and_Slab_Prior__: This is the final report, which contains the overall graphical model and the details of Gibbs sampler.

####Content of BayesianProbTensorDecomp
- basic.cpp: This is adopted from Shuo Yang's eQTL project. It contains some useful basic functions.
- global_var.hpp: This header file defines all the global variables.
- log_likehood.cpp: This file contains all the functions to compute the loglikihood for each distribution.
- main.cpp: This file initialize the global variable and contains the main function.
- sampler.cpp: This file contains all the functions to sample different distribution.
- tensor.cpp: This file contains functions to prepare dataset. 


Last edited by Chuqiao Ren 01/04/2016.
