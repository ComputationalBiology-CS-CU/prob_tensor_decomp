This is the Python version of the program (without sparsity prior).

This folder contains all the scripts developed by Shuo, Chuqiao, Mengqing (timely order).

Currently used scripts:

1. sampler\_mw.py (tensor decomposition Gibbs sampler)
2. sampler\_mw\_mf.py (matrix decomposition Gibbs sampler)
3. FA.py
4. PCA.py
5. plot\_pca.py
6. plot\_loglike.py
7. plot\_fm\_heatmap.py

In terms of gene expression data preprocessing, we will base on the results after the first gene expression processing stage of [this](https://github.com/morrisyoung/eQTL_v6_script) project, and pick up tissue samples from specified tissues (brain; I will manually specify this file), and prepare the input for the training program. The training program needs the following input:

1. Tissue.npy
2. Tissue\_list.npy
3. Individual.npy
4. Individual\_list.npy
5. Gene.npy
6. Gene\_list.npy
7. Tensor\_tissue\_x.npy
8. Tensor\_tissue\_x\_list.npy

and we already have the following hyper prior for the training program:

1. kappa.npy
2. mu.npy
3. precision.npy
4. v.npy

The processing script is named as "pre\_process.py" (extracting the samples, rm null genes, normalization and preparing the input) and "pre\_init.py" (init the factor matrices).


