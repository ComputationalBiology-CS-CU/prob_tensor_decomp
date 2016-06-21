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

The processing script is named as "pre\_process.py" (extracting the samples, rm null genes, normalization and preparing the input) and "pre\_init.py" (init the factor matrices). You can find under "./data\_raw/" what input data the first script needs in order to preprocess the data. Essentially, they are:

1. GTEx\_Data\_20150112\_RNAseq\_RNASeQCv1.1.8\_gene\_rpkm.gct\_1\_genotype
2. list\_tissue.txt (fro "phs000424.v6.pht002743.v6.p1.c1.GTEx\_Sample\_Attributes.GRU.txt\_tissue\_type\_count\_0")
3. phs000424.v6.pht002743.v6.p1.c1.GTEx\_Sample\_Attributes.GRU.txt\_sample\_tissue\_type
4. gene\_tss.txt

In terms of initializing the Individual factor matrix and the Tissue factor matrix, we will do the PCA on Sample x Gene matrix, and then re-construct the Tissue x Individual x Factor tensor from the Sample factor matrix. For this small tensor, for each factor (assuming they are independent in initialization), we will do incomplete PCA with R script "pre\_init.R". We need to download "fX\_tissue\_indiv.txt" data from cluster "./data\_inter/" to local "./result/temp/", and the script will generate "fX\_tissue.txt" and "fX\_indiv.txt" in the same location. Then we can upload these two source data into cluster "./data\_inter/", and then run the second session of "pre\_init.py". By now we should have done all the initialization.


