import numpy as np

if __name__=="__main__":
    real_data = np.loadtxt("../../../data_prepare/data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_normalize_load")
    real_data = real_data.T
    np.save('data/real_data', real_data)
    print real_data.shape