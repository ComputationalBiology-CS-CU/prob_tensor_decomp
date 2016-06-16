import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


if __name__=="__main__":
    n_factor = 40
    print "draw sample loading matrix using seaborn..."

    ind_loading = np.load('result/ind_loading_brain_c22_quantile.npy')
    sample_subtissue = []
    with open("../data_prepare/data_processed/sample_subtissue","r") as f:
        for lines in f:
            lines = lines.strip()
            line = lines.split("\t")
            sample_subtissue.append([line[0], int(line[1]), int(line[2])])

    cnt = 0
    sample_list, gene_list = [], []
    with open("../data_prepare/data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_normalize_c22","r") as f1:
        for lines in f1:
            lines = lines.strip()
            line = lines.split("\t")
            if cnt==0:
                sample_list = line[1:]
            else:
                gene_list.append(line[0])
            cnt += 1

    for ele in sample_subtissue:
        # Draw the heatmap using seaborn
        sns.set(context="paper", font="monospace")
        f, ax = plt.subplots(figsize=(12, 9))
        sns_plot = sns.heatmap(ind_loading[ele[1]:ele[2],],yticklabels=False)
        filename = 'result/quantile_brain_sample_'+ele[0]+'.csv'
        #np.savetxt(filename,ind_loading[ele[1]:ele[2],],delimiter=',')
        outfile = open(filename, "w")
        outfile.write(" , ")
        for i in range(n_factor):
            if i<n_factor-1:
                outfile.write(str(i)+", ")
            else:
                outfile.write(str(i)+"\n")

        for i in range(ele[1], ele[2]):
            sample_id_ext = sample_list[i].split("-")
            sample_id = sample_id_ext[0] + "-" + sample_id_ext[1]
            outfile.write(str(sample_id)+", ")
            for j in range(n_factor):
                if j<n_factor-1:
                    outfile.write(str(ind_loading[i][j])+", ")
                else:
                    outfile.write(str(ind_loading[i][j])+"\n")
        outfile.close()

        fig = sns_plot.get_figure()
        fig.savefig("plot/quantile_brain_"+ele[0]+".jpg")

    gene_loading = np.load('result/gene_loading_brain_c22_quantile.npy')
    outfile = open("result/quantile_gene_c22.csv", "w")
    outfile.write(" , ")
    for i in range(n_factor):
        if i<n_factor-1:
            outfile.write(str(i)+", ")
        else:
            outfile.write(str(i)+"\n")
    for i in range(len(gene_list)):
        outfile.write(gene_list[i] + ", ")
        for j in range(n_factor):
            if j<n_factor-1:
                outfile.write(str(gene_loading[i][j])+", ")
            else:
                outfile.write(str(gene_loading[i][j])+"\n")
    outfile.close()

    sns.set(context="paper", font="monospace")
    f, ax = plt.subplots(figsize=(12, 9))
    sns_plot = sns.heatmap(gene_loading)
    fig = sns_plot.get_figure()
    fig.savefig("plot/quantile_c22_gene.jpg")
