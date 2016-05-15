import six

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


if __name__ == "__main__":
    colors = list(six.iteritems(colors.cnames))

    color = [str(ele[0]) for ele in colors]

    ind = np.load("data/real_individual_brain_c22_quantile.npy")

    d1, d2 = ind[:,1], ind[:,2]
    tissue_list = []
    tissue_name = []
    with open("../data_prepare/data_processed/sample_subtissue","r") as f1:
        for lines in f1:
            lines = lines.strip()
            cols = lines.split('\t')
            tissue_list.append((int(cols[1]), int(cols[2])))
            tissue_name.append(cols[0])

    for i in range(len(tissue_list)):
        #t = np.random.rand(1)
        #print t
        #colors = [t[0] for j in range(len(d1[tissue_list[i][0]:tissue_list[i][1]+1]))]
        #print colors
        #print len(d1[tissue_list[i][0]:tissue_list[i][1]+1])
        #print len(d2[tissue_list[i][0]:tissue_list[i][1]+1])
        ##plt.plot(d1[tissue_list[i][0]:tissue_list[i][1]+1], d2[tissue_list[i][0]:tissue_list[i][1]+1], color=color[i], marker='.', alpha=0.5)
        line1, = plt.plot(d1[tissue_list[i][0]:tissue_list[i][1]+1], d2[tissue_list[i][0]:tissue_list[i][1]+1], 'o', color=color[i*3], alpha=0.5, label=tissue_name[i])
        plt.legend(loc='upper left', shadow=True, fontsize='x-small')


    #for i in range(len(d1)):
    #    plt.annotate(i, (d1[i],d2[i]))
    plt.xlabel("PC 2")
    plt.ylabel("PC 3")
    plt.title("Chromosome 22 PCA First Two Components Analysis")
    plt.show()
