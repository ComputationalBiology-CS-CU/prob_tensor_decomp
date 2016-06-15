#This file is to draw heat map for factor matrix

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

f_b = np.load("./result/gene_result.npy")
f_i = np.load("./result/ind_result.npy")



#DEBUG
print "prepare to draw matrix..."

sns.set(context="paper", font="monospace", font_scale=13)
# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(120, 100))

#DEBUG
print "draw gene using seaborn..."

# Draw the heatmap using seaborn
sns_plot_b = sns.heatmap(f_b)
sns_plot_b.set(xlabel='Factors', ylabel='Gene')
fig_b = sns_plot_b.get_figure()
fig_b.savefig("./result/factor_gene_new.png")

#DEBUG
print "prepare to draw matrix..."

sns.set(context="paper", font="monospace", font_scale=13)
# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(120, 100))

#DEBUG
print "draw individual using seaborn..."

sns_plot_i = sns.heatmap(f_i)
sns_plot_i.set(xlabel='Factors', ylabel='Individual')
fig_i = sns_plot_i.get_figure()
fig_i.savefig("./result/factor_individual_new.png")
# f.tight_layout()