## build the qualified sample-tissue rep; to be further used later on
## two criterions: 1. sample has genotype information; and 2. sample is in eQTL tissue (sample size>=#size)
## right now I use 100 as the threashold for eQTL tissues (33 etissues left then)



# below from previous project:
'''
## there are four sequential intermediate files:
##	1. xxx_1_genotype (remove the genotype-missing individuals)
##	2. xxx_2_esample (pick up all samples that are defined in etissues)
##	3. xxx_3_gene_1_null (remove null genes among these esamples)
##	4. xxx_3_gene_2_normalize
## and we need to do #2 before #3 and #4

## in this script, I will actually extract esamples and do #3 and #4
'''


## as this is for Tensor project, I will use the brain tissue list to select the samples

## we need the following from this script as the input of the training program:
'''
Tissue.npy
—> Tissue_list.npy
Individual.npy
—> Individual_list.npy
Gene.npy
—> Gene_list.npy
Tensor_tissue_0.npy (0-x)
—> Tensor_tissue_0_list.npy (0-x)
'''









##=====================
##==== libraries
##=====================
import math
import sys
import time
import os
import numpy as np
from scipy import stats
import scipy.stats as st





##=====================
##==== global variables
##=====================
individual_rep = {}		# hashing all the individuals with genotype information
sample_tissue_map = {}		# mapping all the samples into their tissue types
ratio_null = 0.5		# TODO: at least this portion of genes are expressed over the below value are treated as expressed genes
rpkm_min = 0.1			# TODO: see above






##==================
##==== sub-routines
##==================
# get the "xxx-yyy" from "xxx-yyy-zzz-aaa-qqq", which is defined as the individual ID of the GTEx samples
def get_individual_id(s):
	## naively find the second '-'
	id = ''
	count = 0
	for i in range(len(s)):
		if s[i] == '-':
			count += 1

		if count == 2:
			break

		id += s[i]

	return id




# at least ratio_null portion of genes are expressed over rpkm_min will be treated as an expressed gene
# this can be later re-defined according to other rules
def check_null(l):
	# transform the list first of all (from string to float)
	l = map(lambda x: float(x), l)

	count = 0
	for i in range(len(l)):
		if l[i] > rpkm_min:
			count += 1

	ratio = (count * 1.0) / len(l)

	if ratio > ratio_null:
		return 0
	else:
		return 1






if __name__ == '__main__':




	##======================================================================================================
	##==== extracting eQTL samples (the eTissue is defined as sample size >= #size)
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_#size"
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample
	##======================================================================================================
	## get all the samples for tissue count >= filter
	## eQTL_tissue
	file = open("./data_raw/list_tissue.txt", 'r')
	eQTL_tissue = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		tissue = (line.split('\t'))[0]
		eQTL_tissue[tissue] = []
	file.close()
	list_tissue = []
	for tissue in eQTL_tissue:
		list_tissue.append(tissue)
	list_tissue = np.array(list_tissue)
	## NOTE
	np.save("./data_processed/Tissue_list.npy", list_tissue)
	print "# of tissues:",
	print len(list_tissue)


	## sample_list
	file = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r')
	sample_list = (((file.readline()).strip()).split('\t'))[1:]
	file.close()

	## sample_tissue_rep
	file = open("./data_raw/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_sample_tissue_type", 'r')
	sample_tissue_rep = {}
	while 1:
		line = file.readline()[:-1]
		if not line:
			break

		line = line.split('\t')

		if len(line) < 3:
			print line
			continue

		sample = line[0]
		tissue = line[2]

		sample_tissue_rep[sample] = tissue
	file.close()

	## fill in the eQTL_tissue rep
	for sample in sample_list:
		tissue = sample_tissue_rep[sample]
		if tissue in eQTL_tissue:
			eQTL_tissue[tissue].append(sample)

	# save the rep
	file = open("./data_raw/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_rep", 'w')
	for tissue in eQTL_tissue:
		file.write(tissue + '\t')
		for sample in eQTL_tissue[tissue]:
			file.write(sample + '\t')
		file.write('\n')
	file.close()


	##============ process the rpkm matrix to get eQTL samples ==============
	## get the sample_rep first
	sample_rep = {}
	file = open("./data_raw/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_rep", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')[1:]
		for sample in line:
			sample_rep[sample] = 1
	file.close()

	# save all the individuals (in order)
	individual_rep = {}
	for sample in sample_rep:
		individual = get_individual_id(sample)
		if individual not in individual_rep:
			individual_rep[individual] = 1
	list_individual = []
	for individual in individual_rep:
		list_individual.append(individual)
	list_individual = np.arrray(list_individual)
	## NOTE
	np.save("./data_processed/Individual_list.npy", list_individual)
	print "# of total individuals:",
	print len(list_individual)


	file = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r')
	file1 = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample", 'w')

	# filter all the samples again
	index_rep = {}
	line = (file.readline()).strip()
	line = line.split('\t')
	file1.write(line[0] + '\t')
	for i in range(1, len(line)):
		sample = line[i]
		if sample in sample_rep:
			index_rep[i] = 1
			file1.write(sample + '\t')
	file1.write('\n')
	
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		file1.write(line[0] + '\t')
		for i in range(1, len(line)):
			if i in index_rep:
				file1.write(line[i] + '\t')
		file1.write('\n')

	file.close()
	file1.close()



	##======================================================================================================
	##==== remove all the NULL genes as defined (testing for all samples)
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null
	##======================================================================================================
	file = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample", 'r')
	file1 = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'w')
	line = file.readline()
	file1.write(line)

	list_gene = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		# check this gene
		line = line.split('\t')
		if check_null(line[1:]):
			continue
		else:
			list_gene.append(line[0])
			file1.write(line[0] + '\t')
			for i in range(1, len(line)):
				file1.write(line[i] + '\t')
			file1.write('\n')

	file.close()
	file1.close()

	list_gene = np.array(list_gene)
	## NOTE
	np.save("./data_processed/Gene_list.npy", list_gene)
	print "# of genes:",
	print len(list_gene)

	##==== remove temp data
	os.system("rm " + "./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample")






	##=============================================================================================
	##==== normalizing all the samples (here we use Log normalize other than the previous Quantile)
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize
	##=============================================================================================
	## the following several normalization methods:
	normal_quantile = 0
	normal_log = 0
	normal_z = 0
	normal_Gaussian_rank = 1

	if normal_quantile:
		print "quantile normalization..."

		file = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'r')
		file1 = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize", 'w')
		line = file.readline()
		file1.write(line)


		gene_list = []
		description_list = []
		rpkm_rep = {}  # hashing gene to L1
		L2 = []

		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			## sorting; indexing (additive filling L2)
			line = line.split('\t')
			gene = line[0]
			description = line[1]
			gene_list.append(gene)
			description_list.append(description)

			expression = map(lambda x: float(x), line[2:])
			expression = np.array(expression)

			sort = np.argsort(expression)
			## get the ordering list (or the rank list for this gene for all samples)
			sort1 = []
			for i in range(len(sort)):
				sort1.append(0)
			for i in range(len(sort)):	# the current i is the rank
				index = sort[i]
				sort1[index] = i

			rpkm_rep[gene] = sort1

			if len(L2) == 0:
				for pos in sort:
					rpkm = expression[pos]
					L2.append(rpkm)
			else:
				for i in range(len(sort)):
					pos = sort[i]
					rpkm = expression[pos]
					L2[i] += rpkm
		file.close()

		length = len(gene_list)
		for i in range(len(L2)):
			L2[i] = L2[i] * 1.0 / length


		for i in range(len(gene_list)):
			## two lists:
			## L1 (value as the re-mapped positions of all original elements; each gene has one such list)
			## L2 (containing the normalized/averaged value for each index/position)
			gene = gene_list[i]
			description = description_list[i]
			L1 = rpkm_rep[gene]

			file1.write(gene + '\t' + description + '\t')

			for index in L1:
				value = L2[index]
				file1.write(str(value) + '\t')

			file1.write('\n')

		file1.close()

		##==== remove temp data
		os.system("rm " + "./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null")


	if normal_log:
		print "log normalization..."

		file = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'r')
		file1 = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize", 'w')
		line = file.readline()
		file1.write(line)

		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			gene = line[0]
			file1.write(gene + '\t')

			rpkm_list = map(lambda x: float(x), line[1:])
			for i in range(len(rpkm_list)):
				rpkm = rpkm_list[i]
				rpkm = math.log(rpkm + 0.1)	# NOTE: here is the rule of transformation: shifted logarithm
				file1.write(str(rpkm) + '\t')

			file1.write('\n')

		file.close()
		file1.close()

		##==== remove temp data
		os.system("rm " + "./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null")


	if normal_z:
		print "z normalization..."

		file = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'r')
		file1 = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize_z", 'w')
		line = file.readline()
		file1.write(line)

		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			gene = line[0]
			file1.write(gene + '\t')

			rpkm_list = map(lambda x: float(x), line[1:])
			rpkm_list = np.array(rpkm_list)
			rpkm_list = stats.zscore(rpkm_list)

			for i in range(len(rpkm_list)):
				rpkm = rpkm_list[i]
				file1.write(str(rpkm) + '\t')
			file1.write('\n')

		file.close()
		file1.close()

		##==== remove temp data
		os.system("rm " + "./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null")


	if normal_Gaussian_rank:
		print "Gaussian rank normalization..."

		file = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'r')
		file1 = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize_Gaussian_rank", 'w')
		line = file.readline()
		file1.write(line)

		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			gene = line[0]
			file1.write(gene + '\t')

			rpkm_list = map(lambda x: float(x), line[1:])
			rpkm_list = np.array(rpkm_list)
			rpkm_list_rank = rpkm_list.argsort()
			rpkm_list_rank_1 = map(lambda x: x+1, rpkm_list_rank)
			for rank in rpkm_list_rank_1:
				zscore = st.norm.ppf(rank*1.0/(len(rpkm_list_rank_1)+1))
				file1.write(str(zscore) + '\t')
			file1.write('\n')

		file.close()
		file1.close()

		##==== remove temp data
		os.system("rm " + "./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null")


	





	##======================================================================================================
	##==== separating all the esamples into their tissues
	##======================================================================================================
	# list_sample_all, Data
	file = open("./data_raw/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_2_normalize_Gaussian_rank", 'r')
	line = (file.readline()).strip()
	line = line.split('\t')[1:]
	list_sample_all = line
	print "total # of samples:",
	print len(list_sample_all)

	Data = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break
		line = line.split('\t')[1:]
		rpkm_list = map(lambda x: float(x), line)
		Data.append(rpkm_list)
	file.close()
	Data = Data.T
	np.save("./data_processed/Data", Data)
	print "shape of Data matrix:",
	print Data.shape



	list_tissue = np.load("./data_processed/Tissue_list.npy")
	#sample_tissue_rep	# we have this
	for i in range(len(list_tissue)):
		tissue = list_tissue[i]

		list_sample = []
		Y = []

		for j in range(len(list_sample_all)):
			sample = list_sample_all[j]
			if sample_tissue_rep[sample] == tissue:		# take j column for all genes
				list_sample.append(sample)
				Y.append(Data[j])
		list_sample = np.array(list_sample)
		Y = np.array(Y)
		np.save("./data_processed/Tensor_tissue_" + str(i) + "_list.npy", list_sample)
		np.save("./data_processed/Tensor_tissue_" + str(i) + ".npy", Y)

		print tissue,
		print list_sample.shape,
		print Y.shape





