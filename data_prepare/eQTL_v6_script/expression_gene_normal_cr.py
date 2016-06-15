## build the qualified sample-tissue rep; to be further used later on
## two criterions: 1. sample has genotype information; and 2. sample is in eQTL tissue (sample size>=#size)
## right now I use 100 as the threashold for eQTL tissues (33 etissues left then)



## there are four sequential intermediate files:
##	1. xxx_1_genotype (remove the genotype-missing individuals)
##	2. xxx_2_esample (pick up all samples that are defined in etissues)
##	3. xxx_3_gene_1_null (remove null genes among these esamples)
##	4. xxx_3_gene_2_normalize
## and we need to do #2 before #3 and #4

## in this script, I will actually extract esamples and do #3 and #4





##=====================
##==== libraries
##=====================
import math
import sys
import time
import os
import numpy as np
import scipy.stats as st






##=====================
##==== global variables
##=====================
individual_rep = {}		# hashing all the individuals with genotype information
sample_tissue_map = {}		# mapping all the samples into their tissue types
filter = 100			# TODO: we can change this to get more or less eTissues
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




def getSubtissue(tname):
	tissueDict = {}
	with open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_0","r") as f:
		for line in f:
			line = line.strip()
			line = line.split("\t")

			word = line[0].split(" ")
			if word[0]==tname:
				tissueDict[line[0]]=[]

	f.close()
	return tissueDict



if __name__ == '__main__':


	##======================================================================================================
	##==== extracting eQTL samples (the eTissue is defined as sample size >= #size)
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_#size"
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample
	##======================================================================================================
	## get all the samples for tissue count >= filter
	## eQTL_tissue

	# mw

	tissue_test = "brain"
	eQTL_tissue = getSubtissue("Brain")
	cid = 22


	## sample_list
	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r')
	sample_list = (((file.readline()).strip()).split('\t'))[1:]
	file.close()

	## sample_tissue_rep
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_sample_tissue_type", 'r')
	sample_tissue_rep = {}
	while 1:
		line = file.readline()[:-1]
		if not line:
			break

		line = line.split('\t')

		if len(line) < 3:
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

	cnt0, cnt1 = 0, 0
	## save the rep
	file1 = open("../data_processed/sample_subtissue","w")
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_" + tissue_test, 'w')
	for tissue in eQTL_tissue:
		file.write(tissue + '\t')
		for sample in eQTL_tissue[tissue]:
			file.write(sample + '\t')
			cnt1 += 1
		file.write('\n')
		file1.write(tissue + "\t" + str(cnt0) + "\t" + str(cnt1) + "\n")
		cnt0 = cnt1+1

	file.close()
	file1.close()


	##============ process the rpkm matrix to get eQTL samples ==============
	## get the sample_rep first
	all_sample_we_need = []
	all_ind_we_need = {}
	sample_rep = {}
	tissue_list = []
	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_sample_" + tissue_test, 'r')
	file1 = open("../data_processed/individual_tissue","w")
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line_ele = line.split('\t')
		tissue_list.append(line_ele[0])
		all_sample_we_need += line_ele[1:]

		for s in line_ele[1:]:
			ind_id = get_individual_id(s)
			if ind_id not in all_ind_we_need:
				all_ind_we_need[ind_id] = []
			all_ind_we_need[ind_id].append(line_ele[0])

		for sample in line_ele[1:]:
			sample_rep[sample] = 1
	file.close()

	for k,v in all_ind_we_need.iteritems():
		file1.write(str(k)+"\t")
		for ele in v[:-1]:
			file1.write(str(ele)+"\t")
		file1.write(str(v[-1])+"\n")

	file1.close()


	sample_ids = []
	sample_id_map = {}
	tensor_data = []
	gene_name = []

	# remove row and col headers
	with open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r') as f:
		i = 0
		for line in f:
			line = line.strip()
			line_ele = line.split("\t")
			if i==0:
				sample_ids = line_ele[1:]
				for j in range(len(line_ele[1:])):
					sample_id_map[line_ele[1:][j]] = j
			else:
				gene_name.append(line_ele[0])
				tensor_data.append([float(ele) for ele in line_ele[1:]])

			i += 1

	f.close()
	tensor_data = np.matrix(tensor_data)

	tensor_reorg = []
	for s in all_sample_we_need:
		col_id = sample_id_map[s]

		#slice
		a = tensor_data[:,col_id]
		b = a.flatten()
		c = b.tolist()
		d = c[0]
		tensor_reorg.append(d)

	tensor_reorg = np.matrix(tensor_reorg)
	tensor_reorg = tensor_reorg.T


	f1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample_brain", 'w')


	f1.write("Name" + "\t")

	for ele in all_sample_we_need[:-1]:
		f1.write(str(ele)+"\t")
	f1.write(str(all_sample_we_need[-1])+"\n")

	for i in range(len(gene_name)):
		f1.write(str(gene_name[i])+"\t")

		data_write = tensor_reorg[i,:]
		a = data_write.tolist()
		b = a[0]

		for d in b[:-1]:
			f1.write(str(d)+"\t")
		f1.write(str(b[-1])+"\n")

	f1.close()




	##======================================================================================================
	##==== remove all the NULL genes as defined (testing for all samples)
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null
	##======================================================================================================
	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_2_esample_brain", 'r')
	file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'w')
	line = file.readline()
	file1.write(line)

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		# check this gene
		line = line.split('\t')
		if check_null(line[1:]):
			continue
		else:
			file1.write(line[0] + '\t')
			for i in range(1, len(line)):
				file1.write(line[i] + '\t')
			file1.write('\n')

	file.close()
	file1.close()


	##=============================================================================================
	##==== Filter genes on some chromosome
	##=============================================================================================
	#file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene", 'w')
	file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_chromosome_c22", 'w')
	gc_map = {}
	with open("../data_processed/gencode.v19.genes.patched_contigs.gtf_gene", "r") as gc_map_file:
		for line in gc_map_file:
			ele = line.split('\t')
			gc_map[ele[3].strip()] = ele[0]

	with open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_1_null", 'r') as gene_data:
		lineid = 0
		for line in gene_data:
			line = line.strip()
			ele = line.split('\t')
			if lineid==0 or gc_map[ele[0]]==str(cid):
				for i in range(0, len(ele)):
					file1.write(ele[i] + '\t')
				file1.write('\n')
			lineid += 1

	file1.close()
	gene_data.close()
	gc_map_file.close()

	'''

	##=============================================================================================
	##==== transform the data into gaussian
	##=============================================================================================

	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_chromosome_c22", 'r')
	file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_gaussian_c22", 'w')
	line = file.readline()
	file1.write(line)

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		gene_ind = []
		for i in range(1, len(line)):
			gene_ind.append(float(line[i]))

		gene_ind = np.array(gene_ind)
		gene_rank = gene_ind.argsort()
		gene_rank_1 = [ele+1 for ele in gene_rank]
		file1.write(gene + "\t")
		for ele in gene_rank_1:
			zscore = st.norm.ppf(ele*1.0/(len(line)+1))
			file1.write(str(zscore) + "\t")
		file1.write("\n")

	file.close()
	file1.close()
	'''


	#===========================================gaussian ends===================================================#


	##=============================================================================================
	##==== transform the data into quantile
	##=============================================================================================


	print "quantile normalization..."

	#file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_chromosome_c"+str(cid), 'r')
	#file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_quantile_c"+str(cid), 'w')
	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_chromosome_c22", 'r')
	file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_quantile_c22", 'w')
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



	#===========================================quantile ends===================================================#


	f1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_quantile_load_c22", 'w')
	with open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_3_gene_quantile_c22", 'r') as f2:
		lineid = 0
		for line in f2:
			if lineid==0:
				lineid += 1
				continue
			line = line.strip()
			ele = line.split('\t')
			for i in range(1, len(ele)):
				f1.write(ele[i] + '\t')
			f1.write('\n')
			lineid += 1

	f2.close()
	f1.close()
