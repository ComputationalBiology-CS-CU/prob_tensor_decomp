## build the qualified sample-tissue rep; to be further used later on
## two criterions: 1. sample has genotype information; and 2. sample is in eQTL tissue (sample size>=#size)
## right now I use 100 as the threashold for eQTL tissues (33 etissues left then)



## there are four sequential intermediate files:
##	1. xxx_1_genotype (remove the genotype-missing individuals)
##	2. xxx_2_esample (pick up all samples that are defined in etissues)
##	3. xxx_3_gene_1_null (remove null genes among these esamples)
##	4. xxx_3_gene_2_normalize
## and we need to do #2 analysis before we do #2, #3 and #4

## in this script, I will do the sample distribution analysis, and practically extract samples for etissues later on




##=====================
##==== libraries
##=====================
import math






##=====================
##==== global variables
##=====================
individual_rep = {}		# hashing all the individuals with genotype information
sample_tissue_map = {}		# mapping all the samples into their tissue types
filter = 100			# TODO: we can change this to get more or less eTissues






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





if __name__ == '__main__':




	##========================================================================================
	##==== get the sample-tissue mapping for all the samples
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type
	##========================================================================================

	file = open("../data/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt", 'r')
	file1 = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_sample_tissue_type", 'w')
	count = 0
	while 1:
		line = file.readline()[:-1]	# can't strip all "\t\t\t\t", as they are place holder
		count += 1
		if count <= 11:  ## 11
			continue

		if not line:
			break

		line = line.split('\t')
		sample = line[1]
		tissue1 = line[12]
		tissue2 = line[14]

		file1.write(sample + '\t' + tissue1 + '\t' + tissue2 + '\n')

	file.close()
	file1.close()


	# sample_tissue_type is the tissue count for each individual


	##========================================================================================
	##==== remove samples that have no genotype information
	##==== target: GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype
	##========================================================================================
	individuals = {}
	file = open("../data/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_sample_IDs.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		id = get_individual_id(line)
		individual_rep[id] = 1

	file.close()



	file = open("../data/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", 'r')
	file1 = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'w')
	file.readline()
	file.readline()

	# get the effective list for selection
	index_map = {}
	line = (file.readline()).strip()
	line = line.split('\t')
	file1.write(line[0] + '\t')
	for i in range(2, len(line)):
		sample = line[i]
		individual = get_individual_id(sample)
		if individual not in individual_rep:
			continue
		else:
			file1.write(line[i] + '\t')
			index_map[i] = 1
	file1.write('\n')

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		file1.write(line[0] + '\t')
		for i in range(2, len(line)):
			if i not in index_map:
				continue
			else:
				file1.write(line[i] + '\t')
		file1.write('\n')

	file.close()
	file1.close()


	##===================================================================================================
	##==== counting samples in each tissue
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_#size
	##===================================================================================================


	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_sample_tissue_type", 'r')
	sample_tissue_map = {}
	while 1:
		line = file.readline()[:-1]	# can't strip '\t\t\t'
		if not line:
			break

		line = line.split('\t')

		if len(line) < 3:
			print line
			continue

		sample = line[0]
		tissue = line[2]

		sample_tissue_map[sample] = tissue
	file.close()



	# counting
	file = open("../data_processed/GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r')
	sample_list = (((file.readline()).strip()).split('\t'))[1:]
	file.close()

	print "there are",
	print len(sample_list),
	print "different samples from the rpkm file."

	counting = {}
	for sample in sample_list:
		tissue = sample_tissue_map[sample]
		if tissue not in counting:
			counting[tissue] = 1
		else:
			counting[tissue] += 1

	print "they are distributed among",
	print len(counting),
	print "different tissue types."

	file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count", 'w')
	for tissue in counting:
		file.write(tissue + '\t' + str(counting[tissue]) + '\n')
	file.close()


	#==== filtering the counts
	filter_list = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200]
	for i in range(len(filter_list)):
		filter = filter_list[i]
		file = open("../data_processed/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_" + str(filter), 'w')
		count = 0
		for tissue in counting:
			if counting[tissue] >= filter:
				count += 1
				file.write(tissue + '\t' + str(counting[tissue]) + '\n')
		print "# of tissue with sample size >= " + str(filter) + ":",
		print count
		file.close()
