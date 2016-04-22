## some mics scripts to process the source data

# libraries
import sys
import time
import os




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




def get_gene_id(s):
	id = ''
	count = 0
	for i in range(len(s)):
		if s[i] == '"':
			count = i + 1
			break

	for i in range(count, len(s)):
		if s[i] == '"':
			break
		id += s[i]

	return id





if __name__ == "__main__":

	"""
	rep = {}

	n_tissue = 33
	for i in range(n_tissue):
		index_tissue = i + 1

		file = open("/Users/shuoyang/Desktop/Genetics_GeneExpression/GTEx/workbench_v.6/data_processed/expression_by_etissue/tissue_" + str(index_tissue) + ".txt", 'r')

		line = (file.readline()).strip()
		line = line.split('\t')
		for sample in line[1:]:
			individual = get_individual_id(sample)
			if individual not in rep:
				rep[individual] = 1
		file.close()


	print len(rep)
	#print rep
	"""





	##==========================================================
	##==== get the chr and pos of all genes from annotation file
	##==========================================================
	file = open("/Users/shuoyang/Desktop/Genetics_GeneExpression/GTEx/workbench_v.6/data_source/gencode.v19.genes.patched_contigs.gtf", 'r')
	file1 = open("/Users/shuoyang/Desktop/Genetics_GeneExpression/GTEx/workbench_v.6/data_processed/gencode.v19.genes.patched_contigs.gtf_gene", 'w')
	for i in range(5):
		file.readline()

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		type = line[2]
		if type == 'transcript':
			chr = line[0]
			start = line[3]
			end = line[4]
			gene = line[8]
			gene = get_gene_id(gene)
			#print chr
			#print start
			#print end
			#print gene
			file1.write(chr + '\t' + start + '\t' + end + '\t' + gene + '\n')

	file.close()
	file1.close()


