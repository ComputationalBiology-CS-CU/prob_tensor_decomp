## function:
##	split the full vcf file into header file and per-chromosome file



#==== global variable
line = ""



def get_chr(string):
	char = ''
	i = 0
	while 1:
		char += string[i]
		i=i+1
		if string[i] == '\t':
			break
	return char



if __name__ == "__main__":


	print "now parsing the genotype vcf file..."


	##========================================================================================================================
	## func: split the VCF file into header and files of each chromosome
	##========================================================================================================================
	filename = "../GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf"
	file_in = open(filename, 'r')


	##==== get the header
	n_header = 18					# the amount of lines in the header
	file_out = open("./genotype_vcf/header.txt", 'w')
	for i in range(n_header):
		line = file_in.readline()
		file_out.write(line)
	file_out.close()

	##==== parsing each chromosome
	line = file_in.readline()			# compensate for the first line
	for i in range(22):

		filename = "./genotype_vcf/chr" + str(i+1) + ".txt"
		file_out = open(filename, 'w')
		file_out.write(line)

		while 1:
			line = file_in.readline()
			if not line:
				break

			chr = get_chr(line)
			if chr != str(i+1):
				break
			else:
				file_out.write(line)

		file_out.close()

	##==== finish
	file_in.close()



