## function:
## 1. parsing the VCF file, to get the following for each chromosome:
##	chrX.tped
##	chrX.dosage
##	chrX.maf05.exclusion.snplist.txt
## 2. get the individual information: (NOTE: we already have had this from the dataset, thus no need)
##	chrX.tfam



## notes (file format):
'''
*.tped: (white-space delimited file)
	chromosome (1-22, X, Y or 0 if unplaced)
	rs# or snp identifier
	Genetic distance (morgans; 0 as not necessary for us now)
	Base-pair position (bp units)
	genotype list (e.g. 1,2,3,4 or A,C,G,T or anything else except 0 which is, by default, the missing genotype character)

*.tfam: (white-space delimited file)
	Family ID
	Individual ID
	Paternal ID (0)
	Maternal ID (0)
	Sex (1=male; 2=female; other=unknown)
	Phenotype (-9: missing)

*.dosage: (white-space delimited file)
	chromosome (1-22, X, Y or 0 if unplaced)
	rs# or snp identifier
	Base-pair position (bp units)
	genotype dosage list
'''


# libraries
import sys
import time
import os





#==== global variable
line = ""





if __name__ == "__main__":



	print "now creating the required genotype files..."

	##========================================================================================================================
	## func: extract the genotype information:
	##	chrX.tped
	##	chrX.dosage
	##	chrX.maf05.exclusion.snplist.txt
	##========================================================================================================================
	#chr = 1			# TODO: for different chromosomes



	### a simple local test case
	#filename = "../genotype/test_chr22_main.vcf"
	#file_in = open(filename, 'r')
	#filename = "../genotype/test_chr22.tped"
	#file_tped = open(filename, 'w')
	#filename = "../genotype/test_chr22.dosage"
	#file_dosage = open(filename, 'w')
	#filename = "../genotype/test_chr22.elist"
	#file_elist = open(filename, 'w')



	for chr in range(1, 23):


		filename = "./genotype_vcf/chr" + str(chr) + ".txt"
		filename_vcf = "./genotype_vcf/chr" + str(chr) + ".txt"
		file_in = open(filename, 'r')
		filename = "./genotype_processed/chr" + str(chr) + ".tped"
		file_tped = open(filename, 'w')
		filename = "./genotype_processed/chr" + str(chr) + ".dosage"
		file_dosage = open(filename, 'w')
		filename = "./genotype_processed/chr" + str(chr) + ".maf05.exclusion.snplist.txt"
		file_elist = open(filename, 'w')


		while 1:
			line = (file_in.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			chr = line[0]
			pos = line[1]
			snp = line[2]
			allele1 = line[3]
			allele2 = line[4]
			#
			filter = line[6]
			#
			type_list = line[8].split(':')		# retrieve "GT:GL:DS" information
			index_genotype = 0
			index_dosage = 0
			for count in range(len(type_list)):
				type = type_list[count]
				if type == 'GT':
					index_genotype = count
				if type == 'DS':
					index_dosage = count
			genotype_list = line[9:]


			##== tped/dosage
			file_tped.write(chr + ' ')
			file_tped.write(snp + ' ')
			file_tped.write('0' + ' ')		# the genetic position
			file_tped.write(pos + '  ')

			file_dosage.write(chr + ' ')
			file_dosage.write(snp + ' ')
			file_dosage.write(pos + ' ')
			file_dosage.write(allele1 + ' ')
			file_dosage.write(allele2 + ' ')

			for element in genotype_list:
				element = element.split(':')
				pair = element[index_genotype].split('/')
				dosage = element[index_dosage]

				# allele#1
				if pair[0] == '0':
					file_tped.write(allele1 + ' ')
				elif pair[0] == '1':
					file_tped.write(allele2 + ' ')
				else:
					file_tped.write("0" + ' ')

				# allele#2
				if pair[1] == '0':
					file_tped.write(allele1 + ' ')
				elif pair[1] == '1':
					file_tped.write(allele2 + ' ')
				else:
					file_tped.write("0" + ' ')

				file_dosage.write(dosage + ' ')
			
			file_tped.write('\n')
			file_dosage.write('\n')

			##== elist (NOTE: we start from info04;maf01;hwe6, thus we only need maf05)
			if filter == 'maf05':
				file_elist.write(snp + '\n')


		file_in.close()
		file_tped.close()
		file_dosage.close()
		file_elist.close()



		##==== delete vcf source file of this chromosome
		os.system("rm " + filename_vcf)



