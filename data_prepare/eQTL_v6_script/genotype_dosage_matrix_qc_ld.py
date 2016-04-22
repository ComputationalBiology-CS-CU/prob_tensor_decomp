## func: parsing the original dosage data into matrix format (chr x individual), and filtering based on QC and LD results




import os




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


	##========= get 450 individuals list ===========
	file = open("./genotype_processed/chrX.tfam", 'r')
	list_individual = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split(' ')
		individual = get_individual_id(line[0])
		list_individual.append(individual)
	file.close()




	for chr in range(1, 23):

		##==== create chr folder
		if not os.path.isdir("./genotype_450_dosage_matrix_qc/chr" + str(chr)):
			os.mkdir("./genotype_450_dosage_matrix_qc/chr" + str(chr))



		##==== load "./genotype_post_prune/chrX.prune.in", to get all SNPs after QC and LD
		rep_snp = {}	# NOTE: key data
		file = open("./genotype_post_prune/chr" + str(chr) + ".prune.in", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			snp = line
			rep_snp[snp] = 1
		file.close()

		print "working on chr#",
		print chr,
		print "and there are",
		print len(rep_snp),
		print "un-pruned SNPs."



		##==== parser dosage file, to fill data (snp name, snp pos, snp data for all individuals)
		list_snp_name = []
		list_snp_pos = []
		rep_individual_snp = {}
		for individual in list_individual:
			rep_individual_snp[individual] = []

		file = open("./genotype_processed/chr" + str(chr) + ".dosage", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split(' ')
			snp = line[1]
			pos = line[2]
			list_dosage = line[5:]
			
			if snp in rep_snp:
				list_snp_name.append(snp)
				list_snp_pos.append(pos)
				for i in range(len(list_dosage)):
					dosage = list_dosage[i]
					individual = list_individual[i]
					rep_individual_snp[individual].append(dosage)
		file.close()




		##==== create "SNP_info.txt" and "SNP_dosage_individualID.txt"
		file = open("./genotype_450_dosage_matrix_qc/chr" + str(chr) + "/SNP_info.txt", 'w')	# "snp_name snp_pos"
		for i in range(len(list_snp_name)):
			snp = list_snp_name[i]
			pos = list_snp_pos[i]
			file.write(snp + " " + pos + "\n")
		file.close()

		for individual in rep_individual_snp:
			filename = "./genotype_450_dosage_matrix_qc/chr" + str(chr) + "/SNP_dosage_" + individual + ".txt"
			file = open(filename, 'w')
			for i in range(len(rep_individual_snp[individual])):
				dosage = rep_individual_snp[individual][i]
				file.write(dosage + "\n")
			file.close()



