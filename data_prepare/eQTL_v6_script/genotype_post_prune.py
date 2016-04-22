## function: get the xxx.ld file as input, output all strongly associated SNPs (un-pruned) with their target SNPs (pruned)
##		then, reverse the list, get all the strongly associated SNPs, with all their strongly associated pruned SNPs
## algorithm:
##	1. if there are no other associated SNPs (un-pruned SNPs), simply drop this pruned SNP;
##	2. else, pick up the most associated SNP (un-pruned SNP);
##	3. build a hashtable with the associated SNP as the key, and the pruned SNP (maybe several) as the value


chr_num = ''
pruned_rep = {}


if __name__ == "__main__":



	##=========== record all the pruned SNPs ===========
	print "record all the pruned SNPs"
	file = open("plink.prune.out", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		pruned_rep[line] = 1
	file.close()


	##=========== find out the most associated un-pruned SNP (if there is) for all pruned SNPs ===========
	print "find out the most associated un-pruned SNP (if there is) for all pruned SNPs"
	file = open("plink.ld", 'r')
	line = file.readline()  # remove the first unnecessary line

	rep = {}

	count = 0
	while 1:
		line = (file.readline()).strip()
		count += 1
		if not line:
			break

		line = line.split()  # remove the first unnecessary line

		if count == 1:
			chr_num = line[0]

		pruned = line[2]
		associated = line[5]
		r2 = float(line[6])

		if pruned not in rep:
			if associated in pruned_rep:
				rep[pruned] = (associated, -1)
			else:
				rep[pruned] = (associated, r2)
		else:
			if associated in pruned_rep:
				continue
			else:
				if r2 > rep[pruned][1]:
					rep[pruned] = (associated, r2)

	file.close()


	##=========== remove pruned SNPs which are not associated with any un-pruned SNP ===========
	print "remove pruned SNPs which are not associated with any un-pruned SNP"
	remove_list = []
	for pruned in rep:
		if rep[pruned][1] == -1:
			remove_list.append(pruned)
	# TODO we may need this number
	print len(remove_list),
	print "will be removed..."

	for remove in remove_list:
		del rep[remove]


	##=========== reversed hashing ===========
	print "reversed hashing and saving the file"
	rep1 = {}
	for pruned in rep:
		associated = rep[pruned][0]
		r2 = rep[pruned][1]
		if associated in rep1:
			rep1[associated].append((pruned, r2))
		else:
			rep1[associated] = [(pruned, r2)]

	file = open("chr" + chr_num + ".post_prune.txt", 'w')
	for associated in rep1:
		file.write(associated + '\t')
		for pair in rep1[associated]:
			pruned = pair[0]
			r2 = pair[1]
			file.write(pruned + ' ' + str(r2) + '\t')
		file.write('\n')
	file.close()
