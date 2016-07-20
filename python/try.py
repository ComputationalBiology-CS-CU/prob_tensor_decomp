import numpy as np




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



if __name__ == "__main__":



	count = 0
	
	for i in range(13):
		list_sample = np.load("data/Tensor_tissue_" + str(i) + "_list.npy")

		## DEBUG
		count += len(list_sample)

		rep_indiv = {}
		for j in range(len(list_sample)):
			sample = list_sample[j]
			individual = get_individual_id(sample)
			if individual in rep_indiv:
				print i,
				print individual
			else:
				rep_indiv[individual] = 1


	print count



