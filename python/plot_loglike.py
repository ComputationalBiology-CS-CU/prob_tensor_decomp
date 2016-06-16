import matplotlib.pyplot as plt

if __name__=="__main__":
	arr = []
	#with open("result/simu_loglike_20.txt", "r+") as f:
	with open("result/simu_loglike_100.txt", "r+") as f:
		i = 0
		for line in f:
			i += 1
			if i==1:
				continue
			a = float(line[:-1])
			#print a
			arr.append(a)

	plt.plot(arr)
	plt.xlabel("Number of Iterations")
	plt.ylabel("Joint Log Likelihood")
	plt.title("Joint Log Likelihood during Model Training")
	plt.show()


