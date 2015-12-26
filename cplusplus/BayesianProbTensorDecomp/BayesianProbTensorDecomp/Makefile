GSL_HOME=/Users/rachelren/Documents/prob_tensor_decomp/dev_CP/Bayesian_Prob_Tensor_Decomp/gsl-2.1
GSL_LIB=/usr/local/lib

main:
	g++ -c -I $(GSL_HOME)/ main.cpp
	g++ -o main -L$(GSL_LIB) -lgsl -lgslcblas main.o

clean:
	rm *.o main
