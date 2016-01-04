//
// basic.hpp
// BayesianProbTensorDecomp
//
// function: some basic functions that maybe used by some parts of the whole program

#ifndef basic_hpp
#define basic_hpp

#include <utility>
#include <string>


using namespace std;


// trimming functions
void rtrim(char *);
void ltrim(char *);
void trim(char *);


// split the pair of (snp, R^2) and return pair variable
void pair_split(char *, pair<string, float> *);


// transform the string to char sequence
void StrToCharSeq(char *, string);


// transform the sample ID (like "GTEX-R55E-0826-SM-2TC5M") into individual ID (here is the first 9 digits)
string sample_to_individual(string);




#endif

// end of basic.hpp
