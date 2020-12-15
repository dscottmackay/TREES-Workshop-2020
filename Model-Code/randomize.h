//****************************************************
//This is the function definitions for generating
//random numbers - implementation in randomization.cpp
//2/8/2007 (c)Sudeep Samanta
//****************************************************
#ifndef _RANDOMIZE_H
#define _RANDOMIZE_H

using namespace std;

//returns a pseudo random number in (0,1)
//the argument is the number of simulation run
//so that the random number is along a sequence
//if the seed happens to remain same

typedef unsigned short Xsubi_array[3];

//initializes the array for randomization
void initialize(Xsubi_array& xsubi);

double unit_random(Xsubi_array& xsubi);

//returns a pseudo random number between min and max
//from a uniform distribution
double uniform_random(double min, double max, Xsubi_array& xsubi);

//returns a pseudo random number from a normal distribution
//with specified mean and standard deviation
double stdNormal(Xsubi_array& xsubi);
double normal_random(double mean, double std_dev, Xsubi_array& xsubi);

//get p density normal for x with mean and sd
double p_norm(double x, double mu, double sd);

//deals with vectors - containers are passed in to keep memeory management
//at the top level, values passed by reference
//ndim i.i.d. N(0,1) variables
void stdNormal(double* stdNvect, int ndim, Xsubi_array& xsubi);
//ndim correlated normal variables with 0 mean and covariance matrix with lower triangular CholeskyF
void corNormal(double* corNvect, double* stdNvect, double** CholeskyF, int ndim, Xsubi_array& xsubi);
//above with meanvect
void corNormal(double* corNvect, double* stdNvect, double* meanvect, double** CholeskyF, int ndim, Xsubi_array& xsubi);

//get p-value for Gamma distribution with theta, alpha and beta
double p_gamma(double theta, double alpha, double beta);

//given a p-value and parameters and interval (accuracy), get the value of theta
double val_gamma(double accuracy, double p_val, double alpha, double beta);

//get p-value for inverse Gamma distribution with theta, alpha and beta
double p_invgamma(double theta, double alpha, double beta);


#endif
