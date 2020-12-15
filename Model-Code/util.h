//**************************************************************
//definition of basic stat functions & file utilities
//(c)2007 Sudeep Samanta - 1/4/07
//**************************************************************
#ifndef UTIL_H
#define UTIL_H

#include <fstream>
#include <iomanip>
#include <string.h>
#include <math.h>

//some elementary functions
#define max(a,b)        ((a) > (b) ? (a) : (b))
#define min(a,b)        ((a) < (b) ? (a) : (b))

using namespace std;

//***********************
//DEFINITIONS
//***********************
const int MAX_LINE  = 1024;       //max chars on a line
typedef char String[MAX_LINE+1]; //a string of characters

//************************************
//some elementary numerical functions
//***********************************
int big_int(int oneint, int othrint);
int small_int(int oneint, int othrint);
double absdbl(double x); //returns abs value

//***********************************
//FILE HANDLING AND RELATED FUNCTIONS
//***********************************
//opens file for input - returns 0 if fails to open
int open_infile(const String& filename, std::ifstream& in);

//opens file for output
int open_outfile(const String& filename, std::ofstream& out);

//increments a string of digits by '1' to add as extension to a file name
void increment_string(String& thestr);

//returns the number of lines in a file
//designed for use with files like .p - a line may start with digit or ch
int lines(const String &filename);

//returns the number of lines of data in a file
//stops if a line starts with a character other than a digit
//unless it's the first line (header)
int datalines(const String &filename);

//returns the number of columns in a file
int columns(const String &filename);

//******************************
//conversions & stuff
//******************************
//checks if year is a leap year
//1 true, 0 false
int leap_year(int year);

//returns radian conversion of deg
double deg2rad(double deg);

//returns degree conversion of rad
double rad2deg(double rad);

//returns the conversion to K from deg C
double C2K(double deg_C);

//returns the conversion to deg C from K
double K2C(double Kelvin);

//returns the longitude correction in hours at the location lat long
double long_corr(double lat, double longi);

//******************************
//sorts in ascending order
//******************************
//n_items is the size
void shell_sort(double* data_array, int n_items);
void gap_sort(double* data_array, int s, int g, int n);

//***********************************
//STAT & MATRIX FUNCTIONS
//***********************************
double median(double* data_array, int n_items);
double mean(double* data_array, int n_items);
double quantile(double* data_array, double quant, int n_items);

//call if the array is already sorted
double sorted_qntl(double* data_array, double quant, int n_items);

double calc_var(double* data_array, int n_items);
double calc_covar(double* datax, double* datay, int n_items);

//returns the covariance matrix in covmat by reference
//nvar in # of variables (col), nval is # of values for each var
//vardata is nvar x nval, covmat is nvar x nvar
void calc_covmat(double** vardata, double** covmat, int nvar, int nval);

//calculates the cholesky factor (lower tringular) by decomposition
//returns by ref in choleskyf (initialized before) 
void calc_cholesky(double** covmat, double** choleskyf, int ndim);

//multiplies a square matrix with a scalar
void scale_sqmatrix(double** sqmat, double scalef, int ndim);

//copies a vector into another
void copy_vect(double* origv, double* copyv, int ndim);

//TO BE REMOVED
//retains a history of past selected parameters
//no_par rows, history+1 cols, the sd of each row is updated
//by this function and kept at col index history for use in mcmc
void calc_stdev(double** par_history, int no_par, int history);

#endif

