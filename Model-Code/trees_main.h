//**********************************************************************************
//added dependence on Ca - March, 2007
//Multivariate normal proposal distribution, trace option removed - 2/12/07
//LINEAR D, lfscl, Q, T, and swc - 10/10/06
//parameter values turn constraints off & on - 2/12/07
//adding additional data col swp in kPa (negative)
//made more compact, uses canopy D for gs and above canopy D for pm and gHa - 271004
//(c)2007 - Sudeep Samanta
//Added carbon cycling (c)2007, 2008 - Scott Mackay
//Added plant water balance model (c) 2010, 2011 David Roberts & Scott Mackay
//Added canopy phenology (c) 2013, 2014 Phil Savoy & Scott Mackay
//**********************************************************************************
#ifndef _MAIN_H
#define _MAIN_H

#include "simulator2.h"
#include "time_keeper.h"
#include "data.h"
#include "data_store.h"
#include "parameter.h"
#include "randomize.h"
#include "state_store.h"

using namespace std;

//calculate p(y|theta)p(theta) - using normal errors with equal variance
//called from simdensity
double logpostdens(const Data_Frame& df, double sd, double& Dtheta);

double logpostdens(const Data_Frame& df1, const Data_Frame& df2, double sd1, double sd2, 
						double sd1_weight, double sd2_scalar, double& Dtheta);

//this function calls the deterministic simulation functions for data and given params
//returns the p(y|theta)p(theta) value
double simdensity(Data_Store& inputdata, State_Store& state, Parameter* parameters, 
						trees_params& treesParams, Data_Frame& Ecframe, 
						Data_Frame& NEEframe, double& Dtheta, int n_agg);

double simdensity(double weight, Data_Store& inputdata, State_Store& state, Parameter* parameters, 
						trees_params& treesParams, Data_Frame& Ecframe, 
						Data_Frame& NEEframe, double& Dtheta, int n_agg);

double simdensityEc(Data_Store& inputdata, State_Store& state, Parameter* parameters, 
						Data_Frame& Ecframe, trees_params& treesParams, 
						Data_Frame& NEEframe, double& Dtheta);

double simdensityNEE(Data_Store& inputdata, State_Store& state, Parameter* parameters, 
						Data_Frame& Ecframe, trees_params& treesParams, 
						Data_Frame& NEEframe, double& Dtheta);

//generate candidate parameter values returned in candval
//the generated values are correlated using Cholesky factor cholf (updated) with mean @ curpval
//values checked against minlim and maxlim to keep in admissible range
void gen_proppval(double* candpval, double* curpval, double* tempval, double* minlim, 
					double* maxlim, int rp, double** cholf, Xsubi_array& xsubi);

//returns the number of values in candpval out of range, called by gen_proppval
int check_range(double* candpval, double* minlim, double* maxlim, int rp);

//copies parameter values generated from proposal distribution into parameter set
//candp has total np parameters, rp is # of randomized parameters in candval
void copy_proppval(Parameter* candp, double* candpval, int np, int rp);

//copy np parameters from one parameter array origp into another
void copy_params(Parameter* origp, Parameter* copyp, int np);

#endif
