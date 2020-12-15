//***************************************************************************
//model parameters - stores name, fixed or random, current, min & max values
//distributions done outside tracking the chain - use for MCMC only
//(c)2007 Sudeep Samanta
//***************************************************************************
#ifndef _PARAMETER_H
#define _PARAMETER_H
#include <fstream>
#include "randomize.h"
#include "util.h"

using namespace std;

enum Distribution {FIXED, RANDOM};
const double BAD_VALUE = -99999.9;
//converts a String to a Distribution type
Distribution string_to_dist(String); 

class Parameter
{
	public:
		Parameter();
		Parameter(String& nmstr, double val);
		Parameter(String& nmstr, Distribution dist, double val, double minval, double maxval);
		Parameter(const Parameter &orig);
		~Parameter();
//functions to set members
		void set_name(const String& nmstr);
		void set_dtype(Distribution dist);
		void set_value(double val);
		void set_min(double minval);
		void set_max(double maxval);
		void set_param(const String& nmstr, double val); //for fixed
		void set_param(const String& nmstr, double val, double minval, double maxval);
//functions to access current members
		void get_name(String& nmstr) const;
		Distribution get_dtype() const;
		double get_value() const;
		double get_min() const;
		double get_max() const;
//overloaded operators
		Parameter & operator =(const Parameter &rhs);
		friend std::istream& operator>>(std::istream&, Parameter&);
//prints the value of a parameter (name also on screen)
		friend std::ostream& operator<<(std::ostream&, const Parameter&);
	private:
		String name;
		Distribution dtype;
		double value;
		double min;
		double max;
};

//reads parameters into an array and returns the # of parameters read
//from file. needs the size of parray in np
int read_params(const String& param_file, int np, Parameter* parray);

#endif

