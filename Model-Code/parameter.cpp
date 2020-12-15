//************************************************************
//implementation file for parameter.h - (c)2007 Sudeep Samanta
//************************************************************
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include "randomize.h"
#include "parameter.h"


//converts string to distribution
Distribution string_to_dist(String d_string)
{
	switch (d_string[0])
	{
	case ('r'):
		return RANDOM;
		break;
	case ('R'):
		return RANDOM;
		break;
        case ('m'):
                return RANDOM;
                break;
	default:
		return FIXED;
	}
}

//constructors
Parameter::Parameter()
{
	strcpy(name, "NONE");
	dtype = FIXED;
	value = min = max = BAD_VALUE;
}

Parameter::Parameter(String& nmstr, double val)
{
	strcpy(name, nmstr);
	set_dtype(FIXED);
	value = min = max = val;
}

Parameter::Parameter(String& nmstr, Distribution dist, double val, double minval, double maxval)
{
	strcpy(name, nmstr);
	set_dtype(dist);
	if (minval == maxval)
	{
		set_dtype(FIXED);
	}
	if (dtype == FIXED)
	{
		value = min = max = val;
	}
	else
	{
		if (minval > maxval)
		{
			max = minval;
			min = maxval;
		}
		else
		{
			min = minval;
			max = maxval;
		}
		if (val < min)
		{
			value = min;
		}
		else if (val > max)
		{
			value = max;
		}
		else
		{
			value = val;
		}
	}
}

Parameter::Parameter(const Parameter &orig)
{
	orig.get_name(name);
	set_dtype(orig.get_dtype());
	value = orig.get_value();
	min = orig.get_min();
	max = orig.get_max();
}

Parameter::~Parameter()
{;}

//functions to set members
void Parameter::set_name(const String& nmstr)
{
	strcpy(name, nmstr);
}

void Parameter::set_dtype(Distribution dist)
{
	dtype = dist;
}

void Parameter::set_value(double val)
{
	value = val;
}

void Parameter::set_min(double minval)
{
	min = minval;
}

void Parameter::set_max(double maxval)
{
	max = maxval;
}

void Parameter::set_param(const String& nmstr, double val)
{
	set_name(nmstr);
	dtype = FIXED;
	set_value(val);
}

void Parameter::set_param(const String& nmstr, double val, double minval, double maxval)
{
	set_name(nmstr);
	dtype = RANDOM;
        if (minval == maxval)
	{
                set_dtype(FIXED);
	}
        if (dtype == FIXED)
	{
                value = min = max = val;
	}
        else
        {
                if (minval > maxval)
                {
                        max = minval;
                        min = maxval;
                }
                else
                {
                        min = minval;
                        max = maxval;
                }
                if (val < min)
		{
                        value = min;
		}
                else if (val > max)
		{
                        value = max;
		}
                else
		{
                        value = val;
		}
        }
}

//functions to access current members
void Parameter::get_name(String& nmstr) const
{
	strcpy(nmstr, name);
}

Distribution Parameter::get_dtype() const
{
	return dtype;
}

double Parameter::get_value() const
{
	return value;
}

double Parameter::get_min() const
{
	return min;
}

double Parameter::get_max() const
{
	return max;
}

Parameter &Parameter::operator =(const Parameter &rhs)
{
	if (this != &rhs)
        {
		rhs.get_name(name);
        	set_dtype(rhs.get_dtype());
        	value = rhs.get_value();
		min = rhs.get_min();
		max = rhs.get_max();
	}
	return *this;
}

//input from keyboard or file
istream& operator>>(istream& in, Parameter& param)
{
	char c = ' '; //to check if it is a space, initialized as ' '
	char diststr[MAX_LINE];
	String nmstr;
	Distribution dist;
	double val, minval, maxval;
//read the distribution
	if (&in == &cin)
	{
		cout << "Fixed<f or a value> or Random<r> parameter? "; 
	}
	while (c <= 32 || c > 126)//not a meaningful character
	{
		in.get(c);
	}
	in.putback(c); //c is no longer ' ' or meaningless
//if a value
	if((c >= '0' && c <= '9')||(c == '+' || c == '-'))
	{
		in >> val;
		param.set_dtype(FIXED);
	}
	else
	{
		c = ' '; //again initialize
		in >> diststr;
		if(&in == &cin)
			in.ignore(MAX_LINE, '\n');
		dist = string_to_dist(diststr);
		param.set_dtype(dist);
		switch(param.dtype)
		{
			case RANDOM:
				if (&in == &cin)
				{
					cout << "Starting value: ";
				}
				in >> val;
				if (&in == &cin)
				{
					cout << "Permissible range (min max): ";
				}
				in >> minval >> maxval;
				break;
			case FIXED:
				if (&in == &cin)
				{
					cout << "Value: ";
				}
				in >> val;
				break;
			default:
				break;
		}
	}
        if (&in == &cin)
        {
                cout << "Parameter name: ";
                in.get(nmstr,MAX_LINE);
                in.ignore(MAX_LINE, '\n');
        }
        else
        {
                in.get(c);
                if (c != '\n')
                {
                        while (c == ' '||c == '\t')//ignore spaces
			{
                                in.get(c);
			}
                        in.putback(c); //not white space
                        in.get(nmstr,MAX_LINE);
                }
        }
	switch(param.dtype)
	{
		case RANDOM:
			param.set_param(nmstr, val, minval, maxval);
			break;
		case FIXED:
			param.set_param(nmstr, val);
			break;
		default:
			strcpy(nmstr, "NONE");
			val = BAD_VALUE;
			param.set_param(nmstr, val);
			break;
	}
        in.ignore(MAX_LINE, '\n');
        return in;
}

ostream& operator<<(ostream& out, const Parameter& param)
{
        if (&out == &cout)
	{
		switch(param.dtype)
		{
			case RANDOM:
				out << param.name << " " << param.value << " "
			    	    << param.min << " " << param.max;
				break;
			case FIXED:
				out << param.name << " " << param.value;
				break;
			default:
				out << "unknown distribution!?";
				break;
		}
	}
        else
	{
                out << param.value;
	}
        return out;
}

int read_params(const String& param_file, int np, Parameter* parray)
{
	int i;
	ifstream in;
	if (open_infile(param_file, in) == 0)
	{
		return 0;
	}
	else
	{
		i = 0;
		while (i < np && (in.peek() != EOF))
		{
			in >> parray[i];
			i++;
		}
		in.close();
		return i;
	}
}

