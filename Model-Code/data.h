//******************************************************************
//class definitions for a data vector and associated methods
//added agreement and efficiency functions - 7704
//******************************************************************
#ifndef _DATA_H
#define _DATA_H
#include <assert.h>
#include <math.h>
#include <iostream>
#include "util.h"

using namespace std;

const double INVALID = -999.0;	//check for invalid data
 
//*********************** class Data_Point *************************
//a simple class for storing an independent and a dependent variable
 
class Data_point
{
        public:
                Data_point();
                Data_point(double x, double y);
                ~Data_point();
 
                //data i/o
                void set_x(double x);
                void set_y(double y);
                void set_xy(double x, double y);
 
                double get_x() const;
                double get_y() const;
 
                //overloaded operators (compares the _x, independent variable)
		//used for sorting
                bool operator < (const Data_point& rhs) const;
                bool operator > (const Data_point& rhs) const;
                bool operator <= (const Data_point& rhs) const;
                bool operator >= (const Data_point& rhs) const;

		//assignment operator
		Data_point operator =(const Data_point& rhs);

                //output method
                friend std::ostream& operator << (std::ostream& out, const Data_point& dp);
 
                //compares _y
                bool is_bigger(double threshold) const;
 
        private:
                double _x;      //independent variable
                double _y;      //dependent variable
};

//*********************** class Regr_Res  ***************************
//a class of objects for storing a regression result
//the flag indicates if the regression is upto date
class Regr_Res
{
        public:
                //default constructor
                Regr_Res();
                Regr_Res(int flag, double slope, double intercept, double r_square);
                ~Regr_Res();

		//input method
		void set_values(int flag, double slope, double intercept, double r_square);
		void set_flag(int flag);

		//operator for assignment
		Regr_Res operator =(const Regr_Res &rhs);

		//access methods
		int check_flag() const;
		double get_slope() const;
                double get_intercept() const;
                double get_r_square() const;

		//output method
                friend std::ostream& operator << (std::ostream& out, const Regr_Res& res);
 
	private:
		int    _flag;
                double _slope;
                double _intercept;
                double _r_square;
};

//*********************** class Data_Frame **************************
//holds a series of data points in an array, methods for sorting
//construct "envelopes" etc. as required by Brent Ewers
class Data_Frame
{
        public:
                Data_Frame(); //default
                Data_Frame(int points); //with a specified # of points
		Data_Frame(const Data_Frame &rhs);
		~Data_Frame();

		//sets all data to INVALID
		void clear();
		//used for disposing a Data_Frame cleanly
		void dispose();
		//assignment operator
                Data_Frame operator =(const Data_Frame &rhs);

		friend Data_Frame operator +(const Data_Frame &arg1, const Data_Frame &arg2);

		//output method
		friend std::ostream& operator << (std::ostream& out, const Data_Frame& df);

                //data input methods
		void set_value(int index, double x_val, double y_val); //input at index
		void set_x_value(int index, double x_val);
		void set_y_value(int index, double y_val);

		//data access methods
		Regr_Res get_result() const; //returns the regression results
		int how_many() const; //returns the # of points in the frame
		double get_x(int index) const; //gets values @ index
                double get_y(int index) const;

		//summary stats calculation methods
		//remove invalid data before calling these
		double sum_x() const;	//sum of independent variable
		double sum_y() const;	//sum of dependent variable
		double sum_xsq() const;	//sum of square of independent variable
                double sum_ysq() const;	//sum of square of dependent variable
		double sum_xy() const;	//sum of xy product
		double mean_x() const;	//x mean
                double mean_y() const;	//y mean
                double sd_x() const;  //x standard deviation
                double sd_y() const;  //y standard deviation

		//calculation of errors
		double max_abs_err() const;	//maximum |x-y|
		double mae() const;		//mean abs error
		double rmse() const;		//root mean square error
		double efficiency() const;	//Coefficient of Efficiency - 7704
		double agreement() const;	//Index of Agreement - 7704

                //returns a Data_Frame with points where _x values >= x_value
                Data_Frame chop_x_left(double x_value);
                void ln_transform();    //does a natural log transform of _x
		void linear_reg();	//does a linear regression stores results in _result
		void shell_sort();	//shell sort ascending _x
		void remove_invalid_data(); //removes invalid data by throwing out -999 values

		//constructs an envelope required by Brent from an existing Data_Frame
		//w_int - width of the interval
		//sd_scale - construct the threshold by mean + sd_scale*std_dev
		//at least min_pts data points has to be in the envelope
		//uses in_interval and over_top methods
		Data_Frame make_envelope(double w_int, double sd_scale, int min_pts);

	private:
		//data members
		int _points; //# of points in the data set
		Regr_Res _result; //result of linear regression
		Data_point* _data; //array of Data_points

		//methods used in shell_sort()
		void shell_sort(Data_point* dat, int pts);
		void gap_sort(Data_point* dat, int s, int g, int n);

		//MAKE SURE THAT ONLY A SORTED Data_Frame is used
		//for the following methods

                //returns the # of data points
                //within [min,max) interval for _x
                //if max is the max value of _x, interval is [min, max]
		//the index corresponding to the min value is passed by ref
		//returns -1 if min higher than the highest value
		//or max lower than the lowest value or invalid interval
                int in_interval(double min, double max, int& start_index);

		//constructs a new Data_Frame with given number of data points
		//starting at start_index
		Data_Frame segment(int start_index, int number);
 
                //constructs a new Data_Frame with the data points
                //that have _y greater than top
                Data_Frame over_top(double top);

};

#endif
