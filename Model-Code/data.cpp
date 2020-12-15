//******************************************************
//data.cpp - implementation of classes defined in data.h
//******************************************************
#include <iostream>
#include <fstream>
#include "data.h"
#include "util.h"

 
//*********************** class Data_Point *************************
Data_point::Data_point():_x(INVALID), _y(INVALID)
{ }
 
Data_point::Data_point(double x, double y)
        :_x(x), _y(y)
{ }
 
Data_point::~Data_point()
{
        _x = _y = 0.0;
}

void Data_point::set_x(double x)
{
        _x = x;
}
 
void Data_point::set_y(double y)
{
        _y = y;
}
 
void Data_point::set_xy(double x, double y)
{
        _x = x;
        _y = y;
}
 
double Data_point::get_x() const
{
        return _x;
}
 
double Data_point::get_y() const
{
        return _y;
}


bool Data_point::operator < (const Data_point& rhs) const
{
        return _x < rhs._x;
}
 
bool Data_point::operator > (const Data_point& rhs) const
{
        return _x > rhs._x;
}

bool Data_point::operator <= (const Data_point& rhs) const
{
        return _x <= rhs._x;
}
 
bool Data_point::operator >= (const Data_point& rhs) const
{
        return _x >= rhs._x;
}

Data_point Data_point::operator =(const Data_point& rhs)
{
	if (this != &rhs)
        {
		_x = rhs._x;
		_y = rhs._y;
	}
	return *this;
}

ostream& operator << (ostream& out, const Data_point& dp)
{
        if (&out == &cout)
	{
                out << "x = "  << setw(14) << setprecision(10) << dp._x
		    << " y = " << setw(14) << setprecision(10) << dp._y << endl;
	}
	else
	{
		out << setw(14) << setprecision(10) << dp._x
                    << setw(14) << setprecision(10) << dp._y << endl;
	}
	return out;
}

bool Data_point::is_bigger(double threshold) const
{
        return _y > threshold;
}

//*********************** class Regr_Res  ***************************
Regr_Res::Regr_Res() : _flag(0), _slope(0.0), _intercept(0.0), _r_square(0.0)
{ }

Regr_Res::Regr_Res(int flag,
                   double slope,
                   double intercept,
                   double  r_square)
        :_flag(flag), _slope(slope), _intercept(intercept), _r_square(r_square)
{ }

Regr_Res::~Regr_Res()
{ }

void Regr_Res::set_values(int flag, double slope, double intercept, double r_square)
{
	_flag      = flag;
	_slope     = slope;
	_intercept = intercept;
	_r_square  = r_square;
}

void Regr_Res::set_flag(int flag)
{
	_flag      = flag;
        _slope     = 0.0;
        _intercept = 0.0;
        _r_square  = 0.0;
}

Regr_Res Regr_Res::operator =(const Regr_Res &rhs)
{
	if (this != &rhs)
	{
        	_flag      = rhs._flag;
        	_slope     = rhs._slope;
        	_intercept = rhs._intercept;
        	_r_square  = rhs._r_square;
	}
        return *this;
}

int Regr_Res::check_flag() const
{
	if (_flag > 0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

double Regr_Res::get_slope() const
{
	return _slope;
}
 
double Regr_Res::get_intercept() const
{
        return _intercept;
}

double Regr_Res::get_r_square() const
{
	return _r_square;
}

ostream& operator << (ostream& out, const Regr_Res& res)
{
        if (&out == &cout)
	{
                out << "slope = " << res._slope
                    << " intercept = " << res._intercept
                    << " R square = " << res._r_square;
	}
        else
	{
                out << setw(10) << setprecision(6) << res._slope << '\t'
                    << setw(10) << setprecision(6) << res._intercept << '\t'
                    << setw(10) << setprecision(6) << res._r_square;
//        out << endl;
	}
        return out;
}

//*********************** class Data_Frame **************************
Data_Frame::Data_Frame()
{
	_points = 0;
	_data   = NULL;
}

Data_Frame::Data_Frame(int points)
{
	_points = points;
	_data   = new Data_point[points];
}

Data_Frame::Data_Frame(const Data_Frame &rhs)
{
	_points = rhs._points;
	_result = rhs._result;
	_data   = new Data_point[_points];
	int i;
	for (i = 0; i < _points; i++)
	{
		_data[i] = rhs._data[i];
	}
}
	
Data_Frame::~Data_Frame()
{
	dispose();
}


void Data_Frame::dispose()
{
	delete [] _data;
	_data   = NULL;
	_points = 0;
}

void Data_Frame::clear()
{
	for (int i = 0; i < _points; i++)
	{
                        _data[i].set_xy(INVALID, INVALID);
	}
}

Data_Frame Data_Frame::operator =(const Data_Frame &rhs)
{
        if (this != &rhs)
        {
                dispose();
                _points = rhs._points;
                _result = rhs._result;
                _data   = new Data_point[_points];
                int i;
                for(i = 0; i < _points; i++)
                        _data[i] = rhs._data[i];
        }
        return *this;
}

Data_Frame operator +(const Data_Frame &arg1, const Data_Frame &arg2)
{
	int new_pts;	//# of points in new object
	int i, j;	//counters
	new_pts = arg1._points + arg2._points;
	Data_Frame df_new(new_pts); //the new Data_Frame

	if (new_pts != 0) //data in at least one of the arguments
	{
		if (arg1._points == 0) //all the data in arg2
		{
			df_new = arg2;
		}
		else if (arg2._points == 0) //all the data in arg1
		{
			df_new = arg1;
		}
		else	//data in both arguments
		{
			//get the data from arg1
			j = 0;
			for (i = 0; i < arg1._points; i++)
			{
				df_new.set_value(j, arg1._data[i].get_x(), arg1._data[i].get_y());
				j++;
			}
                        //get the data from arg2
                        for (i = 0; i < arg2._points; i++)
                        {
                                df_new.set_value(j, arg2._data[i].get_x(), arg2._data[i].get_y());
                                j++;
                        }
		}//end else
	}//end if
	return df_new;
} 

ostream& operator << (ostream& out, const Data_Frame& df)
{
	int i;
	for (i = 0; i < df._points; i++)
	{
		out << df._data[i];
	}
	return out;
}

void Data_Frame::set_value(int index, double x_val, double y_val)
{
	if (index < _points)	//index within bounds
	{
		_data[index].set_xy(x_val, y_val);
		_result.set_flag(0);	//invalid regression results due to input
	}
	else
	{
		cout << "ERROR! index out of bounds." << endl;
	}
}

void Data_Frame::set_x_value(int index, double x_val)
{
        if (index < _points)     //index within bounds
        {
                _data[index].set_x(x_val);
                _result.set_flag(0);    //invalid regression results due to input
        }
        else
	{
                cout << "ERROR! index out of bounds." << endl;
	}
}

void Data_Frame::set_y_value(int index, double y_val)
{
        if (index < _points)     //index within bounds
        {
                _data[index].set_y(y_val);
                _result.set_flag(0);    //invalid regression results due to input
        }
        else
	{
                cout << "ERROR! index out of bounds." << endl;
	}
}

Regr_Res Data_Frame::get_result() const
{
	return _result;
}

int Data_Frame::how_many() const
{
	return _points;
}

double Data_Frame::get_x(int index) const
{
	if (index >= _points)
	{
		return INVALID;
	}
	else
	{
		return _data[index].get_x();
	}
}

double Data_Frame::get_y(int index) const
{
        if (index >= _points)
	{
                return INVALID;
	}
        else
	{
                return _data[index].get_y();
	}
}

double Data_Frame::sum_x() const
{
        double sum = 0.0; //sum of the numbers

	//sum the values 
	for (int i = 0; i < _points; i++)
	{
		sum += _data[i].get_x();
	}
	return sum;
}

double Data_Frame::sum_y() const
{
        double sum = 0.0; //sum of the numbers
 
        //sum the values 
        for (int i = 0; i < _points; i++)
	{
                sum += _data[i].get_y();
 	}
        return sum;
}

double Data_Frame::sum_xsq() const
{
        double sum = 0.0; //sum of the numbers
	double tempdat;   //temporary variable
 
        //sum the values 
        for (int i = 0; i < _points; i++)
	{
		tempdat = _data[i].get_x();
                sum += tempdat*tempdat;
	}
        return sum;
}

double Data_Frame::sum_ysq() const
{
        double sum = 0.0; //sum of the numbers
        double tempdat;   //temporary variable
 
        //sum the values 
        for (int i = 0; i < _points; i++)
        {
                tempdat = _data[i].get_y();
                sum += tempdat*tempdat;
        }
        return sum;
}

double Data_Frame::sum_xy() const
{
        double sum = 0.0; //sum of the numbers
 
        //sum the values 
        for (int i = 0; i < _points; i++)
	{
                sum += _data[i].get_x()*_data[i].get_y();
 	}
        return sum;
}

double Data_Frame::mean_x() const
{
	return sum_x()/(double)_points;
}

double Data_Frame::mean_y() const
{
        return sum_y()/(double)_points;
}

double Data_Frame::sd_x() const
{
        double sum_diff = 0.0; //sum of (Xi - Xbar)^2
        double data_mean = mean_x(); //data mean
        double curr_diff; //current value of (Xi - Xbar)
 
        for (int i = 0; i < _points; i++)
        {
                curr_diff = _data[i].get_x() - data_mean;
                curr_diff = curr_diff*curr_diff;//square of the difference
                sum_diff += curr_diff;
        }
 
        return sqrt(sum_diff/(double)(_points - 1));
}

double Data_Frame::sd_y() const
{
        double sum_diff = 0.0; //sum of (Yi - Ybar)^2
        double data_mean = mean_y(); //data mean
        double curr_diff; //current value of (Yi - Ybar)
 
        for (int i = 0; i < _points; i++)
        {
                curr_diff = _data[i].get_y() - data_mean;
                curr_diff = curr_diff*curr_diff;//square of the difference
                sum_diff += curr_diff;
        }
 
        return sqrt(sum_diff/(double)(_points - 1));
}

double Data_Frame::max_abs_err() const
{
	double error, new_err;
	error = 0.0;
	for (int i = 0; i < _points; i++)
        {
		new_err = absdbl(_data[i].get_x() - _data[i].get_y());
		if (new_err  > error)
		{
			error = new_err;
		}
	}
	return error;
}

double Data_Frame::mae() const
{
        double sumae;
        sumae = 0.0;
        for (int i = 0; i < _points; i++)
	{
                sumae += absdbl(_data[i].get_x() - _data[i].get_y());
	}

        return sumae/(double)_points;
}

double Data_Frame::rmse() const
{
        double sumse;
        sumse = 0.0;
        for (int i = 0; i < _points; i++)
	{
                sumse += pow((_data[i].get_x() - _data[i].get_y()), 2.0);
	}
 
        return sqrt(sumse/(double)_points);
}

//Coefficient of Efficiency - 7704
double Data_Frame::efficiency() const
{
	double ssqoipi, ssqoiobar, cofeff;
	double obsi;
	double obar = mean_x();
	ssqoipi = ssqoiobar = 0.0;
	for (int i = 0; i < _points; i++)
	{
		obsi = _data[i].get_x();
		ssqoipi += pow((obsi-_data[i].get_y()), 2.0);
		ssqoiobar += pow((obsi-obar), 2.0);
	}
	cofeff = 1.0-(ssqoipi/ssqoiobar);
	return cofeff;
}

//Index of Agreement - 7704
double Data_Frame::agreement() const
{
	double ssqoipi, ssqdenomterm, denomterm, iofagr;
	double obsi, predi;
	double obar = mean_x();
	ssqoipi = ssqdenomterm = 0.0;
	for (int i = 0; i < _points; i++)
        {
		obsi = _data[i].get_x();
		predi = _data[i].get_y();
		ssqoipi += pow((obsi - predi), 2.0);
		denomterm = pow(absdbl(predi-obar) + absdbl(obsi-obar),2.0);
		ssqdenomterm += denomterm;
	}
	iofagr = 1.0-(ssqoipi/ssqdenomterm);
	return iofagr;
} 

void Data_Frame::ln_transform()
{
        for (int i = 0; i < _points; i++)
	{
		_data[i].set_x(log(_data[i].get_x()));//set _x at i to log
	}
}

Data_Frame Data_Frame::chop_x_left(double x_value)
{
	int count;
	count = 0;
	for (int i = 0; i < _points; i++)
	{
		if (_data[i].get_x() >= x_value)
		{
			count++;
		}
	}
	Data_Frame df(count);
	count = 0;
	for (int i = 0; i < _points; i++)
	{
                if (_data[i].get_x() >= x_value)
		{
			df.set_value(count, _data[i].get_x(), _data[i].get_y());
			count++;
		}
	}
	return df;
}
			

void Data_Frame::linear_reg()
{
	int flag; //flag value
	double slope, intercept, r_sq;
	double temp1, temp2, temp3; //temporary variables

	//do it only the current regression is invalid
	if (!_result.check_flag())
	{
		//numerator term for slope using least squares fit
		temp1 = sum_xy() - sum_x()*sum_y()/(double)_points;

		//denominator term for slope
		temp3 = sum_x();
		temp2 = sum_xsq() - temp3*temp3/(double)_points;
		slope = temp1/temp2;

		//calculate intercept
		intercept = mean_y() - slope*mean_x();

		//to get R square the predictions have to be constructed
		//first step is to calculate the SSY
		temp1 = 0.0; //SSY will be in temp1
		temp2 = mean_y();
		for (int i = 0; i < _points; i++)
		{
			temp3  = _data[i].get_y() - temp2; //Yi - Ybar
			temp1 += temp3*temp3;
		}

		//then get the SSE, sum sq error
		//and put it in temp2
		temp2 = 0.0;
		for (int i = 0; i < _points; i++)
                {
			temp3 = _data[i].get_y() - (intercept + slope*_data[i].get_x());
			temp2 += temp3*temp3;
                }

		r_sq = (temp1 - temp2)/temp1;

		//put these values in results
		flag = 1; //validate flag
		_result.set_values(flag, slope, intercept, r_sq);
	}//end if
}

//****************** methods used for sorting ***************************
void Data_Frame::shell_sort()
{
	shell_sort(_data, _points);
}

void Data_Frame::shell_sort(Data_point* dat, int pts)
{
	int gap = pts/3 + 1;

	while (gap > 1)
	{
		for (int sublist = 0; sublist < gap; sublist++)
		{
			gap_sort(dat, sublist, gap, pts);
			gap = gap/3 + 1;
		}
	}
	gap_sort(dat, 0, 1, pts);
}

void Data_Frame::gap_sort(Data_point* dat, int s, int g, int n)
{
	int unsorted, loc;
        Data_point nextitem;
        for (unsorted = s+g; unsorted < n; unsorted += g)
        {
                nextitem = dat[unsorted];
                loc = unsorted;
 
                for (; loc > s && dat[loc - g] > nextitem; loc -=g)
		{
                        dat[loc] = dat[loc-g];
		}
                dat[loc] = nextitem;
        }
}

//method for removing invalid data
void Data_Frame::remove_invalid_data()
{
	int j;//counter
	int valid_points = 0;
	for (int i = 0; i < _points; i++)
	{
		if (_data[i].get_x() != INVALID && _data[i].get_y() != INVALID)//invalid data is INVALID
		{
			valid_points++;
		}
	}
	Data_Frame temp_df(valid_points);
	j = 0; //initialize
	for (int i = 0; i < _points && j < valid_points; i++)
        {
		if (_data[i].get_x() != INVALID && _data[i].get_y() != INVALID)
		{
			temp_df.set_value(j, _data[i].get_x(), _data[i].get_y());
			j++;
		}
	}//end for
	*this = temp_df;
}
	
//****************** methods used to construct envelopes **************
int Data_Frame::in_interval(double min, double max, int& start_index)
{
	int hi_index; //index corresponding to max value
	double temp_value;

        //switch values if needed
        if (min > max)
        {
                temp_value = min;
                min = max;
                max = temp_value;
        }

	start_index = 0; //initialize

	//check if valid interval
	if (min == max || max < _data[0].get_x() || min > _data[_points - 1].get_x())
	{
		return -1;
	}
	else
	{
        	//check for the start_index
		while (_data[start_index].get_x() < min)
		{
			start_index++;
		}
		//to get a closed interval at the end
		if (max == _data[_points - 1].get_x())
		{
			hi_index = _points;
		}
		else	//find the hi index
		{
			hi_index = start_index; //initialize
			while (hi_index < _points && _data[hi_index].get_x() < max)
			{
				hi_index++;
			}
		}//end else
		return (hi_index - start_index);
	}//end else
}

Data_Frame Data_Frame::segment(int start_index, int number)
{
	int i = 0;
	Data_Frame df(number);
	while (i < number && i + start_index < _points)
	{
		df.set_value(i, _data[i + start_index].get_x(), _data[i + start_index].get_y());
		i++;
	}
	return df;
}

Data_Frame Data_Frame::over_top(double top)
{
	int size = 0; //size of the resulting object
	int j;
	//determine the size
	for (int i = 0; i < _points; i++)
	{
		if (_data[i].get_y() >= top)
		{
			size++;
		}
	}

	Data_Frame df(size);
	if (size != 0)
	{
		j = 0;
		//put the values
		for (int i = 0; i < _points; i++)
        	{
			if (_data[i].get_y() >= top)
			{
				df.set_value(j, _data[i].get_x(), _data[i].get_y());
				j++;
			}
		}
	}//end if
	return df;
}

Data_Frame Data_Frame::make_envelope(double w_int, double sd_scale, int min_pts)
{
	double min_x, max_x;     //min and max in the object
        double lo_lim, up_lim;   //lower and upper limits of current interval
        int new_index, new_points; //start index and # of points in interval
	double lo_bound;	   //lower _y bound for envelope in current interval
	Data_Frame envelope;     //object to hold envelope
	Data_Frame temp_df;      //temporary data frame for an interval

	//first sort the Data_Frame for the functions to work correctly
	shell_sort();

	min_x = _data[0].get_x();  //at index 0 after sorting
	max_x = _data[_points - 1].get_x();  //at last index after sorting

	//in this loop consider each interval
	//construct an envelope corresponding to the data within that
	//and add this new data to the existing envelope
	lo_lim = min_x;
	new_index = 0;
	while (lo_lim < max_x && (_points - new_index > min_pts)) //until out of range
	{
		up_lim = lo_lim + w_int; //set new upper limit

		//calculate the number of points in the interval
		//and the starting index for it
		new_points = in_interval(lo_lim, up_lim, new_index);
		if (new_points >= min_pts) //otherwise no point! (sorry!)
		{
			//extract this segment into a seperate Data_Frame
			temp_df = segment(new_index, new_points);

			//get the lower bound _y value for the envelope in this interval
			//this is calculated as mean of _y + sd_scale*st_dev of _y
			lo_bound = temp_df.mean_y() + sd_scale*temp_df.sd_y();

			//construct the envelope within this interval
			temp_df = temp_df.over_top(lo_bound);

			if (temp_df._points >= min_pts) //if more than required #
			{
				envelope = envelope + temp_df; //add it to the envelope
			}
		}
		lo_lim = up_lim; //shift the interval
	}//end while

	return envelope;
}

