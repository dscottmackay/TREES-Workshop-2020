//****************************************************
//implementation of Time object - SS June, 2002
//****************************************************
#include "time_keeper.h"


Time::Time()
{
	year = 0;
        jday = 0;
        hour = 0;
        min  = 0;
}

Time::Time(const Time& t)
{
	set_time(t);
}
 
Time::~Time()
{ ; }

void Time::validate()
{
	int carry = 0;
        int leap;
 
        carry = min/MINS;
        if(min >= MINS)
	{
                min = min%MINS;
 	}
        hour += carry;
        carry = hour/HOURS;
        if(hour >= HOURS)
	{
                hour  = hour%HOURS;
 	}
        jday += carry;
        leap = leap_year(year);
        while(jday >= DAYS+leap)
        {
                jday -= DAYS+leap;
                year++;
                leap = leap_year(year);
        }
}

void Time::save_interval(int y, int yday, int hr, int m)
{
	//considers one year passed to this as defined in constants.h days
        //only to be used for advancing time
        year = 0;
        jday = y*DAYS + yday;
        hour = hr;
        min  = m;
}

void Time::set_time(int y, int yday, int hr, int m)
{
	year = y;
	jday = yday;
        hour = hr;
        min  = m;
	validate();
}

void Time::set_time(const Time& tm)
{
	set_time(tm.year, tm.jday, tm.hour, tm.min);
}

void Time::get_time(int& y, int& yday, int& hr, int& m) const
{
	//copy data in reference args
	y = year;
	yday = jday;
	hr = hour;
	m = min;
}

int Time::get_day() const
{
	return jday;
}

int Time::get_year() const
{
        return year;
}

void Time::advance(const Time& intvl)
{
	min += intvl.min;
        hour += intvl.hour;
        jday += intvl.jday;
        year += intvl.year;
        validate();
}

Time& Time::operator=(const Time& rhs)
{
	if (this != &rhs)
	{
		set_time(rhs);
	}
	return *this;
}

bool Time::operator>(const Time& rhs) const
{
	if (year != rhs.year) 
	{
		return year > rhs.year;
	}
	else
	{ //same year
		if (jday != rhs.jday) 
		{
			return jday > rhs.jday;
		}
		else //same day
		{
			if (hour != rhs.hour) 
			{
				return hour > rhs.hour;
			}
			else //same hour
			{
				return min > rhs.min;
			}
		}
	}	
}

bool Time::operator<(const Time& rhs) const
{
        if (year != rhs.year) 
	{
		return year < rhs.year;
	}
        else
        { //same year
                if (jday != rhs.jday) 
		{
			return jday < rhs.jday;
		}
                else //same day
                {
                        if(hour != rhs.hour) return hour < rhs.hour;
                        else //same hour
                                return min < rhs.min;
                }
        }
}

bool Time::operator==(const Time& rhs) const
{
	return (year == rhs.year && jday == rhs.jday && hour == rhs.hour && min == rhs.min);
}

bool Time::operator!=(const Time& rhs) const
{
	return (year != rhs.year || jday != rhs.jday || hour != rhs.hour || min != rhs.min);
}


ostream& operator << (ostream& out, const Time& tm)
{
        out << tm.year << ':'
            << tm.jday << ':'
            << tm.hour << ':'
            << tm.min;
        return out;
}

//************* input function to be modified if needed **************
//this is for Doug/Brent's file format as of June 2002
//the time is in the form yyyyddd hh.m
//ddd - year day (jan1 = 1), hh - 24hr format, m = 0 (top) or 0.3 (bottom)
istream& read_time(istream& in, Time& tm)
{
	int yd;
        double hm;
        int y, j, h, m;
        in >> yd >> hm;
 
        y = yd/1000; //throwing out the last 3 digits
        j = yd%1000; //the rest is the day
        hm = hm*100; //last 2 digits = mins
        h = (int)hm/100;
        m = (int)hm%100;
 
        tm.set_time(y, j, h, m);
        return in;
}

//from time stamp output - yyyy:ddd:hh:mm (not fixed width)
istream& read_time_stamp(istream& in, Time& tm)
{
	int y, j, h, m;
	char tempch;
	in >> y >> tempch >> j >> tempch >> h >> tempch >> m;
	tm.set_time(y, j, h, m);
        return in;
}

 
//from Brent's measured temp data
//the format is mm.dd   hh.mm   yday.fractionday
//no year info so year is in argument
istream& read_time_tdat(istream& in, Time& tm, int year)
{
        int j, h, m;
	double temp;

	in >> temp >> temp; //ignore mm.dd, keep hh.mm
	temp = temp*100; //last 2 digits = mins
	h = (int)temp/100;
	if (h == 24) 
	{
		h = 0; //as the first hr of day is 24 in data
	}
	m = (int)temp%100;
        in >> temp;
	j = (int)temp;
	tm.set_time(year, j, h, m);
        return in;
}

int diff_in_min(const Time& t1, const Time& t2)
{
        Time tg, tl;
        int day_diff, hour_diff, min_diff, ymid;
        if (t1 == t2)
	{
                min_diff = 0;
	}
        else
        {
                if (t1 > t2)
                {
                        tg = t1;
                        tl = t2;
                }
                else
                {
                        tg = t2;
                        tl = t1;
                }
                ymid = tl.year;
		day_diff = 0;
                while (ymid<tg.year)
                {
                        day_diff += (DAYS + leap_year(ymid)); //add days for intervening years
                        ymid++;
                }
		day_diff += (tg.jday - tl.jday);
                hour_diff = day_diff*HOURS + tg.hour - tl.hour;
                min_diff = hour_diff*MINS+tg.min-tl.min;  //whole hours
        }
        return min_diff;
}

//********************* class Time_Array *****************
Time_Array::Time_Array()
{
	steps = 0;
	series = NULL;
}

Time_Array::Time_Array(const Time_Array& ta)
{
	assert(copy_array(ta) > 0);
}

Time_Array::~Time_Array()
{
	reset();
}

void Time_Array::reset()
{
	if (series != NULL)
	{
        	delete [] series;
        	series = NULL;
	}
        steps = 0;
}

int Time_Array::allocate(int max)
{
	int flag = 0;
	if (max > 0)
	{
		series = new Time[max];
		if (series != NULL)
		{
			flag = steps = max;
		}
	}
	return flag;
}

int Time_Array::n_steps() const
{
	return steps;
}

int Time_Array::set_time(const Time& t, int st)
{
	int flag = -1;
	if (st >= 0 && st < steps)
	{
		series[st].set_time(t);
		flag = st;
	}
	return flag;
}
	
int Time_Array::copy_array(const Time_Array& ta)
{
	int i, flag = 0;
	int ns = ta.n_steps();
	flag = allocate(ns); //see if ok
	if (steps != ta.n_steps())
	{
		flag = 0;
	}
	else
	{
		for (i=0; i<steps; i++)
		{
			series[i] = ta.series[i];
		}
	}
	return flag;
}

Time Time_Array::get_time(int step) const
{
	Time new_time;
	if (step < steps)
	{
		new_time.set_time(series[step]);
	}
	else
	{
		new_time.set_time(0,0,0,0);
	}
	return new_time;
}

istream& read_step(istream& in, Time_Array& ta, int step)
{
	if (step < ta.steps)
	{
		read_time(in, ta.series[step]);
	}
	return in;
}

ostream& print_step(ostream& out, const Time_Array& ta, int step)
{
	if (step < ta.steps)
	{
		out << ta.series[step];
	}
	return out;
}


