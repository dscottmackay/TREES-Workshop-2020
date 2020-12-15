//*************************************************
//defines a class of objects to keep time 
//in the format year, julian day, hour and minute
//*************************************************
#ifndef _TREESTIME_H
#define _TREESTIME_H

#include <assert.h>
#include "constants.h"
#include "util.h"

using namespace std;

class Time
{
	public:
		Time();
		Time(const Time& t);
		~Time();

        //data access
		//saving w/o validation, y converted to days
		void save_interval(int y, int yday, int hr, int m);
		//validating the values
		void set_time(int y, int yday, int hr, int m);
		void set_time(const Time& tm);

		//to get back the differnt components, it is done by reference
		void get_time(int& y, int& yday, int& hr, int& m) const;
		int get_day() const;	//returns the current day
		int get_year() const;    //returns the current year

		//advancing with a time interval, year should be 0
		void advance(const Time& intvl);

	//overloaded operators
		Time& operator=(const Time& rhs);
                bool operator>(const Time& rhs) const; //later than rhs
                bool operator<(const Time& rhs) const; //earlier than rhs
                bool operator==(const Time& rhs) const; //same time
                bool operator!=(const Time& rhs) const;

	//for output to stream
                friend std::ostream& operator << (std::ostream& out, const Time& tm);
	//input functions - to be modified/added to if needed
		//ddd - year day (jan1 = 1), hh - 24hr format, m = 0 (top) or 0.3 (bottom)
		friend std::istream& read_time(std::istream& in, Time& tm);

		//from time stamp output - yyyy:ddd:hh:mm (not fixed width)
		friend std::istream& read_time_stamp(std::istream& in, Time& tm);

		//the format is mm.dd	hh.mm	yday.fractionday
		//no year info so year is in argument
		friend std::istream& read_time_tdat(std::istream& in, Time& tm, int year);

//added Jan 04
                //returns the difference between two times in min
                //+ve if t2 is later, -ve if t1 is later
                friend int diff_in_min(const Time& t1, const Time& t2);

	private:
		int year;
                int jday;       //counting 1 for Jan1
                int hour;       //24 hr format
                int min;
                void validate(); //takes care of spillovers
};

class Time_Array
{
	public:
                Time_Array();
		Time_Array(const Time_Array& ta);
                ~Time_Array();

		void reset();
		int allocate(int max); //returns # of steps, 0 if fails

		//returns the # time steps in it
		int n_steps() const;
		//puts in the Time t at step st (counting from 0)
		//returns st or -1
		int set_time(const Time& t, int st);
		//copies the time series, returns 0 if fails, steps otherwise
		int copy_array(const Time_Array& ta);
		//returns the Time @ step
		Time get_time(int step) const;

		friend std::istream& read_step(std::istream& in, Time_Array& ta, int step);
                friend std::ostream& print_step(std::ostream& out, const Time_Array& ta, int step);
	private:
		Time* series;
		int steps;
};

#endif
