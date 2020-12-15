//**********************************************************************
//this is the basic data structure for TREES model
//the data is contained in a 2D array of doubles
//the first index corresponds to the index of variable in the dictionary
//the second index corresponds to the time step
//this data is linked to the dictionary
//******************** SS 2002 ***************************
#ifndef _DATA_STORE_H
#define _DATA_STORE_H
#include "time_keeper.h"
#include "dictionary.h"

using namespace std;

const double INVALID_DATA = -998;
class Data_Store
{
	public:
		Data_Store();
		Data_Store(const Data_Store& dl);
		~Data_Store();
         //allocates the 2D array vlaues based on the dictionary and time series
         //assumes that the series and dictionary have been allocated already
		void allocate();
	//throws out the existing values and clears
		void clear_values();

//memory allocation methods - allocation depends on the dictionary
		int get_dictionary(const Var_Dictionary& vd);
	//set number of time steps
		int set_steps(int mstep);
	//allocation from scratch
		void allocate(const Var_Dictionary& vd, int mstep);

//inserting values in a Data_Store
	//input time @ step st
		int set_time(const Time& t, int st);
	//input a value for a variable @ step st; returns 0 if fails, 1 OK
	//at index ind
		int put_val_at(double val, int ind, int st);
	//for var name vname
		int put_val_for(double val, const Var_name& vname, int st);
        //copy the whole col of variable vname from orig 
                int insert_col(Var_name& vname, const Data_Store& orig);

//data verification & access
	//checks if any variable @ the_step has invalid data
	//the_step is index in time_series
	//returns # of invalid variables @ the_step, -1 if invalid step#
		int check_step(int the_step) const;
	//returns the time at step
		Time get_time(int step) const;
	//get the value of a variable at time step st
	//by index in dictionary
		double get_val_at(int ind, int st) const; //by index
	//by var name
		double get_val_for(const Var_name& vname, int st) const;

		//puts in variables from start to end from the const argument orig
		//returns 0 if fails, 1 OK
		int carve(const Data_Store& orig, int start, int end);

//returns the number of steps available
                int no_of_steps() const;
//returns the number of variables available
                int no_of_vars() const;
 
		//write out the step
		friend std::ostream& print_step(std::ostream& out, const Data_Store& dl, int step);

	private:
		Var_Dictionary dictionary; //dictionary for variable reference
		Time_Array time_series;    //the time series
		double** values;           //data values for variable @ each time step
};

//*****************************************************************************
//a function to get data from a file into a Data_Store object declared in the 
//calling function. Allocation is done here to suit the data in file. 
//returns the number data lines in file obs_file.
//imports data from file into a data store object. Modify for time format.
//This one handles jday (yyyyddd) time (hr.min), then any number of variables
//neglects first 2 strings on first line as they are not variable names.
//*****************************************************************************
int get_obs(String& obs_file, Data_Store& in_data);

//used to retrieve shoot and root module parameters for the sperry model
double getData(std::ifstream&);

#endif
