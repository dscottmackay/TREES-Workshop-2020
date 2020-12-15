//**********************************************************************
//this structure is for state variables in the TREES model
//the data is contained in a 1D array of doubles
//data are indexed using the dictionary
//***************************** DSM 2009 *******************************
#ifndef _STATE_STORE_H
#define _STATE_STORE_H
#include "dictionary.h"

using namespace std;

const double INVALID_STATE = -999.0;
class State_Store
{
	public:
		State_Store();
		State_Store(const State_Store& sl);
		~State_Store();
	//allocates a 1D array based on the dictionary
	//assumes that the dictionary has been allocated already
		void allocate();
	//throws out the existing values and clears
		void clear_values();

//memory allocation methods - allocation depends on the dictionary
		int get_dictionary(const Var_Dictionary& vd);
	//allocation from scratch
		void allocate(const Var_Dictionary& vd);

//inserting values in a State_Store
	//set a value for a state variable; returns 0 if fails, 1 OK
	//at index ind
		int set_val_at(double val, int ind);
	//for var name vname
		int set_val_for(double val, const Var_name& vname);
	//update a value for a state variable; returns 0 if fails, 1 OK
	//at index ind
		int update_val_at(double delta_val, int ind);
	//for var name vname
		int update_val_for(double delta_val, const Var_name& vname);

//accessing values in a State_Store
	//get the value of a state variable
	//by index in dictionary
		double get_val_at(int ind) const; //by index
	//by var name
		double get_val_for(const Var_name& vname) const;

//returns the number of variables available
	int no_of_vars() const;

	//write out the current state
	friend std::ostream& print_state(std::ostream& out, const State_Store& sl);

	private:
		Var_Dictionary dictionary; //dictionary for variable reference
		double* values;		   //data values for variable
};

#endif

