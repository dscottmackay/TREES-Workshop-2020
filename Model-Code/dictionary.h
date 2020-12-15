//**********************************************************
//a class of objects designed to hold a collection of 
//variable names (string of characters), so that these
//variables can be referred by their names and the file format
//for input becomes less critical
//hope to use it as the basis for a flexible modeling system
//***************** SS 2002 *******************************
#ifndef _DICTIONARY_H
#define _DICTIONARY_H

#include <assert.h>
#include "util.h"

using namespace std;

//defines a type for storing variable name for data
const int MAX_VAR_NAME_LENGTH = 31;
typedef char Var_name[MAX_VAR_NAME_LENGTH+1];
const int INITIAL = 20;  //initially made for these many
const int INCREMENT = 5; //increment in capacity each time the limit is reached

//implemented as a dynamic array
class Var_Dictionary
{
	public:
		Var_Dictionary(int can = INITIAL);
		Var_Dictionary(const Var_Dictionary& dorig);
		~Var_Dictionary();

	//method for input
		int insert(const Var_name& vname); //insert a variable in list
	//returns the index at which the data is stored, -1 if not there in list
		int index(const Var_name& vname) const;
	//returns the name of variable by ref @ at, returns 0 if invalid location
		int var_at(Var_name& var, int at) const;
	//returns the # of variables currently available
		int n_vars() const;

        //appends contents of dictionary dorig, returns 0 if fails, 1 if OK
                int append_dictionary(const Var_Dictionary& dorig);
        //appends between indices start to end in dorig (inclusive)
                int append_dictionary(const Var_Dictionary& dorig, int start, int end);
	//copies contents of dictionary dorig, returns 0 if fails, 1 if OK
		int copy_dictionary(const Var_Dictionary& dorig);
	//copies between indices start to end, validates and copies as many available
                int copy_dictionary(const Var_Dictionary& dorig, int start, int end);

	//printing the dictionary
		friend std::ostream& print_dictionary(std::ostream& out, const Var_Dictionary& vd);

	private:
		int capacity; //available capacity
		int current;  //currently valid # of vars
		Var_name* list;

	//methods for memory management
		int space_left();  //returns 0 if no space left
		int create_space(int extra = INCREMENT); //1 if successful
		void reset();
//		void initialize();
};

//a method for data input from file
//constructs the variable names from the first line in the input stream
std::istream& read_dictionary(std::istream& in, Var_Dictionary& vd);

#endif
