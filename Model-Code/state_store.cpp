//***********************************
//implementation of State_Store object
//******************** DSM 2009 *****

#include "state_store.h"


State_Store::State_Store()
{
	values = NULL;
}

State_Store::State_Store(const State_Store& sl)
{
//copy the dictionary
	int flag;
        flag = dictionary.copy_dictionary(sl.dictionary);
        assert(flag != 0);
//copy values
        allocate();
        for (int i = 0; i < dictionary.n_vars(); i++)
	{
        	values[i] = sl.values[i];
	}
}

State_Store::~State_Store()
{
	clear_values();
}

void State_Store::allocate()
{
        int nv = dictionary.n_vars(); //to be allocated for so many vars
        assert(nv != 0);
        clear_values(); //release what is already there
//now allocate
        values = new double[nv];
        assert(values != NULL);
}

void State_Store::clear_values()
{
        if(values != NULL)
        {
		delete[] values;
		values = NULL;
	}
}


int State_Store::get_dictionary(const Var_Dictionary& vd)
{
        return dictionary.copy_dictionary(vd);
}

void State_Store::allocate(const Var_Dictionary& vd)
{
//make the dictionary
        assert(get_dictionary(vd) != 0);
//allocate the data array
        allocate();
}

int State_Store::set_val_at(double val, int ind)
{
        int flag = 0;
        if (ind < dictionary.n_vars() && ind >= 0) //make sure it is a valid index
        {
                values[ind] = val;
                flag = 1;
        }
        return flag;
}

int State_Store::set_val_for(double val, const Var_name& vname)
{
        int reqind = dictionary.index(vname);
        return set_val_at(val, reqind);
}

int State_Store::update_val_at(double delta_val, int ind)
{
        int flag = 0;
//make sure it is a valid index
        if(ind < dictionary.n_vars() && ind >= 0) 
        {
                values[ind] += delta_val;
                flag = 1;
        }
        return flag;
}

int State_Store::update_val_for(double delta_val, const Var_name& vname)
{
        int reqind = dictionary.index(vname);
        return update_val_at(delta_val, reqind);
}

double State_Store::get_val_at(int ind) const
{
        double value = INVALID_STATE;
//make sure it is a valid index
        if(ind < dictionary.n_vars() && ind >= 0) 
	{
                value = values[ind];
	}
        return value;
}

double State_Store::get_val_for(const Var_name& vname) const
{
        int reqind = dictionary.index(vname);
        return get_val_at(reqind);
}

int State_Store::no_of_vars() const
{
        return dictionary.n_vars();
}

ostream& print_state(ostream& out, const State_Store& sl)
{
        for(int i = 0; i < sl.dictionary.n_vars(); i++)
	{
        	out << '\t' << sl.values[i];
	}
        return out;
}

