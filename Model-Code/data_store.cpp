//***********************************
//implementation of Data_Store object
//******************** SS 2002 ******

#include "data_store.h"


Data_Store::Data_Store()
{
	values = NULL;
}

Data_Store::Data_Store(const Data_Store& dl)
{
//copy the dictionary
	int flag;
	flag = dictionary.copy_dictionary(dl.dictionary);
	assert(flag != 0);
	flag = time_series.copy_array(dl.time_series);
	assert(flag != 0);
//copy values
	allocate();
	for (int i = 0; i < dictionary.n_vars(); i++)
	{
		for (int j = 0; j < time_series.n_steps(); j++)
		{
			values[i][j] = dl.values[i][j];
		}
	}
}
	
Data_Store::~Data_Store()
{
	clear_values();
}

void Data_Store::allocate()
{
//to be allocated for so many vars
	int nv = dictionary.n_vars(); 
	int tsteps = time_series.n_steps();
	assert(nv != 0 && tsteps != 0);
//release what is already there
	clear_values(); 
//now allocate
	values = new double*[nv];
	assert(values != NULL);
	for (int i = 0; i < nv; i++)
	{
		values[i] = new double[tsteps];
		assert(values[i] != NULL);
	}
}

void Data_Store::clear_values()
{
        if (values != NULL)
        {
                int n = dictionary.n_vars();
                for (int i = 0; i < n; i++)
                {
                        if (values[i] != NULL)
                        {
                                delete[] values[i];
                                values[i] = NULL;
                        }
                }
                delete[] values;
                values = NULL;
        }
}


int Data_Store::get_dictionary(const Var_Dictionary& vd)
{
	return dictionary.copy_dictionary(vd);
}

int Data_Store::set_steps(int mstep)
{
	return time_series.allocate(mstep);
}

void Data_Store::allocate(const Var_Dictionary& vd, int mstep)
{
//make the dictionary
	assert(get_dictionary(vd) != 0);
//allocate the correct steps in time series
	assert(set_steps(mstep) != 0);
//allocate the data array
	allocate();
}

int Data_Store::set_time(const Time& t, int st)
{
	return time_series.set_time(t,st);
}

int Data_Store::put_val_at(double val, int ind, int st)
{
	int flag = 0;
//make sure it is a valid index
	if (ind < dictionary.n_vars() && st < time_series.n_steps() && ind >= 0 && st >= 0) 
	{
		values[ind][st] = val;
		flag =1;
	}
	return flag;
}

int Data_Store::put_val_for(double val, const Var_name& vname, int st)
{
	int reqind = dictionary.index(vname);
	return put_val_at(val, reqind, st);
}

int Data_Store::insert_col(Var_name& vname, const Data_Store& orig)
{
        int flag = 0;
        int tocopy, in_here, in_orig;
//the lower number of steps is copied, hopefully both same
        tocopy = small_int(orig.time_series.n_steps(), time_series.n_steps());
 
//get the required indices
        in_here = dictionary.index(vname);
        in_orig = orig.dictionary.index(vname);
        if (in_here >= 0 && in_orig >= 0)
        {
                for (int i = 0; i < tocopy; i++)
		{
                        values[in_here][i] = orig.values[in_orig][i];
		}
                flag = 1;
        }
        return flag;
}

int Data_Store::check_step(int the_step) const
{
	int validity = 0;
//max the_step < no_of_steps
	if (the_step >= no_of_steps()) 
	{
		validity = -1;
	}
	else
	{
		for (int i = 0; i < no_of_vars(); i++)
		{
			if (get_val_at(i, the_step) == INVALID_DATA)
			{
				validity++;
			}
		}
	}
	return validity;
}

Time Data_Store::get_time(int step) const
{
	Time cur_time = time_series.get_time(step);
	return cur_time;
}

double Data_Store::get_val_at(int ind, int st) const
{
	double value = INVALID_DATA;
//make sure it is a valid index
	if (ind < dictionary.n_vars() && st < time_series.n_steps() && ind >= 0 && st >= 0) 
	{
		value = values[ind][st];
	}
	return value;
}

double Data_Store::get_val_for(const Var_name& vname, int st) const
{
	int reqind = dictionary.index(vname);
	return get_val_at(reqind, st);
}

int Data_Store::carve(const Data_Store& orig, int start, int end)
{
	int flag = 1;
	Var_name vn;

//copy the relevant part of the dictionary
	flag = flag*dictionary.copy_dictionary(orig.dictionary, start, end);
//copy the time series
	if (flag != 0)
	{
		flag = flag*time_series.copy_array(orig.time_series);
	}
	if (flag != 0)
	{
		allocate();
//copy varaiable values
		for (int i = 0; i < dictionary.n_vars(); i++)
		{
//get the variable name
			dictionary.var_at(vn, i); 
//copy the values
			flag = flag*insert_col(vn, orig);
		}
	}
	return flag;
}//end carve

int Data_Store::no_of_steps() const
{
	return time_series.n_steps();
}

int Data_Store::no_of_vars() const
{
	return dictionary.n_vars();
}

ostream& print_step(ostream& out, const Data_Store& dl, int step)
{
	if (step < dl.time_series.n_steps())
	{
		print_step(out, dl.time_series, step);
		for (int i = 0; i < dl.dictionary.n_vars(); i++)
		{
			out << '\t' << dl.values[i][step];
		}
	}
	return out;
}

int get_obs(String& obs_file, Data_Store& in_data)
{
        ifstream fin;
        Var_Dictionary vd;
        Time time;
        int n_steps;    //# of time steps in file
        double temp_val;        //temporary
        String temp_string;
//find n_steps
        n_steps = datalines(obs_file);
        if (n_steps <= 0)
	{
                return 0;       //no data or file problem
	}
        else
        {
                open_infile(obs_file, fin);
//get the first 2 strings out (not vars)
                for (int i = 0; i < 2; i++)
		{
                        fin >> temp_string;
		}
//get the dictionary
                read_dictionary(fin, vd);
//allocate memeory
                in_data.allocate(vd, n_steps);  
//read data
                for (int i = 0; i < n_steps; i++)
                {
//refer time_keeper.h
                        read_time(fin, time); 
                        in_data.set_time(time, i);
                        for (int j = 0; j < vd.n_vars(); j++)
                        {
                                fin >> temp_val;
                                in_data.put_val_at(temp_val,j,i);
                        }
                        fin.ignore(MAX_LINE, '\n');
                }
                fin.close();
                return n_steps;
        }//end else
}//end get_obs

/*
 * used to retrieve shoot and root module parameters for the sperry model
 */
double getData(ifstream &in)
{
    	double data;

    	while( in.peek() < 45 || in.peek() > 57 ) 
	{
        	in.ignore(256,' ');
    	}

    	in >> data;
    	return data;
}


