//*****************************************************
//implementation of dictionary.h
//***************** SS 2002 ***************************
#include "dictionary.h"


Var_Dictionary::Var_Dictionary(int can)
:capacity(can), current(0), list(NULL)
{
	assert(capacity > 0);
	list = new Var_name[capacity];
	assert(list != NULL);
}

Var_Dictionary::Var_Dictionary(const Var_Dictionary& dorig)
{
	assert(copy_dictionary(dorig) != 0);
}

Var_Dictionary::~Var_Dictionary()
{
	reset();
}

void Var_Dictionary::reset()
{
	if (list != NULL)
	{
		delete[] list;
        	list = NULL;
	}
        capacity = 0;
        current = 0;
}

int Var_Dictionary::space_left()
{
	if (current < capacity)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int Var_Dictionary::create_space(int extra)
{
	int flag = 1; //check
	int inc_space = extra + capacity;

	Var_name* nlist; //temporary storage
	nlist = new Var_name[current];
	if (nlist == NULL)
	{
		flag = 0;
	}
	else
	{
	//copy the existing ones
		for (int i = 0; i < current; i++)
		{
			strcpy(nlist[i], list[i]);
		}
	//allocate extra memory
		delete[] list;
		list = new Var_name[inc_space];
		if (list == NULL)
		{
                	flag = 0; //failed
		}
		else
		{
			capacity = inc_space;
			//copy back
			for (int i = 0; i < current; i++)
			{
				strcpy(list[i], nlist[i]);
			}
			delete[] nlist;
			nlist = NULL;
			flag = 1;
		}
	}
	return flag;
}
	
	
int Var_Dictionary::insert(const Var_name& vname)
{
	int flag = 1; //check
	if (!space_left() )
	{
		flag = create_space(); //flag becomes 0 if space not created
	}
	if (flag)
	{
		strcpy(list[current], vname);
		current++;
	}
	return flag;
}

int Var_Dictionary::index(const Var_name& vname) const
{
	int i = 0;
	while (i < current && (strcmp(list[i], vname) != 0))
	{
		i++;
	}
	if (i == current) //not found
	{
		i = -1;
	}
	return i;
}

int Var_Dictionary::var_at(Var_name& var, int at) const
{
	int flag = 0;
	if (at < current)
	{
		strcpy(var, list[at]);
		flag = 1;
	}
	return flag;
}

int Var_Dictionary::n_vars() const
{
	return current;
}

int Var_Dictionary::append_dictionary(const Var_Dictionary& dorig, int start, int end)
{
        int temp, flag=0;
        Var_name vname;
	//setting up for start and end
	if (start > end) //in correct order
	{
		temp = start;
		start = end;
		end = temp;
	}
	//proper indices
	end = small_int(end, dorig.current-1);
	start = big_int(0, start);
	//get the variables
        for (int i = start; i <= end; i++)
        {
                flag = dorig.var_at(vname, i);
                flag = insert(vname);
        }
        return flag;
}

int Var_Dictionary::append_dictionary(const Var_Dictionary& dorig)
{
	return append_dictionary(dorig, 0, dorig.current);
}

int Var_Dictionary::copy_dictionary(const Var_Dictionary& dorig)
{
	current = 0; //puts in starting at the first variable
	return append_dictionary(dorig);
}

int Var_Dictionary::copy_dictionary(const Var_Dictionary& dorig, int start, int end)
{
	current = 0; //puts in starting at the first variable
	return append_dictionary(dorig, start, end);
}

ostream& print_dictionary(ostream& out, const Var_Dictionary& vd)
{
	for (int i = 0; i < vd.current; i++)
	{
		out << vd.list[i] << '\t';
	}
	return out;
}

istream& read_dictionary(istream& in, Var_Dictionary& vd)
{
        char next_char;
        Var_name vname;
        in >> next_char;
        while (next_char != '\n') //not end of line
        {
                if (next_char == ' ' || next_char == '\t')
		{
                        next_char = in.get();
		}
                else
                {
                        in.unget(); //put it back
                        in >> vname;
                        //put it in dictionary
                        vd.insert(vname);
                        next_char = in.get();
                }
        }
        return in;
}

