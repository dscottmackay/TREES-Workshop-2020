//***************************************************************
//implementation for util.h - SS
//(c)2007 Sudeep Samanta - 1/4/07
//****************************************************************
#include "util.h"
#include <iostream>


int big_int(int oneint, int othrint)
{
	if(oneint > othrint)
	{
		return (oneint);
	}
	else
	{
		return (othrint);
	}
}

int small_int(int oneint, int othrint)
{
        if (oneint < othrint)
	{
        	return (oneint);
	}
        else
	{
                return (othrint);
	}
}

double absdbl(double x)
{
        if (x < 0.0)
	{
                x *= -1.0;
	}
        return (x);
}

//returns 0 if fails to open
int open_infile(const String& filename, ifstream& in)
{
        in.open(filename, ios::in);
        if(in.fail())
        {
            	std::cerr << "ERROR!! Unable to open file " << filename << endl;
                return (0);
        }
        else
	{
                return (1);
	}
}//end open_infile
 
//opens file for output
int open_outfile(const String& filename, ofstream& out)
{
        out.open(filename, ios::out);
        if(out.fail())
        {
            	std::cerr << "ERROR!! Unable to open file " << filename << endl;
                return (0);
        }
        else
	{
                return (1);
	}
}//end open_outfile

void increment_string(String& thestr)
{
        enum progress_flag { PROCEED, STOP, ERROR} decision;
        decision = PROCEED;

        //a string used to append a 1
        String unit_string;
        strcpy(unit_string, "1");

        int length = strlen(thestr);

        int i = length - 1;

        //not an empty string and do it at least once
        while (( i >= 0 ) && (decision == PROCEED))
        {
                //if the last digit is 0 to 8
                //then increase it by 1 and we are done
                if((thestr[i] >= '0') && (thestr[i] < '9'))
                {
                        thestr[i] = thestr[i]+1;
                        decision = STOP;
                }
                //if the digit is 9 then set it to 0
                //and check if the digit to the left is < 9
                else if(thestr[i] == '9')
                {
                        thestr[i] = '0';
                        decision = PROCEED;
                        i--;
                }
                else
		{
                        decision = ERROR;
		}
        }//end while

        if(decision == PROCEED)
        //all previous digits were '9' and now they are set to 0
        //so append this string of 0 s to a 1
        {
                strcat(unit_string, thestr);
                strcpy(thestr,unit_string);
                decision = STOP;
        }

}//end increment_string

int lines(const String &filename)
{
        int i = 0; //counter
        ifstream in;
 
        if(open_infile(filename, in) == 0)
	{
                return (-1);
	}
        else
        {
                while(in.peek() != EOF)
                {
                        in.ignore(512, '\n'); //flush the line
                        i++;
                }
                in.close();
                return (i); //first line is header
        }
}

int datalines(const String &filename)
{
        int i = 0; //counter
        ifstream in;
 
        if(open_infile(filename, in) == 0)
	{
                return (-1);
	}
        else
        {
                if(!isdigit(in.peek()))
		{
                        in.ignore(MAX_LINE, '\n'); //ignore first line if not data
		}
                while(in.peek() != EOF && isdigit(in.peek()))
                {
                        in.ignore(512, '\n'); //flush the line
                        i++;
                }
                in.close();
                return (i); //first line is header
        }
}

//returns the number of columns in a file
int columns(const String &filename)
{
        int cols = 0;
        char next_char;
        String tstr;
        ifstream in;
        open_infile(filename, in);
//added to address problem with interpreting columns with parameter descriptions as header
        if (!isdigit(in.peek()))
	{
        	in.ignore(MAX_LINE, '\n'); //ignore first line if not data
 	}
        in >> next_char;
        while(next_char != '\n') //not end of line
        {
                if (next_char == ' ' || next_char == '\t')
		{
                        next_char = in.get();
		}
                else
                {
                        in.unget(); //put it back
                        in >> tstr;
                        //put it in dictionary
                        cols++;
                        next_char = in.get();
                }
        }
        in.close();
        return (cols);
}

int leap_year(int year)
{
        if (year%4 != 0) 
	{
		return (0);
	}
        else if (year%100 == 0 && year%400 != 0) 
	{
		return (0);
	}
        else 
	{
		return (1);
	}
}

double deg2rad(double deg)
{
	double rad;

	rad = M_PI*deg/180.0;

	return (rad);
}

double rad2deg(double rad)
{
	double deg;

	deg = 180.0*rad/M_PI;

	return (deg);
}

double C2K(double deg_C)
{
	double deg_K;

	deg_K = deg_C + 273.15;

	return (deg_K);
}

double K2C(double Kelvin)
{
	double deg_C;

	deg_C = Kelvin - 273.15;

	return (deg_C);
}

double long_corr(double lat, double longi)
{
	int cen_meridian;
	double correction;
	lat = absdbl(lat);
	longi = absdbl(longi);
//calculate the central meridian
	cen_meridian = 15.0*(int)((longi + 7.5)/15.0); //to round off
	correction = ((double)cen_meridian - longi)/15.0;

	return (correction);
}

//ascending
void shell_sort(double* data_array, int n_items)
{
        int gap = n_items/3 + 1;
        int sublist;
 
        while(gap > 1)
        {
                for(sublist = 0; sublist < gap; sublist++)
                {
                        gap_sort(data_array, sublist, gap, n_items);
                        gap = gap/3 + 1;
                }
        }
        gap_sort(data_array, 0, 1, n_items);
}
 
void gap_sort(double* data_array, int s, int g, int n)
{
        int unsorted, loc;
        double nextitem;

        for (unsorted = s+g; unsorted < n; unsorted += g)
        {
                nextitem = data_array[unsorted];
                loc = unsorted;
 
                for(; loc > s && data_array[loc - g] > nextitem; loc -=g)
		{
                        data_array[loc] = data_array[loc-g];
		}
                data_array[loc] = nextitem;
        }
}
 
double median(double* data_array, int n_items)
{
        double the_median;
	double* temp_array;

	temp_array = new double[n_items];
	for(int i = 0; i < n_items; i++)
	{
		temp_array[i] = data_array[i];
	}
	shell_sort(temp_array, n_items); //to call w/o sorting externally
        if(n_items > 1)
        {
                if(n_items%2 != 0)
		{
                        the_median = temp_array[n_items/2];
		}
                else
		{
                        the_median = (temp_array[n_items/2]
                                        + temp_array[(n_items/2)-1])/2.0;
		}
        }
        else
	{
                the_median = temp_array[0];
	}
	delete [] temp_array;
	temp_array = NULL; 

        return (the_median);
}

double mean(double* data_array, int n_items)
{
	double count, mean, sum;

	sum = 0.0;
	count = 0.0;
	for (int i = 0; i < n_items; i++)
	{
		count += 1.0;
		sum += data_array[i];
	}
	mean = sum/count;
	//return sum/(double)n_items;
	return (mean);
}

double quantile(double* data_array, double quant, int n_items)
{
        double the_quantile;
        double* temp_array;
        int i;
        temp_array = new double[n_items];
        for(i = 0; i < n_items; i++)
	{
                temp_array[i] = data_array[i];
 	}
        shell_sort(temp_array, n_items); //to call w/o sorting externally

//get the index at the quantile
	if (quant >= 1.0)
	{
		i = n_items - 1;
	}
	else
	{
		i = (int)(absdbl(quant)*n_items);
	}
//cout << " " << quant << " " << i << endl;	
        the_quantile = temp_array[i];
 
        delete [] temp_array;
        temp_array = NULL;

        return (the_quantile);
}

//if the array is already sorted - do not repeatedly sort if it is to be thrown anyway
double sorted_qntl(double* data_array, double quant, int n_items)
{
        double the_quantile;
        int i;
 
        //get the index at the quantile
        if (quant >= 1.0)
	{
                i = n_items - 1; //max
	}
	else if (quant == 0.0)
	{
		i = 0;	//min
	}
        else
	{
                i = (int)(absdbl(quant)*n_items);
	}
        the_quantile = data_array[i];

        return (the_quantile);
}


double calc_var(double* data_array, int n_items)
{
	double variance, sumsq, samp_mean;
	sumsq = 0.0;
	samp_mean = mean(data_array, n_items);
	for(int i = 0;i < n_items;i++)
	{
		sumsq += (data_array[i]-samp_mean)*(data_array[i]-samp_mean);
	}
	variance = sumsq/(double)(n_items - 1);

	return (variance);
}

double calc_covar(double* datax, double* datay, int n_items)
{
	double covariance, sumsq, meanx, meany;
	sumsq = 0.0;
	meanx = mean(datax, n_items);
	meany = mean(datay, n_items);
	for(int i=0;i < n_items;i++)
	{
		sumsq += (datax[i]-meanx)*(datay[i]-meany);
	}
	covariance = sumsq/(double)(n_items - 1);

	return (covariance);
}

void calc_covmat(double** vardata, double** covmat, int nvar, int nval)
{
	for (int i = 0; i < nvar; i++)
	{
		for (int j = 0; j < nvar; j++)
		{
			if (i == j)
			{
				covmat[i][j] = calc_var(vardata[i], nval);
			}
			else
			{
				covmat[i][j] = calc_covar(vardata[i], vardata[j], nval);
			}
		}
	}
}

void calc_cholesky(double** covmat, double** choleskyf, int ndim)
{
	double tdbl;

	for (int i = 0; i < ndim; i++)
	{
		for (int j = 0; j < ndim; j++)
		{
			if (i == j)
			{
				tdbl = 0.0;
				if (i > 0)
				{
					for (int k = 0; k <= (i-1); k++)
					{
						tdbl += pow(choleskyf[i][k],2.0);
					}
				}
				choleskyf[i][j]=sqrt(covmat[i][j]-tdbl);
			}
			else if (i > j)
			{
				tdbl = 0.0;
				if (j > 0)
				{
					for (int k = 0; k <= (j-1); k++)
					{
						tdbl += choleskyf[i][k]*choleskyf[j][k];
					}
				}
				choleskyf[i][j] = (covmat[i][j]-tdbl)/choleskyf[j][j];
			}
			else
			{
				choleskyf[i][j] = 0.0;
			}
		}
	}
}

void scale_sqmatrix(double** sqmat, double scalef, int ndim)
{
	for (int i = 0; i < ndim; i++)
	{
		for (int j=0; j < ndim; j++)
		{
			sqmat[i][j] = scalef*sqmat[i][j];
		}
	}
}

void copy_vect(double* origv, double* copyv, int ndim)
{
	for (int i = 0; i < ndim; i++)
	{
		copyv[i] = origv[i];
	}
}

//TO BE REMOVED
void calc_stdev(double** par_history, int no_par, int history)
{
	int i, j;
	double Mean, sum, sumsq;
	for(i=0; i<no_par; i++)
	{
		sum = 0.0;
		for(j=0; j<history; j++)
			sum += par_history[i][j];
		Mean = sum/(double)history;
		sumsq = 0.0;
		for(j=0; j<history; j++)
			sumsq += pow((par_history[i][j] - Mean), 2.0);
		par_history[i][history] = sqrt(sumsq/(double)(history-1));
	}
}

