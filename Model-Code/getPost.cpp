//**********************************************************
//(c)2007 Sudeep Samanta - 3/14/07
//write out posteriors (samples and averages) from a set of MCMC chains
//**********************************************************
#include <iostream>
#include "util.h"



int main()
{
	String setname, chainfile, postfile, optfile, simstr;
	String tstr;
	String* pnames;
	ifstream fin;
	ofstream fout;
	int n_chn, n_iter, n_burn, n_rest, n_sampl, n_var;
	int i, j, k, tint, subgap;
	double** postarray;
	double** psidotj;
	double* psiddot;
//user input
	cout << "Input MCMC set name: ";
	//a specific chain in a set of MCMC is setname#.chn, "#.chn" is auto added
	//all chains should be identical in terms of iterations and variables
	//1st col in these chains should be the iteration # starting at 1
	cin >> setname;
	cout << "Number of chains in set: ";
	cin >> n_chn;
	cout << "Subsampling interval: ";
	cin >> subgap;
//get numbers from 1st chain
	strcpy(chainfile, setname);
	strcpy(simstr,"0");
	increment_string(simstr);
	strcat(chainfile, simstr);
	strcat(chainfile, ".chn");
	strcpy(postfile, setname);
	strcat(postfile, ".postd");
	strcpy(optfile, setname);
	strcat(optfile, ".postp");

	n_iter = datalines(chainfile);
cout << "Number of lines = " << n_iter << endl;
	n_burn = n_iter/2; //ignore 1st half
	n_rest = n_iter - n_burn;
	n_sampl = n_chn*(n_rest/subgap);
	n_var = columns(chainfile)-1; //1st col is iteration#
cout << "Number of variables = " << n_var << endl;
//var names
        pnames = new String[n_var];
//allocate posterior array
        postarray = new double*[n_var];
        for(i=0; i<n_var; i++)
                postarray[i] = new double[n_sampl];
//get variable names
        tint = open_infile(chainfile,fin);
        if(tint == 0)
                return 0;
        fin >> tstr;
        for(i=0;i<n_var;i++)
                fin >> pnames[i];
        fin.close();
//data from files
	for(j=0;j<n_chn;j++)
	{
		tint = open_infile(chainfile,fin);
        	if(tint == 0)
                	return 0;
		fin.ignore(MAX_LINE, '\n'); //header
		for(k=0;k<n_burn;k++)
			fin.ignore(MAX_LINE, '\n');//burn in
		for(k=0;k<n_rest;k++)
		{
			if(k%subgap ==0)
			{
			
				fin >> tint; //iteration #
				for(i=0;i<n_var;i++)
					fin >> postarray[i][(j*n_rest+k)/subgap];
			}
			fin.ignore(MAX_LINE, '\n');
		}
		fin.close();
		cout << "Done " << chainfile << endl;
		increment_string(simstr);
		strcpy(chainfile, setname);
        	strcat(chainfile, simstr);
        	strcat(chainfile, ".chn");
	}
//write out posteriors
        open_outfile(postfile, fout);
//header
        for(i=0;i<n_var;i++)
        {
                fout << pnames[i];
                if(i < n_var - 1)
                        fout << '\t';
                else
                        fout << endl;
        }
//samples
	for(j=0;j<n_sampl;j++)
	{
		for(i=0;i<n_var;i++)
		{
			fout << postarray[i][j];
			if(i < n_var - 1)
				fout << '\t';
			else
				fout << endl;
		}
	}
	fout.close();
	cout << "finished writing " << postfile << endl;
//summary stats
	open_outfile(optfile, fout);
	fout << "Mean values" << endl;
	for(i=0;i<n_var;i++)
		fout << mean(postarray[i], n_sampl) << " " << pnames[i] << endl;
	fout << "\nMedian values" << endl;
	for(i=0;i<n_var;i++)
		fout << median(postarray[i], n_sampl) << " " << pnames[i] << endl;
	fout << "\n95% posterior interval" << endl;
	for(i=0;i<n_var;i++)
	{
		fout << "[" << quantile(postarray[i],0.025,n_sampl) << ","
			<< quantile(postarray[i],0.975,n_sampl) << "]" << pnames[i] << endl;
	}
//output
	fout.close();
	cout << "finished writing " << optfile << endl;
	return 1;
}
