//**********************************************************
//(c)2007 Sudeep Samanta - 2/9/07
//calculate Rhat for the estimands from multiple MCMC chains
//**********************************************************
#include <iostream>
#include "util.h"

int main()
{
	String setname, chainfile, rhatfile, simstr;
	String tstr;
	String* pnames;
	ifstream fin;
	ofstream fout;
	int n_chn, n_iter, n_burn, n_val, n_var;
	int i, j, k, tint;
	double** vardata;
	double** psidotj;
	double** sjsq;
	double* psiddot;
	double* Bvar;
	double* Wvar;
	double* rtRhat;
//user input
	cout << "Input MCMC set name: ";
	//a specific chain in a set of MCMC is setname#.chn, "#.chn" is auto added
	//all chains should be identical in terms of iterations and variables
	//1st col in these chains should be the iteration # starting at 1
	cin >> setname;
	cout << "Number of chains in set: ";
	cin >> n_chn;
//get numbers from 1st chain
	strcpy(chainfile, setname);
	strcpy(simstr,"0");
	increment_string(simstr);
	strcat(chainfile, simstr);
	strcat(chainfile, ".chn");
	strcpy(rhatfile, setname);
	strcat(rhatfile, ".Rhat");
	n_iter = datalines(chainfile);
	n_burn = n_iter/2; //ignore 1st half
	n_val = n_iter - n_burn;
	n_var = columns(chainfile)-1; //1st col is iteration#
//var names
        pnames = new String[n_var];
//data from one chain
        vardata = new double*[n_var];
        for(i=0; i<n_var; i++)
                vardata[i] = new double[n_val];
	psidotj = new double*[n_var];
	for(i=0;i<n_var;i++)
		psidotj[i] = new double[n_chn];
	sjsq = new double*[n_var];
	for(i=0;i<n_var;i++)
		sjsq[i] = new double[n_chn];
//data from all chains
        psiddot = new double[n_var];
        Bvar = new double[n_var];
        Wvar = new double[n_var];
        rtRhat = new double[n_var];
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
		for(k=0;k<n_val;k++)
		{
			fin >> tint; //iteration #
			for(i=0;i<n_var;i++)
				fin >> vardata[i][k];
			fin.ignore(MAX_LINE, '\n');
		}
		for(i=0;i<n_var;i++)
		{
			psidotj[i][j] = mean(vardata[i], n_val);
			sjsq[i][j] = calc_var(vardata[i], n_val);
//cout << psidotj[i][j] << " " << sjsq[i][j] << endl;
		}
		fin.close();
		cout << "Done " << chainfile << endl;
		increment_string(simstr);
		strcpy(chainfile, setname);
        	strcat(chainfile, simstr);
        	strcat(chainfile, ".chn");
	}
//calculations
	for(i=0;i<n_var;i++)
	{
		psiddot[i] = 0.0;
		Bvar[i] = 0.0;
		Wvar[i] = 0.0;
		for(j=0;j<n_chn;j++)
			psiddot[i] += psidotj[i][j];
		psiddot[i] = psiddot[i]/(double)n_chn;
		for(j=0;j<n_chn;j++)
		{
			Bvar[i] += (psidotj[i][j]-psiddot[i])*(psidotj[i][j]-psiddot[i]);
			Wvar[i] += sjsq[i][j];
		}
		Bvar[i] = Bvar[i]*(double)n_val/(double)(n_chn - 1);
		Wvar[i] = Wvar[i]/(double)(n_chn);
		rtRhat[i] = sqrt(((double)(n_val-1)+(Bvar[i]/Wvar[i]))/(double)n_val);
	}
//output
	open_outfile(rhatfile, fout);
	fout << "Parameter\t";
	for(j=1;j<=n_chn;j++)
		fout << "psi." << j << '\t' << "s" << j << "^2\t";
	fout << "psi..\tB\tW\tRhat^0.5" << endl;
	for(i=0;i<n_var;i++)
	{
		fout << pnames[i] << '\t';
		for(j=0;j<n_chn;j++)
			fout << psidotj[i][j] << '\t' << sjsq[i][j] << '\t';
		fout << psiddot[i] << '\t' << Bvar[i] << '\t' 
			<< Wvar[i] << '\t' << rtRhat[i] << endl;
	}
	fout.close();
	return 1;
}
