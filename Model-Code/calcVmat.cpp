//*************************************
//(c)2007 Sudeep Samanta - 2/7/07
//calculate covariance matrix from chains
//*************************************
#include <iostream>
#include "util.h"
int main()
{
	String filename, chainfile, covarfile;
	ifstream fin;
	ofstream fout;
	int n_iter, n_burn, n_val, n_var;
	int i, j, tint;
	double** vardata;
	double** covmat;
	double** cholf;
	cout << "Input chain filename: "; //dont give .chn ext
	cin >> filename;
	strcpy(chainfile, filename);
	strcat(chainfile, ".chn");
	strcpy(covarfile, filename);
	strcat(covarfile, ".cov");

	n_iter = datalines(chainfile);
	n_burn = n_iter/2; //ignore 1st half
	n_val = n_iter - n_burn;
	n_var = columns(chainfile);
	n_var = n_var - 2; //1st col is iteration#, last is Dtheta

	vardata = new double*[n_var];
	for(i=0; i<n_var; i++)
		vardata[i] = new double[n_val];
	covmat = new double*[n_var];
        for(i=0; i<n_var; i++)
                covmat[i] = new double[n_var];
	cholf = new double*[n_var];
        for(i=0; i<n_var; i++)
                cholf[i] = new double[n_var];
	tint = open_infile(chainfile,fin);
	if(tint == 0)
		return 0;
	for(i=0;i<=n_burn;i++)
		fin.ignore(MAX_LINE, '\n');
	for(i=0;i<n_val;i++)
	{
		for(j=0;j<=n_var;j++)
		{
			if(j==0)
				fin >> tint;
			else
				fin >> vardata[j-1][i];
		}
		fin.ignore(MAX_LINE, '\n');
	}
	calc_covmat(vardata, covmat, n_var, n_val);
	calc_cholesky(covmat, cholf, n_var);
	open_outfile(covarfile, fout);
fout << "covariance matrix:\n";
        for(i=0; i<n_var; i++)
        {
                for(j=0;j<n_var;j++)
                {
                        fout << covmat[i][j];
                        if(j<n_var-1)
                                fout << ' ';
                }
                fout << endl;
        }
fout << "\nCholesky factor:\n";
	for(i=0; i<n_var; i++)
	{
		for(j=0;j<n_var;j++)
		{
			fout << cholf[i][j];
			if(j<n_var-1)
				fout << ' ';
		}
		fout << endl;
	}
	return 1;
}
