//************************************************************************************************
//for creating central posterior intervals from posterior param distributions generated with mcEcl
//added a number nperp for each param in postp - nov 04
//posterior parameter file has one line header and all parameter values for one iteration on one line.
//(c)2007 - Sudeep Samanta
//****************************************************************************************************
#include <iostream>
#include "simulator2.h"
#include "time_keeper.h"
#include "data.h"
#include "data_store.h"
#include "parameter.h"
#include "randomize.h"
#include "state_store.h"
 
int main()
{
//declarations
	//for hard coded checks
	int reqp = 21; //for MR
	int reqobs = 12;

	String pfile, param_file_modules, dfile, pdfile, outfile, outfile2; //file names
	ifstream in;
	ofstream out, out2;

	int np, mp, checkp; //total # and # of mc params in .p file
	int npostp, ipostp, npostsim; //# of params, iterations in .pd, number of simulated values per step
	int nperp;	//# of sims for each parameter
	double inpval;	//param value from file
	double detEc, errEc, randEc, sd1; //deterministic, error term, together and sd param for Ec
	double detNEE, errNEE, randNEE, sd2; //deterministic, error term, together and sd param for NEE
	int nsteps; //# of steps in data file
	int i, j, k, m;
	struct trees_params treesParams;

	Data_Store indata;
	State_Store state;
	Parameter* currp;
	double** postEc;	//big sim set, nperp*ipostp
	double** postNEE;	//big sim set, nperp*ipostp

	Xsubi_array xsubi, xsubi2; 
	initialize(xsubi);
	initialize(xsubi2);

//For Sperry
        int smodules;
        int rmodules;
        double al[10];
        double dslat[10];
        double dsax[10];
        double ar[10];
        double drlat[10];
        double drax[10];


//get input
	cout << "Observation File Name: ";
	cin >> dfile;
	nsteps = get_obs(dfile, indata);
        if(nsteps == 0)
        {
                cout << "\nINVALID OR EMPTY FILE " << dfile << endl;
                return 0;
        }
        if(indata.no_of_vars() < reqobs)
        {
                cout << "\nINSUFFICIENT # OF OBSERVED VARIABLES\n" << reqobs
                        << " NEEDED " << indata.no_of_vars() << " AVAILABLE" << endl;
                return 0;
        }
//some decl. here to allocate correctly
	Data_Frame Ecframe(nsteps), NEEframe(nsteps); //to hold obs Ec and NEE in x, simulated in y
        Var_name ecname, neename;
        //strcpy(ecname, "simET");
        strcpy(ecname, "Ec");
        strcpy(neename, "NEEobs");
        for(i=0; i<nsteps; i++)
	{
                Ecframe.set_x_value(i, indata.get_val_for(ecname, i)); //OK
                NEEframe.set_x_value(i, indata.get_val_for(neename, i)); //OK
	}
        cout << dfile << " checked\n";
        cout << "\nName of parameter file: ";
        cin >> pfile;
        np = lines(pfile);
        currp = new Parameter[np];
        checkp = read_params(pfile, np, currp);
        if(checkp < reqp)
        {
                cout << "\nINSUFFICIENT # OF PARAMETERS\n"
                        << reqp << " NEEDED "
                        << checkp << " IN " << pfile << endl;
                return 0;
        }
        cout << endl << pfile << " OK\n";
	mp = 0;	//these would have to be obtained from .pd file
        for(i=0;i<np;i++)
                if(currp[i].get_dtype() == RANDOM)
                        mp++;

//Plant Shoot and Root Module Parameter Read In
        cout << "Shoot and Root Modules File Name: ";
        cin >> param_file_modules;

        ifstream pmodules;
        pmodules.open(param_file_modules);
        if(pmodules.is_open() == true)
        {
                smodules = getData(pmodules);
                treesParams.smodules = smodules;

                for(i=1; i<=smodules; i++)
                {
                    treesParams.al[i] = getData(pmodules);
                    treesParams.dslat[i] = getData(pmodules);
                    treesParams.dsax[i] = getData(pmodules);
                }

                rmodules = getData(pmodules);
                treesParams.rmodules = rmodules;

                for(i=1; i<=rmodules; i++)
                {
                    treesParams.ar[i + smodules + 1] = getData(pmodules);
                    treesParams.drlat[i + smodules + 1] = getData(pmodules);
                    treesParams.drax[i + smodules + 1] = getData(pmodules);
                }
        }


	cout << "\nFile name for posterior distribution of parameters: ";
	cin >> pdfile;
	ipostp = datalines(pdfile);
	npostp = columns(pdfile);
	if(mp > npostp)
	{
		cout << "\nINSUFFICIENT # OF PARAMETERS in " << pdfile << endl;
		cout << mp << " required " << npostp << " available\n";
		return 0;
	}
	cout << endl << pdfile << " OK\n";
	cout << "How many simulations for each parameter set in posterior file? ";
	cin >> nperp;

//set up the matrix with simulated values of transpiration
	postEc = new double*[nsteps];
	postNEE = new double*[nsteps];
	npostsim = ipostp*nperp;
	for(i=0; i<nsteps; i++)
	{
		postEc[i] = new double[npostsim]; //nov04
		postNEE[i] = new double[npostsim]; //aug09
	}
	cout << "\nOutput file name for ET: ";
	cin >> outfile;
	cout << outfile << endl;
	cout << "\nOutput file name for NEE: ";
	cin >> outfile2;
	cout << outfile2 << endl;
 
//prepare for simulations
	open_infile(pdfile, in);
	if(!isdigit(in.peek()))
		in.ignore(MAX_LINE, '\n'); //ignore first line if not data
	if(open_outfile(outfile, out)==0)
	{
		cout << "ERROR opening output file for ET " << outfile << endl;
		return 0;
	}
	if(open_outfile(outfile2, out2)==0)
	{
		cout << "ERROR opening output file for NEE " << outfile2 << endl;
		return 0;
	}
	out << "Time\tobserved\tmean\tmedian\tresidual\tmin\tp0.5%\tp2.5%\tp5%\tp12.5%\tp25%\tp75%\tp87.5%\tp95%\tp97.5%\tp99.5%\tmax" << endl;
	out2 << "Time\tobserved\tmean\tmedian\tresidual\tmin\tp0.5%\tp2.5%\tp5%\tp12.5%\tp25%\tp75%\tp87.5%\tp95%\tp97.5%\tp99.5%\tmax" << endl;
//run sims
	for(i=0; i<ipostp; i++)
	{
		//get parameters
		for(j=0; j<np; j++)
		{
			if(currp[j].get_dtype() == RANDOM)
			{
				in >> inpval;
				currp[j].set_value(inpval);
			}
		}
		in.ignore(MAX_LINE, '\n'); //ignore rest

		//runsim and add error
		simulator(indata, state, currp, treesParams, Ecframe, NEEframe, sd1, sd2);
		for(k=0; k<nsteps; k++)
		{
			detEc = Ecframe.get_y(k);
			for(m=0; m<nperp; m++) //nov04
			{
				randEc = -1.0; //to get in loop at least once
				while (randEc < 0)
					randEc = normal_random(detEc, sd1, xsubi);
				postEc[k][i*nperp+m] = randEc;	//get sim results in matrix
			}
			detNEE = NEEframe.get_y(k);
			for(m=0; m<nperp; m++) //nov04
			{
				randNEE = normal_random(detNEE, sd2, xsubi2);
				postNEE[k][i*nperp+m] = randNEE;       //get sim results in matrix
			}
		}
	}
//output results
	for(k=0; k<nsteps; k++)
	{
		shell_sort(postEc[k], npostsim);
		out << indata.get_time(k) << '\t'
		    << Ecframe.get_x(k) << '\t'
                    << mean(postEc[k], npostsim) << '\t'
                    << median(postEc[k], npostsim) << '\t'
                    << Ecframe.get_x(k) - mean(postEc[k], npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.0, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.005, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.025, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.05, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.125, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.25, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.75, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.875, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.95, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.975, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 0.995, npostsim) << '\t'
                    << sorted_qntl(postEc[k], 1.0, npostsim) << '\n';

		shell_sort(postNEE[k], npostsim);
		out2 << indata.get_time(k) << '\t'
		    << NEEframe.get_x(k) << '\t'
                    << mean(postNEE[k], npostsim) << '\t'
                    << median(postNEE[k], npostsim) << '\t'
                    << NEEframe.get_x(k) - mean(postNEE[k], npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.0, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.005, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.025, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.05, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.125, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.25, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.75, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.875, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.95, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.975, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 0.995, npostsim) << '\t'
                    << sorted_qntl(postNEE[k], 1.0, npostsim) << '\n';
	}
//delete postEc and postNEE
	in.close();
	out.close();
	for(i=nsteps; i>0; i--)
	{
		delete [] postEc[i-1];
		postEc[i-1] = NULL;
		delete [] postNEE[i-1];
		postNEE[i-1] = NULL;
	}
	delete [] postEc;
	postEc = NULL;
	delete [] postNEE;
	postNEE = NULL;

	return 1;
}
