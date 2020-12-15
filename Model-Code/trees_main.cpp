//
//trees_main.cpp
//
//main(): First function called when TREES in run
//
//Allows TREES to run in either Bayesian statistical mode when one or more parameters is preceded by an 'm',
//    or otherwise as a single, deterministic simulation
//
//Bayesian mode implements simulations described in Samanta et al 2007, 2008 Water Resources Research
//    and Mackay et al 2012 Journal of Hydrology
//
//Plant hydraulics implementation described in Mackay et al 2015 Water Resources Research
//
//Canopy photosynthesis implementation described in Loranty et al 2010 WRR and JGR-BGC
//
//
//
//added dependence on Ca - March, 2007
//Multivariate normal proposal distribution, trace option removed - 2/12/07
//LINEAR D, lfscl, Q, T, and swc - 10/10/06
//parameter values turn constraints off & on - 2/12/07
//adding additional data col swp in kPa (negative)
//made more compact, uses canopy D for gs and above canopy D for pm and gHa - 271004
//(c)2007 - Sudeep Samanta
//Added carbon cycling (c)2007, 2008 - Scott Mackay
//Added plant water balance model (c) 2010, 2011 David Roberts & Scott Mackay
//Added canopy phenology (c) 2013, 2014 Phil Savoy & Scott Mackay
//
//

#include <iostream>
#include "trees_main.h"

int main()
{
	enum SizeFlag {SMALL, MEDIUM, LARGE};
	int i, j;
	String pfile, param_file_modules, dfile, cfile, sfile, hfile, lfile, efile, vfile; //filenames
	ifstream fin;
	ofstream chnout, effout, fluxout, hydrout, lfareaout;
	Data_Store indata;
	State_Store state;
	int cursim, total_sims;
	double p_jump, cur_dens, can_dens, curDtheta, canDtheta;
	double cur_dens_Ec, cur_dens_NEE;
	//double cur_dens_Ec, cur_dens_NEE, cur_dens_Theta, can_dens_Ec, can_dens_NEE, can_dens_Theta;
	double weight;
	int np, rp, n_agg, nsteps, checkp, use_alternating_chains;
	int nsamples, isamplsml, isamplmed, isampllrg, updtg, ipastupdt, lastupdt, n_accept, acc_pastupdt;
	double mineff, maxeff, cureff;
	String pname, tempstr;
	SizeFlag samplsize;
	Parameter* currp;
	Parameter* candp;
	double* curpval;
	double* candpval;
	double* minpval;
	double* maxpval;
	double* tmpval;
	double** psamples_sml;
	double** psamples_med;
	double** psamples_lrg;
	double** covmat;
	double** cholskf;
	Xsubi_array xsubi; //for random# gen
	initialize(xsubi);

//For hydraulic model
	struct trees_params treesParams;
        int smodules;
        int rmodules;

	cout << "Input data file: ";
	cin >> dfile;
	nsteps = get_obs(dfile, indata);
	cout << "Number of time steps to aggregate for MCMC evaluations (suggest 1, 2, 4, 6, 8, 12):";
	cin >> n_agg;
//to hold obs in x, simulated in y
	Data_Frame Ecframe(nsteps/n_agg), NEEframe(nsteps/n_agg), Thetaframe(nsteps/n_agg); 
	Var_name ecname, neename;
	//strcpy(ecname, "simET");
	strcpy(ecname, "ET");
	//strcpy(theta30name, "theta30");
	strcpy(neename, "NEEobs");
/*
        for(i=0; i<nsteps; i++)
	{
                Ecframe.set_x_value(i, indata.get_val_for(ecname, i)); //OK
                NEEframe.set_x_value(i, indata.get_val_for(neename, i)); //OK
                //Thetaframe.set_x_value(i, indata.get_val_for(theta30name, i));
	}
*/
        cout << dfile << " contains " << indata.no_of_vars() << " variables " << 
							nsteps << " data points" << endl;
	cout << n_agg << " steps will be aggregated to produce " << nsteps/n_agg << 
							" data points for evaluation" << endl;
	cout << "Parameter file: ";
	cin >> pfile;
	np = lines(pfile);
	currp = new Parameter[np];
	checkp = read_params(pfile, np, currp);
	cout << pfile << " OK\n";
	rp = 0;
	for(i=0;i<np;i++)
	{
		if(currp[i].get_dtype() == RANDOM)
		{
			rp++;
		}
	}

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

//Parameter File Validation - **no validation done for shoot and root module read in (add later)

	if(rp>0)//MCMC
	{
                candp = new Parameter[np];
		copy_params(currp, candp, np);
		curpval = new double[rp];
		candpval = new double[rp];
		minpval = new double[rp];
		maxpval = new double[rp];
		tmpval = new double[rp];
		i=j=0;
		while(i<np && j<rp)
		{
			if(currp[i].get_dtype() == RANDOM)
			{
				curpval[j]=candpval[j]=currp[i].get_value();
				minpval[j]=currp[i].get_min();
				maxpval[j]=currp[i].get_max();
				j++;
			}
			i++;
		}
		if(i!=np && j!=rp)
		{
			cout << "parameter mismatch!\n";
			return 0;
		}
		cout << "Output file: ";
		cin >> cfile;
		cout << cfile << endl;
		strcpy(efile, cfile);
		strcat(cfile, ".chn");
		strcat(efile, ".eff");
		cout << "# of iterations to run: ";
		cin >> total_sims;
		cout << total_sims << endl;
		cout << "Initial covariance matrix file: ";
		cin >> vfile;
		cout << vfile << endl;
		cout << "Minimum target efficiency: ";
		cin >> mineff;
		cout << mineff << endl;
		cout << "Maximum target efficiency: ";
		cin >> maxeff;
		cout << maxeff << endl;
		cout << "# of samples to use for covariance updating: ";
		cin >> nsamples;
		cout << nsamples << endl;
		cout << "covariance update interval: ";
		cin >> updtg;
		cout << updtg << endl;
		cout << "use alternating chains for ET and NEE (1=yes; 0=no): ";
		cin >> use_alternating_chains;
		cout << use_alternating_chains << endl;
//end user input
		psamples_sml = new double*[rp];
		for (i = 0; i <rp; i++)
		{
			psamples_sml[i]=new double[nsamples];
		}
                psamples_med = new double*[rp];
                for (i = 0; i < rp; i++)
		{
                        psamples_med[i]=new double[2*nsamples];
		}
                psamples_lrg = new double*[rp];
                for (i = 0; i < rp; i++)
		{
                        psamples_lrg[i]=new double[4*nsamples];
		}
		covmat = new double*[rp];
		for (i = 0; i < rp; i++)
		{
			covmat[i]=new double[rp];
		}
		cholskf = new double*[rp];
		for (i = 0; i < rp; i++)
		{
			cholskf[i]=new double[rp];
		}
		open_infile(vfile, fin);
		for ( i = 0; i < rp; i++)
		{
			for (j = 0; j < rp; j++)
			{
				fin >> covmat[i][j];
			}
			fin.ignore(MAX_LINE, '\n');
		}
		fin.close();
		calc_cholesky(covmat, cholskf, rp);
		open_outfile(cfile, chnout);
		open_outfile(efile, effout);
		chnout << "iteration\t";
		for (i = 0; i < np; i++)
		{
			if (currp[i].get_dtype() == RANDOM)
			{
				currp[i].get_name(pname);
				chnout << pname << '\t';
			}
		}
		chnout << "Dtheta\n";
		effout << "input data file: " << dfile << endl;
		effout << "parameter file: " << pfile << endl;
		effout << "initial covariance matrix: " << vfile << endl;
		effout << "output file: " << cfile << endl;

        	cout << endl;
        	cout << ">>> TREES, v. 3.1.4, stochastic parameter mode <<< \n\n";
        	cout << endl;

		cout << "Running MCMC - " << total_sims << " iterations\n";
		if (use_alternating_chains == 1)
		{
			weight = 1.0; //start with Ecframe in the likelihood function
			cur_dens = simdensity(weight, indata, state, currp, treesParams, 
								Ecframe, NEEframe, curDtheta, n_agg);
			cur_dens_Ec = cur_dens;
			weight = 0.0; //start with NEEframe in the likelihood function
			cur_dens = simdensity(weight, indata, state, currp, treesParams, 
								Ecframe, NEEframe, curDtheta, n_agg);
			cur_dens_NEE = cur_dens;
		}
		else
		{
			weight = -1.0;
			cur_dens = simdensity(weight, indata, state, currp, treesParams, 
								Ecframe, NEEframe, curDtheta, n_agg);
		}
		n_accept = acc_pastupdt = 0;
		ipastupdt = lastupdt = 0;
		samplsize = SMALL;
		for (cursim = 1; cursim <= total_sims; cursim++) //MCMC run
		{
			ipastupdt++;
			chnout << cursim << '\t';
			for (i = 0; i < rp; i++)
			{
				chnout << curpval[i] << '\t';
			}
			chnout << curDtheta << endl;
			gen_proppval(candpval, curpval, tmpval, minpval, maxpval, rp, cholskf, xsubi);
			copy_proppval(candp, candpval, np, rp);
			p_jump = unit_random(xsubi);
			can_dens = simdensity(weight, indata, state, candp, treesParams, 
								Ecframe, NEEframe, canDtheta, n_agg);
			if (use_alternating_chains == 1)
			{
				if (weight == 1.0)
				{
					cur_dens = cur_dens_Ec;
				}
				else
				{
					cur_dens = cur_dens_NEE;
				}
			}
			if (can_dens > cur_dens || p_jump <= exp(can_dens - cur_dens))
			{
				//accept
				cur_dens = can_dens;
				if (use_alternating_chains == 1)
				{
					if (weight == 1.0)
					{
						cur_dens_Ec = can_dens; 
					}
					else
					{
						cur_dens_NEE = can_dens;
					}
				}
				copy_vect(candpval, curpval, rp);
				copy_proppval(currp, curpval, np, rp);
				curDtheta = canDtheta;
				n_accept++;
				acc_pastupdt++;
				isamplsml = n_accept%nsamples;
				isamplmed = n_accept%(2*nsamples);
				isampllrg = n_accept%(4*nsamples);
				for (i = 0; i < rp; i++)
				{
					psamples_sml[i][isamplsml] = curpval[i];
					psamples_med[i][isamplmed] = curpval[i];
					psamples_lrg[i][isampllrg] = curpval[i];
				}
			}
//changed this - 3/6/07 - look at mcEcl.cpp.1psampl for the previous version
			if (ipastupdt >= updtg && n_accept > nsamples)
			{
				cureff = (double)acc_pastupdt/(double)(cursim - lastupdt);
				if (cursim <= total_sims/2)
				{
					if (samplsize == SMALL)
					{
						if (n_accept > 6*nsamples)
						{
							calc_covmat(psamples_med, covmat, rp, 2*nsamples);
							calc_cholesky(covmat, cholskf, rp);
							samplsize = MEDIUM;
							effout << cursim << '\t' << weight << '\t' << 
									cureff << "\tto medium" << endl;
						}
						else
						{
							calc_covmat(psamples_sml, covmat, rp, nsamples);
							calc_cholesky(covmat, cholskf, rp);
							effout << cursim << '\t' << weight << '\t' << 
									cureff << "\tupdated" << endl;
						}
					}
					else if (samplsize == MEDIUM)
					{
						if (n_accept > 10*nsamples)
                                        	{
                                                	calc_covmat(psamples_lrg, covmat, rp, 4*nsamples);
                                                	calc_cholesky(covmat, cholskf, rp);
                                                	samplsize = LARGE;
							effout << cursim << '\t' << weight << '\t' << 
									cureff << "\tto large" << endl;
                                        	}
						else
						{
							calc_covmat(psamples_med, covmat, rp, 2*nsamples);
							calc_cholesky(covmat, cholskf, rp);
							effout << cursim << '\t' << weight << '\t' << 
									cureff << "\tupdated" << endl;
						}
					}
					else
					{
						if (cureff > maxeff || cureff < mineff)
						{
							calc_covmat(psamples_lrg, covmat, rp, 4*nsamples);
							calc_cholesky(covmat, cholskf, rp);
							strcpy(tempstr, "\tupdated");
						}
						effout << cursim << '\t' << weight << '\t' << 
									cureff << tempstr << endl;
					}

				}
				else
				{
					effout << cursim << '\t' << weight << '\t' << cureff << endl;
				}
				ipastupdt = cursim%updtg;
				lastupdt = cursim;
				acc_pastupdt = 0;

//switch likelihood function each time covariance is updated
				if (use_alternating_chains == 1)
				{
					if (weight == 1.0)
					{
						weight = 0.0;
						cur_dens = cur_dens_NEE;
					}
					else
					{
						weight = 1.0;
						cur_dens = cur_dens_Ec;
					}
				}
			}
		}//end for(cursim)
		effout << "---------------- Final Covarince Matrix -----------------\n";
		for (i = 0; i < rp; i++)
		{
			for (j = 0; j < rp; j++)
			{
				effout << covmat[i][j];
				if (j < rp-1)
				{
					effout << '\t';
				}
			}
			effout << endl;
		}
		effout.close();
		chnout.close();
		for (i = rp-1; i >= 0; i--)
		{
			delete [] psamples_sml[i];
			psamples_sml[i] = NULL;
                        delete [] psamples_med[i];
                        psamples_med[i] = NULL;
                        delete [] psamples_lrg[i];
                        psamples_lrg[i] = NULL;
			delete [] covmat[i];
			covmat[i] = NULL;
			delete [] cholskf[i];
			cholskf[i] = NULL;
		}
		delete [] psamples_sml;
		delete [] psamples_med;
		delete [] psamples_lrg;
		delete [] covmat;
		delete [] cholskf;
		delete [] candp;
		delete [] curpval; 
		delete [] candpval;
		delete [] minpval;
		delete [] maxpval;
		delete [] tmpval;
	}//end if(rp>0)
	else if (rp == 0)//simulation with constant parameters
	{
                cout << "Output file: ";
                cin >> sfile;
                cout << sfile << endl;
		strcpy(hfile, sfile);
        strcpy(lfile, sfile);
                strcat(sfile, ".sim");
                strcat(hfile, ".hyd");
                strcat(lfile, ".leaf");
		open_outfile(sfile, fluxout);
		open_outfile(hfile, hydrout);
        open_outfile(lfile, lfareaout);
		fluxout << "ti\tsimET\tWPlant_K\t";
		fluxout << "Soil_Psi\tLeaf_Psi\tPsi_Crit\tEcrit\tEc\t";
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "RhizFlux" << i << "\t";
                }
		fluxout << "Gs\tLAI\tSLA\tliveLAI\tRmaint\t";
		fluxout << "Rgrowth\treproduction\tleafNSC\tstemNSC\trootNSC\tchloroStarch\tchloroSugar\twaterStress\tlitterH2O\t";
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "theta" << i << "\t";
                }
		fluxout << "thetaRoot\t";
		fluxout << "Can_Evap\tSnowpack\tSnowEdef\tVcmax25\tVcmax_sun\tVcmax_shd\tJmax25\t";
		fluxout << "J_sun\tJ_shd\tAsun\tAshd\tLsun\tLshd\tTsun\tTshd\tDsun\tDshd\tCi_sun\t";
		fluxout << "Ci_shd\tPARsun\tPARshd\tgs_sun\tgs_shd\tNEE\tNPP\tR_total\tR_ag\tR_bg\t";
		fluxout << "Rd_sun\tRd_shd\tCsapwood\t";
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "FibRootC" << i << "\t";
                }
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "FineRootC" << i << "\t";
                }
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "TotRootC" << i << "\t";
                }
/*
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "CoarseRoot" << i << "\t";
                }
*/
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "FineRootCN" << i << "\t";
                }
		fluxout << "LeafCN\t";

//Croot1\tCroot2\tCroot3\tCroot4\tCroot5\tCroot6\tCroot7\t";
//		fluxout << "Croot8\tCroot9\tCroot10\tRhizCl\tRhizNl\t";
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "humusC" << i << "\t";
                }
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "RhizCl" << i << "\t";
                }
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "RhizNl" << i << "\t";
                }
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "AAexudateC" << i << "\t";
                }
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "SugarExudateC" << i << "\t";
                }
		//fluxout << "RhizCl\tRhizNl\t";
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "MicrobC" << i << "\t";
                }
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "MicrobN" << i << "\t";
                }
/*
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "Nleach" << i << "\t";
                }
*/
		fluxout << "RhizN-\tRhizN+\tPlantN\tPlantNstat\tRLA\t";
		for (i = 0; i < treesParams.rmodules; i++)
                {
                        fluxout << "ar" << i << "\t";
                }
		fluxout << endl;

		hydrout << "ti\tlatStemK\t";
		for (i = 0; i < treesParams.rmodules; i++)
                {
                	hydrout << "latRootK" << i << "\t";
                }
		hydrout << "StemAxialYm\tStemLatYm\t";
		for (i = 0; i < treesParams.rmodules; i++)
		{
			hydrout << "RootAxialYm" << i << "\t";
		}
		for (i = 0; i < treesParams.rmodules; i++)
		{
			hydrout << "RootLatYm" << i << "\t";
		}
		hydrout << "StemAxialKm\tStemLatKm\t";
		for (i = 0; i < treesParams.rmodules; i++)
		{
			hydrout << "RootAxialKm" << i << "\t";
		}
		for (i = 0; i < treesParams.rmodules; i++)
		{
			hydrout << "RootLatKm" << i << "\t";
		}
		hydrout << "StemAxial_b\tStemLat_b\t";
		for (i = 0; i < treesParams.rmodules; i++)
		{
			hydrout << "RootAxial_b" << i << "\t";
		}
		for (i = 0; i < treesParams.rmodules; i++)
		{
			hydrout << "RootLat_b" << i << "\t";
		}
		hydrout << "StemAxial_c\tStemLat_c\t";
		for (i = 0; i < treesParams.rmodules; i++)
		{
			hydrout << "RootAxial_c" << i << "\t";
		}
		for (i = 0; i < treesParams.rmodules; i++)
		{
			hydrout << "RootLat_c" << i << "\t";
		}
		hydrout << "\n";
        
        
        lfareaout << "ti\t";
        for (i = 1; i < 500; i++)
        {
            lfareaout << "Area_Leaf_" << i << "\t";
        }
        lfareaout << "\n";

		simulator(indata, state, currp, treesParams, Ecframe, NEEframe, fluxout, hydrout, lfareaout);
		fluxout.close();
		hydrout.close();
        lfareaout.close();
		cout << "\nSimulated values in " << sfile << endl;
		cout << "\nHydraulic diagnostics in " << hfile << endl;
	}//end if(rp==0)
	delete [] currp;
	return 1;
}

void copy_params(Parameter* origp, Parameter* copyp, int np)
{
	for (int i = 0; i < np; i++)
	{
		copyp[i] = origp[i];
	}
}

double logpostdens(const Data_Frame& df, double sd, double& Dtheta)
{
        int points, validPoints;
        double d_val, sum_sq, xval, yval;
 
        sum_sq = 0.0;
        validPoints = points = df.how_many();
        for (int i = 0; i < points; i++)
	{
		xval = df.get_x(i);
		yval = df.get_y(i);
		if (xval == INVALID || yval == INVALID)
		{
			validPoints--;
		}
		else
		{
                	sum_sq += pow((xval-yval),2.0);
		}
	}
        d_val = -((points+2.0)*log(sd) + sum_sq/(2.0*pow(sd,2.0)));
	Dtheta = 2.0*(double)points*log(sd) + sum_sq/pow(sd,2.0);
	return d_val;
}

//Modified 07/14/09 DSM
//Uses two weighted likelihood functions 
// - assumes sd1_weight is between 0 and 1
double logpostdens(const Data_Frame& df1, const Data_Frame& df2, double sd1, double sd2, 
					double sd1_weight, double sd2_scalar, double& Dtheta)
{
        int points1, points2, validPoints1, validPoints2;
        double d_val, d_val1, d_val2, sum_sq1, sum_sq2, xval1, yval1, xval2, yval2, Dtheta1, Dtheta2;
//log posterior density for first data frame
        sum_sq1 = 0.0;
        validPoints1 = points1 = df1.how_many();
        for (int i = 0; i < points1; i++)
        {
                xval1 = df1.get_x(i);
                yval1 = df1.get_y(i);
                if (xval1 == INVALID || yval1 == INVALID)
		{
                        validPoints1--;
		}
                else
		{
                        sum_sq1 += pow((xval1-yval1),2.0);
		}
        }
        d_val1 = -((validPoints1+2.0)*log(sd1) + sum_sq1/(2.0*pow(sd1,2.0)));
        Dtheta1 = 2.0*(double)validPoints1*log(sd1) + sum_sq1/pow(sd1,2.0);

//log posterior density for second data frame
        sum_sq2 = 0.0;
        validPoints2 = points2 = df2.how_many();
        for (int i = 0; i < points2; i++)
        {
                xval2 = df2.get_x(i);
                yval2 = df2.get_y(i);
                if (xval2 == INVALID || yval2 == INVALID)
		{
                        validPoints2--;
		}
                else
		{
			xval2 = xval2*sd2_scalar;
			yval2 = yval2*sd2_scalar;
                        sum_sq2 += pow((xval2-yval2),2.0);
		}
        }
        d_val2 = -((validPoints2+2.0)*log(sd2) + sum_sq2/(2.0*pow(sd2,2.0)));
        Dtheta2 = 2.0*(double)validPoints2*log(sd2) + sum_sq2/pow(sd2,2.0);

//weighted log posterior density using both data frames
	d_val = (sd1_weight*d_val1 + (1.0-sd1_weight)*d_val2);
	Dtheta = (sd1_weight*Dtheta1 + (1.0-sd1_weight)*Dtheta2);
        return d_val;
}


double simdensity(double weight, Data_Store& inputdata, State_Store& state, Parameter* parameters, 
	trees_params& treesParams, Data_Frame& Ecframe, Data_Frame& NEEframe, double& Dtheta, int n_agg)
{
	double sd1, sd2, sd1_weight, sd2_scalar, lpd;
	simulator(inputdata, state, parameters, treesParams, Ecframe, NEEframe, sd1, sd2, sd1_weight, n_agg);
	sd2_scalar = 3.0E-6; //scale magnitude of error on NEE to match order of magnitude for ET
	sd2 = sd2*sd2_scalar;
	if (weight == -1.0) //use weighting from .p file
	{
		lpd = logpostdens(Ecframe, NEEframe, sd1, sd2, sd1_weight, sd2_scalar, Dtheta);
	}
	else
	{
		lpd = logpostdens(Ecframe, NEEframe, sd1, sd2, weight, sd2_scalar, Dtheta);
	}
	return lpd;
}

double simdensity(Data_Store& inputdata, State_Store& state, Parameter* parameters, 
	trees_params& treesParams, Data_Frame& Ecframe, Data_Frame& NEEframe, double& Dtheta, int n_agg)
{
	double sd1, sd2, sd1_weight, sd2_scalar, lpd;
	simulator(inputdata, state, parameters, treesParams, Ecframe, NEEframe, sd1, sd2, sd1_weight, n_agg);
	sd2_scalar = 3.0E-6; //scale magnitude of error on NEE to match order of magnitude for ET
	sd2 = sd2*sd2_scalar;
	lpd = logpostdens(Ecframe, NEEframe, sd1, sd2, sd1_weight, sd2_scalar, Dtheta);
	return lpd;
}

void gen_proppval(double* candpval, double* curpval, double* tempval, double* minlim, 
				double* maxlim, int rp, double** cholf, Xsubi_array& xsubi)
{
	int range_chk;
	range_chk = 1;
	while(range_chk != 0)
	{
		corNormal(candpval, tempval, curpval, cholf, rp, xsubi);
		range_chk = check_range(candpval, minlim, maxlim, rp);
	}
}

int check_range(double* candpval, double* minlim, double* maxlim, int rp)
{
	int out_range;
	out_range = 0;
	for (int i = 0; i < rp; i++)
	{
		if (candpval[i] < minlim[i] || candpval[i] > maxlim[i])
		{
			out_range++;
		}
	}
	return out_range;
}

void copy_proppval(Parameter* candp, double* candpval, int np, int rp)
{
	int i, j;
	i = j = 0;
	while (i < np && j < rp)
	{
		if(candp[i].get_dtype() == RANDOM)
		{
			candp[i].set_value(candpval[j]);
			j++;
		}
		i++;
	}
}


