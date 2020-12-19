//
//hydraulicModel.cpp
//
//DoHydraulicModel(): 
//    Implements John Sperry's 1998 Plant water balance model, Plant Cell and Environment
//
//    Original model was written in Visual Basic and attached to Excel; 
//        some ideosynchratic elements of VB remain
//    Conversion from VB to C++ by Aga Shirazi; implementation as part of TREES by David Roberts, 
//        and full functionality by D. Scott Mackay
//
//Some components of the original model were turned off so that it could be controlled by TREES
//
//Components related to memory of cavitation differ from their original implementation
//
//Modifications made to this code can produce unexpected behavior, and so modify with caution
//
//
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "simulator2.h"


/*
*	-	For setup(), nodeTyp should be changes from strign to String
*	-	Incremental values of i, j, k need to be locatlized.
*	-	convert/uncomment   nodetype = "lateral shoot";
*/

               
void DoHydraulicModel(bool silent,
			HydraulicModel *hydraulicModel, 
			int callHydraulicModel,
			bool callRewet,
			double Ysoil[],
                	double  ev, //don't confuse this with e its e_raw!
			double& Kl,
			double& leafpsi,
			double& Ecrit,
			int &HydraulicModelFailCond,
			double& PsiCrit,
			int hydraulicModel_flag,
			trees_params& treesParams,
			HydraulicModel::outputStruct &ecrit_k_psi, 
			HydraulicModel::outputStruct &k_p_e, 
			double ecritOut[], double pc[], double klpred[],
                    	double ppredOut[], double epredOut[], double nodeFail[], 
			char nodeTyp[], double ***psi, double **rflux, 
			double soilpsiavg[], double nts[], double evap[], double kroot[],
                    	double axr[], double latr[], double kshoot[], 
			double dslat[], double dsax[], double drlat[],
                    	double drax[], double l[], double ksat[][MD], 
			double bsat[][MD], double ccsat[][MD],
                    	double newb[][MD], double newc[][MD], 
			double ip[][ULAT], double b[][ULAT], double b1[][ULAT],
                    	double n[][ULAT], double n1[][ULAT], double r[][ULAT], 
			double vu[][ULAT], double vl[][ULAT],
                    	double cf[][ULAT], double pe[][ULAT], double ws_[][ULAT], 
			double ks[][ULAT], double p[][ULAT],
                    	double dpp[][ULAT], double jl[][ULAT], 
			double wnu[][ULAT], double wu[][ULAT], double wnl[][ULAT],
                    	double wl[][ULAT], double cpu[][ULAT], double cpl[][ULAT], 
			double cpp[][ULAT], double ku[][ULAT],
                    	double kl[][ULAT], double f[][ULAT], 
			double cc[][2], double psinode[][NMAX], double psimin[][NMAX],
                    	double jmatrix[][NMAX], double jmatrix2[][NMAX], 
			double percent[], double rssoil[], double rs[],
                    	double dp[], double ff[], int col[], int row[], int indxx[], 
			double subtract[], double pressure[],
                    	double plc[], double plcweib[], double rsquare[], 
			double soilpsi[], double al[], double ar[],
                    	double saturatedKs[][4], double &ptarg, double &e, 
			double &rr, double &rzw, double &rd, double &einc,
                    	int &SoilRhizElements, double &dt, double &gmd, double &gsd, 
			double &bkd, double &fs, double &fc,
                    	double &Fabs, double &cpplant, double &saxlat, double &raxlat, 
			double &rfract, double &kla,
                    	double &aral, double &latot, double &ratot, 
			int &shootElements, int &ktotal, double &soilpsimin,
                    	int &ShootModules, int &rootElements, 
			int &RootModules, int &totmodules, double &axShoot_b,
                    	double &axShoot_c, double &latShoot_b, double &latShoot_c, 
			double &axRoot_b, double &axRoot_c,
                    	double &latRoot_b, double &latRoot_c, 
			double &plco, double &pdecrement, double &sumy1,
                    	double &sumy2, double &sumprod, double &cincrement, 
			double &soilvol, double &rlateral,
                    	double &rLat_base, double &slateral, double &sLat_base, 
			double &raxial, double &pleafave,
                    	int &tnode, double fac, double &modelledKL)
{
	double Ypd[ULAT], maxPredawn;

	for (int i = 0; i < treesParams.rmodules; i++)
	{
		//maxPredawn = 0.5*(Ysoil[i]+treesParams.pd_at_sat_kl);
		maxPredawn = 0.2*Ysoil[i]+0.8*treesParams.pd_at_sat_kl;
		Ypd[i] = min(Ysoil[i], maxPredawn);
	}
/*
	for (int i = 0; i < treesParams.rmodules; i++)
	{
		maxPredawn = 0.8*treesParams.pd_at_sat_kl;
		Ypd[i] = min(Ysoil[i], maxPredawn);
	}
*/

// zero transpiration not handled gracefully, and so we use 0.1% of saturated transpiration as a minimum
// DSM May 2012
		
	if (ev < 0.001*treesParams.e_at_saturated_kl)
	{
		ev = 0.001*treesParams.e_at_saturated_kl;
	}		
	
// a hydraulicModel object no longer needs to take in any params since we did away with the read in
// DER8/16/07
		
       HydraulicModelFailCond = 0;
 	       
       hydraulicModel->solveWater(silent, callRewet, Ypd, ev, Kl,
			       	leafpsi, Ecrit, HydraulicModelFailCond,
			       	PsiCrit, hydraulicModel_flag,
			       	callHydraulicModel, treesParams,
			       	ecrit_k_psi, k_p_e, ecritOut, pc, klpred, ppredOut, epredOut,
                               	nodeFail, nodeTyp, psi, rflux, soilpsiavg, nts, evap, kroot,
                               	axr, latr, kshoot, dslat, dsax, drlat, drax, l, ksat, bsat,
                               	ccsat, newb, newc, ip, b, b1, n, n1, r, vu, vl, cf, pe, ws_, ks,
                               	p, dpp, jl, wnu, wu, wnl, wl, cpu, cpl, cpp, ku, kl, f, cc, psinode,
                               	psimin, jmatrix, jmatrix2, percent, rssoil, rs, dp, ff, col, row,
                               	indxx, subtract, pressure, plc, plcweib, rsquare, soilpsi, al, ar,
                               	saturatedKs, ptarg, e, rr, rzw, rd, einc, SoilRhizElements, dt, gmd,
                               	gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract, kla, aral,
                               	latot, ratot, shootElements, ktotal, soilpsimin, ShootModules, 
				rootElements, RootModules, totmodules, axShoot_b, axShoot_c, 
				latShoot_b, latShoot_c, axRoot_b, axRoot_c, latRoot_b, latRoot_c, 
				plco, pdecrement, sumy1, sumy2, sumprod,
                               	cincrement, soilvol, rlateral, rLat_base, slateral, sLat_base, raxial,
                               	pleafave, tnode, fac, modelledKL);
				
}


void free3DArray(double ***array, int x, int y) 
{
    	for (int i = 0; i < x; ++i) 
	{
        	for (int j = 0; j < y; ++j)
		{
            		free(array[i][j]);
		}
        	free(array[i]);
    	}
    	free(array);
}

void free2DArray(double **array, int x) 
{
    	for (int i = 0; i < x; ++i) 
	{
        	free(array[i]);
    	}
    	free(array);
}

double*** get3DArray(int x, int y, int z) 
{
    	double ***ret;
    	ret = (double ***)malloc(sizeof(*ret) * x);
    	if (ret == NULL)
	{
        	cerr << "ERROR\n";
	}
    	for (int i = 0; i < x; ++i) 
	{
        	ret[i] = (double **)calloc(y, sizeof(**ret));
        	if (ret[i] == NULL)
		{
            		cerr << "ERROR WITH malloc";
		}
        	for (int j = 0; j < y; ++j) {
            		ret[i][j] = (double *)calloc(z, sizeof(***ret));
            		if (ret[i][j] == NULL)
			{
                		cerr << "ERROR WITH MALLOC\n";
			}
        	}
    	}
    	return ret;
}

double** get2DArray(int x, int y) 
{
    	double **ret;
    	ret = (double **)calloc(x, sizeof(*ret));
    	if (ret == NULL)
	{
        	cerr << "ERROR\n";
	}
    	for (int i = 0; i < x; ++i) 
	{
        	ret[i] = (double *)calloc(y, sizeof(**ret));
        	if (ret[i] == NULL)
		{
            		cerr << "ERROR WITH malloc";
		}
    	}
    	return ret;
}


HydraulicModel::HydraulicModel()
{

}


//The primary purpose of this function is to update pressures and Weibull parameters at tne end of 
//  a period of drought. The function here is modified to allow either refilling of plant modules
//  or to retain full memory of the drought, using the "useRefilling" parameter
//  [ DSM January 2012 update ]

void HydraulicModel::readin(bool silent, bool callRewet, 
				double soilpsi[ULAT], double soilpsiavg[1], double ks[][ULAT],
            			double pressure[NKINC], double plc[NKINC], double **l_raw,
            			double ksat[MD][MD], double psimin[NMAX][NMAX], double bsat[MD][MD],
            			double ccsat[MD][MD], double newc[MD][MD], double plcweib[NKINC],
            			double rsquare[NKINC], double b[MD][ULAT], double newb[MD][MD],
            			double p[MD][ULAT], double md_raw[], double e_raw[],
            			double la_raw[], double arAl[], double cc[MD][2],
            			const int ShootModules, const int totmodules, const int RootModules,
            			const int SoilRhizElements, const int rootElements, 
				double &soilpsimin, const int t, double &e, 
				double &latot, double &aral,
            			double &ratot, double &ptarg,
	    			trees_params& treesParams) 
{

    	int i = 1, k=0, ktotal, j=0;
     	double plco=0.0, pdecrement=0.0, cincrement=0.0, sumy1=0.0, sumy2=0.0, sumprod=0.0;
	double totArea, delta;
	double xylemScalar = treesParams.xylemScalar;

       	ptarg = md_raw[i];
       	e = e_raw[i];
       	//latot = la_raw[i];
       	aral = arAl[i];
       	//ratot = aral * latot;        //sets total ra

	for (j = ShootModules+2; j < ShootModules+RootModules+2; j++)
	{
		soilpsi[j] = l_raw[i][j-2];
	}

//Average soilpsi for i root modules with root area proportions defined in param_mod input file
	soilpsiavg[t] = 0.0;
/*
	for (i = ShootModules+2; i < ShootModules+RootModules+2; i++)
	{
		soilpsiavg[t] += treesParams.ar[i]*soilpsi[i];
	}
*/
	totArea = 0.0;
	for (i = ShootModules+2; i < ShootModules+RootModules+2; i++)
	{
		if (psimin[i][1] != -8.0)
		{
			totArea += treesParams.ar[i];
		}
	}
	
	for (i = ShootModules+2; i < ShootModules+RootModules+2; i++)
	{
		if (psimin[i][1] != -8.0)
		{
			soilpsiavg[t] += treesParams.ar[i]*soilpsi[i]/totArea;
		}
	}

	delta = soilpsiavg[t] - soilpsimin;
	if (delta > 0.005 || treesParams.updatedHydraulics == true)
	{
		callRewet = true; //update Weibull parameters, not necessarry re-fill xylem
	}

	treesParams.useRefilling = false;

//A rewet involves resetting Ks and Weibull parameters only if useRefilling is true
//Now done when there is soil water recharge or once BB attack occurs and callRewet is true
//if there's no rewetting, then keep weibull curves and ksat's
	if (callRewet == false) 
	{ 
		if (soilpsiavg[t] < soilpsimin)
		{
            		soilpsimin = soilpsiavg[t]; // minimum since last wetting event
		}
      	}
//if drought is relieved and there is ability to re-fill, then ksat and weibulls need resetting
	else 
	{
		if (treesParams.forceRefilling == true)
		{
			treesParams.useRefilling = true; //ability to re-fill
		}
         	soilpsimin = soilpsiavg[t];

//Here is where a resetting of pressures should take place if wetting up --> re-filling
		if (treesParams.useRefilling == true)
		{
			for (i = 1; i <= totmodules+1; i++)
			{
				if ( i != (ShootModules+1))  // skip over "hidden" module
				{
					if (i <= ShootModules)
					{
						if (treesParams.allowLeafRefilling == true)
						{
							psimin[i][1] = soilpsiavg[t];
						}
						psimin[i][0] = soilpsiavg[t];
					}
					else if (i > (ShootModules+1))
					{
						psimin[i][1] = soilpsi[i];
						psimin[i][0] = soilpsi[i];
					}
				}
			}
//for large changes in soil water potential, also adjust precision of E increment
			if (delta > 0.5)
			{
				treesParams.E_inc *= 4.0*delta;
			}
		}

		if (!silent)
		{
			if (xylemScalar < 1.0)
			{
				cout << endl << "-> Predawn: updating Weibulls with xylem loss = " << 
						100.0*(1.0-xylemScalar) << "%" << endl;
				if (xylemScalar == 0.96 || xylemScalar == 0.97 || xylemScalar == 0.98 || xylemScalar == 0.99)
					xylemScalar = 1.0;
			}
			else if (xylemScalar > 1.0)
			{
				cout << endl << "-> Predawn: updating Weibulls with xylem loss = " << 
						100.0*(xylemScalar-1.0) << "%" << endl;
			}
			else
			{
				//cout << endl << "-> Predawn: updating Weibulls with soil Psi change = " <<
				//		 delta << endl;
			}
		}

                pdecrement = 0.01;
                cincrement = 0.01;
          	for (i = 1; i <= totmodules+1; i++) 
		{	
//resets data for lateral modules. i is plant element, k pressure increment, j c iteration
             		if (i != (ShootModules+1))  // skip over "hidden" module
			{
//lateral shoot(s) infected with bluestain fungi
				if (i <= ShootModules) 
				{
                  			ks[i][1] = xylemScalar*ksat[i][1] * 
						exp(-pow((-psimin[i][1] / bsat[i][1]), ccsat[i][1]));
//plc at new predawn
                  			plco = (1.0 - xylemScalar*exp(-pow((-psimin[i][1] / 
								bsat[i][1]), ccsat[i][1]))) * 100.0; 
				}
				else
				{
                  			ks[i][1] = xylemScalar*ksat[i][1] * 
						exp(-pow((-psimin[i][1] / bsat[i][1]), ccsat[i][1]));
//plc at new predawn
                  			plco = (1.0 - xylemScalar*exp( -pow((-psimin[i][1] / 
								bsat[i][1]), ccsat[i][1]))) * 100.0; 
				}
                  		k = 0;
                  		pressure[k] = 0.0;

                  		do 
				{ //loop calculates new plc's based on new predawn and finds new b
                      			k++;
                      			pressure[k] = pressure[k-1] - pdecrement;
                      			if ( pressure[k] >= psimin[i][1] ) 
					{
                          			plc[k] = plco;
                      			}
					else
					{
//lateral shoot(s) infected with bluestain fungi
						if (i <= ShootModules) 
						{
                          				plc[k] = (1.0-xylemScalar*exp(-pow((-pressure[k] /
								 bsat[i][1]), ccsat[i][1]))) * 100.0;
						}
						else
						{
                          				plc[k] = (1.0-xylemScalar*exp(-pow((-pressure[k] /
								 bsat[i][1]), ccsat[i][1]))) * 100.0;
						}
				   	}
                          		plc[k] = (1.0 - (100.0 - plc[k]) / (100.0 - plco)) * 100.0;
                          		if ( plc[k] < 66.0 ) 
					{
                              			newb[i][1] = fabs(pressure[k]);
                          		}
                 		} while ( plc[k] <= 99.50 && k < NKINC-1);

                 		ktotal = k - 1;
                 		newc[i][1] = ccsat[i][1];
                 		j = 0;
                 		rsquare[j] = 0.0;

                 		do 
				{ //loop finds new c
                      			j++;
                      			newc[i][1] += cincrement;
                      			sumy1 = 0.0;
                      			sumy2 = 0.0;
                      			sumprod = 0.0;
                      			for (k = 1; k <= ktotal; k++) 
					{
						if (i <= ShootModules) 
						{
//lateral shoot(s) infected with bluestain fungi
                           				plcweib[k] = (1.0 - exp(-pow((-pressure[k] / 
								newb[i][1]), newc[i][1]) )) * 100.0;
						}
						else
						{
                           				plcweib[k] = (1.0 - exp(-pow((-pressure[k] / 
								newb[i][1]), newc[i][1]) )) * 100.0;
						}
                           			sumy1 += pow(plcweib[k], 2.0);
                           			sumy2 += pow(plc[k], 2.0);
                           			sumprod += plcweib[k] * plc[k];
                      			}
                      			rsquare[j] = pow(sumprod, 2.0) / (sumy1 * sumy2);
                  		} while ( rsquare[j] >= rsquare[j-1] );

                  		newc[i][1] -= cincrement;
//assigns new b's and c's to plant
                  		b[i][1] = newb[i][1];    
                  		cc[i][1] = newc[i][1];
              		}	// end of if
          	}	// end of for

//resets data for axial modules. i is plant element, k pressure increment, j c iteration
        	for ( i = 1; i <= totmodules+1; i++) 
		{
           		if ( i != ShootModules + 1 ) 
			{
				if (i <= ShootModules) 
//axial shoot(s) infected with bluestain fungi
				{
                			ks[i][0] = xylemScalar*ksat[i][0] * exp(-pow((-psimin[i][0] / 
								bsat[i][0]), ccsat[i][0]) );
//plc at new predawn
                			plco = (1.0 - xylemScalar*exp(-pow((-psimin[i][0] / 
								bsat[i][0]), ccsat[i][0]))) * 100.0; 
				}
				else
				{
                			ks[i][0] = xylemScalar*ksat[i][0] * exp(-pow((-psimin[i][0] / 
								bsat[i][0]), ccsat[i][0]) );
//plc at new predawn
                			plco = (1.0 - xylemScalar*exp(-pow((-psimin[i][0] / 
								bsat[i][0]), ccsat[i][0]))) * 100.0; 
				}
                		k = 0;
                		pressure[k] = 0.0;
                		do 
				{   //loop calculates new plc's based on new predawn and finds new b
                     			k++;
                     			pressure[k] -= pdecrement;
                     			if ( pressure[k] >= psimin[i][0] ) 
					{
                        			plc[k] = plco;
                     			}
					else
					{
//lateral shoot(s) infected with bluestain fungi
						if (i <= ShootModules) 
						{
                          				plc[k] = (1.0-xylemScalar*exp(-pow((-pressure[k] /
								 bsat[i][0]), ccsat[i][0])) ) * 100.0;
						}
						else
						{
                          				plc[k] = (1.0-xylemScalar*exp(-pow((-pressure[k] /
								bsat[i][0]), ccsat[i][0])) ) * 100.0;
						}
                     			}
                     			plc[k] = (1.0 - (100.0 - plc[k]) / (100.0 - plco)) * 100.0;
                     			if ( plc[k] < 66.0 ) 
					{
                          			newb[i][0] = fabs(pressure[k]);
                     			}
                 		} while ( plc[k] <= 99.50 && k < NKINC-1);

                 		ktotal = k - 1;
                 		newc[i][0] = ccsat[i][0];
				j = 0;
                 		rsquare[j] = 0.0;

                 		do 
				{  //loop finds new c
                      			j++;
                      			newc[i][0] += cincrement;
                      			sumy1 = 0.0;
                      			sumy2 = 0.0;
                      			sumprod = 0.0;
                      			for ( k = 1; k <= ktotal; k++) 
					{
//lateral shoot(s) infected with bluestain fungi
						if (i <= ShootModules) 
						{
                           				plcweib[k] = (1.0 - exp(-pow((-pressure[k] / 
								newb[i][0]), newc[i][0])) ) * 100.0;
						}
						else
						{
                           				plcweib[k] = (1.0 - exp(-pow((-pressure[k] / 
								newb[i][0]), newc[i][0])) ) * 100.0;
						}
                           			sumy1 += pow(plcweib[k], 2.0);
                           			sumy2 += pow(plc[k], 2.0);
                           			sumprod += plcweib[k] * plc[k];
                      			}
                      			rsquare[j] = pow(sumprod, 2.0) / (sumy1 * sumy2);
                 		} while ( rsquare[j] >= rsquare[j-1] );

                 		newc[i][0] -= cincrement;
//assigns new b's and c's to plant
                 		b[i][0] = newb[i][0];        
                 		cc[i][0] = newc[i][0];
                       }	// end of if
              	}	// end of for
	}	//	end of else

//assign initial pressures to plant
    	for( i = 1; i <= totmodules+1; i++) 
	{
      		if (i <= ShootModules+1) 
		{
          		p[i][0] = soilpsiavg[t];
          		p[i][1] = soilpsiavg[t];
      		}
		else
		{
          		p[i][0] = soilpsi[i];
          		p[i][1] = soilpsi[i];
      		}
    	}

//assign pressures to soil elements
    	for (i = ShootModules+2; i <= totmodules+1; i++) 
	{
        	for ( j = 2; j <= 1+SoilRhizElements; j++) 
		{
             		p[i][j] = soilpsi[i];
        	}
    	}
} // end of readin()


void HydraulicModel::solveWater(bool silent, bool callRewet, double Ypd[],
				double&  ev, double&  Kl, double&  leafpsi, double&  Ecrit,
				int &HydraulicModelFailCond,
				double &PsiCrit, int hydraulicModel_flag, int callHydraulicModel,
				trees_params& treesParams,
				outputStruct &ecrit_k_psi, outputStruct &k_p_e, 
				double ecritOut[], double pc[], double klpred[],
                    		double ppredOut[], double epredOut[], double nodeFail[], 
				char nodeTyp[], double ***psi, double **rflux, 
				double soilpsiavg[], double nts[], double evap[], double kroot[],
                    		double axr[], double latr[], double kshoot[], 
				double dslat[], double dsax[], double drlat[],
                    		double drax[], double l[], double ksat[][MD], 
				double bsat[][MD], double ccsat[][MD],
                    		double newb[][MD], double newc[][MD], 
				double ip[][ULAT], double b[][ULAT], double b1[][ULAT],
                    		double n[][ULAT], double n1[][ULAT], double r[][ULAT], 
				double vu[][ULAT], double vl[][ULAT],
                    		double cf[][ULAT], double pe[][ULAT], 
				double ws_[][ULAT], double ks[][ULAT], double p[][ULAT],
                    		double dpp[][ULAT], double jl[][ULAT], 
				double wnu[][ULAT], double wu[][ULAT], double wnl[][ULAT],
                    		double wl[][ULAT], double cpu[][ULAT], double cpl[][ULAT], 
				double cpp[][ULAT], double ku[][ULAT],
                    		double kl[][ULAT], double f[][ULAT], 
				double cc[][2], double psinode[][NMAX], double psimin[][NMAX],
                    		double jmatrix[][NMAX], double jmatrix2[][NMAX], 
				double percent[], double rssoil[], double rs[],
                    		double dp[], double ff[], int col[], int row[], 
				int indxx[], double subtract[], double pressure[],
                    		double plc[], double plcweib[], double rsquare[], 
				double soilpsi[], double al[], double ar[],
                    		double saturatedKs[][4], double &ptarg, double &e, 
				double &rr, double &rzw, double &rd, double &einc,
                    		int &SoilRhizElements, double &dt, 
				double &gmd, double &gsd, double &bkd, double &fs, double &fc,
                    		double &Fabs, double &cpplant, double &saxlat, double &raxlat, 
				double &rfract, double &kla,
                    		double &aral, double &latot, double &ratot, 
				int &shootElements, int &ktotal, double &soilpsimin,
                    		int &ShootModules, int &rootElements, int &RootModules, 
				int &totmodules, double &axShoot_b,
                    		double &axShoot_c, double &latShoot_b, double &latShoot_c, 
				double &axRoot_b, double &axRoot_c,
                    		double &latRoot_b, double &latRoot_c, 
				double &plco, double &pdecrement, double &sumy1,
                    		double &sumy2, double &sumprod, double &cincrement, 
				double &soilvol, double &rlateral,
                    		double &rLat_base, double &slateral, 
				double &sLat_base, double &raxial, double &pleafave,
                    		int &tnode, double fac, double &modelledKL)
{

    	double final_pe, final_ks, final_b, final_ws, final_soilvol, final_rd;
    	int nodefailure = 0, nits, ii;
    	int i, j, t = 1, h, k, m, check, kkk; //index for number of lines in data set
    	bool error = false;
   
// input array containers
    	double *md_raw, *e_raw, *la_raw, *arAl, **l_raw;
    	double base1, base2, exponent, result1, result2;
    	double oldpress, sumpress, se, aa, nn, phiu, gammp, gamma; 
	double phil, lim, abv, deltap, Sum, epred, ppred, pcrit, x, ecrit; 

    	ecritOut = new double[MD], pc = new double[MD],  klpred = new double[MD];
    	ppredOut = new double[MD], epredOut = new double[MD], nodeFail = new double[MD];
    	md_raw = new double[MD], e_raw = new double[MD], la_raw = new double[MD], arAl = new double[MD];
    	l_raw = get2DArray(5,MD);
    
/*
    	preSetup( md_raw, treesParams.Md,  e_raw, ev,  la_raw, 
			treesParams.Al, arAl, treesParams.Ar_Al, l_raw, Ypd);
*/
    	preSetup( md_raw, treesParams.Md,  e_raw, ev,  la_raw, 
			treesParams.lai, arAl, treesParams.Ar_Al, l_raw, Ypd);

//Original HydraulicModel code looped through time steps here
//Here it retained memory of drought effects on k and Weibulls
//To retain this memory setup needs to be moved into the TREES end

        readin(silent, callRewet, soilpsi,  soilpsiavg,  ks, pressure,  plc,  l_raw, ksat,  
			psimin,  bsat, ccsat, newc, plcweib, rsquare,  b,  newb, p,  
			md_raw,  e_raw, la_raw,  arAl,  cc,ShootModules,
               		totmodules,   RootModules, SoilRhizElements, rootElements, soilpsimin, t, e,
               		latot,  aral,ratot,  ptarg, treesParams);

	free2DArray(l_raw, 5);

        initialize(r, pe, b, b1, n, n1, ks, ws_, rssoil, cf, vu, vl, saturatedKs, cpp, wnu, wu, wnl,
                   	wl, ip, p, al, l, SoilRhizElements, ShootModules, totmodules, latot, 
			RootModules, rr, rd, Fabs, cpplant, gsd, fc, fs, bkd, gmd, soilvol );

        nits = 0;	 //sets # iterations to 0
        oldpress = 0.0;
        sumpress = 0.0;
        ii = 0;
	
//***************************************************************************************************
//instead of incrementing using "EINC", i have set the initial value as E calculated by PM from TREES
//***************************************************************************************************
	evap[ii] = ev;

        do 
	{ //this loop increments Q
            	h = 0;
            	error = false;
            	do 
		{ //this loop finds steady-state at each Q
                	error = false;
                	oldpress = sumpress;
                	nits = 0;
                	do 
			{ //this loop finds new potentials and contents at each time step
		    		error = false;
                    		se = 0.0; //se is sum of node mass balance
                    		for (  i = 1; i <= ShootModules; i++ ) 
		    		{
                        		jl[i][2] = evap[ii] * al[i]; //sets evaporation rate
                    		}
//makes sure outermost soil node is constant at psisoil
                    		for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{ 
                        		p[i][SoilRhizElements + 1] = soilpsi[i];
                    		}
//sets upper and lower elements for lateral elements of shoot
                    		for (i = 1; i <= ShootModules; i++ ) 
				{ 
                        		if ( p[i][1] > ip[i][1] ) 
					{
                            			wnu[i][1] = ws_[i][1] * vu[i][1];
                            			cpu[i][1] = 0.01 * vu[i][1] / dt;
                        		}
                        		else 
					{
//content of upper half
                            			wnu[i][1] = ws_[i][1] * vu[i][1] - 
							( ( ip[i][1] - p[i][1] ) * 0.5 * cpp[i][1] ); 
//dw/dp *1/dt for jacobian
                            			cpu[i][1] = 0.5 * cpp[i][1] / dt; 
                        		}

                        		if ( p[i][0] > ip[i][1] ) 
					{
                            			wnl[i][1] = ws_[i][1] * vl[i][1];
//strictly][this should be zero][but set to arbitrary small #?
                            			cpl[i][1] = 0.01 * vl[i][1] / dt; 
                        		}
                        		else 
					{
                            			wnl[i][1] = ws_[i][1] * vl[i][1] - 
							( ( ip[i][1] - p[i][0] ) * 0.5 * cpp[i][1] ); 
// content][lower half of plant element
//dw/dp *1/dt for jacobian
                            			cpl[i][1] = 0.5 * cpp[i][1] / dt; 
                        		}
                    		}

                    		for (  i = 1; i <= ShootModules; i++ ) 
				{ //sets upper and lower elements for axial elements of shoot
                        		if ( p[i][0] > ip[i][0] ) 
					{
                            			wnu[i][0] = ws_[i][0] * vu[i][0];
                            			cpu[i][0] = 0.01 * vu[i][0] / dt;
                        		}
                        		else 
					{
//content of upper half
                            			wnu[i][0] = ws_[i][0] * vu[i][0] - 
							( ( ip[i][0] - p[i][0] ) * 0.5 * cpp[i][0] ); 
//dw/dp *1/dt for jacobian
                            			cpu[i][0] = 0.5 * cpp[i][0] / dt; 
                        		}

                        		if ( p[i + 1][0] > ip[i][0] ) 
					{
                            			wnl[i][0] = ws_[i][0] * vl[i][0];
//[Strictly][this should be zero][but set to arbitrary small #?
                            			cpl[i][0] = 0.01 * vl[i][0] / dt; 
                        		}
                        		else 
					{
                            			wnl[i][0] = ws_[i][0] * vl[i][0] - 
							( ( ip[i][0] - p[i + 1][0] ) * 0.5 * cpp[i][0] );
// content][lower half of plant element
//dw/dp *1/dt for jacobian
                            			cpl[i][0] = 0.5 * cpp[i][0] / dt; 
                        		}
                    		}

                    		for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{ //sets upper and lower elements for axial elements of root
                        		if ( p[i - 1][0] > ip[i][0] ) 
					{
                            			wnu[i][0] = ws_[i][0] * vu[i][0];
                            			cpu[i][0] = 0.01 * vu[i][0] / dt;
                        		}
                        		else 
					{
//content of upper half
                            			wnu[i][0] = ws_[i][0] * vu[i][0] - 
							( ( ip[i][0] - p[i - 1][0] ) * 0.5 * cpp[i][0] ); 
                            			cpu[i][0] = 0.5 * cpp[i][0] / dt; //dw/dp *1/dt for jacobian
                        		}

                        		if ( p[i][0] > ip[i][0] ) 
					{
                            			wnl[i][0] = ws_[i][0] * vl[i][0];
//strictly][this should be zero][but set to arbitrary small #?
                            			cpl[i][0] = 0.01 * vl[i][0] / dt; 
                        		}
                        		else 
					{
                            			wnl[i][0] = ws_[i][0] * vl[i][0] - 
							( ( ip[i][0] - p[i][0] ) * 0.5 * cpp[i][0] );
// content][lower half of plant element
//dw/dp *1/dt for jacobian
                            			cpl[i][0] = 0.5 * cpp[i][0] / dt; 
                        		}
                    		}

                    		for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{ //sets upper and lower elements for lateral elements of root
                        		if ( p[i][0] > ip[i][1] ) 
					{ //remember: pressures are nodal][other values are elements!!
                            			wnu[i][1] = ws_[i][1] * vu[i][1];
                            			cpu[i][1] = 0.01 * vu[i][1] / dt;
                        		}
                        		else 
					{
//content of upper half
                            			wnu[i][1] = ws_[i][1] * vu[i][1] - 
							( ( ip[i][1] - p[i][0] ) * 0.5 * cpp[i][1] ); 
//dw/dp *1/dt for jacobian
                            			cpu[i][1] = 0.5 * cpp[i][1] / dt; 
                        		}

                        		if ( p[i][1] > ip[i][1] ) 
					{
                            			wnl[i][1] = ws_[i][1] * vl[i][1];
//strictly][this should be zero][but set to arbitrary small #?
                            			cpl[i][1] = 0.01 * vl[i][1] / dt; 
                        		}
                        		else 
					{
                            			wnl[i][1] = ws_[i][1] * vl[i][1] - 
							( ( ip[i][1] - p[i][1] ) * 0.5 * cpp[i][1] );
// content][lower half of plant element
//dw/dp *1/dt for jacobian
                            			cpl[i][1] = 0.5 * cpp[i][1] / dt; 
                        		}
                    		}

                    		for (i = 1; i <= ShootModules; i++ ) 
				{ //get upper and lower k's for axial shoot elements
                        		base1 = -p[i + 1][0] / b[i][0];
                        		base2 = -p[i][0] / b[i][0];
                        		exponent = cc[i][0];

                        		if( !verifyPowCondition(base1,exponent) || 
								!verifyPowCondition(base2,exponent) ) 
					{
                            			error = true;
                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
			    			//cout << "failure 1 " << endl;
		//			    	nodetype = "axial shoot";
                            			//cout << "from 4d\t" << i <<"\n";
                            			goto innerWhile1;
                        		}

                        		result1 = pow( base1, exponent);
                        		result2 = pow( base2, exponent);

                        		if ( !verifyExpCondition(result1) || 
								!verifyExpCondition(result2) ) 
					{
                            			//HydraulicModelFailCond = 1;
			    			nodefailure = i;
			    			//cout << "failure 2 " << endl;
//			    			nodetype = "axial shoot";
                            			error = true;
                            			//cout << "from 4c\n";
                            			goto innerWhile1;
                        		}

//weibull function for lower half of element i
                        		kl[i][0] = ks[i][0] * exp ( -result1 );  
//weibull function for upper half of element i
                        		ku[i][0] = ks[i][0] * exp ( -result2 );  

//assigns a value for incomplete gamma function
                        		aa = 1.0 / cc[i][0];    
//n value for upper half in preparation for pgamma function
                        		nn = result2;	      

                        		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
					{
                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
//			    			nodetype = "lateral shoot";
			    			//cout << "failure 3" << endl;
                            			error = true;
                            			//cout << "from 4b\n";
                            			goto innerWhile1;
                        		}

//phi for upper half
                        		phiu = ( 1.0 - gammp ) * gamma * ks[i][0] * b[i][0] / cc[i][0]; 
                        		nn = result1; //n value for lower half

//calls pgamma again with lower half n value
                        		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
					{
                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
//			    			nodetype = "lateral shoot";
			    			//cout << "failure 4" << endl;
                            			error = true;
                            			//cout << "from 4a\n";
                            			goto innerWhile1;
                        		}

//phi for lower half
                        		phil = ( 1.0 - gammp ) * gamma * ks[i][0] * b[i][0] / cc[i][0]; 
//flux of element i--CF converts conductivity to conductance
                        		jl[i][0] = cf[i][0] * ( phil - phiu ); 
                    		}

                    		for (  i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{ //get upper and lower k's for axial rhoot elements
                        		base1 = -p[i][0] / b[i][0];
                        		base2 = -p[i - 1][0] / b[i][0];
                        		exponent = cc[i][0];

                        		if( !verifyPowCondition(base1,exponent) || 
							!verifyPowCondition(base2,exponent) ) 
					{
                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
//			    			nodetype = "axial root";
			    			//cout << "failure 5" << endl;
                            			error = true;
                            			//cout << "from 3d\n";
                            			goto innerWhile1;
                        		}

                        		result1 = pow(base1, exponent);
                        		result2 = pow(base2, exponent);

                        		if ( !verifyExpCondition(result1) || 
							!verifyExpCondition(result2) ) 
					{
                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
//			    			nodetype = "axial root";
			    			//cout << "failure 6" << endl;
                            			error = true;
                            			//cout << "from 3c\n";
                            			goto innerWhile1;
                        		}

//weibull function for lower half of element i
                        		kl[i][0] = ks[i][0] * exp( -result1 ); 
//weibull function for upper half of element i
                        		ku[i][0] = ks[i][0] * exp( -result2 ); 

//n value for upper half in preparation for pgamma function
                        		nn = result2; 
//assigns a value for incomplete gamma function
                        		aa = 1.0 / cc[i][0]; 

                        		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
					{

                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
//			    			nodetype = "lateral shoot";
			    			//cout << "failure 7" << endl;
                            			error = true;
                            			//cout << "from 3b\n";
                            			goto innerWhile1;
                        		}

//phi for upper half
                        		phiu = ( 1.0 - gammp ) * gamma * ks[i][0] * b[i][0] / cc[i][0]; 
                        		nn = result1; //n value for lower half

                        		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
					{
                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
//			    			nodetype = "lateral shoot";
			    			//cout << "failure 8" << endl;
                            			error = true;
                            			//cout << "from 3a\n";
                            			goto innerWhile1;
                        		}

//phi for lower half
                        		phil = ( 1.0 - gammp ) * gamma * ks[i][0] * b[i][0] / cc[i][0]; 
//flux of element i--CF converts conductivity to conductance
                        		jl[i][0] = cf[i][0] * ( phil - phiu ); 
                    		}

                    		for (i = 1; i <= ShootModules; i++ ) 
				{ //get upper and lower k's for lateral shoot elements
                        		base1 = -p[i][0] / b[i][1];
                        		base2 = -p[i][1] / b[i][1];
                        		exponent = cc[i][1];

                        		if( !verifyPowCondition(base1,exponent) || 
								!verifyPowCondition(base2,exponent) ) 
					{
                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
//			    			nodetype = "lateral shoot";
			    			//cout << "failure 9" << endl;
                            			error = true;
                            			//cout << "from 2d\n";
                            			goto innerWhile1;
                        		}

                        		result1 = pow(base1,exponent);
                        		result2 = pow(base2,exponent);

                        		if ( !verifyExpCondition(result1) || 
								!verifyExpCondition(result2) ) 
					{
                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
//			    			nodetype = "lateral shoot";
			    			//cout << "failure 10" << endl;
                            			error = true;
                            			//cout << "from 2c\n";
                            			goto innerWhile1;
                        		}

//weibull function for lower half of element i
                        		kl[i][1] = ks[i][1] * exp( -result1 ); 
//weibull function for upper half of element i
                        		ku[i][1] = ks[i][1] * exp( -result2 ); 
//n value for upper half in preparation for pgamma function
                        		nn = result2; 
//assigns a value for incomplete gamma function
                        		aa = 1.0 / cc[i][1]; 

                        		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
					{
                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
//			    			nodetype = "lateral shoot";
			    			//cout << "failure 11" << endl;
                            			error = true;
                            			//cout << "from 2b\n";
                            			goto innerWhile1;
                        		}

//phi for upper half
                        		phiu = ( 1.0 - gammp ) * gamma * ks[i][1] * b[i][1] / cc[i][1]; 
//n value for lower half
                        		nn = result1; 
//calls pgamma again with lower half n value
                        		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
					{
                            			nodefailure = i;
			    			//HydraulicModelFailCond = 1;
//			    			nodetype = "lateral shoot";
			    			//cout << "failure 12" << endl;
                            			error = true;
                           			//cout << "from 2a\n";
                            			goto innerWhile1;
                        		}

//phi for lower half
                        		phil = ( 1.0 - gammp ) * gamma * ks[i][1] * b[i][1] / cc[i][1]; 
//flux of element i--CF converts conductivity to conductance
                        		jl[i][1] = cf[i][1] * ( phil - phiu ); 
                    		}

                        	for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
				{ //get upper and lower k's for lateral root elements
                            		base1 = -p[i][1] / b[i][1];
                            		base2 = -p[i][0] / b[i][1];
                            		exponent = cc[i][1];

                            		if( !verifyPowCondition(base1,exponent) || 
								!verifyPowCondition(base2,exponent) ) 
					{
                                		nodefailure = i;
						//HydraulicModelFailCond = 1;
//						nodetype = "lateral shoot";
						//cout << "failure 13" << endl;
                                		//cout << "from 1d\n";
                                		goto innerWhile1;
                            		}

                            		result1 = pow( base1, exponent );
                            		result2 = pow( base2, exponent );

                            		if ( !verifyExpCondition(result1) || 
								!verifyExpCondition(result2) ) 
					{
                                		nodefailure = i;
						//HydraulicModelFailCond = 1;
    //						nodetype = "lateral shoot";
						//cout << "failure 14" << endl;
                                		error = true;
                                		//cout << "from 1c\n";
                                		goto innerWhile1;
                            		}

//weibull function for lower half of element i
                            		kl[i][1] = ks[i][1] * exp( -result1 ); 
//weibull function for upper half of element i
                            		ku[i][1] = ks[i][1] * exp( -result2 ); 

//n value for upper half in preparation for pgamma function
                            		nn = result2; 
//assigns a value for incomplete gamma function
                            		aa = 1.0 / cc[i][1]; 

                            		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
					{
                                		nodefailure = i;
			        		//HydraulicModelFailCond = 1;
    //						nodetype = "lateral shoot";
						//cout << "failure 15" << endl;
                                		error = true;
                                		//cout << "from 1b\n";
                                		goto innerWhile1;
                            		}

//phi for upper half
                            		phiu = ( 1.0 - gammp ) * gamma * ks[i][1] * b[i][1] / cc[i][1]; 
//n value for lower half
                            		nn = result1; 

//calls pgamma again with lower half n value
                            		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
					{
                                		//HydraulicModelFailCond = 1;
						//cout << "failure 16" << endl;
						nodefailure = i;
    //						nodetype = "lateral shoot";
                                		error = true;
                                		//cout << "from 1a\n";
                                		goto innerWhile1;
                            		}

//phi for lower half
                            		phil = ( 1.0 - gammp ) * gamma * ks[i][1] * b[i][1] / cc[i][1]; 
//flux of element i--CF converts conductivity to conductance
                            		jl[i][1] = cf[i][1] * ( phil - phiu ); 
                        	}

                        	for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{ //sets upper and lower element properties for soil
                            		for (j = 2; j <= SoilRhizElements + 1; j++ ) 
					{
                                		if ( p[i][j] < pe[i][j] ) 
						{
                                    			wnl[i][j] = ws_[i][j] * vl[i][j] * 
								pow( ( pe[i][j] / p[i][j] ), b1[i][j] );
// content lower half of element [eqn 5.9][campbell]
//dw/dp*1/dt for jacobian][lower half of element
                                    			cpl[i][j] = -wnl[i][j] / (b[i][j] * p[i][j] * dt); 
//k[i] biased to lower
                                    			kl[i][j] = ks[i][j] * 
								pow( ( pe[i][j] / p[i][j] ), n[i][j] ); 
// node p [eqn. 6.14][campbell]
//a hybrid between phi[i] and phi[i+1][eqn. 8.19 of campbell]
                                    			phil = kl[i][j] * p[i][j] / n1[i][j]; 
                                		}
                                		else 
						{
                                    			wnl[i][j] = ws_[i][j] * vl[i][j];
//I think this is an arbitrary small number because
// dw/dp would otherwise be zero
                                    			cpl[i][j] = 0.01 * vl[i][j] / dt; 
                                    			kl[i][j] = ks[i][j];
//This is integral for p > pe
                                    			phil = ks[i][j] * ( pe[i][j] * 
									n[i][j] / n1[i][j] + p[i][j] ); 
                                		}

                                		if ( p[i][j - 1] < pe[i][j] ) 
						{
                                    			wnu[i][j] = ws_[i][j] * vu[i][j] * 
								pow( ( pe[i][j] / p[i][j - 1]), b1[i][j]);
// content of upper half of element
//dw/dp*1/dt 'dw/dp*1/dt for
                                    			cpu[i][j] = -wnu[i][j] / ( b[i][j] * 
								p[i][j - 1] * dt ); 
// jacobian][upper half of element
//k[i]][no bias
                                    			ku[i][j] = ks[i][j] * pow( ( pe[i][j] / 
								p[i][j - 1] ), n[i][j] ); 
//phi[i]][no bias
                                    			phiu = ku[i][j] * p[i][j - 1] / n1[i][j]; 
                                		}
                                		else 
						{
                                    			wnu[i][j] = ws_[i][j] * vu[i][j];
//arbitrary small number to get model going at p above pe
                                    			cpu[i][j] = 0.01 * vu[i][j] / dt; 
                                    			ku[i][j] = ks[i][j];
//integral for p>pe
                                    			phiu = ku[i][j] * ( pe[i][j] * n[i][j] / 
								n1[i][j] + p[i][j - 1] ); 
                                		}
//flux of element i--CF converts conductivity to conductance
                                		jl[i][j] = cf[i][j] * ( phil - phiu ); 
                            		}
                        	}

                        	for (  i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{
                            		rflux[ii][i] = jl[i][2] / ar[i];
                        	}

//tie up loose ends of network for mass balance
                        	jl[ totmodules + 2][0] = 0.0;
                        	jl[0][0] = 0.0;
                        	wnu[ totmodules + 2][0] = 0.0;
                        	wu[ totmodules + 2][0] = 0.0;
                        	wnl[0][0] = 0.0;
                        	wl[0][0] = 0.0;

                        	for (  i = 1; i <= ShootModules; i++ ) 
				{
                            		wnl[i][2] = 0.0;
                            		wl[i][2] = 0.0;
                        	}

                        	for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{
                            		wnu[i][SoilRhizElements + 2] = 0.0;
                            		wu[i][SoilRhizElements + 2] = 0.0;
                            		jl[i][SoilRhizElements + 2] = 0.0;
                        	}

                        	cpl[0][0] = 0.0;
                        	cf[0][0] = 0.0;
                        	cpu[ totmodules + 2][0] = 0.0;
                        	cf[ totmodules + 2][0] = 0.0;

                        	for (i = 1; i <= ShootModules; i++ ) 
				{ //mass balance for shoot lateral nodes
                            		f[i][1] = ( wnu[i][1] - wu[i][1] + wnl[i][2] - 
							wl[i][2] ) / dt + jl[i][2] - jl[i][1];
                        	}

                        	for (i = 1; i <= (ShootModules + 1); i++ ) 
				{ //mass balance for shoot axial nodes
                            		if ( i == ShootModules + 1 ) 
					{
                                		f[i][0] = ( wnu[i + 1][0] - wu[i + 1][0] + 
								wnl[i - 1][0] - wl[i - 1][0] ) / 
								dt + jl[i - 1][0] - jl[i + 1][0];
                            		}
                            		else 
					{
                                		f[i][0] = ( wnu[i][0] - wu[i][0] + wnl[i][1] - 
								wl[i][1] + wnl[i - 1][0] - wl[i - 1][0]) /
								 dt + jl[i][1] + jl[i - 1][0] - jl[i][0];
                            		}
                        	}

                        	for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{ //mass balance for axial root nodes
                            		f[i][0] = ( wnu[i][1] - wu[i][1] + wnu[i + 1][0] - 
								wu[i + 1][0] + wnl[i][0] - wl[i][0] ) / 
								dt + jl[i][0] - jl[i][1] - jl[i + 1][0];
                        	}

                        	for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{ //mass balance for lateral root and soil nodes
                            		for (j = 1; j <= SoilRhizElements + 1; j++ ) 
					{
                                		f[i][j] = ( wnu[i][j + 1] - wu[i][j + 1] + 
								wnl[i][j] - wl[i][j] ) / dt + 
								jl[i][j] - jl[i][j + 1];
                                		if ( j == SoilRhizElements + 1 ) 
						{
                                        		f[i][j] = 0.0;
                                		}
                            		}
                        	}

                        	//next sum the se's
                        	se = 0.0;
                        	for (i = 1; i <= (totmodules + 1); i++ ) 
				{
                            		se += f[i][0];
                        	}
                        	for (i = 1; i <= ShootModules; i++ ) 
				{
                            		se += f[i][1];
                        	}
                        	for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{
                            		for (j = 1; j <= SoilRhizElements + 1; j++ ) 
					{
                                    		se += f[i][j];
                            		}
                        	}
                        	tnode = ShootModules + totmodules + 1 + ( SoilRhizElements + 1 ) * 
								RootModules; //total# of NODES
                        	k = 1;

                        	for (i = 1; i <= (ShootModules + 1); i++ ) 
				{ //reset f's to ff's [for gauss] in order of equations in matrix
                            		if ( i == ShootModules + 1 ) 
					{
                                		ff[k] = f[i][0];
                                		k++;
                            		}
                            		else 
					{
                                		for (j = 0; j <= 1; j++ ) 
						{
                                    			ff[k] = f[i][j];
                                    			k++;
                                		}
                            		}
                        	}

                        	for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{
                            		for (  j = 0; j <= SoilRhizElements + 1; j++ ) 
					{
                                		ff[k] = f[i][j];
                                		k++;
                            		}
                        	}

                        	for (i = 1; i <= tnode; i++ ) 
				{ //zeros the jacobian matrix
                            		for (  j = 1; j <= tnode; j++ ) 
					{
                                		jmatrix[i][j] = 0.0;
                            		}
                        	}

                        	//calculate and load shoot lateral nodes
                        	i = 2;
                        	j = 2;
                        	for (k = 1; k <= ShootModules; k++ ) 
				{
                            		jmatrix[i][j] = ku[k][1] * cf[k][1] + cpu[k][1]; //dF[k1]/dY[k1]
                            		jmatrix[i][j - 1] = -kl[k][1] * cf[k][1]; //dF[k,1]/dY[k,0]
                            		i = i + 2;
                            		j = j + 2;
                        	}

                        	//load shoot axial nodes
                        	i = 1;
                        	j = 1;
                        	for ( k = 1; k <= ShootModules; k++ ) 
				{
//dF[k0]/dY[k0]
                            		jmatrix[i][j] = ku[k][0] * cf[k][0] + kl[k][1] * cf[k][1] + 
								kl[k - 1][0] * cf[k - 1][0] + cpu[k][0] + 
								cpl[k - 1][0] + cpl[k][1]; 
//dF[k,0]/dY[k-1,0]
                            		jmatrix[i][j - 2] = -ku[k - 1][0] * cf[k - 1][0]; 
                            		jmatrix[i][j + 1] = -ku[k][1] * cf[k][1]; //dF[k,0]/dY[k,1]
                            		jmatrix[i][j + 2] = -kl[k][0] * cf[k][0]; //dF[k,0]/dY[k+1,0]
                            		i = i + 2;
                            		j = j + 2;
                        	}

                        	//load root collar node
                        	k = ShootModules + 1;
//dF[k,0]/dY[k,0]
                        	jmatrix[i][j] = kl[k - 1][0] * cf[k - 1][0] + ku[k + 1][0] * 
							cf[k + 1][0] + cpu[k + 1][0] + cpl[k - 1][0]; 
                        	jmatrix[i][j - 2] = -ku[k - 1][0] * cf[k - 1][0]; //dF[k,0]/dY[k-1,0]
                        	jmatrix[i][j + 1] = -kl[k + 1][0] * cf[k + 1][0]; //dF[k,0]/dY[k+1,0]

                        	//load axial root nodes
                        	i++;
                        	j++;

                        	for (k = ShootModules + 2; k <= (totmodules + 1); k++ ) 
				{
                            		jmatrix[i][j] = kl[k][0] * cf[k][0] + ku[k + 1][0] * 
							cf[k + 1][0] + ku[k][1] * cf[k][1] + cpl[k][0] + 
							cpu[k + 1][0] + cpu[k][1]; //dF[k,0]/dY[k,0]
                            		if ( k == ShootModules + 2 ) 
					{
                                		jmatrix[i][j - 1] = -ku[k][0] * cf[k][0];
                            		}
                            		else 
					{
                                		jmatrix[i][j - SoilRhizElements - 2] = -ku[k][0] * 
									cf[k][0];//dF[k,0]/dY[k-1,0]
                            		}
                            		jmatrix[i][j + 1] = -kl[k][1] * cf[k][1]; //dF[k,0]/dY[k,1]
                            		jmatrix[i][j + SoilRhizElements + 2] = -kl[k + 1][0] * 
									cf[k + 1][0]; //dF[k,0]/dY[k+1,0]
                            		i = i + SoilRhizElements + 2;
                            		j = j + SoilRhizElements + 2;
                        	}

                        	//load lateral root and soil nodes
                        	i = 2 * ShootModules + 3;
                        	j = 2 * ShootModules + 3;
                        	for (k = ShootModules + 2; k <= (totmodules + 1); k++ ) 
				{
                            		for (m = 1; m <= SoilRhizElements + 1; m++ ) 
					{
                                		jmatrix[i][j] = kl[k][m] * cf[k][m] + ku[k][m + 1] * 
									cf[k][m + 1] + cpl[k][m] + 
									cpu[k][m + 1]; //dF[k,m]/dY[k,m]
//dF[k,m]/dY[k,m-1]
                                		jmatrix[i][j - 1] = -ku[k][m] * cf[k][m]; 
//dF[k,m]/dY[k,m+1]
                                		jmatrix[i][j + 1] = -kl[k][m + 1] * cf[k][m + 1]; 
                                		i++;
                                		j++;
                            		}
                            		i++;
                            		j++;
                        	}

                        	gauss(ShootModules, SoilRhizElements, tnode, jmatrix, col, ff, 
										indxx, subtract, dp);

                        	k = 1;
                        	for (i = 1; i <= ShootModules; i++ ) 
				{
                            		dpp[i][0] = dp[k];
                            		dpp[i][1] = dp[k + 1];
                            		k = k + 2;
                        	}

                        	dpp[ ShootModules + 1][0] = dp[k];
                        	k++;

                        	for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{
                            		dpp[i][0] = dp[k];
                            		for (m = 1; m <= (SoilRhizElements + 1); m++ ) 
					{
                                		k++;
                                		dpp[i][m] = dp[k];
                            		}
                            		k++;
                        	}

                        	//assign new p's
                        	for (i = 1; i <= (totmodules + 1); i++ ) 
				{
                            		lim = 0.8 * fabs( p[i][0] );
                            		abv = fabs( dpp[i][0] );

                            		if ( ( abv > lim ) && ( p[i][0] < -0.1 ) ) 
					{
                                		dpp[i][0] *= lim / abv;
                            		}
                            		p[i][0] -= dpp[i][0];
                        	}

                        	for (i = 1; i <= ShootModules; i++ ) 
				{
                            		lim = 0.8 * fabs( p[i][1] );
                            		abv = fabs( dpp[i][1] );
                            		if ( ( abv > lim ) && ( p[i][1] < -0.1 ) ) 
					{
                                		dpp[i][1] *= lim / abv;
                            		}
                            		p[i][1] -= dpp[i][1];
                        	}

                        	for (i = ShootModules + 2; i <= (totmodules + 1); i++ ) 
				{
                            		for ( j = 1; j <= SoilRhizElements + 1; j++ ) 
					{
                                		lim = 0.8 * fabs( p[i][j] );
                                		abv = fabs( dpp[i][j] );
                                		if ( ( abv > lim ) && ( p[i][j] < -0.1 ) ) 
						{
                                    			dpp[i][j] *= lim / abv;
                                		}
                                		p[i][j] -= dpp[i][j];
                            		}
                        	}
                        	nits = nits + 1;
                   	} while ( ( se >= IM ) && ( nits <= MAXITS ) ); //convergence or divergence

    			innerWhile1: cout <<"";
	

                    	if ( error == true || nits>MAXITS )
			{
                      		//cout << "exiting innerWhile2\t"<<nits<<"\n";
                        	goto innerWhile2;
                    	}

                    	//add up all pressures for steady-state test
                    	sumpress = 0.0;
                    	for ( i = 1; i <= totmodules + 1; i++) 
			{
                        	sumpress += fabs(p[i][0] );
                    	}

                    	for ( i = 1; i <= ShootModules; i++) 
			{
                        	sumpress += fabs(p[i][1]);
                    	}

                    	for ( i = ShootModules + 2; i <= totmodules + 1; i++) 
			{
                        	for ( j = 1; j <= SoilRhizElements + 1; j++) 
				{
                            		sumpress += fabs(p[i][j]);
                        	}
                    	}

                    	deltap = fabs(sumpress - oldpress);

//reset water contents of all ELEMENTS
                    	for ( i = 1; i <= totmodules + 1; i++) 
			{
                        	if ( i != ShootModules + 1 ) 
				{
                            		wl[i][0] = wnl[i][0];
                            		wu[i][0] = wnu[i][0];
                            		wl[i][1] = wnl[i][1];
                            		wu[i][1] = wnu[i][1];
                        	}
                    	}

                    	for (i = ShootModules + 2; i <= totmodules + 1; i++) 
			{
                        	for ( j = 2; j <= SoilRhizElements + 1; j++) 
				{
                            		wl[i][j] = wnl[i][j];
                            		wu[i][j] = wnu[i][j];
                        	}
                    	}
                    	h++;
            	} while ( (deltap >= 0.001) && (h <= MAXH) );

    		innerWhile2: cout <<"";
            	if ( (nits > MAXITS) || (h > MAXH) || error) 
		{
			//cout << "Exit outer\t"<<ssss<<"\t"<<h<<" \n";
                	goto outerWhile;
            	}

            	for ( i = 1; i <= totmodules + 1 ; i++)  
		{    //ii is the e-increment index
                	psi[ii][i][0] = p[i][0];
                	psi[ii][i][1] = p[i][1];
            	}

//This is when the final pressures for the rhizosphere nodes are outputted.  DER 6/25/08
            	for ( i = ShootModules+2; i <= totmodules+1; i++) 
		{
                	for ( j = 2; j <= SoilRhizElements+1; j++) 
			{
                    		psi[ii][i][j] = p[i][j];
                	}
            	}

            	//calculate the predicted e to match ptarg
            	//first calculate the average p in canopy
            	Sum = 0.0;
            	for ( i = 1; i <= ShootModules; i++) 
		{
                	Sum += psi[ii][i][1];
            	}

            	pleafave = Sum / ShootModules;
	    	//cout << "pleafave " << pleafave << endl;
	    	//cout << "ptarg " << ptarg << endl;
            	//if (pleafave == ptarg) {
//WHY WAS THIS SET TO e == evap[i]?
            	if (pleafave > ptarg) 
		{
//	    	if(e == evap[ii])
	    	//{
                	epred = evap[ii];
		
                	if (e == 0.0) 
			{    //if no e measurements][get pressures from p data
                    		for ( i = ShootModules; i <= totmodules+1; i++) 
				{
                        		psinode[i][0] = psi[ii][i][0];
                        		kkk = ii;
                        		if ( i != ShootModules+1) 
					{
                            			psinode[i][1] = psi[ii][i][1];
                        		}
                    		}
                	}
            	}
	    	//if(hydraulicModel_flag == 0)
	    	//{
	    	//	leafpsi = pleafave;
	    	//	hydraulicModel_flag = 1;
	    	//}
	    
	    	if (e == evap[ii]) 
	    	{
                	ppred = pleafave;
			leafpsi = pleafave;
                	for ( i = ShootModules; i <= totmodules + 1; i++) 
			{  //get pressures from p data
                    		psinode[i][0] = psi[ii][i][0];
                    		kkk = ii;
                    		if (i != ShootModules+1 ) 
				{
                        		psinode[i][1] = psi[ii][i][1];
                    		}
                	}
            	}
        	evap[ii + 1] = evap[ii] + treesParams.E_inc;
        	ii++;
	    
        } while (true && ii < NREC-2);  
//this loop continues until there is failure - note when ii==NREG  segmentation violation

        outerWhile: cout <<"";
	error = false;
        ii--;
        check = 0;
        if ( (e == 0.0) && (ptarg==0.0) ) 
	{     //no predictions or pressures
            	check = 1;
	    	//cout << "no pred" << endl;
            	epred = 0.0;
            	ppred = 0.0;
            	for ( i = ShootModules; i <= totmodules + 1; i++) 
		{
                	psinode[i][0] = 0.0;
                	if (i != ShootModules + 1 ) 
			{
                    		psinode[i][1] = 0.0;
                	}
            	}
        }
        if( ii == -1) 
	{
            	nodeFail[t] = 0;
        }
        else 
	{
            	//go on to record ecrits][pcrits][and failure][if any
            	ecrit = evap[ii]; //sets ecrit
	    	Ecrit = ecrit;
            	ecritOut[t] = ecrit;
            	for ( i = 1; i <= ShootModules; i++) 
	    	{   //sets pcrit to the lowest canopy p
                	if (i == 1)
			{
                    		pcrit = psi[ii][i][1];
				if (isnan(pcrit))
				{
					 pcrit = psi[ii][i][1]= ptarg;
				}
                	}
			else
			{
                    		if (psi[ii][i][1] < psi[ii][i - 1][1] && !isnan(psi[ii][i][1])) 
		    		{
                            		pcrit = psi[ii][i][1];
                    		}
                	}
            	}
            	pc[t] = pcrit;
            	if (ptarg != 0.0) 
	    	{
                	epredOut[t] = epred;
                	if (ptarg < pcrit) 
	        	{
		    		cout << "fail " << endl;
		    		cout << "ptarg " << ptarg << endl;
		    		//cout << "pcrit " << pcrit << endl;
		    		HydraulicModelFailCond = 1; //notifies Trees of failure
                    		nodeFail[t] = nodefailure;
    //		    		nodeTyp[t] = nodetype;
                    		check = 1;
		    		cin.get();
                	}
            	}
            	if (e != 0.0) 
	    	{
                	ppredOut[t] = ppred;
                	if (e > ecrit) 
			{
		    		cout << "fail " << endl;
		    		cout << "e " << e << endl;
		    		cout << "ecrit " << ecrit << endl;
		    		HydraulicModelFailCond = 1; //notifies Trees of failure
                    		nodeFail[t] = nodefailure;
    //		    		nodeTyp[t] = nodetype;
                    		check = 1;
                    		cin.get();
			}
        	}
	
		PsiCrit = pcrit;

        	if(ii == -1)
		{
			HydraulicModelFailCond = 1;
		}
		else
		{
        		if (check == 0) 
			{    //if no failure and if data is available][get avg. pressures
	    			HydraulicModelFailCond = 0; //no failure continue.
            			klpred[t] = evap[kkk] / (soilpsiavg[t] - ppred);

            			for ( i = 1; i<= totmodules+1; i++ ) 
				{   //calculate average pressures
                			if (i != ShootModules + 1) 
					{
                    				psinode[i][1] = 0.5 * (psinode[i][1] + psinode[i][0]);
                    				if ( callHydraulicModel == 1 ) 
						{
                        				psimin[i][1] = psinode[i][1];
                    				}
                			}

                			psinode[i][0] = 0.5 * (psinode[i][0] + psinode[i + 1][0]);
                			if ( callHydraulicModel == 1 ) 
					{
                    				psimin[i][0] = psinode[i][0];
                			}
            			}

            			for ( i = 1; i <= totmodules + 1; i++ ) 
				{ //reset psimin if appropriate
                			if (psinode[i][0] < psimin[i][0]) 
					{
                    				psimin[i][0] = psinode[i][0];
                			}
                			if (i != ShootModules + 1 ) 
					{
                    				if (psinode[i][1] < psimin[i][1]) 
						{
                            				psimin[i][1] = psinode[i][1];
                    				}
                			}
            			}
        		}

        		j = 0;
		
        		for (  i = 1; i <= ShootModules; i++) 
			{
//output conductances at ecrit and at e or p
            			ecrit_k_psi.run[j+i] = t;
            			ecrit_k_psi.module[j+i] = i;
            			ecrit_k_psi.axialK[i+j] = ks[i][0] * exp( -pow( (-psi[ii][i][0] / 
										b[i][0]), cc[i][0]));
	    
            			ecrit_k_psi.latK[i+j] = ks[i][1] * exp( -pow( (-psi[ii][i][1] / 
										b[i][1]), cc[i][1]) );
            			ecrit_k_psi.rhizK[i+j] = psi[295][i][1];

            			k_p_e.run[j+i] = t;
            			k_p_e.module[j+i] = i;
            			if (check == 0) 
				{
                			k_p_e.axialK[i+j] = ks[i][0] * exp( -pow( (-psi[kkk][i][0] / 
										b[i][0]), cc[i][0]) );
                			k_p_e.latK[i+j] = ks[i][1] * exp( -pow( (-psi[kkk][i][1] / 
										b[i][1]), cc[i][1]) );
            			}
            			else 
				{
                			k_p_e.axialK[i+j] = 0.0;
                			k_p_e.latK[j + i] = 0.0;
            			}
        		}
        		for ( i = ShootModules + 2; i <= totmodules + 1; i++) 
			{
            			ecrit_k_psi.module[j+i-1] = i;
            			ecrit_k_psi.axialK[i+j-1] = ks[i][0] * exp( -pow( (-psi[ii][i][0] / 
										b[i][0]), cc[i][0]) );
            			ecrit_k_psi.latK[i+j-1] = ks[i][1] * exp( -pow( (-psi[ii][i][1] / 
										b[i][1]), cc[i][1]) );

            			k_p_e.module[j+i-1] = i;
            			if (check == 0) 
				{
                			k_p_e.axialK[i+j-1] = ks[i][0] * exp( -pow( (-psi[kkk][i][0] / 
										b[i][0]), cc[i][0]) );
                			k_p_e.latK[i+j-1] = ks[i][1] * exp( -pow( (-psi[kkk][i][1] / 
										b[i][1]), cc[i][1]) );
            			}
            			else 
				{
                			k_p_e.axialK[i+j-1] = 0.0;
                			k_p_e.latK[i+j-1] = 0.0;
            			}
        		}

        		for ( i = ShootModules + 2; i <= totmodules + 1; i++) 
			{
            			ecrit_k_psi.rhizK[i+j-1] = (ecrit * latot) / 
					(psi[ii][i][SoilRhizElements + 1] -psi[ii][i][1]); //rhizosphere k
            			ecrit_k_psi.rhizFlux[i+j-1] = rflux[ii][i];
	    
            			if (check == 0) 
				{
                			k_p_e.rhizFlux[i+j-1] = rflux[kkk][i];
                			k_p_e.rhizK[i+j-1] = (evap[kkk] * latot) / 
						(psi[kkk][i][SoilRhizElements + 1] -psi[kkk][i][1]);
            			}
            			else 
				{
                			k_p_e.rhizK[i+j-1] = 0.0;
                			k_p_e.rhizFlux[i+j-1] = 0.0;
            			}
        		}

    		}
//Final whole-plant hydraulic conductance on a per unit leaf basis
//Note: klpred is mathematically correct for whole plant K per unit leaf area
//      However, klpred will be negative and appears to have a bug that causes a numerical problem
//      Strictly, we don't need the water balance model at night, but in practice we always call it
//      Suggestion is to use root and shoot K components at night
// Note: This assumes you have 1 shoot and 3 root modules defined in your model
// DSM July 2012
/* OLD
    		bool isNegRhizFlux = false;
    		double rhizK1, rhizK2, rhizK3, alr1, alr2, alr3;
    		double shootK, root1K, root2K, root3K, rootK, latK;
    		double latKL, latKR1, latKR2, latKR3, axKL, axKR1, axKR2, axKR3;
    		latKL = k_p_e.latK[1];
    		latKR1 = k_p_e.latK[2];
    		latKR2 = k_p_e.latK[3];
    		latKR3 = k_p_e.latK[4];
    		axKL = k_p_e.axialK[1];
    		axKR1 = k_p_e.axialK[2];
    		axKR2 = k_p_e.axialK[3];
    		axKR3 = k_p_e.axialK[4];
    		rhizK1 = k_p_e.rhizK[2]/latot;
    		rhizK2 = k_p_e.rhizK[3]/latot;
    		rhizK3 = k_p_e.rhizK[4]/latot;
    		alr1 = ar[3]/ratot;
    		alr2 = ar[4]/ratot;
    		alr3 = ar[5]/ratot;
    		shootK = 1/(1/latKL + 1/axKL);
    		root3K = 1/(1/axKR3 + 1/latKR3);
    		root2K = 1/(1/(root3K+latKR2)+1/axKR2);
    		rootK = 1/(1/(root2K+latKR1)+1/axKR1);
*/

//NEW
    		Kl = ev / (soilpsiavg[t] - ppred); //output k at the leaf.
/*
    		Kl = 0.0;
    		if (soilpsiavg[t] > (ppred+0.01))
    			Kl = ev / (soilpsiavg[t] - ppred); //output k at the leaf.

    		if (Kl > ev/treesParams.e_at_saturated_kl*treesParams.saturated_kl_for_whole_plant)
			Kl = ev/treesParams.e_at_saturated_kl*treesParams.saturated_kl_for_whole_plant;
*/

/*
//determine if at least rhizosphere flux is negative (flux from root to soil)
    		for (i = 0; i < treesParams.rmodules; i++)
    		{
			if (k_p_e.rhizK[i+2] < 0.0)
			{
	    			isNegRhizFlux = true;
			}
    		}

//if there is negative rhizosphere flux
    		if (Kl == 0.0 && isNegRhizFlux)
    		{
    			latKL = k_p_e.latK[1];
    			axKL = k_p_e.axialK[1];
    			shootK = 1.0/(1.0/latKL + 1.0/axKL);
    			latKR1 = k_p_e.latK[treesParams.rmodules+1];
    			axKR1 = k_p_e.axialK[treesParams.rmodules+1];
        		rhizK1 = k_p_e.rhizK[treesParams.rmodules+1];
			if (rhizK1 < 0.0)
			{
    				rootK = 1.0/(1.0/axKR1 + 1.0/latKR1);
			}
			else
			{
    				rootK = 1.0/(1.0/axKR1 + 1.0/latKR1 + 1.0/rhizK1);
			}
    			for (i = treesParams.rmodules-1; i > 0; i--)
    			{
            			latKR2 = k_p_e.latK[i+1];
            			axKR2 = k_p_e.axialK[i+1];
            			rhizK2 = k_p_e.rhizK[i+1];
	    			if (rhizK2 < 0.0)
				{
    	    				root2K = 1.0/(1.0/(rootK+latKR2)+1.0/axKR2);
				}
	    			else
				{
    	    				root2K = 1.0/(1.0/(rootK+latKR2)+1.0/axKR2+1.0/rhizK2);
				}
            			rootK = root2K;
    			}
			Kl = 1.0/(1.0/rootK+1.0/shootK);
    		}
    		else
    		{
			Kl = klpred[t];
    		}
*/
	


/*
//    		if (klpred[t] <= 0)
    		if (rhizK1 < 0 || rhizK2 < 0 || rhizK3 < 0)
    		{
    			Kl = 1/(1/rootK+1/shootK);
			//Kl = 1/(1/latKL+1/axKL+1/axKR1+1/axKR2+1/axKR3+alr1/latKR1+alr2/latKR2+alr3/latKR3+max(0,alr1/rhizK1)+max(0,alr2/rhizK2)+max(0,alr3/rhizK3));
    		}
    		else
    		{
			Kl = klpred[t];
    		}

    		//Kl = klpred[t];
    		if (ev < treesParams.e_at_saturated_kl && Kl > treesParams.saturated_kl_for_whole_plant)
    		{
			//Kl = treesParams.saturated_kl_for_whole_plant*ev/treesParams.e_at_saturated_kl;
			Kl = treesParams.saturated_kl_for_whole_plant;
    		}
*/


    		//Kl = 1/(latot/(k_p_e.latK[1])+latot/(k_p_e.axialK[1])+latot/(k_p_e.latK[2])+latot/(k_p_e.axialK[2]));
   		// Kl = latot/(latot/(k_p_e.latK[1])+latot/(k_p_e.axialK[1])+latot/(k_p_e.latK[2])+latot/(k_p_e.axialK[2]));


    		final_pe = pe[totmodules][2] * 1000.0;	//soil air entry pressure, kPa
    		final_ks = ks[totmodules][2] * 0.001;	//soil sat. k, mol/s m
    		final_b = b[totmodules][2];			//soil b value
    		final_ws = ws_[totmodules][2] * 0.001;	//soil sat. content, mol/m3
    		final_soilvol = soilvol;		        //total volume in m3
    		final_rd = rd;    
    	}
}


void HydraulicModel::print( double final_pe, double final_ks, double final_b,
                double final_ws, double final_soilvol, double final_rd,
                outputStruct k_p_e, outputStruct ecrit_k_psi,
                double epredOut[], double ecritOut[],
                double ppredOut[], double saturatedKs[][4], double klpred[],
                double md_raw[], double la_raw[], double arAl[], double e_raw[],
                double pc[] ) 
{         //outputs data to input sheet at end of simulation run


        ofstream outData("final.out");
        outData << "Ye,kPa " << final_pe << "\nks,mol/s_MPa_m " << final_ks << "\nb " << final_b
                        << "\nWs,moles/m3 " << final_ws << "\ntotal_volume,m3 " << final_soilvol << 
                        "\nroot_density_(m/m3) "<< final_rd << endl;
        outData.close();

        outData.open("k_at_p_e.out");
        outData << "axial_k's lat_k's rhiz_k rhiz_flux\n";
        for( int i = 1; i < 209; i++ )
	{
                outData << k_p_e.axialK[i] << " " << k_p_e.latK[i] << " " << k_p_e.rhizK[i] << 
						" " << k_p_e.rhizFlux[i]
                 << endl;
        }
        outData.close();

        outData.open("ecrit.out");
        outData << "axial_k's lat_k's P's,rhiz_k rhiz_flux\n";
        for( int i = 1; i < 209; i++ )
	{
            	outData << ecrit_k_psi.axialK[i] << " " << ecrit_k_psi.latK[i] << " " << 
					ecrit_k_psi.rhizK[i] << " " << ecrit_k_psi.rhizFlux[i] << endl;
        }
        outData.close();

        outData.open("out_data.out");
        outData << "MD klpred Al Ar:Al E E_pred Ecrit ppred pc\n";
        for( int i = 1; i < 5; i++ ) 
	{
	  	outData << md_raw[i] <<" "<< klpred[i] <<" "<< la_raw[i] << " " << arAl[i] << 
					" " << e_raw[i] << " " << epredOut[i] << " " << 
					ecritOut[i] << " " << ppredOut[i] << " " << pc[i] << endl;
        }
        outData.close();

        outData.open("saturated_ks.out");
        outData << "axial_k's lat_k's rhiz_k\n";
        for( int i=0; i<ULAT; i++)
	{
            	outData << saturatedKs[i][1] <<" "<< saturatedKs[i][2] <<" "<< saturatedKs[i][3] << endl;
        }
        outData.close();
}


//preSetup - among other things this is where we feed the soil PSI values in from TREES
//This is called for every time step
void HydraulicModel::preSetup( double md_raw[], double Md, double e_raw[], 
				double ev, double la_raw[], double Al,
                      		double arAl[], double Ar_Al, double **l_raw, double Ypd[])
{

    	int i = 1;
    	int j = 1;

	while (j < 20)
	{
        	l_raw[i][j] = Ypd[j-1];
		j++;
	}

        md_raw[i] = Md;
        la_raw[i] = Al;
        arAl[i] = Ar_Al;
        e_raw[i] = ev;
}

//Setup is called at the beginning of simulation only
void HydraulicModel::setup( 	bool reset, outputStruct &ecrit_k_psi, outputStruct &k_p_e, 
				double ecritOut[], double pc[], double klpred[],
		    		double ppredOut[], double epredOut[], double nodeFail[], 
				char nodeTyp[], double ***psi, double **rflux, 
				double soilpsiavg[], double nts[], double evap[], double kroot[],
		    		double axr[], double latr[], double kshoot[], 
				double dslat[], double dsax[], double drlat[],
		    		double drax[], double l[], double ksat[][MD], 
				double bsat[][MD], double ccsat[][MD],
 		    		double newb[][MD], double newc[][MD], double ip[][ULAT], 
				double b[][ULAT], double b1[][ULAT], 
		    		double n[][ULAT], double n1[][ULAT], double r[][ULAT], 
				double vu[][ULAT], double vl[][ULAT],
             	    		double cf[][ULAT], double pe[][ULAT], 
				double ws_[][ULAT], double ks[][ULAT], double p[][ULAT],
		    		double dpp[][ULAT], double jl[][ULAT], 
				double wnu[][ULAT], double wu[][ULAT], double wnl[][ULAT],
		    		double wl[][ULAT], double cpu[][ULAT], double cpl[][ULAT], 
				double cpp[][ULAT], double ku[][ULAT],
		    		double kl[][ULAT], double f[][ULAT], double cc[][2], 
				double psinode[][NMAX], double psimin[][NMAX],
             	    		double jmatrix[][NMAX], double jmatrix2[][NMAX], 
				double percent[], double rssoil[], double rs[], 
		    		double dp[], double ff[], int col[], int row[], int indxx[], 
				double subtract[], double pressure[],
		    		double plc[], double plcweib[], double rsquare[], 
				double soilpsi[], double al[], double ar[],
		    		double saturatedKs[][4], double &ptarg, double &e, 
				double &rr, double &rzw, double &rd, double &einc,
		    		int &SoilRhizElements, double &dt, double &gmd, double &gsd, 
				double &bkd, double &fs, double &fc,	
		    		double &Fabs, double &cpplant, double &saxlat, double &raxlat, 
				double &rfract, double &kla,
		    		double &aral, double &latot, double &ratot, 
				int &shootElements, int &ktotal, double &soilpsimin,
             	    		int &ShootModules, int &rootElements, int &RootModules, 
				int &totmodules, double &axShoot_b,
              	    		double &axShoot_c, double &latShoot_b, double &latShoot_c, 
				double &axRoot_b, double &axRoot_c,
		    		double &latRoot_b, double &latRoot_c, 
				double &plco, double &pdecrement, double &sumy1,
		    		double &sumy2, double &sumprod, double &cincrement, 
				double &soilvol, double &rlateral,
             	    		double &rLat_base, double &slateral, 
				double &sLat_base, double &raxial, double &pleafave,
		    		int &tnode, double fac, double &modelledKL, 
				int &HydraulicModelFailCond, trees_params treesParams)
{

//assign the params out of the hydraulicModel_param struct to variables
        double midday_at_sat_kl = treesParams.midday_at_sat_kl;
        double e_at_saturated_kl = treesParams.e_at_saturated_kl;
        double rhizosphere_width = treesParams.rhizosphere_width;
        int soilshells = treesParams.soilshells;
        int time_steps = 1800.0; //seconds
        double GMP = treesParams.GMP;
        double GSD = treesParams.GSD;
        double BD = treesParams.BD;
        double silt_fraction = treesParams.silt_fraction;
        double clay_fraction = treesParams.clay_fraction;
        double frac_absorbing_length = treesParams.frac_absorbing_length;
        double Capacitance = treesParams.Capacitance;
        double axK_latKl_shoot_modules = treesParams.axK_latKl_shoot_modules;
        double axKr_latKr_root_modules = treesParams.axKr_latKr_root_modules;
        double per_total_R_in_root_system = treesParams.per_total_R_in_root_system;
        double saturated_kl_for_whole_plant = treesParams.saturated_kl_for_whole_plant;
        double aral_at_sat_kl = treesParams.aral_at_sat_kl;
        double lai_at_sat_kl = treesParams.lai;
        double pd_at_sat_kl = treesParams.pd_at_sat_kl;
        axShoot_b = treesParams.ax_Shoot_b_value;
        axShoot_c = treesParams.ax_Shoot_c_value;
        latShoot_b = treesParams.lat_Shoot_b_value;
        latShoot_c = treesParams.lat_Shoot_c_value;
        axRoot_b = treesParams.ax_Root_b_value;
        axRoot_c = treesParams.ax_Root_c_value;
        latRoot_b = treesParams.lat_Root_b_value;
        latRoot_c = treesParams.lat_Root_c_value;
        double initial_conductivity_root = treesParams.initial_conductivity_root;
        double decrement_root = treesParams.decrement_root;
        double initial_conductivity_shoot = treesParams.initial_conductivity_shoot;
        double decrement_shoot = treesParams.decrement_shoot;
	double psimin_in[NMAX][NMAX];

	if (reset == false)
	{
    		for ( int j = 0; j < NMAX; j++) 
		{
          		for ( int y = 0; y < NMAX; y++) 
			{
				psimin_in[j][y] = psimin[j][y];
			}
		}
	}

	if (reset == true)
	{
    		for( int i = 0; i < 209; i++ )
		{
          		ecrit_k_psi.axialK[i] = 0.0;
          		ecrit_k_psi.latK[i] = 0.0;
          		ecrit_k_psi.module[i] = 0.0;
          		ecrit_k_psi.rhizFlux[i] = 0.0;
          		ecrit_k_psi.rhizK[i] = 0.0;
          		ecrit_k_psi.run[i] = 0.0;
          		k_p_e.axialK[i] = 0.0;
          		k_p_e.latK[i] = 0.0;
          		k_p_e.module[i] = 0.0;
          		k_p_e.rhizFlux[i] = 0.0;
          		k_p_e.rhizK[i] = 0.0;
          		k_p_e.run[i] = 0.0;
    		}

    		for (int j = 0; j < MD; j++ )
		{
          		ecritOut[j] = 0.0;
          		pc[j] = 0.0;
          		klpred[j] = 0.0;
          		ppredOut[j] = 0.0;
          		epredOut[j] = 0.0;
          		nodeTyp[j] = ' ';
          		nodeFail[j] = 0;
    		}

    		for (int z = 0; z < NREC; z++) 
		{
          		for (int j = 0; j < MD; j++) 
			{
                  		for (int y = 0; y < ULAT; y++) 
				{
                          		psi[z][j][y] = 0.0;
                  		}
          		}

          		for (int j = 0; j < NMAX; j++) 
			{
                  		rflux[z][j] = 0.0;
          		}
          		soilpsiavg[z] = 0.0; //NOW GETTING SOIL MOISTURE INITIALIZED FROM TREES
          		nts[z] = 0;
          		evap[z] = 0.0;
    		}

    		for ( int j = 0; j < MD; j++) 
		{
          		kroot[j]=0.0; 
	  		axr[j]=0.0; 
	  		latr[j]=0.0; 
	  		kshoot[j]=0.0; 
	  		dslat[j]=0.0; 
	  		dsax[j]=0.0; 
	  		drlat[j]=0.0; 
	  		drax[j]=0.0; 
	  		l[j]=0.0;
          		for ( int y = 0; y < MD; y++) 
			{
                  		ksat[j][y]=0.0; 
		  		bsat[j][y]=0.0; 
		  		ccsat[j][y]=0.0; 
		  		newb[j][y]=0.0; 
		  		newc[j][y]=0.0;
          		}

          		for ( int y = 0; y < ULAT; y++) 
			{
                		ip[j][y]=0.0; 
				b[j][y]=0.0; 
				b1[j][y]=0.0; 
				n[j][y]=0.0; 
				n1[j][y]=0.0; 
				r[j][y]=0.0; 
				vu[j][y]=0.0; 
				vl[j][y]=0.0;
                		cf[j][y]=0.0; 
				pe[j][y]=0.0; 
				ws_[j][y]=0.0; 
				ks[j][y]=0.0; 
				p[j][y]=0.0; 
				dpp[j][y]=0.0;
                		jl[j][y]=0.0; 
				wnu[j][y]=0.0; 
				wu[j][y]=0.0; 
				wnl[j][y]=0.0;
                		wl[j][y]=0.0; 
				cpu[j][y]=0.0; 
				cpl[j][y]=0.0;
                		cpp[j][y]=0.0; 
				ku[j][y]=0.0; 
				kl[j][y]=0.0;
                		f[j][y] = 0.0;
          		}
          		for ( int y = 0; y < 2; y++) 
			{
                  		cc[j][y] = 0.0;
          		}
    		}

    		for ( int j = 0; j < NMAX; j++) 
		{
          		for ( int y = 0; y < NMAX; y++) 
			{
                		psinode[j][y] = 0.0; 
				psimin[j][y] = 0.0;
                		jmatrix2[j][y] = 0.0;
                		jmatrix[j][y] = 0.0;
          		}
          		percent[j] = 0.0; 
	  		rssoil[j] = 0.0; 
	  		rs[j] = 0.0; 
	  		dp[j] = 0.0; 
	  		ff[j] = 0.0; 
	  		col[j] = 0; 
	  		row[j] = 0; 
	  		indxx[j] = 0; 
	  		subtract[j] = 0.0;
    		}

    		for ( int j = 0; j < NKINC; j++)
		{
          		pressure[j] = 0.0; 
 	  		plc[j] = 0.0; 
	  		plcweib[j] = 0.0; 
	  		rsquare[j] = 0.0;
    		}

    		for ( int j = 0; j < ULAT; j++) 
		{
          		soilpsi[j] = 0.0; 
	  		al[j] = 0.0; 
	  		ar[j] = 0.0;
    		}

    		for ( int j = 0; j < MD; j++)
		{
          		for (int y = 0; y < 4; y++)
			{
                  		saturatedKs[j][y] = 0.0;
          		}
    		}

	}//end if reset

    	int t = 1;
      	ptarg = midday_at_sat_kl;

      	e = e_at_saturated_kl;
//This modification allows for reasonable scaling of maximum hydraulic conductance
//  for LAI values from less than one to greater than one
//  It slightly exaggerates Kmax for very tiny plants with low LAI
//  It scales Kmax using whole plant saturation transpiration
	e *= max(sqrt(treesParams.lai+0.00000001),treesParams.lai);
	//e *= treesParams.lai;

    	rr = 0.0001;                           //root radius in m
      	rzw = rhizosphere_width;

    	rzw = rzw * 0.001;                    //rhizosphere width, converted to m
    	rd = pow( ( 1.0 / rzw ), 2.0 ) / M_PI;      //root length density

    	SoilRhizElements = soilshells;
    	dt = (double) time_steps;
    	gmd = GMP;
    	gsd = GSD;
    	bkd = BD;
    	fs = silt_fraction;
    	fc = clay_fraction;

    	Fabs = frac_absorbing_length;
    	cpplant = Capacitance;
    	saxlat = axK_latKl_shoot_modules;
    	raxlat = axKr_latKr_root_modules;
    	rfract = per_total_R_in_root_system;
    	kla = saturated_kl_for_whole_plant;
    	aral = aral_at_sat_kl;
    	//latot = lai_at_sat_kl;

//ratot now computed in bgc 
    	//ratot = aral * latot;                 //sets total ralatot = Cells(3, 22)
    	rfract = rfract / 100.0;              //fraction of hydraulic resistance in root system
    	cpplant = cpplant * 1000.0;             //converts moles to mmol
//assign weibulls
    	int i = 0;

    	shootElements = treesParams.smodules;
    
    	for(i=1; i<=shootElements; i++) 
	{
    		al[i] = treesParams.al[i];
    		dslat[i] = treesParams.dslat[i];
    		dsax[i] = treesParams.dsax[i];
    	}

    	rootElements = treesParams.rmodules;
    
    	for(i=1+shootElements+1; i<=(rootElements+shootElements+1); i++) 
	{
        	ar[i] = treesParams.ar[i];
        	drlat[i] = treesParams.drlat[i];
        	drax[i] = treesParams.drax[i];
    	}

    	for (i=1; i<=shootElements; i++) 
	{
        	//al[i] = al[i] * latot;             //leaf area in each canopy module
        	al[i] = al[i] * 1.0;             //leaf area in each canopy module
    	}

    	ShootModules = shootElements;             //shoot modules

/*
    	for (i=1; i<=rootElements; i++) 
	{
        	ar[i + ShootModules + 1] *= ratot;
    	}
*/

    	for (i=1; i<=rootElements; i++) 
	{
        	l[i + ShootModules + 1] = (ratot*ar[i + ShootModules + 1] / (2.0 * M_PI * rr));
    	}

    	RootModules = rootElements;                    //number of root modules
    	totmodules = ShootModules + RootModules;

	if (reset == true)
	{
    		bsat[1][0] = axShoot_b;
    		bsat[1][1] = latShoot_b;
    		ccsat[1][0] = axShoot_c;
    		ccsat[1][1] = latShoot_c;
	}
    	b[1][1] = bsat[1][1];
    	cc[1][1] = ccsat[1][1];
    	b[1][0] = bsat[1][0];
    	cc[1][0] = ccsat[1][0];

    	for ( i = ShootModules+2; i <= totmodules+1; i++) 
	{
		if (reset == true)
		{
        		bsat[i][0] = axRoot_b;
        		bsat[i][1] = latRoot_b;
        		ccsat[i][0] = axRoot_c;
        		ccsat[i][1] = latRoot_c;
		}
        	b[i][1] = bsat[i][1];
        	cc[i][1] = ccsat[i][1];
        	b[i][0] = bsat[i][0];
        	cc[i][0] = ccsat[i][0];
    	}

//***************************************************************************************************
//		IMPORTANT:	SHOULD LOOP TILL Cells(22, 3 + i) = Empty. BUT IMPLEMENTION WAS
//		INCORRECT IN XLS FILE.... SoilElements HAS BEEN IMPLEMENTED CAUSE IN XLS
//		THE LOOP RUNNS FOR 3 ITERATIONS..... NEED TO CONFIRM
//***************************************************************************************************

    	soilpsiavg[t] = pd_at_sat_kl;  //NOW GETTING SOIL MOISTURE INITIALIZED FROM TREES
    	soilpsimin = soilpsiavg[t];
    	for (i = 1; i <= rootElements; i++)
	{
        	soilpsi[i + ShootModules + 1] = soilpsiavg[t]; //assign soil avg to root modules
    	}

//assign initial pressures to plant
    	for ( i = 1; i <= totmodules + 1; i++ ) 
	{
        	if ( i <= ShootModules + 1 ) 
		{
            		p[i][0] = soilpsiavg[t];
            		p[i][1] = soilpsiavg[t];
        	}
        	if ( i > ShootModules + 1 ) 
		{
            		p[i][0] = soilpsi[i];
            		p[i][1] = soilpsi[i];
        	}
    	}

//assign pressures to soil elements
    	for ( i = ShootModules + 2; i <= totmodules + 1; i++ ) 
	{
        	for (  int j = 2; j <= 1 + SoilRhizElements; j++ ) 
		{
            		p[i][j] = soilpsi[i];
        	}
    	}
    	ksolve(  r, pe, b, b1, n, n1, ks, ws_, rssoil, cf, vu, vl, saturatedKs, cpp,
              	wnu, wu, wnl, wl, ip, p,  al, l, rlateral, rLat_base, slateral, sLat_base,  ratot,
              	raxial,  raxlat,  ShootModules, totmodules, ar, axr, drax, latr, drlat,
              	ksat, soilpsiavg,  rfract, ptarg,  saxlat,  latot, dsax,  dslat, pleafave,
              	psimin, jl, soilpsi, cpu, cpl, f, ff, jmatrix, dpp, dp,   col,   indxx, subtract,
              	RootModules,  rr,  rd, Fabs,  cpplant,  gmd, gsd,  fc,  fs,  bkd,
              	tnode,  e,  dt,  fac,  SoilRhizElements, cc, kl, ku, rflux, modelledKL, 
		soilvol, initial_conductivity_root, decrement_root, initial_conductivity_shoot, 
		decrement_shoot, HydraulicModelFailCond);

//we want to retain a memory of the plant's hydraulic status
//to be restored after new saturated Ks have been computed
	if (reset == false)
	{
    		for ( i = 1; i <= totmodules+1; i++) 
		{
        		psimin[i][0] = psimin_in[i][0];
        		psimin[i][1] = psimin_in[i][1];
    		}
	}
	else
	{
    		for ( i = 1; i <= totmodules+1; i++) 
		{
        		psimin[i][0] = p[i][0];
        		psimin[i][1] = p[i][1];
    		}
	}

        pdecrement = 0.01;
        cincrement = 0.01;

    	for ( i = 1; i <= totmodules + 1; i++ ) 
	{ //resets data for lateral modules. i is plant element][k pressure increment][j c iteration
        	if ( i != ShootModules + 1 ) 
		{
//plc at new predawn
            		plco = (1.0 - exp(-pow((-psimin[i][1] / bsat[i][1]), ccsat[i][1]))) * 100.0; 
            		int k = 0;
            		pressure[k] = 0.0;
//loop calculates new plc//s based on new predawn and finds new b
            		do 
			{
                		k++;
                		pressure[k] = pressure[k - 1] - pdecrement;
                		if ( pressure[k] >= psimin[i][1] )
				{
                    			plc[k] = plco;
                		}
				else
				{
                    			plc[k] = (1.0 - exp( -pow( (-pressure[k] / bsat[i][1]), 
										ccsat[i][1] ))) * 100.0;
                		}
                		plc[k] = ( 1.0 - ( 100.0 - plc[k] ) / ( 100.0 - plco ) ) * 100.0;
                		if ( plc[k] < 66.0 ) 
				{
                    			newb[i][1] = fabs( pressure[k] );
                		}
            		} while ( plc[k] <= 99.50 );

            		ktotal = k - 1;
            		newc[i][1] = ccsat[i][1];
            		int j = 0;
            		rsquare[j] = 0.0;

//loop finds new c
            		do 
			{
                		j++;
                		newc[i][1] += cincrement;
                		sumy1 = 0.0;
                		sumy2 = 0.0;
                		sumprod = 0.0;

                		for ( k = 1; k <= ktotal; k++ ) 
				{
                    			plcweib[k] = (1.0 - exp(-pow( (-pressure[k] / newb[i][1]), 
									newc[i][1]) )) * 100.0;
                    			sumy1 += pow( plcweib[k], 2.0 );
                    			sumy2 += pow( plc[k], 2.0 );
                    			sumprod += plcweib[k] * plc[k];
                		}
                		rsquare[j] = pow( sumprod, 2.0 ) / ( sumy1 * sumy2 );
            		} while ( rsquare[j] >= rsquare[j - 1] );

            		newc[i][1] -= cincrement;
            		b[i][1] = newb[i][1]; //assigns new b's and c's to plant
            		cc[i][1] = newc[i][1];
        	}	// end of if
    	}	// end of for

//  resets data for axial modules. i is plant element][k pressure increment][j c iteration
	pdecrement = 0.01;
        cincrement = 0.01;
    	for (  i = 1; i <= totmodules + 1; i++ ) 
	{
        	if ( i != ShootModules + 1 ) 
		{
            		plco = ( 1.0 - exp( ( -pow( (-psimin[i][0] / 
					bsat[i][0]), ccsat[i][0]) ))) * 100.0; //plc at new predawn
            		int k = 0;
            		pressure[k] = 0.0;

//  loop calculates new plc//s based on new predawn and finds new b
            		do 
			{
                		k++;
                		pressure[k] = pressure[k - 1] - pdecrement;
                		if ( pressure[k] >= psimin[i][0] ) 
				{
                    			plc[k] = plco;
                		}
				else
				{
                    			plc[k] = ( 1.0 - exp( (-pow( (-pressure[k] / 
							bsat[i][0]), ccsat[i][0]) ))) * 100.0;
                		}
                		plc[k] = ( 1.0 - ( 100.0 - plc[k] ) / ( 100.0 - plco ) ) * 100.0;
                		if ( plc[k] < 66.0 ) 
				{
                    			newb[i][0] = fabs( pressure[k] );
                		}
            		} while ( plc[k] <= 99.50 );

            		ktotal = k - 1;
            		newc[i][0] = ccsat[i][0];
            		int j = 0;
            		rsquare[j] = 0.0;

// loop finds new c
            		do 
			{
                		j++;
                		newc[i][0] += cincrement;
                		sumy1 = 0.0;
                		sumy2 = 0.0;
                		sumprod = 0.0;
                		for ( k = 1; k <= ktotal; k++ ) 
				{
                    			plcweib[k] = ( 1.0 - exp( (-pow( (-pressure[k] / 
								newb[i][0]), newc[i][0]) ))) * 100.0;
                    			sumy1 += pow( plcweib[k], 2.0 );
                    			sumy2 += pow( plc[k], 2.0 );
                    			sumprod += plcweib[k] * plc[k];
                		}
                		rsquare[j] = pow( sumprod, 2.0 ) / ( sumy1 * sumy2 );
            		} while ( rsquare[j] >= rsquare[j - 1] );

            		newc[i][0] -= cincrement;
            		b[i][0] = newb[i][0];   // assigns new b's and c's to plant
            		cc[i][0] = newc[i][0];
        	}	// end of if
    	}	// end of for
}

void HydraulicModel::ksolve( double r[][ULAT], double pe[][ULAT], double b[][ULAT], double b1[][ULAT],
            		     double n[][ULAT], double n1[][ULAT], double ks[][ULAT], double ws_[][ULAT],
            		     double rssoil[], double cf[][ULAT], double vu[][ULAT], double vl[][ULAT],
            		     double saturatedKs[][4], double cpp[][ULAT], double wnu[][ULAT], 
			     double wu[][ULAT], double wnl[][ULAT], double wl[][ULAT],
            		     double ip[][ULAT], double p[][ULAT], double al[],
            		     double l[], double &rlateral, double &rLat_base,
            		     double slateral, double sLat_base, double ratot,
            		     double &raxial, double raxlat, const int ShootModules,
            		     const int totmodules, double ar[], double axr[],
            		     double drax[], double latr[], double drlat[],
            		     double ksat[][MD], double soilpsiavg[], double rfract,
            		     double &ptarg, double saxlat, double latot,
            		     double dsax[], double dslat[], double &pleafave,
            		     double psimin[][NMAX], double jl[][ULAT], double soilpsi[], 
			     double cpu[][ULAT], double cpl[][ULAT], double f[][ULAT], 
			     double ff[NMAX], double jmatrix[][NMAX], double dpp[][ULAT], 
			     double dp[NMAX], int col[], int indxx[], double subtract[],
            		     const int RootModules, double rr, double rd,
            		     double Fabs, double cpplant, double gmd,
            		     double gsd, double fc, double fs, double bkd,
            		     int &tnode, double e, double dt, double fac, 
			     const int SoilRhizElements, double cc[][2], double kl[][ULAT], 
			     double ku[][ULAT], double **rflux, double &modelledKL, 
			     double &soilvol, double initial_conductivity_root, double decrement_root, 
            		     double initial_conductivity_shoot, double decrement_shoot, 
			     int &HydraulicModelFailCond) 
{
    	int i=0, j=0, k=0, t=1;
    	double increment=0.0, Sum=0.0, peter=0.0, pRoot=0.0, saxial=0.0, pLeaf=0.0;

//we need to confirm that these variable are getting the correct values.
//added by DR 8/3/07
    	//rlateral = initial_conductivity_root/latot;
    	rlateral = initial_conductivity_root/min(1.0, latot);
    	rLat_base = decrement_root/min(1.0, latot);
    	//slateral = initial_conductivity_shoot/latot;
    	slateral = initial_conductivity_shoot/min(1.0, latot);
    	sLat_base = decrement_shoot/min(1.0, latot);

    	for ( i=1; i <= ShootModules; i++) 
	{
        	ks[i][0] = (double) 10000000.0;
        	ks[i][1] = (double) 10000000.0;
    	}
//solve for root element conductances given distances and areas
    	do 
	{
        	rlateral -= increment;  // decrement k's until pressuures match actual
        	raxial = rlateral * raxlat;

        	for ( i = ShootModules+1; i <= totmodules; i++) 
		{
            		Sum = ar[i+1] * ratot;  // note that ar's are offset by 1 with 2D k array!
            		for ( k = i+1; k <= totmodules; k++) 
			{
                		Sum += ar[k+1] * ratot; // tallies proximal root areas
            		}
// converts conductivities to conductances (resistances)
            		axr[i] = 1.0 / (raxial * Sum * (1.0 / drax[i+1]));
           		latr[i] = 1.0 / (rlateral * ar[i+1] * ratot * (1.0 / drlat[i+1]));
        	}

// put root conductances in 2D array
        	j=1; 
		k=0;
        	for ( i = ShootModules+1; i<=totmodules; i++) 
		{
            		ks[i+1][0] = 1.0 / axr[i];
            		ksat[i+1][0] = ks[i+1][0];
            		ks[i+1][1] = 1.0 / latr[i];
            		ksat[i+1][1] = ks[i+1][1];
            		saturatedKs[j][k] = (double) i+1.0;
            		k++;
            		saturatedKs[j][k] = ks[i+1][0];
            		k++;
            		saturatedKs[j][k] = ks[i+1][1];
            		j++;
            		k = 0;
        	}

        	j = 1;

        	solvek( r, pe, b, b1, n, n1, ks, ws_, rssoil, cf, vu, vl, saturatedKs,
                	cpp, wnu, wu, wnl, wl, ip, p, al, l, jl, soilpsi, cpu, cpl, cc,
                	kl, ku, rflux, ar, f, ff, jmatrix, dpp, dp, col, indxx, subtract,
                	latot, RootModules, rr, rd, Fabs, cpplant, gmd, gsd, fc, fs, bkd,
                	tnode, pleafave, e, dt, fac, SoilRhizElements, ShootModules, 
			totmodules, soilvol, HydraulicModelFailCond );

        	pRoot = p[ShootModules+1][0];
        	peter = p[ShootModules+1][0] - (soilpsiavg[t] - (rfract * (soilpsiavg[t] - ptarg)));

        	if ( peter > 0.3 ) 
		{
            		increment = rLat_base;
        	}
        	else if ( peter <= 0.3 && peter >= 0.1 ) 
		{
            		increment = rLat_base / 4.0;
        	}
        	else if ( peter < 0.1 ) 
		{
            		//increment = rLat_base / ULAT;
            		increment = rLat_base / 20.0;
        	}

    	} while ( p[ShootModules+1][0] >= ( soilpsiavg[t] - (rfract * (soilpsiavg[t] - ptarg)) )  ); 
//loop until root collar pressure is correct

//repeat process for shoot
//solves for shoot elements from leaf areas and ratios

    	do 
	{
        	slateral = slateral - increment;
        	saxial = slateral * saxlat;
        	for ( i = 1;  i <= ShootModules; i++) 
		{
            		Sum = al[i] * latot;
            		for ( k = i-1; k>=1; k--) 
			{
                		Sum += al[k] * latot;   // tallies distal leaf area
            		}
// converts conductivities to conductances (resistances)
            		axr[i] = 1.0 / (saxial * Sum * (1.0 / dsax[i]));
            		latr[i] = 1.0 / (slateral * al[i] * latot * (1.0 / dslat[i]));
        	}

        	j=0; 
		k=0;
// put into 2D array
        	for ( i = 1; i <= ShootModules; i++) 
		{
            		ks[i][1] = 1.0 / latr[i];
            		ksat[i][1] = ks[i][1];
            		ks[i][0] = 1.0 / axr[i];
            		ksat[i][0] = ks[i][0];
            		saturatedKs[j][k] = (double) i;
            		k++;
            		saturatedKs[j][k] = ks[i][0];
            		k++;
            		saturatedKs[j][k] = ks[i][1];
            		j++;
        	}

        	solvek( r, pe, b, b1, n, n1, ks, ws_, rssoil, cf, vu, vl, saturatedKs, cpp, 
			wnu, wu, wnl, wl, ip, p,  al,  l, jl, soilpsi, cpu, cpl, cc, kl, ku, 
			rflux, ar, f, ff, jmatrix, dpp, dp, col, indxx, subtract, latot,  
			RootModules,  rr,  rd, Fabs,  cpplant,  gmd, gsd,  fc,  fs,  bkd,
                	tnode, pleafave, e,  dt,  fac,  SoilRhizElements, ShootModules,  
			totmodules, soilvol, HydraulicModelFailCond );

        	pLeaf = pleafave;
        	peter = pleafave - ptarg;
        	if ( peter > 0.3 ) 
		{
            		increment = sLat_base;
        	}
        	else if ( peter <= 0.3 && peter >= 0.1 ) 
		{
            		increment = sLat_base / 5.0;
        	}
        	else if ( peter < 0.1 ) 
		{
            		increment = sLat_base / 10.0;
        	}

    	} while ( pleafave >= ptarg );    // loop until pressures match

/*
    	for ( i = 1; i <= totmodules+1; i++) 
	{
        	psimin[i][0] = p[i][0];
        	psimin[i][1] = p[i][1];
    	}
*/
    	modelledKL = e / (soilpsiavg[t] - pleafave);     // should match measured klsat
}


void HydraulicModel::solvek( double r[][ULAT], double pe[][ULAT], double b[][ULAT], 
				double b1[][ULAT], double n[][ULAT],
              			double n1[][ULAT], double ks[][ULAT], double ws_[][ULAT], 
				double rssoil[], double cf[][ULAT],
              			double vu[][ULAT], double vl[][ULAT], double saturatedKs[][4],
              			double cpp[][ULAT], double wnu[][ULAT], double wu[][ULAT], 
				double wnl[][ULAT], double wl[][ULAT], double ip[][ULAT], 
				double p[][ULAT], double al[], double l[],
              			double jl[][ULAT], double soilpsi[], double cpu[][ULAT], 
				double cpl[][ULAT], double cc[][2], double kl[][ULAT], 
				double ku[][ULAT], double **rflux, double ar[MD], 
				double f[][ULAT], double ff[NMAX], double jmatrix[][NMAX],
              			double dpp[][ULAT], double dp[NMAX], int col[], 
				int indxx[], double subtract[], double latot, 
				const int RootModules, double rr, double rd,
              			double Fabs, double cpplant, double gmd,
              			double gsd, double fc, double fs, double bkd,
              			int &tnode, double &pleafave,
              			double e, double dt, double fac, 
				const int SoilRhizElements, const int ShootModules, 
				const int totmodules, double &soilvol, int &HydraulicModelFailCond ) 
{
    	int ii=0, i=0, j=0, k=0, m=0, h=0;
    	double se=0.0, phiu=0.0, phil=0.0, nn=0.0, aa=0.0, gammp=0.0, gamma=0.0, x=0.0; 
    	double abv=0.0, lim=0.0, deltap=0.0, Sum=0.0, oldpress=0.0, sumpress=0.0;

    	double base1, base2, exponent, result1, result2;

    	int nits;

    	bool exitError = false;
    	initialize( r,  pe,  b,  b1, n,  n1,  ks,  ws_, rssoil,  cf,  vu,  vl, 
			saturatedKs,  cpp,  wnu,  wu, wnl,  wl,  ip, p, al, l,  
			SoilRhizElements,  ShootModules,  totmodules, latot,  RootModules,
			rr, rd, Fabs, cpplant,  gsd,  fc,  fs, bkd, gmd, soilvol );

    	nits = 0;		//sets # iterations to 0
    	oldpress = 0.0;
    	sumpress = 0.0;
    	ii = 0;

    	do 
	{
        	oldpress = sumpress;
        	nits = 0;

// this loop finds new potentials and contents at each time step
        	do 
		{
            		se = 0.0;
            		for (i = 1; i <= ShootModules; i++ ) 
			{
                		jl[i][2] = e * al[i]; // sets evaporation rate
            		}

// makes sure outermost soil node is constant at psisoil
            		for (i = ShootModules + 2; i <= totmodules + 1; i++ )
			{
                		p[i][SoilRhizElements + 1] = soilpsi[i];
            		}

            		for (i = 1; i <= ShootModules; i++ ) 
			{ //sets upper and lower elements for lateral elements of shoot
                		if ( p[i][1] > ip[i][1] ) 
				{
                    			wnu[i][1] = ws_[i][1] * vu[i][1];
                    			cpu[i][1] = 0.01 * vu[i][1] / dt;
                		}
                		else 
				{
//content of upper half
                    			wnu[i][1] = ws_[i][1] * vu[i][1] - ( ( ip[i][1] - p[i][1] ) * 
								0.5 * cpp[i][1] );
                    			cpu[i][1] = 0.5 * cpp[i][1] / dt; //dw/dp *1/dt for jacobian
                		}

                		if ( p[i][0] > ip[i][1] ) 
				{
                    			wnl[i][1] = ws_[i][1] * vl[i][1];
// strictly, this should be zero, but set to arbitrary small #?
                    			cpl[i][1] = 0.01 * vl[i][1] / dt; 
                		}
                		else 
				{
//fabs content, lower half of plant element
                    			wnl[i][1] = ws_[i][1] * vl[i][1] - ( ( ip[i][1] - p[i][0] ) * 
								0.5 * cpp[i][1] ); 
                    			cpl[i][1] = 0.5 * cpp[i][1] / dt; // dw/dp *1/dt for jacobian
                		}
            		}

//sets upper and lower elements for axial elements of shoot
            		for (i = 1; i <= ShootModules; i++ ) 
			{
                		if ( p[i][0] > ip[i][0] ) 
				{
                    			wnu[i][0] = ws_[i][0] * vu[i][0];
                    			cpu[i][0] = 0.01 * vu[i][0] / dt;
                		}
                		else 
				{
                    			wnu[i][0] = ws_[i][0] * vu[i][0] - ( ( ip[i][0] - p[i][0] ) * 
								0.5 * cpp[i][0] ); //content of upper half
                    			cpu[i][0] = 0.5 * cpp[i][0] / dt; //dw/dp *1/dt for jacobian
                		}

                		if ( p[i + 1][0] > ip[i][0] ) 
				{
                    			wnl[i][0] = ws_[i][0] * vl[i][0];
//strictly, this should be zero, but set to arbitrary small #?
                    			cpl[i][0] = 0.01 * vl[i][0] / dt; 
                		}
                		else 
				{
                    			wnl[i][0] = ws_[i][0] * vl[i][0] - ( ( ip[i][0] - p[i + 1][0] ) * 
								0.5 * cpp[i][0] ); //fabs
// content, lower half of plant element
                    			cpl[i][0] = 0.5 * cpp[i][0] / dt; //dw/dp *1/dt for jacobian
                		}
            		}

// sets upper and lower elements for axial elements of root
            		for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{
                		if ( p[i - 1][0] > ip[i][0] ) 
				{
                    			wnu[i][0] = ws_[i][0] * vu[i][0];
                    			cpu[i][0] = 0.01 * vu[i][0] / dt;
                		}
                		else 
				{
                    			wnu[i][0] = ws_[i][0] * vu[i][0] - ( ( ip[i][0] - p[i - 1][0] ) * 
								0.5 * cpp[i][0] ); //content of upper half
                    			cpu[i][0] = 0.5 * cpp[i][0] / dt; //dw/dp *1/dt for jacobian
                		}

                		if ( p[i][0] > ip[i][0] ) 
				{
                    			wnl[i][0] = ws_[i][0] * vl[i][0];
//strictly, this should be zero, but set to arbitrary small #?
                    			cpl[i][0] = 0.01 * vl[i][0] / dt; 
                		}
                		else 
				{
//fabs content, lower half of plant element
                    			wnl[i][0] = ws_[i][0] * vl[i][0] - ( ( ip[i][0] - p[i][0] ) * 
								0.5 * cpp[i][0] ); 
                    			cpl[i][0] = 0.5 * cpp[i][0] / dt; //dw/dp *1/dt for jacobian
                		}
            		}

            		for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{ //sets upper and lower elements for lateral elements of root
                		if ( p[i][0] > ip[i][1] ) 
				{ //remember: pressures are nodal, other values are elements!!
                    			wnu[i][1] = ws_[i][1] * vu[i][1];
                    			cpu[i][1] = 0.01 * vu[i][1] / dt;
                		}
                		else 
				{
                    			wnu[i][1] = ws_[i][1] * vu[i][1] - ( ( ip[i][1] - p[i][0] ) * 
								0.5 * cpp[i][1] ); //content of upper half
                    			cpu[i][1] = 0.5 * cpp[i][1] / dt; //dw/dp *1/dt for jacobian
                		}

                		if ( p[i][1] > ip[i][1] ) 
				{
                    			wnl[i][1] = ws_[i][1] * vl[i][1];
//strictly, this should be zero, but set to arbitrary small #?
                    			cpl[i][1] = 0.01 * vl[i][1] / dt; 
                		}
                		else 
				{
//fabs content, lower half of plant element
                    			wnl[i][1] = ws_[i][1] * vl[i][1] - ( ( ip[i][1] - p[i][1] ) * 
								0.5 * cpp[i][1] ); 
                    			cpl[i][1] = 0.5 * cpp[i][1] / dt; //dw/dp *1/dt for jacobian
                		}
            		}

            		for (i = 1; i <= ShootModules; i++ ) 
			{ //get upper and lower k's for axial shoot elements
                		base1 = -p[i + 1][0]/b[i][0];
                		base2 = -p[i][0]/b[i][0];
                		exponent = cc[i][0];

                		if( !verifyPowCondition(base1,exponent) || 
							!verifyPowCondition(base2,exponent) ) 
				{
                        		//cout << "ERROR pow AT 1a\n";
                        		exitError = true;
                        		goto exitInner;
                		}

                		result1 = pow(base1,exponent);
                		result2 = pow(base2,exponent);

                		if ( !verifyExpCondition(result1) || !verifyExpCondition(result2) ) 
				{
                        		//cout << "ERROR pow AT 1b\n";
                        		exitError = true;
                        		goto exitInner;
                		}

//weibull function for lower half of element i
                		kl[i][0] = ks[i][0] * exp( -result1 ); 
//weibull function for upper half of element i
                		ku[i][0] = ks[i][0] * exp( -result2 ); 

                		nn = result2;
//assigns a value for incomplete gamma function
                		aa = 1.0 / cc[i][0]; 

//calls function pgamma
                		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

//phi for upper half
                		phiu = (1.0 - gammp) * gamma * ks[i][0] * b[i][0] / cc[i][0]; 
                		nn = result1;

                		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

// pgamma again with lower half n value
//phi for lower half
                		phil = ( 1.0 - gammp ) * gamma * ks[i][0] * b[i][0] / cc[i][0]; 
//flux of element i--CF converts conductivity to conductance
                		jl[i][0] = cf[i][0] * ( phil - phiu ); 
            		}

//get upper and lower k's for axial root elements
            		for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{
                		base1 = -p[i][0]/b[i][0];
                		base2 = -p[i-1][0]/b[i][0];
                		exponent = cc[i][0];

                		if( !verifyPowCondition(base1,exponent) || 
							!verifyPowCondition(base2,exponent) )  
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

                		result1 = pow(base1,exponent);
                		result2 = pow(base2,exponent);

                		if( !verifyExpCondition(result1)  || !verifyExpCondition(result2) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

                		kl[i][0] = ks[i][0] * exp( -result1 );
                		ku[i][0] = ks[i][0] * exp( -result2 );
                		nn = result2;
//assigns a value for incomplete gamma function
                		aa = 1.0 / cc[i][0]; 

                		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

//phi for upper half
                		phiu = ( 1.0 - gammp ) * gamma * ks[i][0] * b[i][0] / cc[i][0]; 
                		nn = result1; //n value for lower half

                		if ( pgamma(aa, nn, fac, x, gamma, gammp) )
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

// pgamma again with lower half n value
//phi for lower half
                		phil = ( 1.0 - gammp ) * gamma * ks[i][0] * b[i][0] / cc[i][0]; 
//flux of element i--CF converts conductivity to conductance
                		jl[i][0] = cf[i][0] * ( phil - phiu ); 
            		}

// get upper and lower k's for lateral shoot elements
            		for (i = 1; i <= ShootModules; i++ ) 
			{
                		base1 = -p[i][0]/b[i][1];
                		base2 = -p[i][1]/b[i][1];
                		exponent = cc[i][1];

                		if ( !verifyPowCondition(base1,exponent) || 
							!verifyPowCondition(base2,exponent) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

                		result1 = pow(base1,exponent);
                		result2 = pow(base2,exponent);

                		if ( !verifyExpCondition(result1) || !verifyExpCondition(result2) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

                		kl[i][1] = ks[i][1] * exp( -result1 );
                		ku[i][1] = ks[i][1] * exp( -result2 );
                		nn = result2;
//assigns a value for incomplete gamma function
                		aa = 1.0 / cc[i][1]; 

                		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

//phi for upper half
                		phiu = ( 1.0 - gammp ) * gamma * ks[i][1] * b[i][1] / cc[i][1]; 
                		nn = result1;

                		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

//phi for lower half
                		phil = ( 1.0 - gammp ) * gamma * ks[i][1] * b[i][1] / cc[i][1]; 
//flux of element i--CF converts conductivity to conductance
                		jl[i][1] = cf[i][1] * ( phil - phiu ); 
            		}

//get upper and lower k's for lateral root elements
            		for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{
                		base1 = -p[i][1]/b[i][1];
                		base2 = -p[i][0]/b[i][1];
                		exponent = cc[i][1];

                		if ( !verifyPowCondition(base1,exponent) || 
							!verifyPowCondition(base2,exponent) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

                		result1 = pow(base1,exponent);
                		result2 = pow(base2,exponent);

                		if ( !verifyExpCondition(result1) || !verifyExpCondition(result2) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

                		kl[i][1] = ks[i][1] * exp( -result1 );
                		ku[i][1] = ks[i][1] * exp( -result2 );
                		nn = result2;
                		aa = 1.0 / cc[i][1]; //assigns a value for incomplete gamma function

                		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

//phi for upper half
                		phiu = ( 1.0 - gammp ) * gamma * ks[i][1] * b[i][1] / cc[i][1]; 
                		nn = result1;

                		if ( pgamma(aa, nn, fac, x, gamma, gammp) ) 
				{
                        		exitError = true;
					HydraulicModelFailCond = 1;
                        		goto exitInner;
                		}

//phi for lower half
                		phil = ( 1.0 - gammp ) * gamma * ks[i][1] * b[i][1] / cc[i][1]; 
//flux of element i--CF converts conductivity to conductance
                		jl[i][1] = cf[i][1] * ( phil - phiu ); 
            		}

// sets upper and lower element properties for soil
            		for ( i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{
                		for ( j = 2; j <= SoilRhizElements + 1; j++ ) 
				{
                    			if ( p[i][j] < pe[i][j] ) 
					{
//absolute content lower half of element [eqn 5.9][campbell]
                        			wnl[i][j] = ws_[i][j] * vl[i][j] * 
								pow( (pe[i][j]/p[i][j]), b1[i][j] ); 	
//dw/dp*1/dt for jacobian][lower half of element
                        			cpl[i][j] = -wnl[i][j] / (b[i][j] * p[i][j] * dt);	 
//k[i] biased to lower node p [eqn. 6.14][campbell]
                        			kl[i][j] = ks[i][j] * pow( (pe[i][j] / p[i][j]), n[i][j]);
//a hybrid between phi[i] and phi[i+1][eqn. 8.19 of campbell]
                        			phil = kl[i][j] * p[i][j] / n1[i][j];
                    			}
                    			else 
					{
                        			wnl[i][j] = ws_[i][j] * vl[i][j];
//I think this is an arbitrary small number because dw/dp would otherwise be zero
                        			cpl[i][j] = 0.01 * vl[i][j] / dt; 
                        			kl[i][j] = ks[i][j];
// This is integral for p > pe
                        			phil = ks[i][j] * (pe[i][j] * n[i][j] / 
									n1[i][j] + p[i][j]); 
                    			}
                    			if ( p[i][j-1] < pe[i][j] ) 
					{
//absolute content of upper half of element
                        			wnu[i][j] = ws_[i][j] * vu[i][j] * 
								pow( (pe[i][j] / p[i][j - 1]), b1[i][j] );
//dw/dp*1/dt 'dw/dp*1/dt for jacobian][upper half of element
                        			cpu[i][j] = -wnu[i][j] / (b[i][j] * p[i][j - 1] * dt);
//k[i]][no bias
                        			ku[i][j] = ks[i][j] * pow( (pe[i][j] / 
								p[i][j - 1]), n[i][j] ); 
//phi[i]][no bias
                        			phiu = ku[i][j] * p[i][j - 1] / n1[i][j];
                    			}
                    			else 
					{
                        			wnu[i][j] = ws_[i][j] * vu[i][j];
//arbitrary small number to get model going at p above pe
                        			cpu[i][j] = 0.01 * vu[i][j] / dt; 
                        			ku[i][j] = ks[i][j];
//integral for p>pe
                        			phiu = ku[i][j] * (pe[i][j] * n[i][j] / 
								n1[i][j] + p[i][j - 1]);
                    			}
//flux of element i--CF converts conductivity to conductance
                    			jl[i][j] = cf[i][j] * (phil - phiu);
                		}
            		}

            		for (i = ShootModules + 2; i <= totmodules + 1; i++ )
			{
                 		rflux[ii][i] = jl[i][2] / ar[i];
            		}

// tie up loose ends of network for mass balance
            		jl[totmodules + 2][0] = 0.0;
            		jl[0][0] = 0.0;
            		wnu[totmodules + 2][0] = 0.0;
            		wu[totmodules + 2][0] = 0.0;
            		wnl[0][0] = 0.0;
            		wl[0][0] = 0.0;

            		for (i = 1; i <= ShootModules; i++ ) 
			{
                		wnl[i][2] = 0.0;
                		wl[i][2] = 0.0;
            		}

            		for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{
                		wnu[i][SoilRhizElements + 2] = 0.0;
                		wu[i][SoilRhizElements + 2] = 0.0;
                		jl[i][SoilRhizElements + 2] = 0.0;
            		}

            		cpl[0][0] = 0.0;
            		cf[0][0] = 0.0;
            		cpu[totmodules + 2][0] = 0.0;
            		cf[totmodules + 2][0] = 0.0;

// mass balance for shoot lateral nodes
            		for (i = 1; i <= ShootModules; i++ ) 
			{
                		f[i][1] = ( wnu[i][1] - wu[i][1] + wnl[i][2] - wl[i][2] ) / 
									dt + jl[i][2] - jl[i][1];
            		}

// mass balance for shoot axial nodes
            		for (i = 1; i <= ShootModules + 1; i++ ) 
			{
                		if ( i == ShootModules + 1 ) 
				{
                    			f[i][0] = ( wnu[i + 1][0] - wu[i + 1][0] + 
							wnl[i - 1][0] - wl[i - 1][0] ) / 
							dt + jl[i - 1][0] - jl[i + 1][0];
                		}
                		else 
				{
                    			f[i][0] = ( wnu[i][0] - wu[i][0] + wnl[i][1] - 
							wl[i][1] + wnl[i - 1][0] - wl[i - 1][0] ) / 
							dt + jl[i][1] + jl[i - 1][0] - jl[i][0];
                		}
            		}

            		for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{ // mass balance for axial root nodes
                		f[i][0] = ( wnu[i][1] - wu[i][1] + wnu[i + 1][0] - 
							wu[i + 1][0] + wnl[i][0] - wl[i][0] ) / 
							dt + jl[i][0] - jl[i][1] - jl[i + 1][0];
            		}

            		for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{ // mass balance for lateral root and soil nodes
                		for (j = 1; j <= SoilRhizElements + 1; j++ ) 
				{
                    			f[i][j] = ( wnu[i][j + 1] - wu[i][j + 1] + 
							wnl[i][j] - wl[i][j] ) / 
							dt + jl[i][j] - jl[i][j + 1];
                    			if ( j == SoilRhizElements + 1 ) 
					{
                        			f[i][j] = 0.0;
                    			}
                		}
            		}

// next sum the se's
            		se = 0.0;
            		for (i = 1; i <= totmodules + 1; i++ ) 
			{
                		se += f[i][0];
            		}

            		for (i = 1; i <= ShootModules; i++ )
			{
                		se += f[i][1];
            		}

            		for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{
                		for (j = 1; j <= SoilRhizElements + 1; j++ ) 
				{
                    			se += f[i][j];
                		}
            		}

//total # of NODES
            		tnode = ShootModules + totmodules + 1 + ( SoilRhizElements + 1 ) * RootModules; 
            		k = 1;

// reset f//s to ff's (for gauss) in order of equations in matrix
            		for (i = 1; i <= ShootModules + 1; i++ ) 
			{
                		if ( i == ShootModules + 1 ) 
				{
                    			ff[k] = f[i][0];
                    			k++;
                		}
                		else 
				{
                    			for (j = 0; j <= 1; j++ ) 
					{
                        			ff[k] = f[i][j];
                        			k++;
                    			}
                		}
            		}

            		for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{
                		for (j = 0; j <= SoilRhizElements + 1; j++ ) 
				{
                    			ff[k] = f[i][j];
                    			k++;
                		}
            		}

//zeros the jacobian matrix
            		for (i = 1; i <= tnode; i++ ) 
			{
                		for (j = 1; j <= tnode; j++ ) 
				{
                    			jmatrix[i][j] = 0.0;
                		}
            		}

            		i = 2;
            		j = 2;
//calculate and load shoot lateral nodes
            		for ( k=1; k <= ShootModules; k++ ) 
			{
                		jmatrix[i][j] = ku[k][1] * cf[k][1] + cpu[k][1]; //dF[k1]/dY[k1]
                		jmatrix[i][j - 1] = -kl[k][1] * cf[k][1]; //dF[k][1]/dY[k][0]
                		i = i + 2;
                		j = j + 2;
            		}

//load shoot axial nodes
            		i = 1;
            		j = 1;
            		for ( k=1; k <= ShootModules; k++ ) 
			{
                		jmatrix[i][j] = ku[k][0] * cf[k][0] + kl[k][1] * cf[k][1] + 
							kl[k - 1][0] * cf[k - 1][0] + cpu[k][0] + 
							cpl[k - 1][0] + cpl[k][1]; //dF[k0]/dY[k0]
                		jmatrix[i][j - 2] = -ku[k - 1][0] * cf[k - 1][0]; //dF[k][0]/dY[k-1][0]
                		jmatrix[i][j + 1] = -ku[k][1] * cf[k][1]; //dF[k][0]/dY[k][1]
                		jmatrix[i][j + 2] = -kl[k][0] * cf[k][0]; //dF[k][0]/dY[k+1][0]
                		i = i + 2;
                		j = j + 2;
            		}

//load root collar node
            		k = ShootModules + 1;
            		jmatrix[i][j] = kl[k - 1][0] * cf[k - 1][0] + ku[k + 1][0] * cf[k + 1][0] + 
							cpu[k + 1][0] + cpl[k - 1][0]; //dF[k][0]/dY[k][0]
            		jmatrix[i][j - 2] = -ku[k - 1][0] * cf[k - 1][0]; //dF[k][0]/dY[k-1][0]
            		jmatrix[i][j + 1] = -kl[k + 1][0] * cf[k + 1][0]; //dF[k][0]/dY[k+1][0]

//load axial root nodes
            		i++;
            		j++;
            		for ( k = ShootModules + 2; k <= totmodules + 1; k++ ) 
			{
                		jmatrix[i][j] = kl[k][0] * cf[k][0] + ku[k + 1][0] * cf[k + 1][0] + 
							ku[k][1] * cf[k][1] + cpl[k][0] + 
							cpu[k + 1][0] + cpu[k][1]; //dF[k,0]/dY[k,0]
                		if ( k == ShootModules + 2 ) 
				{
                    			jmatrix[i][j - 1] = -ku[k][0] * cf[k][0];
                		}
                		else 
				{
//dF[k,0]/dY[k-1,0]
                    			jmatrix[i][j - SoilRhizElements - 2] = -ku[k][0] * cf[k][0]; 
                		}
                		jmatrix[i][j + 1] = -kl[k][1] * cf[k][1]; //dF[k,0]/dY[k,1]
//dF[k,0]/dY[k+1,0]
                		jmatrix[i][j + SoilRhizElements + 2] = -kl[k + 1][0] * cf[k + 1][0]; 
                		i = i + SoilRhizElements + 2;
                		j = j + SoilRhizElements + 2;
            		}

// load lateral root and soil nodes
            		i = 2 * ShootModules + 3;
            		j = 2 * ShootModules + 3;
            		for ( k = ShootModules + 2; k <= totmodules + 1; k++ ) 
			{
                		for ( m = 1; m <= SoilRhizElements + 1; m++ ) 
				{
//dF[k,m]/dY[k,m]
                    			jmatrix[i][j] = kl[k][m] * cf[k][m] + ku[k][m + 1] * 
								cf[k][m + 1] + cpl[k][m] + cpu[k][m + 1]; 
                    			jmatrix[i][j - 1] = -ku[k][m] * cf[k][m]; //dF[k,m]/dY[k,m-1]
//dF[k,m]/dY[k,m+1]
                    			jmatrix[i][j + 1] = -kl[k][m + 1] * cf[k][m + 1]; 
                    			i++;
                    			j++;
                		}
                		i++;
                		j++;
            		}

            		gauss( ShootModules, SoilRhizElements, tnode, jmatrix, col, ff, 
										indxx, subtract, dp );

            		k = 1;
            		for ( i = 1; i <= ShootModules; i++ ) 
			{
                		dpp[i][0] = dp[k];
                		dpp[i][1] = dp[k + 1];
                		k = k + 2;
            		}

            		dpp[ShootModules + 1][0] = dp[k];
            		k++;

            		for ( i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{
                		dpp[i][0] = dp[k];
                		for ( m = 1; m <= SoilRhizElements + 1; m++ ) 
				{
                    			k++;
                    			dpp[i][m] = dp[k];
                		}
                		k++;
            		}

// assign new p//s
            		for ( i = 1; i <= totmodules + 1; i++ ) 
			{
                		lim = 0.8 * fabs( p[i][0] );
                		abv = fabs( dpp[i][0] );
                    		if ( ( abv > lim ) && ( p[i][0] < -0.1 ) ) 
				{
                        		dpp[i][0] *= lim / abv;
                    		}
                		p[i][0] -= dpp[i][0];
            		}

            		for ( i = 1; i <= ShootModules; i++ ) 
			{
                		lim = 0.8 * fabs( p[i][1] );
                		abv = fabs( dpp[i][1] );
                		if ( ( abv > lim ) && ( p[i][1] < -0.1 ) ) 
				{
                    			dpp[i][1] *= lim / abv;
                		}
                		p[i][1] -= dpp[i][1];
            		}

            		for ( i = ShootModules + 2; i <= totmodules + 1; i++ ) 
			{
                		for ( j = 1; j <= SoilRhizElements + 1; j++ ) 
				{
                    			lim = 0.8 * fabs( p[i][j] );
                    			abv = fabs( dpp[i][j] );
                    			if ( ( abv > lim ) && ( p[i][j] < -0.1 ) ) 
					{
                        			dpp[i][j] *= lim / abv;
                    			}
                    			p[i][j] -= dpp[i][j];
                		}
            		}
            		nits++;

        	} while ( (se >= IM) && (nits <= MAXITS) );

        	exitInner: cout<<"";
        	if ( nits > MAXITS || exitError ) 
		{
            		goto exitOuter;
        	}

//add up all pressures for steady-state test
        	sumpress = 0.0;

        	for ( i = 1; i <= totmodules + 1; i++ ) 
		{
            		sumpress += fabs( p[i][0] );
        	}

        	for ( i = 1; i <= ShootModules; i++ ) 
		{
            		sumpress += fabs( p[i][1] );
        	}

        	for ( i = ShootModules + 2; i <= totmodules + 1; i++ ) 
		{
            		for ( j = 1; j <= SoilRhizElements + 1; j++ ) 
			{
                		sumpress += fabs( p[i][j] );
            		}
        	}

        	deltap = fabs( (sumpress-oldpress) );

//reset water contents of all ELEMENTS
        	for ( i = 1; i <= totmodules + 1; i++ ) 
		{
            		if ( i != ShootModules + 1 ) 
			{
                		wl[i][0] = wnl[i][0];
                		wu[i][0] = wnu[i][0];
                		wl[i][1] = wnl[i][1];
                		wu[i][1] = wnu[i][1];
            		}
        	}

        	for ( i = ShootModules + 2; i <= totmodules + 1; i++ ) 
		{
            		for ( j = 2; j <= SoilRhizElements + 1; j++ ) 
			{
                		wl[i][j] = wnl[i][j];
                		wu[i][j] = wnu[i][j];
            		}
        	}
        	h++;
    	} while ( ( deltap >= 0.001 ) && ( h <= MAXH ) );

    	exitOuter: cout<<"";
//calculate pleafave

    	Sum = 0.0;
    	for (  i = 1; i <= ShootModules; i++ ) 
	{
        	Sum += p[i][1];
    	}

    	pleafave = Sum / (double) ShootModules;
}

void HydraulicModel::initialize( double r[MD][ULAT], double pe[MD][ULAT], 
				double b[MD][ULAT], double b1[MD][ULAT],
                 		double n[MD][ULAT], double n1[MD][ULAT], 
				double ks[MD][ULAT], double ws_[MD][ULAT],
                 		double rssoil[], double cf[MD][ULAT], 
				double vu[MD][ULAT], double vl[MD][ULAT],
                 		double saturatedKs[][4], double cpp[MD][ULAT], 
				double wnu[MD][ULAT], double wu[MD][ULAT],
                 		double wnl[MD][ULAT], double wl[MD][ULAT], 
				double ip[MD][ULAT], double p[MD][ULAT],
                 		double al[], double l[], int SoilRhizElements, 
				int ShootModules, int totmodules,
                 		double latot, int RootModules, double rr, double rd, double Fabs,
                 		double cpplant, double gsd, double fc, double fs,
                 		double bkd, double gmd, double &soilvol )
{
    	int snmax=0, i=0, j=0;
    	double wsplant=0.0;

    	snmax = SoilRhizElements + 1; //outermost soil node
    	for (i = ShootModules+2; i <= totmodules+1; i++ ) 
	{
//radius of outer soil bound based on root density, eqn. 11.10, campbell
        	r[i][snmax] = pow( (1.0/(rd*M_PI*Fabs)), (double)0.5 ); 
        	r[i][1] = rr; //sets to root radius
    	}

   	wsplant = cpplant * latot * 10.0; //plant water in mmol yielded at -10 MPa

//assign ks, capacitance, volumes to soil nodes and elements--next two for loops
    	for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
	{
        	for (j = 2; j <= SoilRhizElements + 1; j++ ) 
		{
            		pe[i][j] = -0.5 * pow( gmd , (double) -0.5 ); //eqn 5.10 in campbell
            		b[i][j] = -2.0 * pe[i][j] + 0.2 * gsd; //eqn 5.11 in campbell
            		b1[i][j] = 1.0 / b[i][j];
            		n[i][j] = 2.0 + 3.0 / b[i][j]; //see p. 56, campbell
            		n1[i][j] = 1.0 - n[i][j];
//eqn 5.12 in campbell
            		pe[i][j] *= pow( (double)( bkd/1.3 ) , (double)( 0.67 * b[i][j] ) ); 
            		pe[i][j] = pe[i][j] * 0.001; //conversion from J/kg to MPa
            		ks[i][j] = 0.004 * pow( ( 1.3 / bkd ), ( 1.3 * b[i][j] ) ) * 
						exp( -6.9 * fc - 3.7 * fs ); //eqn 6.12a in campbell
            		ks[i][j] *= (double) 55555556.0; //conversion from kg/s m3 to mmol/MPa s m
//sat content (m3/m3) based on particle density of 2.65Mg/m3; p. 8, campbell
            		ws_[i][j] = ( 1.0 - ( bkd / 2.65 ) ); 
            		ws_[i][j] *= (double) 55555556.0; //conversion from m3/m3 to mmol/m3
//returns ln transform of r
            		r[i][j] = pow( (double) r[i][snmax] / rr, (double) (j-1) / (double)(snmax-1) ) * rr;
        	}
    	}

    	soilvol = 0.0;
    	for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
	{
        	rssoil[i] = 0.0;
        	for (j = 2; j <= SoilRhizElements + 1; j++ ) 
		{
//conductance factor: cf x soil conductivity = soil element conductance in v/t MPa
            		cf[i][j] = 2.0 * M_PI * l[i] / log( r[i][j] / r[i][j - 1] ); 
//upper half of the volume of subtending element
            		vu[i][j] = l[i] * M_PI * ( pow( r[i][j], 2.0 ) - pow( r[i][j - 1], 2.0 ) ) / 2.0; 
            		vl[i][j] = vu[i][j]; //lower half of the volume of subtending element
            		soilvol += 2.0 * vu[i][j]; //total soil volume in m3
            		rssoil[i] += 1.0 / ( cf[i][j] * ks[i][j] );
        	}
    	}

    	for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
	{
            	saturatedKs[i-2][3] = (double)1.0 / rssoil[i];
    	}

//assign capacitances and volumes to plant nodes
    	for (  i = 1; i <= ShootModules; i++ ) 
	{
//water content of leaves for 0.4 of plant water in foliage
        	ws_[i][1] = wsplant * 0.4 * al[i] / latot; 
//water content of stems, for 0.5 of plant water in stems
        	ws_[i][0] = wsplant * 0.5 / (double) ShootModules; 
//absolute capacitance in mmol/MPa for each plant node
        	cpp[i][1] = cpplant * 0.4 * al[i]; 
        	cpp[i][0] = cpplant * 0.5 * latot / ( (double) ShootModules ); //ditto for axial nodes
    	}

    	for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
	{
//water content of axial root segs, for 0.05 plant water in those roots
        	ws_[i][0] = wsplant * 0.05 / (double) RootModules; 
//water content of absorbing roots, for 0.05 plant water in those roots
        	ws_[i][1] = wsplant * 0.05 / (double) RootModules; 
        	cpp[i][0] = cpplant * 0.05 * latot / (double) RootModules;
        	cpp[i][1] = cpplant * 0.05 * latot / (double) RootModules;
    	}

//initializes plant water contents, makes them absolute (in mmol)
    	for (i = 1; i <= totmodules + 1; i++ ) 
	{
        	if ( i != ShootModules + 1 ) 
		{
            		vu[i][0] = 0.5; //upper half of element water
            		vl[i][0] = vu[i][0]; //lower half of element water
//wnu is "new" absolute water content of upper half of element volume
            		wnu[i][0] = ws_[i][0] * vu[i][0]; 
            		wu[i][0] = wnu[i][0]; //old and new water contents are initially the same
            		wnl[i][0] = ws_[i][0] * vl[i][0];
            		wl[i][0] = wnl[i][0];
            		ip[i][0] = p[i][0]; //sets initial pressure for plant nodes
//cf 1 used for plant because k's are already absolute conductances
            		cf[i][0] = 1.0; 
        	}
        	else 
		{
            		ip[i][0] = p[i][0];
        	}
    	}

//lateral elements
    	for (i = 1; i <= totmodules + 1; i++ ) 
	{
        	if ( i != ShootModules + 1 ) 
		{
            		vu[i][1] = 0.5; //upper half of element water
            		vl[i][1] = vu[i][1]; //lower half of element water
//wnu is "new" absolute water content of upper half of element volume
            		wnu[i][1] = ws_[i][1] * vu[i][1]; 
            		wu[i][1] = wnu[i][1]; //old and new water contents are initially the same
            		wnl[i][1] = ws_[i][1] * vl[i][1];
            		wl[i][1] = wnl[i][1];
            		ip[i][1] = p[i][1]; //sets initial pressure for plant nodes
            		cf[i][1] = 1.0; //cf 1 used for plant because k's are already absolute conductances
        	}
    	}

//initializes soil water contents and makes them absolute (in mmol)
    	for (i = ShootModules + 2; i <= totmodules + 1; i++ ) 
	{
        	for (j = 2; j <= SoilRhizElements + 1; j++ ) 
		{
            		if ( p[i][j - 1] > pe[i][j] ) 
			{
                		wnu[i][j] = ws_[i][j] * vu[i][j];
                		wu[i][j] = wnu[i][j];
            		}
            		else 
			{
                		wnu[i][j] = ws_[i][j] * (pow(pe[i][j] / p[i][j-1], b1[i][j])) * vu[i][j];
                		wu[i][j] = wnu[i][j];
            		}

            		if ( p[i][j] > pe[i][j] ) 
			{
                		wnl[i][j] = ws_[i][j] * vl[i][j];
                		wl[i][j] = wnl[i][j];
            		}
            		else 
			{
                		wnl[i][j] = ws_[i][j] * vl[i][j] * pow(( pe[i][j] / p[i][j]), b1[i][j]);
                		wl[i][j] = wnl[i][j];
            		}
        	}
    	}
}

bool HydraulicModel::verifyPowCondition(double base, double exponent)
{
        return ( (base!=0.0 && exponent>0.0)  );
}

bool HydraulicModel::verifyExpCondition( double exponent)
{
        return ( exponent <= 709.7827 && exponent >= -709.7827 );
}

//Loose with the coding here. Doesn't need a goto...
bool HydraulicModel::pgamma( double aa, double nn, double fac, double &x,
                      		double &gamma, double &gammp) 
{

    	bool error = false;
    	double gamser=0.0, gammcf=0.0, gln=0.0;

    	if ( x < ( aa + 1.0 ) ) 
	{
        	error = gser( aa, nn, gamser, x, gln ); //uses gser (series) to get  gamser
        	if( error ) 
		{
            		cout << "ERROR in pgamma\n";
            		error = true;
            		goto stop;
        	}
        	gammp = gamser;
    	}
    	else 
	{
        	//uses continued fraction to get gln and gamcf (the complement to gammser)
        	gcf( aa, fac, nn, gammcf, gln, x );
        	gammp = 1.0 - gammcf; //takes the complement
    	}
   	gamma = exp( gln );

	stop:
    		return error;
}


bool HydraulicModel::gser(double aa, double nn, double &gamser, double &x, double &gln) 
{
    	double ap=0.0, Sum=0.0, del=0.0, gammln=0.0;
    	int k=0;

    	gammap( aa, gammln, x );
    	gln = gammln;
    	if ( nn < 0.0 ) 
	{
        	gamser = 0.0;
    	}
    	else 
	{
        	ap = aa;
        	Sum = 1.0 / aa;
        	del = Sum;
        	k = 1;
        	do 
		{
            		ap += 1.0;
            		del *= nn / ap;
            		if ( del <= MINDOUBLE || del >= MAXDOUBLE ) 
			{
                    		cout << "FOUND ERR\n";
                    		return true;
            		}
            		Sum += del;
            		k++;
        	} while ( ( k != ITMAX ) && ( del >= Sum*EPS ) );
    	}

    	if ( k == ITMAX ) 
	{
        	gamser = -1000.0;
    	}
    	if ( del < Sum*EPS ) 
	{
        	gamser = Sum * exp( -nn + aa * log( nn ) - gln );
    	}

    	return false;
}

void HydraulicModel::gammap(double aa, double &gammln, double &x) 
{
        double tmp=0.0, ser=0.0;
	double tmp1;

    	x = aa - 1.0;
    	tmp1 = x + 5.5;
    	tmp = ( x + 0.5 ) * log( tmp1 ) - tmp1;
    	ser = 1.0 + 76.18009173 / (x + 1.0) - 86.50532033 / (x + 2.0) + 24.01409822 / (x + 3.0) - 
			1.231739516 / (x + 4.0) + 0.00120858003 / (x + 5.0) - 0.00000536382 / (x + 6.0);
    	gammln = tmp + log( STPP * ser ); //ln(gamma(aa))
}

void HydraulicModel::gcf( double aa, double fac, double nn, double &gammcf,
                  		double &gln, double &x) 
{
     	double a11, b11, a0, b0, gold, an=0.0; 
     	double ana=0.0, g=0.0, gammln=0.0, kk=0.0,anf=0.0;
     	int k=0, ttest=0;

     	a11 = nn;
     	b11 = 1.0;
     	a0 = 1.0;
     	b0 = 0.0;
     	gold = 0.0;

     	gammap( aa, gammln, x );
     	gln = gammln;

     	do 
	{
         	k++;
		kk += 1.0;
        	an = kk;
         	ana = an - aa;
         	a0 = ( a11 + a0 * ana ) * fac;
         	b0 = ( b11 + b0 * ana ) * fac;
         	anf = an * fac;
         	a11 *= nn * a0 + anf;
         	b11 *= nn * b0 + anf;
         	if ( a11 != 0.0 ) 
		{
             		fac = 1.0 / a11;
             		g = b11 * fac;
             		if ( (g-gold)/g < EPS ) 
			{
                 		ttest = 1;
             		}
             		else 
			{
                 		ttest = 0;
             		}
             		gold = g;
         	}
         	else 
		{
             		k = ITMAX;
         	}
     	} while ( k != ITMAX && ttest != 1 );

     	if ( k == ITMAX ) 
	{
         	gammcf = -1000.0;
     	}
     	if ( ttest == 1 ) 
	{
         	gammcf = exp( -nn + aa * log( nn ) - gln ) * g;
     	}
}


void HydraulicModel::gauss( int ShootModules, int SoilRhizElements, int tnode, double jmatrix[][NMAX],
                    		int col[], double ff[], int indxx[], double subtract[], double dp[]) 
{
    	int k=0, i=0, j=0;
    	double pivot=0.0, Sum=0.0;

    	for (k = 1; k <= tnode; k++ ) 
	{
//find column of first coefficient in row k
        	for (j = tnode; j >= 1; j-- ) 
		{
            		if ( jmatrix[k][j] != 0.0 )  
			{           // (jmatrix[k][j] < -.00001) && (jmatrix[k][j] > .00001 )
                  		col[k] = j;
            		}
        	}

 //go down that column to identify pivot elements
        	for (i = k; i <= tnode; i++ ) 
		{
            		if ( jmatrix[i][col[k]] != 0.0 ) 
			{
                		pivot = jmatrix[i][col[k]];
                    		for (j = col[k]; j <= tnode; j++ ) 
				{ //divide row by the pivot
                         		jmatrix[i][j] = jmatrix[i][j] / pivot;
                    		}
                		ff[i] /= pivot;
                		indxx[i] = 1; //memory of pivoting/subtracting rows
            		}
            		else 
			{
                		indxx[i] = 0;
            		}
        	}

        	for (j = col[k]; j <= tnode; j++ ) 
		{
            		subtract[j] = jmatrix[k][j]; //identify subtractors on row k
        	}

        	for (j = col[k]; j <= tnode; j++ ) 
		{
            		for (i = k + 1; i <= tnode; i++ ) 
			{
                  		if ( subtract[j] != 0.0 && indxx[i] != 0 ) 
				{
                      			jmatrix[i][j] -= subtract[j];
                  		}
            		}
        	}

        	for (i = k + 1; i <= tnode; i++ ) 
		{
            		if ( indxx[i] != 0 )
			{
                		ff[i] -= ff[k];
            		}
        	}
    	}

    	i = ShootModules * 2 + 1;
    	do 
	{
        	i = i + SoilRhizElements + 2;
        	dp[i] = 0.0;
        	ff[i] = 0.0;
    	} while (i <= tnode);

    	for (i = tnode-1; i >= 1; i--) 
	{		//backsubstitution to solve for dp's
        	Sum = jmatrix[i][i+1] * dp[i+1]; //sum coefficients and dp's to right of diagonal
        	for (j = i+2; j <= tnode; j++) 
		{
            		Sum += jmatrix[i][j] * dp[j];
        	}
        	dp[i] = ff[i] - Sum;
    	}
}
//end of hydraulicModel
