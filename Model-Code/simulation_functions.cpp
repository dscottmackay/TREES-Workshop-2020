//
//simulation_functions.cpp
//
//simulation_functions(): Implents the set of processes or calls to processes to simulation a single
//                        time step of the TREES model
//
//Author: D. Scott Mackay, from original code by Sudeep Samanta and David Roberts
//
//This function is called once per time step
//
//The iterative convergence of plant hydraulic and canopy variables is done in here
//Output variables are returned as a structure (should be made more object-oriented)
//
//One layer two compartment model
//stability correction included - SS Jun03
//units corrected, some area basis confusion to be clarified
//changed heat cap to MJm-2lai-1C-1 3/4 - not using
//changed to molar specific heats etc. (molar units upto Ec) - 042004
//changed do_pm to get zero transpiration in the radiation component
//for less that 0 radiation - 053104
//Coupled farquhar-katul model added - 05/05 DSM
//NEE calculations - 07/07 DSM
//Soil water dynamics - 05/09 DSM
//Evaporation from multple surface types - 02/10 DSM
//Integration of John Sperry's plant water balance model - 2009-11 DER & DSM
//Improved logic to allow for longterm simulations - 2014-4 DSM
//Improved integration of phenology based on GSI - 2014-10 DSM & PS
//Improved integration of root-soil modules - 2014-10 DSM
//Root phenology - 2015 DSM
//Rhizosphere N cycling - 2016 DSM
//
//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <random>
#include "constants.h"
#include "simulator2.h"

struct sim_out simulation_functions(
				bool silent,
                                int ts,  //time step
                                double Ecrit, //Rate of transpiration above which hydraulic failure occurs
                                double &Thresh,     //
                                State_Store& state, //current state
				double thetaSoil[], //current state of soil moisture
				double tempSoil[], //current state of soil temperature
                                const Time ti,     //current time
                                double u_ref,       //wind speed @ ref height, m/s
                                double t_ref,       //deg C @ ref ht
                                double D_ref,       //vapor pressure deficit @ ref ht, kPa
                                double precip,
                                double Qpar,        //umol m-2 s-1 above canopy, horizontal plane
                                double t_canopy,    //deg C inside canopy
                                double D,
				double p_air,	    //atmospheric pressure kPa
                                double co2_atm,     //PPM atmospheric CO2
                                double t_surface,       //Soil temperature at air-soil interface, deg. C
                                double t_soil,  //Soil temperature, deg. C
                                double t_root,  //Rooting depth average temperature, deg. C
                                double Zw,              //Water table depth, m
                                trees_params &treesParams,
                                vector<double> EcState,
                                vector<double> EcStateInner,
                                vector<double> &GSI_vector,
                                double Ec_,
                                int ModelStatus,
                                vector<double> LeafPsiInner,
                                vector<double> KpInner,
                                vector<int> InnerStep,
				BiogeochemicalCycles& bgc,
				HydraulicModel *hydraulicModel,
				HydraulicModel::outputStruct &ecrit_k_psi, 
				HydraulicModel::outputStruct &k_p_e, 
				double ecritOut[], double pc[], double klpred[],
                    		double ppredOut[], double epredOut[], double nodeFail[], 
				char nodeTyp[], double ***psi,
                    		double **rflux, double soilpsiavg[], double nts[], 
				double evap[], double kroot[],
                    		double axr[], double latr[], double kshoot[], 
				double dslat[], double dsax[], double drlat[],
                    		double drax[], double l[], double ksat[][MD], 
				double bsat[][MD], double ccsat[][MD],
                    		double newb[][MD], double newc[][MD], double ip[][ULAT], 
				double b[][ULAT], double b1[][ULAT],
                    		double n[][ULAT], double n1[][ULAT], double r[][ULAT], 
				double vu[][ULAT], double vl[][ULAT],
                    		double cf[][ULAT], double pe[][ULAT], double ws_[][ULAT], 
				double ks2[][ULAT], double p[][ULAT],
                    		double dpp[][ULAT], double jl[][ULAT], double wnu[][ULAT], 
				double wu[][ULAT], double wnl[][ULAT],
                    		double wl[][ULAT], double cpu[][ULAT], double cpl[][ULAT], 
				double cpp[][ULAT], double ku[][ULAT],
                    		double kl[][ULAT], double f[][ULAT], double cc[][2], 
				double psinode[][NMAX], double psimin[][NMAX],
                    		double jmatrix[][NMAX], double jmatrix2[][NMAX], double percent[], 
				double rssoil[], double rs[],
                    		double dp[], double ff[], int col[], int row[], int indxx[], 
				double subtract[], double pressure[],
                    		double plc[], double plcweib[], double rsquare[], double soilpsi[], 
				double al[], double ar[],
                    		double saturatedKs[][4], double &ptarg, double &e, double &rr, 
				double &rzw, double &rd, double &einc,
                    		int &SoilRhizElements, double &dt, double &gmd, double &gsd, 
				double &bkd, double &fs, double &fc,
                    		double &Fabs, double &cpplant, double &saxlat, double &raxlat, 
				double &rfract, double &kla,
                    		double &aral, double &latot, double &ratot, int &shootElements, 
				int &ktotal, double &soilpsimin,
                    		int &ShootModules, int &rootElements, int &RootModules, 
				int &totmodules, double &axShoot_b,
                    		double &axShoot_c, double &latShoot_b, double &latShoot_c, 
				double &axRoot_b, double &axRoot_c,
                    		double &latRoot_b, double &latRoot_c, double &plco, 
				double &pdecrement, double &sumy1,
                    		double &sumy2, double &sumprod, double &cincrement, 
				double &soilvol, double &rlateral,
                    		double &rLat_base, double &slateral, double &sLat_base, 
				double &raxial, double &pleafave,
                    		int &tnode, double fac, double &modelledKL, 
				int &HydraulicModelFailCond)
{
	bool ExitFunction = false;
	bool callRewet;
	int loopcounter=0;
	int MODELFAIL = 0;
	int callHydraulicModel = ts;
	double pm_Ec=0.0; //check that this is equivalent to Ec_t in this version DSM 
        double Et, Ec_t;    //total Ecanopy (mms-1), ground area
	double adj_phiJ_shd, adj_phiJ_sun;
	double r_sun=1.0; //multiplier in gvc (jarvis func)
        double r_shd=1.0; //multiplier in gvc (jarvis func)
        //double gs1_sun_;//sun leaf - before leaf boundl and photo
        //double gs2_sun_;//sun leaf - after leaf boundl and photo
        //double gs1_shd_;//shd leaf - before leaf boundl and photo
        //double gs2_shd_;//shd leaf - after leaf boundl and photo
        bool Converged = false; //only set to true when e>ecrit and the binary search initiates.
        double PsiCrit = 0.0;
	//double altitude = treesParams.altitude;
        double lat = treesParams.lat;
        double longi = treesParams.longi;
        double z_ref = treesParams.z_ref;
        double lai = treesParams.lai;
        double dead_lai = treesParams.dead_lai;
	double total_lai;
        double canopy_cover = treesParams.canopy_cover;
        double canopy_ht = treesParams.canopy_ht;
	double laiFullCanopyHeight = treesParams.laiFullCanopyHeight;
	double temp_ht = canopy_ht;
        //double fcloud = treesParams.fcloud;
        double l_angle = treesParams.l_angle;
        double canopy_e = treesParams.canopy_e;
        double fPAR_beam = treesParams.fPAR_beam;
        double fPAR_diff = treesParams.fPAR_diff;
        double alpha_PAR = treesParams.alpha_PAR;
        double alpha_NIR = treesParams.alpha_NIR;
        double omega = treesParams.omega;
        double p_crown = treesParams.p_crown;
        double d_factor = treesParams.d_factor;
        double zm_factor = treesParams.zm_factor;
        double zh_factor = treesParams.zh_factor;
        double Rd_mult = treesParams.Rd_mult;
        double Jmax_mult = treesParams.Jmax_mult;
        double thetaJ = treesParams.thetaJ;
        double phiJ_sun = treesParams.phiJ_sun;
        double phiJ_shd = treesParams.phiJ_shd;
        double Nleaf = treesParams.Nleaf;
        double N_fixed_proportion = treesParams.N_fixed_proportion;
        double Nrubisco = treesParams.Nrubisco;
        double Kc25 = treesParams.Kc25;
        double q10Kc = treesParams.q10Kc;
        double Ko25 = treesParams.Ko25;
        double q10Ko = treesParams.q10Ko;
        double act25 = treesParams.act25;
        double q10act = treesParams.q10act;
        double Gsref0 = treesParams.Gsref0;
        double delta = treesParams.delta;
        double theta_opt = treesParams.theta_opt;
        double optimal_soil_T = treesParams.optimal_soil_T;
        double growth_resp_proportion = treesParams.growth_resp_proportion;
        double resp_coef_root = treesParams.resp_coef_root;
        double resp_coef_stem = treesParams.resp_coef_stem;
        double resp_coef_leaf = treesParams.resp_coef_leaf;
        double resp_coefficient = treesParams.resp_coefficient;
        double Csoil = treesParams.Csoil;
        double Clitter = treesParams.Clitter;
        double Croot = treesParams.Croot;
        double Cstem = treesParams.Cstem;
//Get root depths from the plant module parameterization (see below)
	double Droot, cumDroot[ULAT];
	double rootDepth[ULAT];
	int rmodules;
        double litter_capacity = treesParams.litter_capacity;
        double SLA = treesParams.SLA;
        double porosity = treesParams.porosity;
        double ks = treesParams.ks;
        double bubbling_pressure = treesParams.bubbling_pressure;
        double pore_size_index = treesParams.pore_size_index;
        double residual = treesParams.residual;
        double fieldCapacity = treesParams.field_capacity;
        double LWP_spring_minimum = treesParams.LWP_spring_minimum;
        double LWP_stomatal_closure = treesParams.LWP_stomatal_closure;
        int is_bryophyte = treesParams.is_bryophyte;
        double capRiseScalar= treesParams.capRiseScalar;
	double drainScalar = treesParams.drainScalar;
	double melt_Rcoef = treesParams.melt_Rcoef;
	double Soil_Psi, Ysoil[ULAT];
        int year, yday, hour, min;
        double growth_resp_frac, rratio, psn, nsc, rmaint, rgrowth, rdefense, nsc_base;
	double snowpack_E_deficit, snowpack, snowmelt, sublimation, turbulent_melt;
	double sumRhizFlux, Rflux[ULAT];
	double D0, D0_sun, D0_shd;
	double ks_in[MD][MD], newb_in[MD][MD], newc_in[MD][MD], psimin_in[MD][MD];
	double totalRootArea, ar1, ar2, ar3;
	int mm;
	double N_avail_rate[ULAT], N_avail_rate_plant, N_neg_demand, N_pos_demand, N_neg_fract, N_pos, N_neg;
	double mineral_fraction, litter_store, canopy_store;
	double capFringe;
	double UP_neg[ULAT][10], UP_pos[ULAT][10];
	double LE_neg[ULAT][10], LE_pos[ULAT][10];
	double kratio, kratio_vector[ULAT], kratioRoot, nscRatio, rootClassRespRate;
	double CostLoading;
	double totalRootCarbon = 0.0;
	double Ci, Cm, Cc, gc0, gc, Pa, Nmax, Navail, NGsun, NSsun, NGshd, NSshd;

//set equal to 0 prior to determination of hydraulic parameters based on input of e.
//set equal to 1 after determination of hydraulic parameters but prior to determination of ecrit.
        int hydraulicModel_flag = 0;

//**********************************************************************************************
        double Ec_Store=0.0; //this variable holds Ec from PM for the Delta > Thresh test.
        double leafpsi = 0.0;
        double Kl = 0.0;
        int is = 1;
        int Step_State;
        bool NoFailureRecorded;

        if (ts == 1)
        {
                Step_State = 0; //if 0 then hydraulicModel has not yet been called
                NoFailureRecorded = true;
        }
        else
        {
                Step_State = 1;
                NoFailureRecorded = false;
        }
	double Delta;

        for(int iii=0; iii<50; iii++)
        {
            	InnerStep[iii]=0.0;
            	KpInner[iii]=0.0;
            	LeafPsiInner[iii]=0.0;
            	EcStateInner[iii]=0.0;
        }

        double Ec_sun, Ec_shd; //transpiration in the two parts
        double Ec_max; //max transpiration unconstrained by gs
	double A_sun, A_shd;  // assimilation in two parts
//values in the two compartments
        double Rabs_sun, Rabs_shd, Rnet_sun, Rnet_shd, Rth, lai_sun, lai_shd;
        double PAR_sun, PAR_shd; //added 042004
        double gvc_sun, gvc_shd, gv_sun, gv_shd;
//values for both compartments
        double rho_mol, ea, zenith_angle, daylen, t_refK, t_canK, t_soilK;
	double t_nearSurface, t_nearSurfaceK;
        double zeta, psi_m, psi_h;      //stability and diabatic correction factors
        double gr, gHa, gva, gHr, gv_t;
	double svp1, svp2, del_svp;
	//double gvb;  // leaf boundary layer conductance mol m-2 s-1
	//double alpha, l_m; //parameters for calculating leaf boundary layer conductance;
	double Kd;
	double R_leaf;
	double optimal_mineral_rate;
	double atm_e = 0.90;  //atmospheric emissivity
	double G = 0.0;   //ground heat flux
        double Qtotal = Qpar*2.12766*0.235; //par mumolm-2s-1 to total Wm-2
	double delta_Sr;  //change in soil moisture in the rooting zone
	double crownRadius, crownArea;

//structures for Farquhar (temporary)
	struct farqin in;
	struct farqout out;

//structure for simulation output
	struct sim_out simOut;

//values for respiration components
        double NEE, R_het, R_root, R_root_segment, R_bg, R_stem, R_total, A_tot, A_tot_kg;
	double RL, stressedLeafLifeSpan, unstressedLeafLifeSpan, leafCfraction, stemAllocation;

//values for soil evaporation
        double Dsoil, Dsurface, dZ, water_stress, evaporative_fraction; 
	double soil_evaporation1, soil_evaporation2, gv_soil, capRise; 
	double capRise_shallow, capRise_mid, capRise_deep, drain, bypassFlow[ULAT];
//variables for transpiration and canopy evaporation
	double waterStress[ULAT], evaporation, transpiration, canopy_evaporation;
	double can_evap_sun, can_evap_shd, eff_precip;
        double Rnet_soil, soil_e;
//values for wet canopy evaporation
	double interception, interception_capacity, canopy_wetness, canopy_store_max;
//display current status of the larget live LAI and stress on the the canopy
	ti.get_time(year, yday, hour, min);
        if ((hour+min) == 5.0 && treesParams.useHydraulics == true && !silent)
	{
		cout << endl << "-------- ";
		cout << year << ":" << yday;
		cout << " --------" << endl;
		cout << "target LAI = " << treesParams.live_lai;
		cout << "; canopy PLC = " << 100.0*(1.0-k_p_e.latK[1]/ksat[1][1]);
		cout << "; N status (mol N+ mol organic N) = " << bgc.plantNstatus[0] << endl;
	}


//***HERE IS WHERE MODIFIED N CYCLING STARTS ***
//To account for total N demand, assume plant takes up nitrate and ammonium in proportion 
//  to their availability
	bgc.getRhizosphereN(N_neg, N_pos);
	N_neg_fract = N_neg/(N_neg+N_pos);
	N_neg_demand = 0.0;
	N_pos_demand = 0.0;

//*****************************************************************************************
//PHENOLOGY
//This is called if useLeafModule parameter is set to false
//adjust lai based on phenological response
//uses potential LAI from the C allocation
//adjusts LAI based on Jolly et al (2005) GSI approach
//added Feb 2013 DSM and PRS
//integrated with carbon allocation March 2013 DSM
//*****************************************************************************************
	int i, daystep, dayindex;
	double delta_lai, delta_nsc;

	double lai_pot = state.get_val_at(CURRENTTARGETLAI);

	daylen = daylength(ti, lat, longi);

//compute relative lateral stem K
//if relative lateral K is under 50% then scale between 0 and 1
	RL = 1.0;
	if (treesParams.useHydraulics == true)
	{
		RL = k_p_e.latK[1]/ksat[1][1];
		if (RL < 0.5)
		{
			RL *= 2.0;
		}
		else
		{
			RL = 1.0;
		}
	}

//This routine works best if the simulations begins before leaf flush
	if (treesParams.usePhenology == true && treesParams.useLeafModule == false)
	{
		compute_GSI_LAI(ts, hour, min, yday, daylen, RL, D_ref, t_canopy, 
					N_neg_fract, N_neg_demand, N_pos_demand,
					bgc, state, simOut, treesParams, GSI_vector);
	}
//determine LAI dynmics directly from C allocation without imposing phenology
//most useful for an annual plant
	else if (treesParams.usePhenology == false && treesParams.useLeafModule == false)
	{
		delta_lai = treesParams.live_lai - treesParams.lai;
		nsc = bgc.getLeafNSC();
		bgc.updateLeafCarbonNitrogenPools(treesParams, delta_lai, RL,
					N_neg_fract, N_neg_demand, N_pos_demand);
	
		treesParams.lai = bgc.getLeafBiomassCarbon()*treesParams.SLA/10000.0;
		delta_nsc = nsc - bgc.getLeafNSC();
		simOut.leaf_growth = 0.86*delta_nsc/48.0;
		simOut.leaf_growth_respiration = 0.14*delta_nsc /48.0;
		treesParams.live_lai = treesParams.lai;
	}
	lai = treesParams.lai;

// ********************************************************************************************
// LEAF GROWTH MODULE - sequential leaf development
// use for modeling single leaves
// ********************************************************************************************
	if (treesParams.useLeafModule == true)
	{
                // if in stochastic leaf mode, then leaf parameters from gamma distributions
        
        if (treesParams.useLeafGamma == true && ts == 1) // update leaf param arrays with sampled values
                {
                    vector <double> limitvals;
                
                    //K parameter
                    limitvals = bgc.getQuantileVals(treesParams.Kalpha, treesParams.Kbeta, 0.1, 0.9);
                    for(i=0; i<500; ++i)
                    {
                        bgc.Karray[i]= bgc.sampleTruncatedGamma(treesParams.Kalpha, treesParams.Kbeta, limitvals[0], limitvals[1]);
                    }
                    
                    //N parameter
                    limitvals = bgc.getQuantileVals(treesParams.Nalpha, treesParams.Nbeta, 0.1, 0.9);
                    for(i=0; i<500; ++i)
                    {
                        bgc.Narray[i]= bgc.sampleTruncatedGamma(treesParams.Nalpha, treesParams.Nbeta, limitvals[0], limitvals[1]);
                    }
                    //r parameter
                    limitvals = bgc.getQuantileVals(treesParams.ralpha, treesParams.rbeta, 0.1, 0.9);
                    for(i=0; i<500; ++i)
                    {
                        bgc.rarray[i]= bgc.sampleTruncatedGamma(treesParams.ralpha, treesParams.rbeta, limitvals[0], limitvals[1]);
                    }
                }
            	double tm_current, delta_thermTime;
            	state.set_val_at(lai_pot, FORECASTTARGETLAI);

// thermal time tracker, plant
            	//double thermal_time;
            	tm_current= t_canopy;

            	delta_thermTime = max(0.0, (tm_current- treesParams.Tbase))*0.5 ;
            	treesParams.therm_plant += delta_thermTime;
            	bgc.phyllo_tracker[0] += delta_thermTime;
            	//thermal_time = treesParams.therm_plant;

// if at first ts, immediately initiate lf
            	if (ts == 1)
            	{
                	bgc.phyllo_tracker[0] = treesParams.phyllochron;
            	}

// checkpoint: enter repro phase?
            	if (treesParams.therm_plant >= treesParams.floweringTime)
             	{
                	cout << "Enter Reproductive Phase" << endl;
            	} 
// else continue in vegetative phase
		else
            	{
                	if (bgc.phyllo_tracker[0] >= treesParams.phyllochron)
                	{
                    		bgc.Lf_idx[0]++; // initiates new leaf
                    		bgc.putSingleLeafAreaPotential(bgc.Narray[bgc.Lf_idx[0]-1], bgc.Lf_idx[0]-1 );
                    		bgc.phyllo_tracker[0] = 0; // resets to 0

                	}

                	if (bgc.Lf_idx[0] > 0)
                	{
// Loop from youngest lf (Lf_idx) to oldest (k of 0 --> Lf_idx = 1) to prioritize youngest lf
                    		for (int k = (bgc.Lf_idx[0]-1); k >= 0; k-- )
                    		{
// updates any emerged leaf (Lf_idx indexs leaf number which starts at 1)
// co-opted code to get thetaRoot for updateLeaf function
                                	double thetaRootTmp = 0.0; 
                                	rmodules = treesParams.rmodules;
                                	for (i = 0; i < rmodules; i++)
                                	{
                                    		thetaRootTmp += thetaSoil[i]*treesParams.ar[i+3];
                                	}  
                                
                        			bgc.updateLeaf(k, delta_thermTime, treesParams, 
							treesParams.pot_size, RL,
							N_neg_fract, N_neg_demand, N_pos_demand, 
							bgc.Karray, bgc.Narray, bgc.rarray, 
							bgc.Lf_idx[0]-1, thetaRootTmp);

                        		//cout << "leaf " << k << " has potA of " << bgc.getSingleLeafAreaPotential(k) << endl;
                    		}
                	}
		}

// update plant-level variables for this ts
// plant LAI
            	treesParams.lai = bgc.calcLAI(bgc.Lf_idx[0], treesParams, bgc.Karray);

// plant-level SLA
            	if (bgc.Lf_idx[0] > 0 && bgc.getSingleLeafArea(0)> 0)
            	{
                	//treesParams.SLA = bgc.calcPlantSLA(bgc.Lf_idx[0], treesParams);
                	//cout << "plant level SLA at this ts is " << bgc.calcPlantSLA(bgc.Lf_idx[0], treesParams) << "\n";
                	//cout << "plant level SLA at this ts is " << treesParams.SLA << "\n";
            	}

// plant-level leaf biomass C
            	//cout << "leaf biomass c = " << bgc.getLeafBiomassCarbon() << endl;

// plant-level leaf biomass N
            	//cout << "leaf biomass N = " << bgc.getLeafBiomassN() << endl;

// total lf growth respiration
		treesParams.live_lai = treesParams.lai;

	}

	lai = treesParams.lai;


//********************************************************************************************
//SOIL - Update with precipitation and drainage
//get volumetric soil moisture fraction (proportion of soil volume)
//********************************************************************************************

	capFringe = bubbling_pressure * 0.01; //convert from cm to m

//New root parameterization - DSM 07/2013
//NEW
//Keep track of multiple root layers - OCT 2014
	rmodules = treesParams.rmodules;
	i = 0;
	rootDepth[i] = treesParams.drax[i+3];
	cumDroot[i] = rootDepth[i];
	Droot = rootDepth[i];
	for (i = 1; i < rmodules; i++)
	{
		rootDepth[i] = treesParams.drax[i+3];
		cumDroot[i] = cumDroot[i-1] + rootDepth[i];
		Droot += rootDepth[i];
	}

	mineral_fraction = 0.97;

//Force water table to just below the rooting depth
//This overrides the input Zw in the meteorological file
//  -- comment this out if your input Zw is reliable (e.g. from observations, groundwater model, etc)
	if (treesParams.useInputWaterTable == false)
	{
//		Zw = Droot-0.01;
//Add capillary fringe to water table depth to make Zw the top of the capillary fringe
		Zw = Droot+capFringe-0.1*rootDepth[rmodules-1];
//		Zw = Droot+0.01+capFringe;
//		Zw = Droot;
	}

//********************************************************************************************
//INTERCEPTION - calculate how much water is intercepted by the canopy
// DSM May 2010
//********************************************************************************************
	canopy_store_max = lai*treesParams.interception_per_leafArea;
	eff_precip = compute_effective_precipitation(state, treesParams, canopy_store_max, precip, t_ref);

//DRAINAGE: precipitation --> litter layer --> shallow soil layer --> deep soil layer
//
	infiltration(eff_precip, state, rootDepth, thetaSoil, bypassFlow, rmodules,
			litter_capacity, porosity, pore_size_index, fieldCapacity, 
			residual, capFringe, mineral_fraction, drainScalar, ks);

//***************************************************************************************
//CANOPY ELEMENTS
//***************************************************************************************
//calculate values for both compartments
        //p_air = air_pr(altitude); //we're currently reading this value in - DSM 9/28/11
        rho_mol = mol_density(p_air, t_ref);
        zenith_angle = zenith(ti, lat, longi);

//
//DSM October 2019
//Added to grow plant height for aerodynamic conductance calculations
//
	//canopy_ht *= min(lai,laiFullCanopyHeight) / max(0.01,laiFullCanopyHeight) + 0.000001;

	if (treesParams.usePhenology == 0 && treesParams.useLeafModule == 0 && lai < laiFullCanopyHeight)
	{
		canopy_ht *= lai/(laiFullCanopyHeight+0.000001);
	}

//stability calculations
	if (u_ref < 0.3)
	{
		u_ref = 0.3;
	}
	else if (u_ref > 15.0)
	{
		u_ref = 15.0;
	}
        zeta = stability_sucs(t_ref, t_canopy, z_ref, zm_factor, zh_factor, d_factor, 
								u_ref, canopy_ht, p_air);
        psi_m = calc_psim(zeta);
        psi_h = calc_psih(zeta);

//conductance values for both compartments
        gr = calc_gr(t_ref, canopy_e);
        gHa = calc_gHa(u_ref, z_ref, canopy_ht, d_factor, zm_factor, zh_factor, rho_mol, psi_m, psi_h);
        gva = gHa;
        gHr = gHa + gr;

//canopy LAI is concentrated where there is canopy cover, DSM 08/23/08
	//lai = lai /(canopy_cover+0.0000001);

//leaf boundary layer conductance - not currently use; 5/1/06 DSM
//	gvb = calc_gvb(zeta, d_factor, h_canopy, zm_factor, z_ref, u_ref, alpha, lai, d_leaf);
//temporary assumption of a leaf width of 8 cm
//this is worth looking into for CO2 diffusion, but for vapor diffusion we're covered 
//with the aerodynamic with stability corrections
/*
	l_m = pow(6*0.08*0.08*canopy_ht/(3.1417*lai),0.3333);
	alpha = sqrt(0.2*lai*canopy_ht/l_m);
	gvb = calc_gvb(zeta, d_factor, canopy_ht, zm_factor, z_ref, u_ref, alpha, lai, 0.056);
*/
 
//partition radiation and lai into sun and shade compartments
	if (state.get_val_at(REGEN) == 1.0)
	{
		total_lai = lai + dead_lai;
		if (dead_lai > 0.01)
		{
			Kd = cnpy_diff_ext(l_angle, total_lai, omega, p_crown);
        		absorb_rad(ti, Qtotal, dead_lai, lai, zenith_angle, Kd,
                        	l_angle, fPAR_beam, fPAR_diff, alpha_PAR, alpha_NIR, omega, p_crown,
                                PAR_sun, PAR_shd, Rabs_sun, Rabs_shd, lai_sun, lai_shd); //passed by ref
		}
		else
		{
			Kd = cnpy_diff_ext(l_angle, lai, omega, p_crown);
			absorb_rad(ti, Qtotal, lai, zenith_angle, Kd,
                        	l_angle, fPAR_beam, fPAR_diff, alpha_PAR, alpha_NIR, omega, p_crown,
                                PAR_sun, PAR_shd, Rabs_sun, Rabs_shd, lai_sun, lai_shd); //passed by ref
		}
		Rnet_sun = Rabs_sun;
		Rnet_shd = Rabs_shd;
	}
	else
	{
		Kd = cnpy_diff_ext(l_angle, lai, omega, p_crown);
        	absorb_rad(ti, Qtotal, lai, zenith_angle, Kd, l_angle, fPAR_beam, fPAR_diff, 
			   alpha_PAR, alpha_NIR, omega, p_crown, PAR_sun, PAR_shd,
                           Rabs_sun, Rabs_shd, lai_sun, lai_shd); //passed by ref
		Rnet_sun = Rabs_sun;
		Rnet_shd = Rabs_shd;
	}
 
//changed 032004
	//simOut.L_sun = lai_sun * canopy_cover;
	//simOut.L_shd = lai_shd * canopy_cover;
	simOut.L_sun = lai_sun;
	simOut.L_shd = lai_shd;
	//simOut.PAR_sun = PAR_sun * 0.47 / 0.235;
	//simOut.PAR_shd = PAR_shd * 0.47 / 0.235;
	PAR_sun *= 4.2; //PAR conversion from W m-2 to umol m-2 s-1, 2.25 x 10^5 J/mol quanta (C&N, 1998)
	PAR_shd *= 4.2; //PAR conversion from W m-2 to umol m-2 s-1, 2.25 x 10^5 J/mol quanta (C&N, 1998)
	simOut.PAR_sun = PAR_sun;
	simOut.PAR_shd = PAR_shd;

//
//Before first use of root areas (treesParams.ar[]) get updated root areas from root carbon
//
	totalRootArea = bgc.computeRootArea(treesParams);
	for (i = 0; i < rmodules; i++)
	{
		treesParams.ar[i+3] = bgc.computeRootArea(treesParams, i)/totalRootArea;
	}

//
//Soil_Psi is computed using a weighted mean based on proportional root areas by root layer
//
	double totArea = 0.0;
	Soil_Psi = 0.0;
        for (i = 0; i < rmodules; i++)
        {
                if (psimin[i+3][1] != -8.0)
                {
                        totArea += treesParams.ar[i+3];
                }
        }
	for (i = 0; i < rmodules; i++)
	{
		if (psimin[i+3][1] != -8.0)
                {
			Ysoil[i] = water_potential(porosity, bubbling_pressure, pore_size_index, residual, thetaSoil[i]);
			Soil_Psi += Ysoil[i]*treesParams.ar[i+3]/totArea;
		}
	}
	simOut.Soil_Psi = Soil_Psi;

//
//This is legacy code used to disconnect roots from the rhizosphere assuming only 3 layers
//Root areas are based on the leaf area proportions in the param_mod input file
//Note: This code will still do as intended for the first two layers
//
	ar1 = ar2 = ar3 = 0.0;
	if (treesParams.useHydraulics == true)
	{
//disconnect roots from soil
                if (treesParams.xylemScalar == 0.96)
                {
                        ar1 = treesParams.ar[3];
                        ar2 = treesParams.ar[4];
                        ar3 = 0.0;
                }
		else if (treesParams.xylemScalar == 0.97)
		{
			ar1 = ar2 = 0.0;
			ar3 = treesParams.ar[5];
		}
		else if (treesParams.xylemScalar == 0.98)
		{
			ar1 = 0.0;
			ar2 = treesParams.ar[4];
			ar3 = treesParams.ar[5];
		}
		else
		{
			ar1 = treesParams.ar[3];
			ar2 = treesParams.ar[4];
			ar3 = treesParams.ar[5];
		}
	}

//
//calculate water stress using previous time-step knowledge of Ecrit, as represented in E_inc
//Note: if plant water balance model is not used then this function is replaced with 
//      that using theta_gsref_reductio(), a simple hyperbolic function of water potential
//
	if (treesParams.useHydraulics == true)
	{
		water_stress = 200.0*treesParams.E_inc/treesParams.e_at_saturated_kl;
	}
	else
	{
		if (treesParams.useInputStress == true)
		{
			water_stress = treesParams.stressScalar;
		}
		else
		{
//Compute water stress by root layer for use when running non-hydraulic mode
			water_stress = 0.0;
			for (i = 0; i < rmodules; i++)
			{
				waterStress[i] = theta_gsref_reduction(porosity, bubbling_pressure, 
									pore_size_index, residual, 
									LWP_spring_minimum, 
									LWP_stomatal_closure, 
									thetaSoil[i], fieldCapacity);
				water_stress += waterStress[i]*rootDepth[i];
			}
			water_stress /= Droot;
			if (water_stress > 1.0)
			{
				water_stress = 1.0;
			}
		}
	}
	if (water_stress < 0.001)
	{
		water_stress = 0.001;
	}
	simOut.waterStress = water_stress;

//
//CANOPY EVAPORATION
//DSM May 2010
//

	canopy_evaporation = compute_canopy_evaporation(state, canopy_store_max, precip, canopy_wetness, Rnet_sun, Rnet_shd,
							simOut.T_sun, simOut.T_shd, D, p_air, gHr, gva, canopy_cover);
	simOut.canopy_evaporation = canopy_evaporation;

//calculate net absorbed radiation assuming that thermal radiation
//from canopy is partitioned in proportion to lai
        t_refK = C2K(t_ref);
        t_canK = C2K(t_canopy);
        t_nearSurface = 0.5*(t_surface+t_soil); // average soil temperature near surface
        t_soilK = C2K(t_soil);
        t_nearSurfaceK = C2K(t_nearSurface);
        soil_e = 0.96;

	simOut.T_sun = t_canopy;
	simOut.T_shd = t_canopy;

//Do coupled photosynthesis-transpiration
//Do sunlit leaves first
        in.pa = p_air*1000.0;   //Atmospheric pressure kPa -> Pa
        in.co2 = co2_atm;       //Ambient CO2 concentration (ppm)
        in.t = t_canopy;        //Canopy leaf temperature

//******************************************************************************************
//COUPLED WATER AND CARBON FLOW LOOP
//--> Insert HydraulicModel components and convergence loop here
// loop contains everything down to just past PM
//******************************************************************************************

//this forces the model to call the empirical calculation for gvc at start of a new timestep.
        ModelStatus = 0; 

//canopy solver works better if VPD canopy is not too small
	if (D < 0.1)
	{
		D0 = 0.1;
	}
	else
	{
		D0 = D;
	}
	D0_sun = D0_shd = D0;

//This will force the model to skip DoHydraulicModel()
	if (treesParams.useHydraulics == false)
	{
		Converged = true;
	}

	callHydraulicModel = 0;

//Save K's, P's, and Weibull's before entering loop
	for (mm = 1; mm <= totmodules+1; mm++)
	{
		if ( mm != ShootModules+1)  // skip over "hidden" module
		{
			ks_in[mm][1] = ks2[mm][1];
			ks_in[mm][0] = ks2[mm][0];
			psimin_in[mm][1] = psimin[mm][1];
			psimin_in[mm][0] = psimin[mm][0];
			newb_in[mm][1] = newb[mm][1];
			newb_in[mm][0] = newb[mm][0];
			newc_in[mm][1] = newc[mm][1];
			newc_in[mm][0] = newc[mm][0];
		}
	}

//this is where roots are disconnected
	if (ar1 == 0.0)
	{
		psimin_in[3][1] = -8.0;
	}
	if (ar2 == 0.0)
	{
		psimin_in[4][1] = -8.0;
	}
	if (ar3 == 0.0)
	{
		psimin_in[5][1] = -8.0;
	}
	r_sun = r_shd = 1.0;

//This is the main solver - iteratively solve canopy biochemistry and plant water balance
	do //coupled plant water balance and carbon assimilation
	{
// INITIAL GUESS AT STOMATAL CONDUCTANCE - FIRST ITERATION ONLY
            	if(ModelStatus == 0)
            	{
                	is = 0;
//canopy element stomatal conductance in mol s-1
                	gvc_sun = calc_gvc(lai_sun, D0_sun, Gsref0, delta)*water_stress*r_sun;
                	gvc_shd = calc_gvc(lai_shd, D0_shd, Gsref0, delta)*water_stress*r_shd;

			simOut.Gs_est = gvc_sun/lai_sun;
            	}

//Call hydraulic model until convergence or hydraulic failure - subsequent iterations
            	else
            	{
			if (EcStateInner[is] > 0.99*state.get_val_at(ECRIT))
			{
				EcStateInner[is] = 0.99*state.get_val_at(ECRIT);
			}

//Set the call to perform an update to Weibulls and plc at 5 am each day
//This is also done on the first time step
//Modify if you need to update more frequently
//Modified to call a refilling when xylemScalar changes
			ti.get_time(year, yday, hour, min);
			treesParams.useRefilling = false;
			callRewet = false;

			if (ts == 1)
			{
				callRewet = true;
				treesParams.useRefilling = true;
				callHydraulicModel = 0;
			}
			else if (fabs(treesParams.xylemScalar-state.get_val_at(XYLEM)) > 0.0 && 
					fabs(treesParams.xylemScalar-state.get_val_at(XYLEM)) < 0.01 )
			{
				if ((hour+min) == 5.0)
				{
					callRewet = true;
//This should ideally be replaced with treesParams.useRefilling
					callHydraulicModel = 0;
				}
			}
			else if (treesParams.forceRefilling == true && (double) (hour+min) == 3.0)
			{
				callRewet = true;
                                callHydraulicModel = 0;
                                treesParams.useRefilling = true;
                        }
			else if (treesParams.xylemScalar != state.get_val_at(XYLEM) && 
								treesParams.forceRefilling == true)
			{
				callRewet = true;
				callHydraulicModel = 0;
				treesParams.useRefilling = true;
			}
			else
			{
				callRewet = false;
				callHydraulicModel = 0;
				treesParams.useRefilling = false;
			}
//Reset to saved K's P's and Weibulls until the algorithm converges
        		for (mm = 1; mm <= totmodules+1; mm++)
        		{
                		if ( mm != ShootModules+1)  // skip over "hidden" module
                		{
                        		ks2[mm][1] = ks_in[mm][1];
                        		ks2[mm][0] = ks_in[mm][0];
                        		psimin[mm][1] = psimin_in[mm][1];
                        		psimin[mm][0] = psimin_in[mm][0];
                        		newb[mm][1] = newb_in[mm][1];
                        		newb[mm][0] = newb_in[mm][0];
                        		newc[mm][1] = newc_in[mm][1];
                        		newc[mm][0] = newc_in[mm][0];
                		}
        		}
//Call the soil-plant hydraulic model (Sperry et al 1998 PCE)
                    	DoHydraulicModel(silent, hydraulicModel, callHydraulicModel, callRewet,
			        	Ysoil, EcStateInner[is], Kl, leafpsi, Ecrit,
                                	HydraulicModelFailCond, PsiCrit, hydraulicModel_flag, 
					treesParams, ecrit_k_psi, k_p_e, ecritOut, 
					pc, klpred, ppredOut, epredOut, nodeFail, nodeTyp, 
					psi, rflux, soilpsiavg, nts, evap, kroot, axr, latr, 
					kshoot, dslat, dsax, drlat, drax, l, ksat, bsat,
                                        ccsat, newb, newc, ip, b, b1, n, n1, r, vu, vl, cf, pe, 
					ws_, ks2, p, dpp, jl, wnu, wu, wnl, wl, cpu, cpl, cpp, ku, kl, 
					f, cc, psinode, psimin, jmatrix, jmatrix2, percent, rssoil, 
					rs, dp, ff, col, row, indxx, subtract, pressure, plc, plcweib, 
					rsquare, soilpsi, al, ar, saturatedKs, ptarg, 
					e, rr, rzw, rd, einc, SoilRhizElements, dt, gmd,
                                        gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract, kla, 
					aral, latot, ratot, shootElements, ktotal, soilpsimin, 
					ShootModules, rootElements, RootModules, totmodules, 
					axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b,
                                        axRoot_c, latRoot_b, latRoot_c, plco, pdecrement, 
					sumy1, sumy2, sumprod, cincrement, soilvol, rlateral, 
					rLat_base, slateral, sLat_base, raxial,
                                        pleafave, tnode, fac, modelledKL);

//cout << "nodeFail = " << nodeFail[1] << endl;


//at low pre-dawn WP nighttime fluxes start to return zero Kl, PsiCit, or Ecrit
//even when the water balance model does not report a failure
//

//This if{} is suspect. Negative KL can occur at night when net flow 
			state.set_val_at(leafpsi, LPSI);
			if (HydraulicModelFailCond == 0)
			{
				if (isnan(Kl))
				{
					Kl = state.get_val_at(KL);
				}
				else if (Kl < 0.001)
				{
					Kl = state.get_val_at(KL);
				}
				else
				{
					state.set_val_at(Kl, KL);
				}
				if (isnan(leafpsi))
				{
					leafpsi = state.get_val_at(LPSI);
				}
				else if (leafpsi > -0.01) 
//No pressure, use previous state, try a lower stomatal conductance
				{
					leafpsi = state.get_val_at(LPSI);
				}
				else
				{
					state.set_val_at(leafpsi, LPSI);
				}
				if (Ecrit < 0.01)
				{
					Ecrit = state.get_val_at(ECRIT);
				}
				else
				{
					state.set_val_at(Ecrit, ECRIT);
				}
				if (PsiCrit == 0.0 || isnan(PsiCrit))
				{
					PsiCrit = state.get_val_at(PSICRIT);
				}
				else
				{
					state.set_val_at(PsiCrit, PSICRIT);
				}
			}

                    	if(HydraulicModelFailCond == 1)
                    	{
                        //	MODELFAIL = 1;
				if (isnan(Kl))
				{
					Kl = state.get_val_at(KL);
				}
                                else if (Kl < 0.001)
				{
                                        Kl = state.get_val_at(KL);
				}
                                else
				{
                                        state.set_val_at(Kl, KL);
				}
				if (isnan(leafpsi))
				{
					leafpsi = state.get_val_at(LPSI);
				}
                                else if (leafpsi > -0.01)
				{
                                        leafpsi = state.get_val_at(LPSI);
				}
                                else
				{
                                        state.set_val_at(leafpsi, LPSI);
				}
                                if (Ecrit < 0.01)
				{
                                        Ecrit = state.get_val_at(ECRIT);
				}
                                else
				{
                                        state.set_val_at(Ecrit, ECRIT);
				}
				if (PsiCrit == 0.0 || isnan(PsiCrit))
				{
					PsiCrit = state.get_val_at(PSICRIT);
				}
				else
				{
					state.set_val_at(PsiCrit, PSICRIT);
				}

//Plant water balance failure, use previous state, try a lower stomatal conductance
				cout << "*** Hydraulic failure - attempting to recover ***\n";
                    	}
		}

//calculate initial canopy stomatal conductance - hydraulics, no light response
            	if(ModelStatus == 1)
            	{
                	LeafPsiInner[is] = leafpsi;
                	KpInner[is] = Kl;
                	InnerStep[is] = is;

//calculate conductance using a combination of Darcy's Law and Fick's Law. returns in mol s-1
//the sign on the gravitational potential gradient is opposite to tension gradient (down vs up)
                	gvc_sun = calc_gvc(Kl, Soil_Psi, leafpsi, lai, temp_ht, treesParams, 
						pm_Ec, D0_sun, t_canopy, p_air)*lai_sun;
                	gvc_shd = calc_gvc(Kl, Soil_Psi, leafpsi, lai, temp_ht, treesParams, 
						pm_Ec, D0_shd, t_canopy, p_air)*lai_shd;

			temp_ht = 0.0;
            	}

//DSM - continues with TREES original elements, noting that some of these calculations 
//      could be moved outside the loop
//Now calculating gv values just before call to Penman-Monteith function

//Note: The iterative solution includes leaf temperature, photosynthesis, conductance, and hydraulics, 
//      and not extensively tested under a variety of conditions for convergence. 
//	However, convergence is fast at low fluxes.


//Adjust sun leaf temperature based on convective-radiative energy exchange
//Now calculates temperature at the leaf surface
//Modified energy balance for sun leaf element - 
//    	solve temperature, then long-wave emission, then temperature 
//(DSM May 2010)
        	Rnet_sun = Rabs_sun - G;
//Sunlit leaf energy balance has numerical problems at high solar zenith angle (horizon)
		if (lai_sun >= 0.5)
		{
			gv_t = (gvc_sun/lai_sun*gva)/(gvc_sun/lai_sun+gva)/42.0;
        		in.t = t_canopy + (Rnet_sun/lai_sun-lheat_vap(t_canopy)*gv_t*D/p_air)/
					(cp_air*gHr + lheat_vap(t_canopy)*calc_s(t_canopy, p_air)*gv_t);
			t_canK = C2K(in.t);
//emission from canopy (W m-2)
			Rth = canopy_e*sigma*pow(t_canK, 4.0)-atm_e*sigma*pow(t_refK, 4.0); 
        		Rnet_sun = Rnet_sun - Rth*lai_sun;
        		in.t = t_canopy + (Rnet_sun/lai_sun-lheat_vap(t_canopy)*gv_t*D0_sun/p_air)/
					(cp_air*gHr + lheat_vap(t_canopy)*calc_s(t_canopy, p_air)*gv_t);
		}
		else
		{
			in.t = t_canopy;
		}
		simOut.T_sun = in.t;

//Adjust shade leaf temperature based on convective-radiative energy exchange
//Now calculates temperature at leaf surface
        	Rnet_shd = Rabs_shd - G;
		if (lai_shd >= 0.5)
		{
			gv_t = (gvc_shd/lai_shd*gva)/(gvc_shd/lai_shd+gva)/42.0;
        		in.t = t_canopy + (Rnet_shd/lai_shd-lheat_vap(t_canopy)*gv_t*D/p_air)/
					(cp_air*gHr + lheat_vap(t_canopy)*calc_s(t_canopy, p_air)*gv_t);
			t_canK = C2K(in.t);
//emission from canopy (W m-2)
			Rth = canopy_e*sigma*pow(t_canK, 4.0)-atm_e*sigma*pow(t_refK, 4.0); 
        		Rnet_shd = Rnet_shd - Rth*lai_shd;
        		in.t = t_canopy + (Rnet_shd/lai_shd-lheat_vap(t_canopy)*gv_t*D0_shd/p_air)/
					(cp_air*gHr + lheat_vap(t_canopy)*calc_s(t_canopy, p_air)*gv_t);
		}
		else
		{
			in.t = t_canopy;
		}
		simOut.T_shd = in.t;

//additional inputs to photosynthesis
//LEAF N
//N is made available in proportion to its concentration in organic form in the leaf
//and inorganic form in the leaves and fine roots
//
//Allow for fixed (biologically or otherwise made available in stored) N sources
//If this is set to 0 then the dependency of photosynthesis on N is entirely a function of its allocation to leaf biomass
//    if set to 1, then photosynthetic dependency on N is based on bgc.getLeafBiomassN()
		if (N_fixed_proportion < 0.0)
		{
			N_fixed_proportion = 0.0;
		}
		else if (N_fixed_proportion > 1.0)
		{
			N_fixed_proportion = 1.0;
		}

		N_avail_rate_plant = bgc.plantNstatus[0];

//proportion of Nleaf determined by its biomass concentration
//bgc.getLeafBiomassN() is converted from kg ha-1 to kg m-2 ground area
//divide by lai (m2 leaf area m-2 ground area) to get kg m-2 leaf area
		double leafBiomassN, leafBiomassC, leafArea;;
		leafBiomassN = 0.0001 * bgc.getLeafBiomassN();
		leafBiomassC = 0.0001 * bgc.getLeafBiomassCarbon();
		leafArea = leafBiomassC * treesParams.SLA_min;

		Nleaf = (1.0-N_fixed_proportion) * leafBiomassN/max(0.000000001,leafArea);
		Nleaf += N_fixed_proportion*treesParams.Nleaf;
//cout << leafArea << '\t' << leafBiomassN << '\t' << leafBiomassN/max(0.000000001,leafArea) << '\t' << Nleaf*Nrubisco*fnr*treesParams.act25*1.0e6 / 60.0<< endl;
		in.lnc = Nleaf;
//shift more available N to rubisco when N is limiting
		in.flnr = Nrubisco;

//adjust quantum yield based on N status
		//adj_phiJ_sun = phiJ_sun*in.lnc/treesParams.Nleaf;
		//adj_phiJ_shd = phiJ_shd*in.lnc/treesParams.Nleaf;
//adjust quantum yield based on water stress
		adj_phiJ_sun = phiJ_sun*min(1.0, water_stress);
		adj_phiJ_shd = phiJ_shd*min(1.0, water_stress);
adj_phiJ_sun = phiJ_sun;
adj_phiJ_shd = phiJ_shd;

//reset lnc to Nleaf, if you think Vcmax is unaffected by environmental conditions
	//	in.lnc = Nleaf;

//Assume PAR = 0.47 * Absorbed Radiation, 0.235J/umol
        	//in.irad = (Rabs_sun/lai_sun) * 0.47 / 0.235;  
        	in.irad = PAR_sun;  
//Molar conversion for H2O to CO2 and mol/s to m/s - w/o aerodynamic
 		in.g = gvc_sun / lai_sun / 1.6 / 42.0;	
//In case you want to include leaf scale boundary-layer conductance, here you go:
//Molar conversion for H2O to CO2 and mol/m2/s to m/s
// 		in.g = 1.0/(1.0/(gvc_sun/1.6/lai_sun)+1.0/(gvb/lai))/42.0;	

		in.Rd = 0.0;	//temporarily set dark respiration to zero (umol/m2/s)
		in.D = D*1000.0;	//convert VPD from kPa to Pa

//calculate sun leaf photosynthesis, intercellular CO2, and stomatal conductance
		in.t = simOut.T_sun;
		NGsun = NGshd = NSsun = NSshd = 0.0;
		Pa = p_air * 10.0; //atmospheric pressure in bars
		if (treesParams.ps_model == 1 || treesParams.ps_model == 2) //C3 photosynthesis
		{
			if (treesParams.ps_model == 1) //use the original TREES photosynthesis model (pre-2019)
			{
				A_sun = bgc.photosynthesis(treesParams, Jmax_mult, thetaJ, adj_phiJ_sun, in, out);
				A_sun = (1.0-canopy_wetness)*A_sun;	
				gvc_sun = out.g * 1.6 * 42.0; //Molar conversion for CO2 to H2O and m/s to mol/m2/s
			}
			else //use new, experimental C3 model
			{
//available nitrogen in umol N m-2 leaf = kgN ha-1 * ha m-2 ground / (m2 leaf m-2 ground)* g kg-1 / gN mol-1 * umol mol-1
				Navail = bgc.getLeafMineralN() / 10000 / (lai_sun + lai_shd) * 1000.0 / 14.0067 * 1.0e6;
//maximum do novo nitrogen supply to the chloroplast, umol N m-2 30 min-1
				Nmax = min(treesParams.Nmax * 1800.0, Navail) / 30.0;
				gc0 = gvc_sun / lai_sun / 1.6; //molar conversion from H2O to CO2
				bgc.coupledA3_gc(treesParams, in.co2, Pa, in.irad, adj_phiJ_sun, in.t, gc0, Nmax, gc, A_sun, Ci, Cc, NGsun, NSsun);
				gvc_sun = 1.6 * gc; //molar conversion from CO2 to H2O
				out.Ci = Ci;
				out.Rd = 0.0;
			}
		}
		else //C4 photosynthesis
		{
			gc0 = gvc_sun / lai_sun / 1.6; //molar conversion from H2O to CO2
			if (treesParams.isAmphistomatous == true)
			{
				gc0 *= 0.5;
			}
			bgc.coupledA4_gc(treesParams, in.co2, in.irad, phiJ_sun, in.t, gc0, gc, A_sun, Ci, Cm);
			gvc_sun = 1.6 * gc; //molar conversion from CO2 to H2O
			if (treesParams.isAmphistomatous == true)
			{
				gvc_sun *= 2.0;
				gc0 *= 2.0;
				bgc.coupledA4_gc(treesParams, in.co2, in.irad, phiJ_sun, in.t, gc0, gc, A_sun, Ci, Cm);
				
			}
			out.Ci = Ci;
			out.Rd = 0.0;
		}

		if (gvc_sun < 0.00001*Gsref0)
		{
			gvc_sun = 0.00001*Gsref0;
		}
		if (gvc_sun > 10.0*Gsref0*water_stress)
		{
			gvc_sun = 10.0*Gsref0*water_stress;
		}
		gvc_sun = gvc_sun * lai_sun; //convert to mol/s

		simOut.G_sun = gvc_sun;
		simOut.A_sun = A_sun;
		simOut.Ci_sun = out.Ci;
		simOut.Rd_sun = out.Rd;
		simOut.Vcmax25 = out.Vcmax25;
		simOut.Vcmax_sun = out.Vcmax;
		simOut.Jmax25 = out.Jmax25;
		simOut.J_sun = out.J;

//calculate shd leaf photosynthesis, intercellular CO2, and stomatal conductance
//Assume PAR = 0.47 * Absorbed Radiation, 0.235 J/umol
        	//in.irad = (Rabs_shd/lai_shd) * 0.47 / 0.235;  
        	in.irad = PAR_shd;  
//Molar conversion for H2O to CO2 and mol/s to m/s
  		in.g = gvc_shd / lai_shd / 1.6 / 42.0;	 
//In case you want to include boundary-layer conductance, here you go:
//Molar conversion for H2O to CO2 and mol/m2/s to m/s
// 		in.g = 1.0/(1.0/(gvc_shd/1.6/lai_shd)+1.0/(gvb/lai))/42.0;	

//calculate shade leaf photosynthesis, intercellular CO2, and stomatal conductance
		in.t = simOut.T_shd;
		if (treesParams.ps_model == 1 || treesParams.ps_model == 2) //C3 photosynthesis
		{
			if (treesParams.ps_model == 1) //use the original TREES photosynthesis model (pre-2019)
			{
				A_shd = bgc.photosynthesis(treesParams, Jmax_mult, thetaJ, adj_phiJ_shd, in, out);
				A_shd = (1.0-canopy_wetness)*A_shd;	
				gvc_shd = out.g * 1.6 * 42.0; //Molar conversion for CO2 to H2O and m/s to mol/m2/s
			}
			else //use new, experimental C3 model
			{
//available nitrogen in umol N m-2 leaf = kgN ha-1 * ha m-2 ground / (m2 leaf m-2 ground)* g kg-1 / gN mol-1 * umol mol-1
				Navail = bgc.getLeafMineralN() / 10000 / (lai_sun + lai_shd) * 1000.0 / 14.0067 * 1.0e6;
//maximum do novo nitrogen supply to the chloroplast, umol N m-2 30 min-1
				Nmax = min(treesParams.Nmax * 1800.0, Navail) / 30.0;
				gc0 = gvc_shd / lai_shd / 1.6; //molar conversion from H2O to CO2
				bgc.coupledA3_gc(treesParams, in.co2, Pa, in.irad, adj_phiJ_shd, in.t, gc0, Nmax, gc, A_shd, Ci, Cc, NGshd, NSshd);
				gvc_shd = 1.6 * gc; //molar conversion from CO2 to H2O
				out.Ci = Ci;
				out.Rd = 0.0;
			}
		}
		else //C4 photosynthesis
		{
			gc0 = gvc_shd / lai_shd / 1.6; //molar conversion from H2O to CO2
			if (treesParams.isAmphistomatous == true)
			{
				gc0 *= 0.5;
			}
			bgc.coupledA4_gc(treesParams, in.co2, in.irad, phiJ_shd, in.t, gc0, gc, A_shd, Ci, Cm);
			gvc_shd = 1.6 * gc; //molar conversion from CO2 to H2O
			if (treesParams.isAmphistomatous == true)
			{
				gvc_shd *= 2.0;
				gc0 *= 2.0;
				bgc.coupledA4_gc(treesParams, in.co2, in.irad, phiJ_shd, in.t, gc0, gc, A_shd, Ci, Cm);
			}
			out.Ci = Ci;
			out.Rd = 0.0;
		}
		if (gvc_shd < 0.0001*Gsref0)
		{
			gvc_shd = 0.0001*Gsref0;
		}
		if (gvc_shd > 10.0*Gsref0*water_stress)
		{
			gvc_shd = 10.0*Gsref0*water_stress;
		}
		gvc_shd = gvc_shd * lai_shd;	//convert to mol/s

		simOut.G_shd = gvc_shd;
		simOut.Vcmax_shd = out.Vcmax;
		simOut.J_shd = out.J;

//output element stomatal conductance
//
		simOut.gs_sun = gvc_sun/lai_sun;
		simOut.gs_shd = gvc_shd/lai_shd;

//calculate canopy average stomatal conductance, Gs
		simOut.Gs = (simOut.gs_sun*lai_sun+simOut.gs_shd*lai_shd)/lai;

		if (ModelStatus == 0)
		{
			simOut.Gs_est = simOut.Gs;
		}

		simOut.A_shd = A_shd;
		simOut.Ci_shd = out.Ci;
		simOut.Rd_shd = out.Rd;

		simOut.A = (A_sun*lai_sun + A_shd*lai_shd)/lai;

//calculate GPP 
		A_tot = canopy_cover*(A_sun*lai_sun + A_shd*lai_shd);
//cout << ModelStatus << "; A_tot = " << A_tot << "; Gs = " << simOut.Gs << endl;


//TRANSPIRATION
//calculate vapor conductances
//For amphistomatous leaves (stomata on both sides), boundary layer conductance (gva) 
//  applied to both sides above
			
		if (treesParams.isAmphistomatous == true)
		{
			gv_sun = 2.0 * (1.0/(1.0/(0.5*gvc_sun) + 1.0/gva));
			gv_shd = 2.0 * (1.0/(1.0/(0.5*gvc_shd) + 1.0/gva));
		}
		else
		{
			gv_sun = 1.0/(1.0/gvc_sun + 1.0/gva);
			gv_shd = 1.0/(1.0/gvc_shd + 1.0/gva);
		}

//call P-M
		Ec_sun = (1.0-0.5*canopy_wetness)*canopy_cover*
				do_pm(Rnet_sun, simOut.T_sun, 0.5*(D0+D0_sun), p_air, gHr, gv_sun);
		Ec_shd = (1.0-0.5*canopy_wetness)*canopy_cover*
				do_pm(Rnet_shd, simOut.T_shd, 0.5*(D0+D0_shd), p_air, gHr, gv_shd);
		if (ModelStatus == 0)
		{
			Ec_max = do_pm(Rnet_sun, simOut.T_sun, 0.5*(D0+D0_sun), p_air, gHr, gva) + 
				     do_pm(Rnet_shd, simOut.T_shd, 0.5*(D0+D0_shd), p_air, gHr, gva);
			Ec_max *= canopy_cover;
			Ec_max = EcConversion_TtoS(Ec_max,treesParams.lai);
		}

//throwing out negative values because other mechanisms will prevent this - Future mod
		if (Ec_max < 0.0) 
		{
			Ec_max = 0.0;
		}
		if (Ec_sun < 0.0) 
		{
			Ec_sun = 0.0;
		}
		if (Ec_shd < 0.0) 
		{
			Ec_shd = 0.0;
		}
		pm_Ec = Ec_sun + Ec_shd;

//Note: Ec_ goes into HydraulicModel routine where it is multiplied by treesParams.lai
		Ec_ = EcConversion_TtoS(pm_Ec, treesParams.lai);
		is++;

		EcStateInner[is] = Ec_;
        	if(ModelStatus == 0)
        	{
            		//force loop
            		ModelStatus = 1;
            		Delta = Ec_max + 1.0;
        	}
        	else
        	{
//Ec_Store is the E from the previous PM calculation
            		Delta = fabs(EcStateInner[is] - Ec_Store);

//choose one of these statements to adjust how iterations are handled
//convergence rarely takes more than two iterations (initial, full hydraulic)
			//if (Soil_Psi >= 1.1*treesParams.pd_at_sat_kl || D < 1.0 || A_sun < 0.01)
			//if ((Ec_ < 0.99*Ecrit && leafpsi > 0.99*PsiCrit) || D < 0.5 || A_sun < 0.05)
			if ((Ec_ < 0.99*Ecrit && leafpsi > 0.99*PsiCrit))
			{
				Converged = true;
			}
        	}
        	Ec_Store = Ec_;

//This forces the model to converge after a set number of iterations.
//Depending on the parameters, the model does not always converge.
//I rarely see more than 5 iterations and usually only after a wet-up
//Typically, if there is a need for more than 5 iterations, then there will be no convergence
        	loopcounter++;
        	if(loopcounter > treesParams.max_iterations)
        	{
            		Delta = 0.0;
            		Converged = true;
        	}

        	if(ExitFunction == true)
        	{
            		Delta = 0.0;
            		Converged = true;
        	}

//various conditions that signal hydraulic failure
//in this case, we set the Gsres0 scalars (r_sun, r_shd) lower and go again
//to attempt to prevent hydraulic failure
//comment out this block if you don't like this
		if ((is > 1) && (PsiCrit >= treesParams.pd_at_sat_kl || 
					isnan(leafpsi) || 
					Ec_ > Ecrit ||
					HydraulicModelFailCond == 1))
		{
			ModelStatus = 0;
			Converged = false;
			Delta = 0.05*Ec_max;
			r_sun *= 0.95;
			r_shd = r_sun;
			loopcounter--;
			HydraulicModelFailCond = 0;
		}

//Re-compute aerodyanmic conductance to account for change in leaf temperatures
//This calculation reduces the number of loop iterations needed in some situations
//stability calculations
		if (Rabs_sun > 50.0 && lai >= 0.1)
		{
			svp1 = sat_vapor_pressure(t_canopy);
			svp2 = sat_vapor_pressure(simOut.T_sun);
			del_svp = svp2-svp1;
			D0_sun = D0 + del_svp;
			if (D0_sun < 0.1)
			{
				D0_sun = 0.1;
			}
                        svp2 = sat_vapor_pressure(simOut.T_shd);
                        del_svp = svp2-svp1;
                        D0_shd = D0 + del_svp;
			if (D0_shd < 0.1)
			{
				D0_shd = 0.1;
			}
			t_canopy = (simOut.T_sun*lai_sun+simOut.T_shd*lai_shd)/lai;

                	zeta = stability_sucs(t_ref, t_canopy, z_ref, zm_factor, zh_factor, d_factor, 
										u_ref, canopy_ht, p_air);
                	psi_m = calc_psim(zeta);
                	psi_h = calc_psih(zeta);
//conductance values for both compartments
                	gr = calc_gr(t_ref, canopy_e);
                	gHa = calc_gHa(u_ref, z_ref, canopy_ht, d_factor, zm_factor, zh_factor, 
										rho_mol, psi_m, psi_h);
                	gva = gHa;
                	gHr = gHa + gr;
		}
		simOut.D0_sun = D0_sun;
		simOut.D0_shd = D0_shd;
//THIS IS THE CRITERIA REQUIRED FOR THE MODEL TO CONTINUE CALLING THE PLANT HYDRAULICS******//
        } while ((Delta > Ec_max*0.01) && (Converged == false));

//update E increment based on current Ecrit, to use higher resolution as Ecrit declines
	treesParams.E_inc = 0.005*Ecrit;

//****************************************************************************//
//AFTER ALL ITERATIONS OF THE SPERRY MODEL ARE COMPLETE, OUTPUT TO SIMOUT AND
//CONTINUE WITH THE FUNCTION.
//This is mostly a lot of clean-up, generally looking for graceful ways of handling
//  results from the hydraulic model, which can sometimes be messy

        EcState[ts] = EcStateInner[is-1];

//store current state of whole-plant Kleaf
	state.set_val_at(Kl, KL);

//no prediction of transpiration, leaf water potential, or conductivity.
        if(ExitFunction == true)
        {
            simOut.Ecrit = Ecrit;
            simOut.Ec = Ec_;
            simOut.Psi_Crit = PsiCrit;
            simOut.WPlant_K = Kl;
            simOut.Leaf_Psi = leafpsi;
            simOut.OutFlag = 0;
        }
//hydraulic failure
        else if(MODELFAIL == 1)
        {
            simOut.Ecrit = Ecrit;
            simOut.Ec = Ec_;
            simOut.Psi_Crit = PsiCrit;
            simOut.WPlant_K = Kl;
            simOut.Leaf_Psi = leafpsi;
            simOut.OutFlag = 1;
//tree died, so stop simulation, or comment out and TREES will try again with a lower E
	    //exit(1); 
        }
//model did not converge
        else if(loopcounter > 50)
        {
            simOut.Ecrit = Ecrit;
            simOut.Ec = Ec_;
            simOut.Psi_Crit = PsiCrit;
            simOut.WPlant_K = Kl;
            simOut.Leaf_Psi = leafpsi;
            simOut.OutFlag = 2;
        }
//success
        else
        {
            simOut.Ecrit = Ecrit;
            simOut.Ec = Ec_;
            simOut.Psi_Crit = PsiCrit;
            simOut.WPlant_K = Kl;
            simOut.Leaf_Psi = leafpsi;
            simOut.OutFlag = 3;
        }
	Ec_t = canopy_cover*EcConversion_StoT(Ec_, treesParams.lai);
	if (Ec_t < 0.0)
	{
		Ec_t = 0.0;
	}
	simOut.Ec_t = Ec_t;

//note these are relative terms, and proportional to their contributions to transpiration
//negative numbers indicate flow from root to rhizosphere
//we'll calculate the absolute rhizosphere fluxes for each root module later
//
	for (i = 0; i < rmodules; i++)
	{
		if (treesParams.useHydraulics == true)
		{
			simOut.rhizFlux[i] = k_p_e.rhizFlux[i+2]*ar[i+3]*lai*treesParams.canopy_cover;
		}
		else
		{
			simOut.rhizFlux[i] = treesParams.ar[i+3]*lai*waterStress[i]/water_stress;
		}
	}

//we'll write out the updated Weibull parameters, b and c at 1pm each day
	if ((hour+min) == 13.0 && treesParams.useHydraulics == true && !silent)
	{
		cout << "Plant hydraulic status (Ks and Weibulls):" << endl;
		cout << "\tShoot:\t" << "axial Ks=" << ks2[1][0] << ", b=" << b[1][0] << 
				", c=" << cc[1][0] << "; lateral Ks=" << ks2[1][1] << 
				", b=" << b[1][1] << ", c=" << cc[1][1] << endl;
		for (i = 0; i < rmodules; i++)
		{
			cout << "\tRoot(" << i << ")\t" << "axial Ks=" << ks2[i+3][0] << 
				", b=" << b[i+3][0] << ", c=" << cc[i+3][0] << "; lateral Ks=" << 
				ks2[i+3][1] << ", b=" << b[i+3][1] << ", c=" << cc[i+3][1] << endl;
		}
	}

//
//EVAPORATION
//This section addresses evaporation from soil or non-vascular plant (bryophyte) surfaces
//Evaporation rates are controlled by:
//  1) radiation and vapor pressure deficit
//  2) conductance based on Darcy's Law and moisture content
//  3) feedback between surface soil moisture and VPD
//
	soil_evaporation1 = soil_evaporation2 = 0.0;

//calculate net absorbed radiation by soil surface under canopy
//emission from canopy (W m-2)
        Rth = soil_e*sigma*pow(t_nearSurfaceK, 4.0)-canopy_e*sigma*pow(t_canK, 4.0); 

//calculate VPD at the soil surface -- uses ambient atm absolute humidity and soil temperature
	ea = sat_vapor_pressure(t_ref) - D;
	Dsoil = sat_vapor_pressure(t_nearSurface) - ea;
	Dsurface = sat_vapor_pressure(t_surface) - ea;
	ks = ks/100.0/3600.0;  //convert from cm/hour to m/s

//allow evaporation from litter, soil, etc
	snowpack = state.get_val_at(SNOWPACK);
	litter_store = state.get_val_at(LITTER);
	if (snowpack == 0.0) 
	{
		evaporative_fraction = (thetaSoil[0]-residual)/(porosity-residual);
//Assume that albedo ranges from 0.1 (wet soil) to 0.3 (dry soil)
        	Rnet_soil = (1.0-(0.3-0.2*evaporative_fraction)) * (Qtotal - Rabs_sun - Rabs_shd) - Rth;

//Calculate soil water flow conductance to limit the rate of evaporation from the soil surface
		dZ = 0.5*rootDepth[0];
		if (is_bryophyte == 0) //This version is for mineral and/or organic soils
		{
        		gv_soil = soil_cond_rate(dZ, Dsoil, t_nearSurface, residual, thetaSoil[0], 
				fieldCapacity, ks, pore_size_index, capFringe, porosity, mineral_fraction);
			gv_soil = 1.0/(1.0/(gva/(lai*canopy_cover))+1.0/gv_soil);
		}
		else //Bryophytes
		{
//This is for Bryophytes, e.g. Sphagnum moss
//Simply put, resistance is based on deep water supply, 
//  flux is scaled with deep root zone soil moisture - hmmm, crude
			gv_soil = soil_cond_rate(dZ, Dsoil, t_nearSurface, residual, 
					thetaSoil[(int)rmodules-1], fieldCapacity, ks, 
					pore_size_index, capFringe, porosity, mineral_fraction);
			gv_soil = 1.0/(1.0/(gva/(lai*canopy_cover))+1.0/gv_soil);
		}
//calculate evaporation from soil surface
        	soil_evaporation1 = (1.0-litter_store/litter_capacity)*
						do_pm(Rnet_soil, t_surface, Dsurface, p_air, 
							gva/(lai*canopy_cover), gv_soil);
//cout << litter_store/litter_capacity << endl;
//cout << litter_store << '\t' << litter_capacity << '\t' << Rnet_soil << '\t' << t_surface << '\t' << Dsurface << '\t' << gva << '\t' << lai << '\t' << canopy_cover << endl;
		if (soil_evaporation1 < 0.0)
		{
			soil_evaporation1 = 0.0;
		}

//assume albedo of 0.1 for litter
        	Rnet_soil = 0.9 * (Qtotal - Rabs_sun - Rabs_shd) - Rth;
		litter_store = state.get_val_at(LITTER);
		evaporative_fraction = 1.0;
//evaporative_fraction = 1.0;
		if (Dsoil < 0.01)
		{
			Dsoil = 0.01;
		}
		gv_soil = 100.0;
		soil_evaporation2 = evaporative_fraction*do_pm(Rnet_soil, t_surface, Dsoil, 
									p_air, gva, gv_soil);
		if (soil_evaporation2*1.8 > litter_store)
		{
			soil_evaporation2 = litter_store/1.8;
		}
		state.update_val_at(-soil_evaporation2, LITTER);
	}
//independently contributing
	Et = Ec_t + soil_evaporation1 + soil_evaporation2 + canopy_evaporation; 
	simOut.Et = Et;


//use rhizosphere flux from the plant water balance model
//rhizFluxes are relative; convert to mmol s-1 
	sumRhizFlux = 0.0;
	for (i = 0; i < rmodules; i++)
	{
		sumRhizFlux += simOut.rhizFlux[i];
	}
	for (i = 0; i < rmodules; i++)
	{
		if (fabs(simOut.rhizFlux[i]) <= 2.0*sumRhizFlux)
		{
			simOut.rhizFlux[i] = simOut.rhizFlux[i]/sumRhizFlux*simOut.Ec;
		}
		else
		{
			if (simOut.rhizFlux[i] < 0.0)
			{
				simOut.rhizFlux[i] = -2.0*simOut.Ec;
			}
			else
			{
				simOut.rhizFlux[i] = 2.0*simOut.Ec;
			}
		}
//Rflux[] is used in BiogeochemicalCycles::computeRootNitrogenUptake() to compute nitrogen uptake
		//Rflux[i] = simOut.rhizFlux[i] / aral;
		Rflux[i] = simOut.rhizFlux[i];
	}

//update soil moisture with ET losses
	evaporation = soil_evaporation1*1.8; // 1800/1000 -> convert mm s-1 to m 30-min-1

//shallow layer, allow transpiration + evaporation
	i = 0;
	transpiration = EcConversion_StoT(simOut.rhizFlux[i], treesParams.lai)*1.8;
	delta_Sr = -(transpiration+evaporation);
	thetaSoil[i] += delta_Sr/(rootDepth[i]*mineral_fraction);
//all other layers, transpiration only
	for (i = 1; i < rmodules; i++)
	{
//water uptake in each rhizosphere contributing to transpiration
		transpiration = EcConversion_StoT(simOut.rhizFlux[i], treesParams.lai)*1.8;
		if (is_bryophyte == 1) //key byophyte evaporation to deeper layer
		{
			delta_Sr = -(transpiration+evaporation);
		}
		else //deep soil, no evaporation from here
		{
			delta_Sr = -transpiration;
		}
		thetaSoil[i] += delta_Sr/(rootDepth[i]*mineral_fraction);
	}

//
//CAPILLARITY - Modified Gardner(1957) type steady capillary rise
//Update unsaturated soil moisture with ET loss and capillary rise gain
//Capillary rise uses a steady flow governed by water table depth and linearly proportional
//  to soil surface dryness
// 04/09 DSM; 05/10 DSM - logic to handle shallow water tables and bryophytes
//Ec_t is converted from mm s-1 to m 30min-1
//
        
	ks = ks*86400.0; // convert from m s-1 to m day-1
	capRise = capRiseScalar*capillary_rise(Zw, bubbling_pressure, pore_size_index, ks)*
							(1.0-(thetaSoil[0]-residual)/(porosity-residual));

	if (capRise > (Et-soil_evaporation1-soil_evaporation2-canopy_evaporation)*1.8)
	{
		capRise = (Et-soil_evaporation1-soil_evaporation2-canopy_evaporation)*1.8;
	}
	if (capRise < 0.0)
	{
		capRise = 0.0;
	}
	if ((Zw - capFringe) >= cumDroot[0])
	{
		capRise_shallow = capRise * (rootDepth[0])/(Zw-capFringe);
	}
	else
	{
		capRise_shallow = 0.0;
	}
	thetaSoil[0] += capRise_shallow/(rootDepth[0]*mineral_fraction);
	if (thetaSoil[0] > porosity || (Zw-capFringe-cumDroot[0]) < 1.0E-6)
	{
		thetaSoil[0] = (cumDroot[0]-(Zw - capFringe))/rootDepth[0]*porosity +
				((Zw - capFringe) - (cumDroot[0] - rootDepth[0]))/rootDepth[0] * fieldCapacity;
	}
	for (i = 1; i < rmodules-1; i++)
	{
		if ((Zw-capFringe) >= cumDroot[i])
		{
			capRise_mid = capRise * (rootDepth[i])/(Zw - capFringe);
		}
		else
		{
			capRise_mid = 0.0;
		}
		thetaSoil[i] += capRise_mid/(rootDepth[i]*mineral_fraction);
		if (thetaSoil[i] > porosity || (Zw-capFringe-cumDroot[i]) < 1.0E-6)
		{
			thetaSoil[i] = (cumDroot[i]-(Zw - capFringe))/rootDepth[i]*porosity +
				((Zw - capFringe) - (cumDroot[i] - rootDepth[i]))/rootDepth[i] * fieldCapacity;
		}
	}
	if ((Zw-capFringe) >= cumDroot[rmodules-1])
	{
		capRise_deep = capRise * (rootDepth[rmodules-1]) / (Zw - capFringe);
	}
	else
	{
		capRise_deep = 0.0;
	}
	thetaSoil[rmodules-1] += capRise_deep/(rootDepth[rmodules-1]*mineral_fraction);
	if (thetaSoil[rmodules-1] > porosity || (Zw-capFringe-cumDroot[rmodules-1]) < 1.0E-6)
	{
		thetaSoil[rmodules-1] = (cumDroot[rmodules-1]-(Zw - capFringe))/rootDepth[rmodules-1]*porosity +
				((Zw - capFringe) - (cumDroot[rmodules-1] - rootDepth[rmodules-1]))/rootDepth[rmodules-1] * fieldCapacity;
	}

//
//STRUCTURAL AND NON-STRUCTURAL CARBON FLOWS
//
//NSC consumption is done in these steps:
//  1) Maintenance respiration
//  2) Growth allocation
//  3) Defense, if so desired
//  4) Refilling, if permissible
//  5) Redistribution through phloem flow
//

	ti.get_time(year, yday, hour, min);

//
//Updates to glycine and serine exported from C3 photosynthesis
// Glycine: C2H5NO2
// Serine: C3H7NO3
//
	bgc.storeGlycineAndSerine(NGsun, NGshd, NSsun, NSshd, lai_sun, lai_shd);

//convert from umol m-2 leaf area to kgC ha-1 30min-1
// umol m-2 s-1 = kgC ha-1 30min-1 * (1/0.012) mol kg-1 * 10^6 umol mol-1 *(1/10000) ha m-2 *
//    (1/1800) 30min s-1 -> (1/0.012) *10^6 / 10^4 / 1800 = 4.6296
        A_tot_kg = A_tot/4.6296;
	
//total non-structural carbon
	nsc = bgc.getLeafNSC()+bgc.getStemNSC()+bgc.getRootNSC();


//Step 1: Maintenance Respiration
//
//Start with root maintenance respiration, computed for each root segment
//  Includes computation of a stress factor, kratio, based on root hydraulic conductance
//  Modified to handle 10 size classes of roots in each soil-root layer - DSM July 2015
//
	R_root = kratioRoot = 0.0;
	rmaint = 0.0;
	for (i = 0; i < rmodules; i++)
	{
		if (treesParams.useHydraulics == true)
		{
			kratio = k_p_e.latK[2+i]/(ksat[i+3][1]+0.00000001);
		}
		else
		{
			kratio = water_stress;
		}
//when relative hydraulic conductance is < 50%, then scale it between 0 and 1
		if (kratio < 0.5)
		{
			kratio *= 2.0;
		}
		else
		{
			kratio = 1.0;
		}
//reduce respiration in low oxygen
		if (thetaSoil[i] > 0.9*porosity)
		{
			kratio *= sqrt((porosity-thetaSoil[i]+0.00000001)/(0.1*porosity));
		}

		kratio = max(0.01, kratio);

		kratio_vector[i] = kratio;
		kratioRoot += kratio * treesParams.ar[i+3];
		rootClassRespRate = resp_coef_root/48.0; //kg kg-1 30min-1 deg
		for (int k = 0; k < 10; k++)
		{
/*
			if (k < 5)
			{
				rootClassRespRate = resp_coef_root/48.0; //kg kg-1 30min-1 deg
			}
			else
			{
				rootClassRespRate = resp_coef_stem/48.0; //kg kg-1 30min-1 deg
			}
*/
			R_root_segment = bgc.root_respiration_rate(rootClassRespRate, resp_coefficient,
						tempSoil[i], bgc.getRootCarbon(i, k), kratio);
			if (R_root_segment > 0.99*bgc.getRootNSC(i, k))
			{
				R_root_segment = 0.99*bgc.getRootNSC(i, k);
				bgc.updateRootNSC(-R_root_segment, i, k);
			}
			else
			{
				bgc.updateRootNSC(-R_root_segment, i, k);
			}
			R_root += R_root_segment;
			rootClassRespRate *= 0.9;
		}
	}

//Begin maintenance respiration calculations for stem and leaves here
	if (treesParams.useHydraulics == true)
	{
		kratio = k_p_e.latK[1]/(ksat[1][1]+0.00000001);
	}
	else
	{
		kratio = water_stress;
	}
//stress begins at loss of conductivity greater than 50%, and so we scale kratio 0-1 
//    if less than or equal to 0.5
	if (kratio < 0.5)
	{
		kratio *= 2.0;
	}
	else
	{
		kratio = 1.0;
	}
	R_leaf = bgc.leaf_respiration_rate(resp_coef_leaf, resp_coefficient, simOut.T_sun, 
						lai_sun, SLA)/48.0*max(1.0, kratio);
	R_leaf += bgc.leaf_respiration_rate(resp_coef_leaf, resp_coefficient, simOut.T_shd, 
						lai_shd, SLA)/48.0*max(1.0, kratio);

        R_stem = bgc.stem_respiration_rate(resp_coef_stem, resp_coefficient, t_canopy, 
						bgc.getLiveStemCarbon())/48.0*max(1.0, kratio);


//take leaf maintenance respiration from photosynthate, then chloroplast starch, 
//    then from other leaf stores of NSC
//if all sources of NSC are limiting then reduce respiration
//Cannell and Thornley, 2000, Annals of Botany
//   Respiratory costs of phloem loading
//   in the range of 0.5 to 0.8 mol CO2 (mol sucrose)-1
//   using 0.7 mol CO2 (mol sucrose)-1 = 0.7/12 = 0.06 g C respired (g NSC loading)-1
	CostLoading = 0.06*R_leaf;

	if (A_tot_kg > R_leaf)
	{
		A_tot_kg -= R_leaf;
	}
	else if ((bgc.getChloroplastSugar()+A_tot_kg) > R_leaf)
	{
		bgc.putChloroplastSugar(bgc.getChloroplastSugar()+A_tot_kg-R_leaf);
		A_tot_kg = 0.0;
	}
	else if (bgc.getChloroplastStarch() > R_leaf)
	{
		bgc.putChloroplastStarch(bgc.getChloroplastStarch()-R_leaf);
	}
	else if (bgc.getLeafNSC() > (R_leaf+CostLoading))
	{
		bgc.putLeafNSC(bgc.getLeafNSC()-R_leaf-CostLoading);
	}
	else
	{
		rratio = bgc.getLeafNSC()/R_leaf;
                R_leaf *= rratio;
		CostLoading *= rratio;
		bgc.putLeafNSC(bgc.getLeafNSC()-R_leaf-CostLoading);
	}

	psn = bgc.getChloroplastSugar();
	double propToStarch = 0.0;
	if (simOut.A_sun > 0.0)
	{
		//propToStarch = (2.0*simOut.A_sun/simOut.Vcmax_sun)*((24.0-daylen)/24.0);
		propToStarch = (1.0*simOut.A_sun/treesParams.Vcmax25)*((24.0-daylen)/24.0);
	}
	if (A_tot_kg > 0.0)
	{
		psn += (1.0-propToStarch)*A_tot_kg;
		state.set_val_at(psn, PSN);
		bgc.putChloroplastSugar(psn);
		bgc.putChloroplastStarch(bgc.getChloroplastStarch() + propToStarch*A_tot_kg);
	}

//take stem maintenance respiration from stored NSC
	if (bgc.getStemNSC() >= 1.01*R_stem)
	{
		bgc.putStemNSC(bgc.getStemNSC()-R_stem);
	}
	else
	{
		rratio = bgc.getStemNSC()/R_stem;
                R_stem *= rratio;
		bgc.putStemNSC(bgc.getStemNSC()-R_stem);
	}

        rmaint = R_root + R_leaf + R_stem;

//
//How much growth respiration
//additional growth other than associated with leaf expansion;
//this growth will:
//1. determine the target LAI for the next year, but only consume resources during actual leaf expansion
//2. determine growth of stem and root elements, in which case resources are consumed
//
//We are assuming here that growth occurs on days when there is some net positive photosynthesis,
//   and growth in the current days uses stored carbohydrates
//Calculate a potential rate proportional to GPP and hydraulic conductance
//Reduction in growth follows a rate proportional to the square of hydraulic limitation
//
//activate growth when there is chloroplast NSC
//
	growth_resp_frac = growth_resp_proportion*kratio*kratioRoot;
	//rgrowth = (bgc.getChloroplastSugar()+bgc.getChloroplastStarch())/24.0 * growth_resp_frac;
	rgrowth = (bgc.getChloroplastStarch())/24.0 * growth_resp_frac;
	if (treesParams.useLeafModule == 1)
	{
		rgrowth += 0.5*(bgc.getLeafNSC()+bgc.getStemNSC()+bgc.getRootNSC())/48.0*growth_resp_frac;
	}
//
//the following (nscRatio) is used to increase growth allocation and allocation to canopy
//  if there is a surplus of non-structural carbon
//  or decrease it if NSC is in short supply
//
	totalRootCarbon = bgc.getRootCarbon();
	if (treesParams.usePhenology == true)
	{
		nsc_base = (lai_pot-state.get_val_at(CURRENTTARGETLAI))/treesParams.SLA*10000.0 + 
			0.5*(treesParams.leafNSCscalar*bgc.getLeafBiomassCarbon() +
			treesParams.leafNSCscalar*totalRootCarbon + 
			0.4*treesParams.leafNSCscalar*bgc.getLiveStemCarbon());
	}
	else
	{
		nsc_base =  (treesParams.leafNSCscalar*bgc.getLeafBiomassCarbon() +
			treesParams.leafNSCscalar*bgc.getRootCarbon() + 
			0.4*treesParams.leafNSCscalar*bgc.getLiveStemCarbon());
	}
	nsc = bgc.getLeafNSC()+bgc.getStemNSC()+bgc.getRootNSC();
        nscRatio = nsc/(nsc_base+0.001);
	if (treesParams.useLeafModule==1)
	{
		nscRatio = 1.0;
	}
	rgrowth *= nscRatio;
	if (rgrowth < 1.0E-6)
	{
		rgrowth = 0.0;
	}

//Update growth and leaf area and root-to-leaf area ratio if not in dormancy (root T < 5 C)
//What to do about LAI - need to modify lai and Al; lai_at_sat_kl is used only once per simulation
//Assuming leaf biomass can be approximated by 86% (carbon in cellulose) of NSC use
//This computes potential lai at saturated kl, with actual determined by phenology
//NEED - LIVE LAI for canopy conductance since brown canopy LAI absorbs energy, but does not transpire
//
	unstressedLeafLifeSpan = treesParams.leafLifeSpan*365.25*48.0;

//Use plantNstatus to set N_avail_rate_plant so that lack of passive N uptake 
//   reduces growth and allocation to leaf
//
	N_avail_rate_plant = bgc.plantNstatus[0];
	double reproduction0, delta_reproduction;
	if (t_root >= -50.0)
	{
		reproduction0 = bgc.getFruitCarbon();
		bgc.computeLeafAllocation(treesParams, newc, N_avail_rate_plant, kratio, nsc, psn, 
				nscRatio, rgrowth, leafCfraction, lai, yday, stressedLeafLifeSpan, SLA);
		delta_reproduction = bgc.getFruitCarbon() - reproduction0;

//During current growing season mark the highest potential live lai to be used for setting 
//   lai target for next year
		lai_pot = treesParams.live_lai;
		if (yday == 180 && hour == 0 && min == 0)
		{
			state.set_val_at(lai_pot, FORECASTTARGETLAI);
		}
		else if (yday >= 180 && yday < 270)
		{
			state.update_val_at(lai_pot, FORECASTTARGETLAI);
		}
		else if (yday == 270 && hour == 0 && min == 0)
		{
			lai_pot = state.get_val_at(FORECASTTARGETLAI) / 90.0 / 48.0;
			state.set_val_at(lai_pot, FORECASTTARGETLAI);
		}
			
		treesParams.aral_at_sat_kl = treesParams.Ar_Al;

//Update stem C and N and add to the N demand
		bgc.updateStemCarbonNitrogenPools(treesParams, N_avail_rate_plant, kratio, nscRatio, rgrowth, 
			stemAllocation, N_neg_fract, N_neg_demand, N_pos_demand, leafCfraction, 
			t_canopy, treesParams.lai, stressedLeafLifeSpan, treesParams.SLA);

//Update root C and N and add to the N demand
//Update root area to leaf area ratio from rootCarbon stores
		bgc.updateRootCarbonNitrogenPools(treesParams, tempSoil, rgrowth, N_neg_fract, 
			N_neg_demand, N_pos_demand, kratio_vector, stressedLeafLifeSpan, 
			unstressedLeafLifeSpan, N_avail_rate_plant, leafCfraction, stemAllocation);
	}

//if nsc drops below zero then reduce live leaf area
//this is a bit of a kludge since this shouldn't happen
	if (nsc < 0.0)
	{
		treesParams.live_lai = treesParams.live_lai + nsc*treesParams.SLA/10000.0;
		nsc = 0.0;
	}

//add NSC use for defense; ramp up with BB attack if NSC can be utilized
//assumes about 5% NSC used over a 10-day period under non-invasion conditions
//this is poorly specified and so rdefense is currently set to zero
	if (treesParams.useHydraulics == true)
	{
		kratio = k_p_e.latK[1]/(ksat[1][1]+0.00000001);
	}
	else
	{
		kratio = water_stress;
	}
        rdefense = 0.0;
	if (kratio > 0.5)
	{
		rdefense = 0.05/kratio * rmaint;
	}
	else
	{
		rdefense = 0.05*2.0*kratio * rmaint;
	}
	bgc.putLeafNSC(bgc.getLeafNSC()-rdefense);
/*
        if (t_root > 0.0 && nsc > 0.0)
	{
                rdefense = min(1.0,t_root/5.0)*(0.05/30.0/48.0*nsc*treesParams.leafLifeSpan)/max(0.25, treesParams.xylemScalar)*max(0.01, kratio);
		nsc -= rdefense;
	}
        else
                rdefense = 0.0;
*/

//What to do about litter capacity
//Not sure about the rate of change, but litter capacity decline with leaf area index
	treesParams.litter_capacity = treesParams.litter_capacity_init * 
							treesParams.lai/treesParams.lai_at_sat_kl;

//This section emulates a circadian control of chloroplast starch and sugar modulated by daylight length
	CostLoading = 0.06*0.5*psn;
	nsc += 0.5*psn;
	bgc.putLeafNSC(bgc.getLeafNSC()+0.5*psn);
	bgc.putChloroplastSugar(bgc.getChloroplastSugar()-0.5*psn-CostLoading);
	psn -= 0.5*psn-CostLoading;

/*
cout << bgc.getChloroplastStarch() << '\t';
	bgc.putChloroplastStarch(bgc.getChloroplastStarch()-rgrowth/nscRatio);
cout << bgc.getChloroplastStarch() << endl;
*/

	if ((hour+min/100.0) >= (12.0-0.5*daylen) && (hour+min/100.0) < (12.5-0.5*daylen))
        {
		if (psn > 0.0)
                {
			CostLoading = 0.06*psn;
                        nsc += psn-CostLoading;
			bgc.putStemNSC(bgc.getStemNSC()+psn-CostLoading);
			bgc.putLeafNSC(bgc.getLeafNSC()+bgc.getChloroplastStarch());
			bgc.putChloroplastSugar(0.0);
			bgc.putChloroplastStarch(0.0);
                        psn = 0.0;
			state.set_val_at(psn, PSN);
                }

//LOW NSC OR GROWTH -> TURN OFF REFILLING
//Now checking this once per day since rewetting is checked once per day
//Low NSC is less than 2% during cold periods or less 0.5% anytime
		if (treesParams.forceRefilling == true && treesParams.useHydraulics == true)
		{
			if (t_root < 10.0 && kratio < 0.5 && 
				nsc/(Croot+Cstem+treesParams.lai/treesParams.SLA*10000.0)*100.0 < 2.0)
			{
				treesParams.useRefilling = false;
				if (!silent)
				{
					cout << "*** Re-filling turned off ***" << endl;
				}
			}
			if (nsc/(Croot+Cstem+treesParams.lai/treesParams.SLA*10000.0)*100.0 < 0.5)
			{
				treesParams.useRefilling = false;
				if (!silent)
				{
                        		cout << "*** Re-filling turned off ***" << endl;
				}
			}
//If NSC recovers to at least 0.5% of live tissue then restore refilling
			if (nsc/(Croot+Cstem+treesParams.lai/treesParams.SLA*10000.0)*100.0 >= 0.5 && 
								t_root >= 10.0 && kratio >= 0.5)
			{
				treesParams.useRefilling = true;
				if (!silent)
				{
					cout << "*** Re-filling turned on ***" << endl;
				}
			}
			if (kratio < 0.2)
			{
				treesParams.useRefilling = false;
                                if (!silent)
				{
					cout << "*** Re-filling turned off ***" << endl;
				}
                        }
			else if (kratio > 0.5)
			{
				treesParams.useRefilling = true;
                                if (!silent)
				{
					cout << "*** Re-filling turned on ***" << endl;
				}
                        }
		}
	}

//Move NSC and N as a function of concentration gradients and hydraulics
	bgc.computeNSCfluxes(treesParams, kratio, kratio_vector);
	bgc.computeNfluxes(treesParams, kratio, kratio_vector);

	state.set_val_at(psn, PSN);
	state.set_val_at(nsc, NSC);
	simOut.psn = psn;

//switch commenting for the next two linse to write out NSC as a proportion of structural carbon
	//simOut.nsc = nsc/(Croot+treesParams.Csapwood+treesParams.lai/treesParams.SLA*10000)*100;
	simOut.nsc = nsc; //this is straight carbon in kg ha-1

//These terms are C allocation to growth and maintenance respiration, not CO2 production
	simOut.rgrowth = rgrowth;
	simOut.rmaint = rmaint;

//
//***********************************************************************************
//HETEROTROPHIC RESPIRATION
//***********************************************************************************
//
//reset heterotrophic respiration values to zero, as these will be accumulated
//
	for (i = 0; i < rmodules; i++)
	{
		bgc.heterotrophicRespiration[i] = 0.0;
	}
//
//ROOT EXUDATE PRODUCTION
//
	//if ((hour+min/100.0) == 5.0)
	//{
		bgc.computeRootExudates(treesParams, thetaSoil, kratio_vector);
	//}
//
// NITROGEN CYCLE UPDATES
// Compute the effects of plant metabolism on root and rhizosphere nitrogen 
//compute nitrogen uptake 
//
	
	bgc.computeRootNitrogenUptake(UP_neg, UP_pos, treesParams, thetaSoil, Rflux, ratot, 
					fieldCapacity, N_neg_demand, N_pos_demand);

	bgc.computeLeaching(LE_neg, LE_pos, treesParams, thetaSoil, rootDepth, bypassFlow);

//update carbon and nitrogen pools in the rhizosphere
	bgc.updateRhizospherePools(treesParams, thetaSoil, tempSoil, rootDepth, UP_neg, UP_pos, 
					LE_neg, LE_pos);
//
//assumes fast heterotrophic respiration occurs in proximity to roots, slow elsewhere
//therefore, use relative root area for the fast respiration, relative depth for slow
//
	R_het = 0.0;
	for (i = 0; i < rmodules; i++)
	{
		R_het += bgc.heterotrophicRespiration[i];
	}

//Rate functions give a daily rate, which is divided by 48 to get a half-hourly rate
        R_bg =  (R_het + R_root + (1.0-leafCfraction-stemAllocation)*0.14*rgrowth);

//Total respiration
//Convert units from KgC/ha/30 min to umol/m2/s 
// umol m-2 s-1 = kgC ha-1 30min-1 * (1/0.012) mol kg-1 * 10^6 umol mol-1 *(1/10000) ha m-2 * 
//    (1/1800) 30min s-1 -> (1/0.012) *10^6 / 10^4 / 1800 = 4.6296
//Assume 50% of C used for defensive compounds goes to respired CO2
//Includes fraction of leaf growth respiration from leaf expansion
//	R_total = (R_bg + R_stem) / 0.012 * pow(10.0,6.0) / pow(10.0,4.0) / 1800.0;
	R_total = (simOut.leaf_growth_respiration + 0.14*rgrowth +
		   R_bg + R_leaf + R_stem + delta_reproduction * 1.14 +
		   0.5*rdefense) * 4.6296; //convert to umol/m2/s
        NEE = -(A_tot - R_total);
        simOut.NEE = NEE;
	simOut.NPP = A_tot - (simOut.leaf_growth_respiration + 0.14*rgrowth +
				0.14*delta_reproduction +
				R_leaf + R_stem + R_root + 0.5*rdefense) * 4.6296; 

//cout << A_tot << '\t' << simOut.NPP << '\t' << 0.14*rgrowth * 4.6296 << endl;
        simOut.R_total = R_total;
        simOut.R_bg = R_bg*4.6296;
	simOut.R_ag = simOut.R_total - simOut.R_bg;
/*
	simOut.R_ag = (simOut.leaf_growth_respiration + (leafCfraction+stemAllocation)*
						0.14*rgrowth + R_leaf + R_stem + 0.5*rdefense)*4.6296;
*/


//*******************************************************************************************
//DO SNOWMELT HERE
//*******************************************************************************************
	snowpack = state.get_val_at(SNOWPACK);
	snowpack_E_deficit = state.get_val_at(SNOWEDEF);
	snowmelt = 0.0;
	if (snowpack > 0.0)
	{
		if (t_ref > 0.0 && snowpack_E_deficit >= 0.0)
		{
//divide Rnet_soil by 5 to account for snow albedo - temporary, need to evolve
			snowmelt += radMelt(snowpack_E_deficit, Rnet_soil/5.0, t_ref, melt_Rcoef, snowpack);
			if (snowmelt > snowpack)
			{
				snowmelt = snowpack;
				snowpack = 0.0;
				snowpack_E_deficit = 0.0;
			}
			else
			{
				snowpack -= snowmelt;
			}
		}
		else if (t_ref <= 0.0)
		{
			sublimation = sublimate(snowpack_E_deficit, Rnet_soil/5.0, t_ref, 
										melt_Rcoef, snowpack);
			if (sublimation <= snowpack)
			{
				snowpack -= sublimation;
				simOut.Et += sublimation;
			}
			else
			{
				sublimation = snowpack;
				snowpack = 0.0;
				snowpack_E_deficit = 0.0;
				simOut.Et += sublimation;
			}
		}
		if (snowpack_E_deficit >= 0.0 && snowpack > 0.0)
		{
			turbulent_melt = sensibleMelt(snowpack_E_deficit, t_ref, p_air, 
									snowpack, Dsoil, 0.1*u_ref);
			turbulent_melt += latentMelt(snowpack_E_deficit, t_ref, snowpack, Dsoil, 0.1*u_ref);
			if (turbulent_melt <= snowpack)
			{
				snowpack -= turbulent_melt;
				snowmelt += turbulent_melt;
			}
			else
			{
				turbulent_melt = snowpack;
				snowpack = 0.0;
				snowpack_E_deficit = 0.0;
				snowmelt += turbulent_melt;
			}
		}
		state.set_val_at(snowpack, SNOWPACK);
		state.set_val_at(snowpack_E_deficit, SNOWEDEF);
	}
	simOut.snowpack = snowpack;
	simOut.snowpack_E_deficit = snowpack_E_deficit;

//******************************************************************************************
//FINAL DRAINAGE
//add snowmelt water here
//update drainage for second half of the time-step
//drain the shallow soil zone first
//1) drainage fron litter layer
//******************************************************************************************

	ks = ks * 100.0 / 24.0;    // convert from m/day to cm/hr
	infiltration(snowmelt, state, rootDepth, thetaSoil, bypassFlow, rmodules,
			litter_capacity, porosity, pore_size_index, fieldCapacity, 
			residual, capFringe, mineral_fraction, drainScalar, ks);

//check that lowest soil-root theta is no higher than porosity
	if (thetaSoil[rmodules-1] > porosity)
	{
		thetaSoil[rmodules-1] = porosity;
	}
//compute full root zone average theta
	simOut.thetaRoot = 0.0;
	for (i = 0; i < rmodules; i++)
	{
		simOut.thetaRoot += thetaSoil[i]*treesParams.ar[i+3];
	}
//end of simulation functions, simOut contains results to be printed
	return simOut;
}


