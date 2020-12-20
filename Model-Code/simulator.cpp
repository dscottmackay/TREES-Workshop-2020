//
//simulator.cpp
//
//simulator(): Implements the front end to the TREES simulations, including the time loop
//
//There are three implentations of simulator(). The first implementation is called by the Bayesian
//    code in trees_main.cpp. The second implementation is called by trees_main.cpp when just a single
//    deterministic simulation is desired. The third implementation is called by getPIL.cpp, a 
//    post-processing program used with the Bayesian analysis, and not part of the TREES executable
//
//Author: D. Scott Mackay, sometime in 2007
//
//
//This function handles:
//    - initialization of state variables, arrays, etc
//    - time loop
//    - post-simulation bookkeeping
//This function has a polymorphic behavior with:
//    - first simulator() called from Bayesian MCMC loop in trees_main
//    - second simulator() called once from trees_main
//    - third simulator() called from get_PIL for post-MCMC analysis
//
//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "simulator2.h"


//--------------------------------------------------------------------------------------------------
//This is called during MCMC simulations
// Note: Due to memory issue the full hydraulic  model is best not used for this call
//--------------------------------------------------------------------------------------------------
void simulator(Data_Store& in_data,
		 State_Store& state,
		 Parameter* current_params,
		 trees_params treesParams,
		 Data_Frame& Ecframe,
		 Data_Frame& NEEframe,
		 double& sd_err1,
		 double& sd_err2,
		 double& sd_err1_wt,
		 int n_agg)
{
	int n_tsteps = in_data.no_of_steps();  //# of time steps in input data
	int i, j, k, l;
	Time ti;
	double u_ref, t_ref, D_ref, precip, Qpar, t_canopy, D_canopy, CO2_atm, Ecobs;
	double p_atm, t_surface, t_soil, t_root, Zw, NEEobs;
	double xylem_capacity_multiplier;
	double thetaSoil[ULAT], tempSoil[ULAT];
	struct sim_out simOut;
	double aggNEEobs, aggETobs, aggNEEsim, aggETsim;
	int isvalidNEE, isvalidET;
	double validNEEcount, validETcount;
        std::vector<double> GSI_vector(22);
//These variables interface with the HydraulicModel model
        std::vector<double> EcState;
        std::vector<double> EcStateInner(50);
        std::vector<double> LeafPsiInner(50);
        std::vector<double> KpInner(50);
        std::vector<int> InnerStep(50);
        int ModelStatus = 0;
        double Ecrit = 0.0;
        double Thresh = 0.0;
        int ts;

	double ks_in[MD][MD], newb_in[MD][MD], newc_in[MD][MD], psimin_in[MD][MD];

	for (i = 1; i <= 21; ++i)
	{
		GSI_vector[i] = 0.0;
	}

//declare and assign parameters
//Note: position in the .p file is important as name matching is not done
        i = 0;
//Radiation parameters
        treesParams.altitude = current_params[i].get_value(); i++;  //altitude (meters)
        treesParams.lat = current_params[i].get_value(); i++;  //latitude (degrees)
        treesParams.longi = current_params[i].get_value(); i++; //longitude (degrees)
        treesParams.z_ref = current_params[i].get_value(); i++;  //ref height, m
        treesParams.lai = current_params[i].get_value(); i++;  //single sided lai
	treesParams.live_lai = treesParams.lai;
	//treesParams.Al = max(1.0, treesParams.lai);
	treesParams.Al = treesParams.lai;
	treesParams.lai_at_sat_kl= treesParams.Al;
        treesParams.canopy_ht = current_params[i].get_value(); i++; //height of canopy, m
        treesParams.lai_at_full_canopy_height = current_params[i].get_value(); i++; //fraction cloud cover
        treesParams.l_angle = current_params[i].get_value(); i++; //leaf angle distribution
        treesParams.canopy_e = current_params[i].get_value(); i++;  //canopy emissivity
        treesParams.fPAR_beam = current_params[i].get_value(); i++;  //fraction of PAR in beam radiation
        treesParams.fPAR_diff = current_params[i].get_value(); i++;  //fraction of PAR in diffuse radiation
        treesParams.alpha_PAR = current_params[i].get_value(); i++;  //alpha for PAR & NIR typically 0.8 & 0.2
        treesParams.alpha_NIR = current_params[i].get_value(); i++;
        treesParams.omega = current_params[i].get_value(); i++;  //parameter used to adjust ext. coefficient (0 to 1)
        treesParams.p_crown = current_params[i].get_value(); i++;  //parameter used in clumping factor calculation (1 to 3.34)
//Aerodynamic parameters
        treesParams.d_factor = current_params[i].get_value(); i++;  //d = d_factor*canopy_ht
        treesParams.zm_factor = current_params[i].get_value(); i++;  //zm = zm_factor*canopy_ht
        treesParams.zh_factor = current_params[i].get_value(); i++;  //zh = zh_factor*zm
	treesParams.ps_model = (int) current_params[i].get_value(); i++;  //photosynthesis model to select (1,2, or 3)
//Photosynthesis parameters
        treesParams.Rd_mult = current_params[i].get_value(); i++;  //Rd = Rd_mult * Vmax
        treesParams.Jmax_mult = current_params[i].get_value(); i++;  //Jmax = Jmax_mult * Vmax
        treesParams.thetaJ = current_params[i].get_value(); i++;  //J curvature parameter
        treesParams.phiJ_sun = current_params[i].get_value(); i++;  //effective quantum yield of sunlit PSII system, e-/umol
        treesParams.phiJ_shd = current_params[i].get_value(); i++;  //effective quantum yield of shaded PSII system, e-/umol
        treesParams.Nleaf = current_params[i].get_value(); i++;  //leaf N concentration (kg/m2)
        treesParams.N_fixed_proportion = current_params[i].get_value(); i++;  //leaf N concentration (kg/m2)
        treesParams.Nrubisco = current_params[i].get_value(); i++;  //leaf proportion of N in rubisco
        treesParams.Kc25 = current_params[i].get_value(); i++;  //(Pa) MM const carboxylase, 25 deg C
        treesParams.q10Kc = current_params[i].get_value(); i++;  //(DIM) Q_10 for kc
        treesParams.Ko25 = current_params[i].get_value(); i++;  //(Pa) MM const oxygenase, 25 deg C
        treesParams.q10Ko = current_params[i].get_value(); i++;  //(DIM) Q_10 for ko
        treesParams.act25 = current_params[i].get_value(); i++;  //(umol/mgRubisco/min) Rubisco activity
        treesParams.q10act = current_params[i].get_value(); i++;  //(DIM) Q_10 for Rubisco activity
//Added for C4 - DSM August 2019
        treesParams.Vcmax25 = current_params[i].get_value(); i++;  //maximum Rubisco activity at 25 C, umol m-2 s-1
        treesParams.Vpmax25 = current_params[i].get_value(); i++;  //maximum PEP carbolylase activity at 25 C, umol m-2 s-1
        treesParams.Jmax25 = current_params[i].get_value(); i++;  //maximum electron transport rate at 25 C, umol m-2 s-1
        treesParams.gammaStar25 = current_params[i].get_value(); i++;  //compensation point, ubar
        treesParams.Kp25 = current_params[i].get_value(); i++;  //Michaelis constant of PEP carboxylase for CO2 at 25 C, ubar
        treesParams.Vpr = current_params[i].get_value(); i++;  //PEP regeneration rate, umol m-2 s-1
        treesParams.f = current_params[i].get_value(); i++;  //correction for spectral quality of light
        treesParams.x = current_params[i].get_value(); i++;  //partitioning factor of electron transport rate
        treesParams.absorptance = current_params[i].get_value(); i++;  //fraction of irradiance absorbed
        treesParams.E_Vcmax = current_params[i].get_value(); i++;  //activation energy, maximum carboxylation rate, kJ mol-1
        treesParams.E_Vpmax = current_params[i].get_value(); i++;  //activation energy, maximum PEP rate, kJ mol-1
        treesParams.E_Jmax = current_params[i].get_value(); i++;  //activation energy, electron transport, kJ mol-1
        treesParams.E_Kp = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of PEP, kJ mol-1
        treesParams.E_Kc = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of carboxylation, kJ mol-1
        treesParams.E_Ko = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of oxygenation, kJ mol-1
        treesParams.E_Rd = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of mitochondrial respiration, kJ mol-1
        treesParams.E_gammaStar = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of compensation point, kJ mol-1
        treesParams.gm = current_params[i].get_value(); i++;  //mesophyll conductance to CO2, mol m-2 s-1
        treesParams.gbs = current_params[i].get_value(); i++;  //conductance of the bundle sheath, mol m-2 s-1
        treesParams.alphaGmax = current_params[i].get_value(); i++;  //fraction of glycolate carbon diverted to glycine during photorespiration
        treesParams.alphaSmax = current_params[i].get_value(); i++;  //fraction of glycolate carbon diverted to serine during photorespiration
        treesParams.Nmax = current_params[i].get_value(); i++;  //maximum rate of de novo nitrogen supply to the chloroplast, umol N m-2 s-1

//Plant hydraulic Parameters
        treesParams.Gsref0 = current_params[i].get_value(); i++;  //maximum gs, molm-2s-1 (leaf hemisurface area)
        double m = current_params[i].get_value(); i++;  //kPa-1 as D in kPa
	treesParams.isAmphistomatous = (bool) current_params[i].get_value(); i++; //stomata on two sides
        treesParams.Md = current_params[i].get_value(); i++;
        treesParams.midday_at_sat_kl = current_params[i].get_value(); i++;
        treesParams.e_at_saturated_kl= current_params[i].get_value(); i++;
        treesParams.rhizosphere_width= current_params[i].get_value(); i++;
        //treesParams.E_inc= current_params[i].get_value(); i++;
	treesParams.E_inc=0.0001;
        treesParams.soilshells= (int) current_params[i].get_value(); i++;
        treesParams.GMP= current_params[i].get_value(); i++;
        treesParams.GSD= current_params[i].get_value(); i++;
        treesParams.BD= current_params[i].get_value(); i++;
	treesParams.porosity = current_params[i].get_value(); i++;  //
        treesParams.silt_fraction= current_params[i].get_value(); i++;
        treesParams.clay_fraction= current_params[i].get_value(); i++;
	treesParams.residual= current_params[i].get_value(); i++;
        treesParams.frac_absorbing_length= current_params[i].get_value(); i++;
        treesParams.Capacitance= current_params[i].get_value(); i++;
        treesParams.axK_latKl_shoot_modules= current_params[i].get_value(); i++;
        treesParams.axKr_latKr_root_modules= current_params[i].get_value(); i++;
        treesParams.per_total_R_in_root_system= current_params[i].get_value(); i++;
        treesParams.pd_at_sat_kl= current_params[i].get_value(); i++;
//calculate kl from hydraulic parameters, kl=e/(pd-md)
	treesParams.saturated_kl_for_whole_plant = treesParams.e_at_saturated_kl/(treesParams.pd_at_sat_kl-treesParams.midday_at_sat_kl);

        treesParams.ax_Shoot_b_value= current_params[i].get_value(); i++;
        treesParams.ax_Shoot_c_value= current_params[i].get_value(); i++;
        treesParams.lat_Shoot_b_value= current_params[i].get_value(); i++;
        treesParams.lat_Shoot_c_value= current_params[i].get_value(); i++;
        treesParams.ax_Root_b_value= current_params[i].get_value(); i++;
        treesParams.ax_Root_c_value= current_params[i].get_value(); i++;
        treesParams.lat_Root_b_value= current_params[i].get_value(); i++;
        treesParams.lat_Root_c_value= current_params[i].get_value(); i++;
        treesParams.initial_conductivity_root= current_params[i].get_value(); i++;
        treesParams.decrement_root= current_params[i].get_value(); i++;
        treesParams.initial_conductivity_shoot= current_params[i].get_value(); i++;
        treesParams.decrement_shoot= current_params[i].get_value(); i++;
//Biogeochemical cycling parameters
	treesParams.theta_opt = current_params[i].get_value(); i++;  //
	treesParams.optimal_soil_T = current_params[i].get_value(); i++;  //
	treesParams.growth_resp_proportion = current_params[i].get_value(); i++;  //
	treesParams.resp_coef_root = current_params[i].get_value(); i++;  //
	treesParams.resp_coef_stem = current_params[i].get_value(); i++;  //
	treesParams.resp_coef_leaf = current_params[i].get_value(); i++;  //
	treesParams.resp_coefficient = current_params[i].get_value(); i++;  //
	treesParams.EaSx = current_params[i].get_value(); i++;  //
	treesParams.kMsx = current_params[i].get_value(); i++;  //
	treesParams.xASx = current_params[i].get_value(); i++;  //

	treesParams.kd = current_params[i].get_value(); i++;  //
	treesParams.kn = current_params[i].get_value(); i++;  //
	treesParams.kea = current_params[i].get_value(); i++;  //
	treesParams.kes = current_params[i].get_value(); i++;  //
	treesParams.kl = current_params[i].get_value(); i++;  //
	treesParams.kh = current_params[i].get_value(); i++;  //

	treesParams.fr_minCN = current_params[i].get_value(); i++;  //
	treesParams.fr_maxCN = current_params[i].get_value(); i++;  //
	treesParams.leaf_minCN = current_params[i].get_value(); i++;  //
	treesParams.leaf_maxCN = current_params[i].get_value(); i++;  //


	treesParams.Cbelowground = current_params[i].get_value(); i++;  //
	treesParams.Clitter_frac = current_params[i].get_value(); i++;  //
	treesParams.Croot_frac = current_params[i].get_value(); i++;  //
	treesParams.Clitter = treesParams.Clitter_frac * treesParams.Cbelowground;
	treesParams.Croot = treesParams.Croot_frac * treesParams.Cbelowground;
	treesParams.Cstem = current_params[i].get_value(); i++;  //
	treesParams.Csapwood = current_params[i].get_value(); i++;  //
	treesParams.Croot_coarse_frac = current_params[i].get_value(); i++;//
	treesParams.Croot_coarse = treesParams.Croot_coarse_frac * treesParams.Cbelowground;
	treesParams.Csoil = (1.0-treesParams.Clitter_frac - treesParams.Croot_frac - treesParams.Croot_coarse_frac) * treesParams.Cbelowground;
	treesParams.interception_per_leafArea = current_params[i].get_value(); i++;//
	treesParams.litter_capacity = current_params[i].get_value(); i++;//
	treesParams.litter_capacity_init = treesParams.litter_capacity;
	treesParams.theta_deep0 = current_params[i].get_value(); i++;//
	treesParams.theta_mid0 = current_params[i].get_value(); i++;//
	treesParams.theta_shallow0 = current_params[i].get_value(); i++;//
	treesParams.litter_store0 = current_params[i].get_value(); i++;//
	treesParams.SLA = current_params[i].get_value(); i++;  //
	treesParams.SLA_instant = treesParams.SLA;
	treesParams.SRL1 = current_params[i].get_value(); i++;  //
	treesParams.minRootDiam = current_params[i].get_value(); i++;  //
	treesParams.maxRootDiam = current_params[i].get_value(); i++;  //
	treesParams.rootDiamMultiplier = pow(treesParams.maxRootDiam/treesParams.minRootDiam,1.0/9.0);
	treesParams.minRootLifespan = current_params[i].get_value(); i++;  //
	treesParams.LWP_spring_minimum = current_params[i].get_value(); i++;  //
	treesParams.LWP_stomatal_closure = current_params[i].get_value(); i++;  //
	treesParams.is_bryophyte = (int) current_params[i].get_value(); i++;  //
	treesParams.capRiseScalar = current_params[i].get_value(); i++;  //
//How to to multiple precipitation by for experimental drought designs
	double precipReduction = current_params[i].get_value(); i++;  //
	if (precipReduction < 0.0)
	{
		precipReduction = 0.0;
	}
//Threshold for Ec convergence
	treesParams.drainScalar= current_params[i].get_value(); i++;
        treesParams.leafNSCscalar= current_params[i].get_value(); i++;
//What are these for?
        treesParams.usePhenology= (bool) current_params[i].get_value(); i++;
        treesParams.leafLifeSpan= current_params[i].get_value(); i++;
//Maximum number of iterations to achieve convergence of Delta
	treesParams.max_iterations = (int) current_params[i].get_value(); i++;
//Conductance to use for Darcy's Law, 1=WholePlant,2=AxialComponents, 3=Shoot,4=AxialRoot, 5=LaterialRoot
	treesParams.microbiomeScalar = current_params[i].get_value(); i++;
// MCH 07082020
	treesParams.microbialrainrate = current_params[i].get_value(); i++; //
//MCH 23092020
        treesParams.raininAmmonium = current_params[i].get_value(); i++; //
        treesParams.raininNitrate = current_params[i].get_value(); i++; //
        treesParams.raininMineralN = current_params[i].get_value(); i++; //
        treesParams.raininLabileC = current_params[i].get_value(); i++; //
//Parameters for snowpack accumulation and melt
	treesParams.snowpack_water_equivalent = current_params[i].get_value(); i++; //
	treesParams.snowpack_E_deficit_max = current_params[i].get_value(); i++; //
	treesParams.melt_Rcoef = current_params[i].get_value(); i++; //
//Run with full hydraulics or not (1 = yes, 0 = no)
	treesParams.useHydraulics = (bool) current_params[i].get_value(); i++; //
	treesParams.useInputStress = (bool) current_params[i].get_value(); i++; //
	treesParams.useInputWaterTable = (bool) current_params[i].get_value(); i++; //
	treesParams.dayToStopMaizeRefilling = current_params[i].get_value(); i++; //
//Leaf growth module parameters
        treesParams.useLeafModule = (bool) current_params[i].get_value(); i++; //
        treesParams.leafAreaMax = current_params[i].get_value(); i++; //
        treesParams.initialLeafSize = current_params[i].get_value(); i++; //
        treesParams.leafArea_Rate = current_params[i].get_value(); i++; //
        treesParams.dur_LeafExpansion= current_params[i].get_value(); i++; //
        treesParams.SLA_max = current_params[i].get_value(); i++; //
    	treesParams.SLA_min = current_params[i].get_value(); i++; //
        treesParams.leaf_insertAngle = current_params[i].get_value(); i++; //
    	treesParams.leaf_len_to_width = current_params[i].get_value(); i++; //
    	treesParams.proportion_CD = current_params[i].get_value(); i++; //
        treesParams.phyllochron = current_params[i].get_value(); i++; //
        treesParams.floweringTime = current_params[i].get_value(); i++; //
        treesParams.Tbase = current_params[i].get_value(); i++; //
        treesParams.therm_plant = current_params[i].get_value(); i++; //
    	treesParams.projectedArea_init = current_params[i].get_value(); i++; //
    	treesParams.pot_size = current_params[i].get_value(); i++; //
    	treesParams.root_to_shoot = current_params[i].get_value(); i++; //
    	treesParams.leaf_to_stem = current_params[i].get_value(); i++; //
        treesParams.useLeafGamma = (bool) current_params[i].get_value(); i++; //
        treesParams.Kalpha = current_params[i].get_value(); i++; //
        treesParams.Kbeta = current_params[i].get_value(); i++; //
        treesParams.Nalpha = current_params[i].get_value(); i++; //
        treesParams.Nbeta = current_params[i].get_value(); i++; //
        treesParams.ralpha = current_params[i].get_value(); i++; //
        treesParams.rbeta = current_params[i].get_value(); i++; //
//MCMC Parameters
	sd_err1 = current_params[i].get_value(); i++; //
	sd_err2 = current_params[i].get_value(); i++; //
	sd_err1_wt = current_params[i].get_value(); 

	treesParams.delta = m*treesParams.Gsref0; // sensitivity of stomata to VPD, mol m-2 s-1

	//treesParams.live_lai = treesParams.live_lai * (1.0 - 0.15/treesParams.leafLifeSpan);

//Calculate Brooks-Corey ks (cm/hr), bubbling pressure (cm), pore-size distribution index, and residual water content
	double ks, bubbling_pressure, pore_size_index, residual, fieldCapacity;
	double pClay = 100.0*treesParams.clay_fraction;
	double pSand = 100.0*(1.0 - treesParams.silt_fraction - treesParams.clay_fraction);
	soil_texture_parameters(treesParams.porosity, pClay, pSand, ks, bubbling_pressure, pore_size_index, residual);
	treesParams.ks = ks;
	treesParams.bubbling_pressure = bubbling_pressure;
	treesParams.pore_size_index = pore_size_index;
	//treesParams.residual = residual;
	residual = treesParams.residual;

//Initialize state variables here
	Var_Dictionary vd;
	Var_name deep_name = "deep";
	Var_name mid_name = "mid";
	Var_name shallow_name = "shallow";
	Var_name litter_name = "litter";
	Var_name canopy_name = "canopy";
	Var_name lwp_name = "lwp";
	Var_name kl_name = "kl";
	Var_name ecrit_name = "ecrit";
	Var_name psicrit_name = "psicrit";
	Var_name snowpack_name = "snowpack";
	Var_name snowEdef_name = "snowEdef";
	Var_name xylemCap_name = "xylem";
	vd.insert(deep_name);
	vd.insert(mid_name);
	vd.insert(shallow_name);
	vd.insert(litter_name);
	vd.insert(canopy_name);
	vd.insert(lwp_name);
	vd.insert(kl_name);
	vd.insert(ecrit_name);
	vd.insert(psicrit_name);
	vd.insert(snowpack_name);
	vd.insert(snowEdef_name);

	vd.insert(xylemCap_name);

        Var_name regen_name = "regen";
        vd.insert(regen_name);

        Var_name currentTargetLai_name = "currenttargetlai";
        vd.insert(currentTargetLai_name);
        Var_name forecastTargetLai_name = "forecasttargetlai";
        vd.insert(forecastTargetLai_name);

	state.allocate(vd);   //allocate memory 
	state.set_val_at(treesParams.theta_deep0, DEEPSW);
	state.set_val_at(treesParams.theta_mid0, MIDSW);
	state.set_val_at(treesParams.theta_shallow0, SHALLOWSW);
	state.set_val_at(treesParams.litter_store0, LITTER);
	state.set_val_at(0.0, CANOPY);
	state.set_val_at(0.0, PSN);
	state.set_val_at(treesParams.lai/treesParams.SLA*10000.0*((0.14/treesParams.leafLifeSpan)+0.10)+treesParams.Croot*0.10+treesParams.Croot_coarse*0.10 + 0.04*treesParams.Csapwood, NSC);
        state.set_val_at(treesParams.Cbelowground*treesParams.Croot_frac, ROOTC);
	state.set_val_at(treesParams.snowpack_water_equivalent, SNOWPACK);
	state.set_val_at(treesParams.snowpack_E_deficit_max, SNOWEDEF);
	state.set_val_at(treesParams.saturated_kl_for_whole_plant, KL);

	state.set_val_at(treesParams.lai, CURRENTTARGETLAI);
	state.set_val_at(0.0, FORECASTTARGETLAI);


//Calculate field capacity
	fieldCapacity = field_capacity(bubbling_pressure, pore_size_index, residual, treesParams.porosity);
	treesParams.field_capacity = fieldCapacity;

//Initializing state variables for HydraulicModel
        for(j=0; j < n_tsteps; j++)
        {
                EcState.push_back(0.0);
        }
        //set to 50 because model should take nowhere near this long to converge.
        for(j=0; j < 50; j++)
        {
                EcStateInner.push_back(0.0);
                LeafPsiInner.push_back(0.0);
                KpInner.push_back(0.0);
                //RInner.push_back(0.0);

        }
        for(j=0; j < 50; j++)
        {
                InnerStep.push_back(0.0);
        }

//Initialize soil theta in each layer (up to ULAT layers)
	for (j = 0; j < treesParams.rmodules; j++)
	{
		if (j == 0 || j == 1)
		{
			thetaSoil[j] = treesParams.theta_shallow0;
		}
		else if (j == 2 || j == 3)
		{
			thetaSoil[j] = treesParams.theta_mid0;
		}
		else
		{
			thetaSoil[j] = treesParams.theta_deep0;
		}
	}


//define a BiogeochemicalCycles object
	BiogeochemicalCycles bgc(treesParams);
//Define the canopy closure based on crown diameter
	treesParams.canopy_cover = bgc.computeCanopyCover(treesParams);
//Initialize root carbon stores in each layer (up to ULAT layers)
	double rootArea = bgc.computeRootArea(treesParams);
	bgc.computeLateralRootWeibulls(treesParams);

	treesParams.Ar_Al = rootArea / treesParams.lai / treesParams.canopy_cover;
	treesParams.aral_at_sat_kl = treesParams.Ar_Al;
	treesParams.Ar_Al_init = treesParams.Ar_Al;

//define a HydraulicModel object
	HydraulicModel *hydraulicModel = new HydraulicModel();

//define parameters for setup()
	bool reset = true;
        HydraulicModel::outputStruct ecrit_k_psi, k_p_e;
        double ecritOut[MD], pc[MD], klpred[MD], ppredOut[MD], epredOut[MD], nodeFail[MD];
        char *nodeTyp = new char[MD];
        double ***psi, **rflux;
        if (treesParams.useHydraulics == true)
        {
                psi = get3DArray(NREC,MD,ULAT);
                rflux = get2DArray(NREC,NMAX);
        }
        double soilpsiavg[NREC],  nts[NREC], evap[NREC], kroot[MD], axr[MD], latr[MD];
        double kshoot[MD], dslat[MD], dsax[MD], drlat[MD],drax[MD], ll[MD];
        double ksat[MD][MD], bsat[MD][MD], ccsat[MD][MD], newb[MD][MD], newc[MD][MD];
        double ip[MD][ULAT], b[MD][ULAT], b1[MD][ULAT], n[MD][ULAT], n1[MD][ULAT];
        double r[MD][ULAT], vu[MD][ULAT], vl[MD][ULAT], cf[MD][ULAT], pe[MD][ULAT], ws_[MD][ULAT];
        double ks2[MD][ULAT], p[MD][ULAT], dpp[MD][ULAT], jl[MD][ULAT], wnu[MD][ULAT], wu[MD][ULAT];
        double wnl[MD][ULAT], wl[MD][ULAT], cpu[MD][ULAT], cpl[MD][ULAT], cpp[MD][ULAT], ku[MD][ULAT];
        double kl[MD][ULAT], f[MD][ULAT], cc[MD][2], psinode[NMAX][NMAX], psimin[NMAX][NMAX];
        double jmatrix[NMAX][NMAX], jmatrix2[NMAX][NMAX];
        double percent[NMAX], rssoil[NMAX], rs[NMAX], dp[NMAX], ff[NMAX];
        int col[NMAX], row[NMAX],  indxx[NMAX];
        double subtract[NMAX], pressure[NKINC], plc[NKINC], plcweib[NKINC];
        double rsquare[NKINC], soilpsi[ULAT], al[MD], ar[MD], saturatedKs[MD][4];
        double ptarg, e=0.0, rr, rzw, rd, einc;
        int SoilRhizElements;
        double dt, gmd, gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract;
        double kla, aral, latot, ratot;
        int shootElements = treesParams.smodules, ktotal;
        double soilpsimin = 0.0;
        int ShootModules, rootElements = treesParams.rmodules, RootModules, totmodules;
        double axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b, axRoot_c, latRoot_b, latRoot_c;
        double plco, pdecrement, sumy1, sumy2, sumprod, cincrement, soilvol=0.0;
        double rlateral, rLat_base, slateral, sLat_base, raxial, pleafave;
        int tnode;
        double fac = 0.0, modelledKL;
        int HydraulicModelFailCond=0;

//call setup() to initialize HydraulicModel objects
	if (treesParams.useHydraulics == true)
	{
    		hydraulicModel->setup(reset, ecrit_k_psi, k_p_e, ecritOut, pc, klpred, ppredOut, epredOut,
         	nodeFail, nodeTyp, psi, rflux, soilpsiavg, nts, evap, kroot,
         	axr, latr, kshoot, dslat, dsax, drlat, drax, ll, ksat, bsat,
         	ccsat, newb, newc, ip, b, b1, n, n1, r, vu, vl, cf, pe, ws_, ks2,
         	p, dpp, jl, wnu, wu, wnl, wl, cpu, cpl, cpp, ku, kl, f, cc, psinode,
         	psimin, jmatrix, jmatrix2, percent, rssoil, rs, dp, ff, col, row,
         	indxx, subtract, pressure, plc, plcweib, rsquare, soilpsi, al, ar,
         	saturatedKs, ptarg, e, rr, rzw, rd, einc, SoilRhizElements, dt, gmd,
         	gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract, kla, aral,
         	latot, ratot, shootElements, ktotal, soilpsimin, ShootModules, rootElements,
         	RootModules, totmodules, axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b,
         	axRoot_c, latRoot_b, latRoot_c, plco, pdecrement, sumy1, sumy2, sumprod,
         	cincrement, soilvol, rlateral, rLat_base, slateral, sLat_base, raxial,
         	pleafave, tnode, fac, modelledKL, HydraulicModelFailCond, treesParams);
	}

	state.set_val_at(treesParams.pd_at_sat_kl, LPSI);
	state.set_val_at(treesParams.saturated_kl_for_whole_plant, KL);
	state.set_val_at(treesParams.e_at_saturated_kl, ECRIT);

//time loop
	for (k = 0; k < n_tsteps/n_agg; k++)
	{
	    aggNEEobs = aggNEEsim = aggETobs = aggETsim = 0.0; //set aggregate variables to zero
	    isvalidNEE = isvalidET = (int) INVALID; //keep track of valid aggregatess
	    validNEEcount = validETcount = 0.0; //number of elements used to calculate the aggregate time average flux
	    for (l = 0; l < n_agg; l++) //test aggrregation to hourly
	    {
		i = k*n_agg+l;
		ts = i+1; //HydraulicModel initialization check, ts=0 initializes Gs
		if(in_data.check_step(i) == 0) //valid step w/ valid data
                {
		//get input data
			ti = in_data.get_time(i);
			j=0;
			u_ref = in_data.get_val_at(j,i); j++;
                	t_ref = in_data.get_val_at(j,i); j++;
                	D_ref = in_data.get_val_at(j,i); j++;
                	precip = in_data.get_val_at(j,i); j++;
			precip *= 0.001; //convert mm to m
			precip *= precipReduction;
                	Qpar = in_data.get_val_at(j,i); j++;
                	t_canopy = in_data.get_val_at(j,i); j++;
                	D_canopy = in_data.get_val_at(j,i); j++;
                	p_atm = in_data.get_val_at(j,i); j++;
                	CO2_atm = in_data.get_val_at(j,i); j++;
                	t_surface = in_data.get_val_at(j,i); j++;
                	t_soil = in_data.get_val_at(j,i); j++;
                	t_root = in_data.get_val_at(j,i); j++;
                	Zw = in_data.get_val_at(j,i); j++;
                	xylem_capacity_multiplier = in_data.get_val_at(j,i); j++;
			treesParams.xylemScalar = xylem_capacity_multiplier;
			if (treesParams.useInputStress == true)
			{
                		treesParams.stressScalar = in_data.get_val_at(j,i); j++;
			}

                	//ignore = in_data.get_val_at(j,i); j++;
                	NEEobs = in_data.get_val_at(j,i); j++;
                	Ecobs = in_data.get_val_at(j,i);

			bool silent = true;

//Initialize soil temp in each layer (up to ULAT layers)
        		for (j = 0; j < treesParams.rmodules; j++)
        		{
                		if (j == 0)
				{
                        		tempSoil[j] = 0.5*(t_surface+t_soil);
				}
                		else if (j == 1)
				{
                        		tempSoil[j] = t_soil;
				}
				else
				{
					tempSoil[j] = t_root;
				}
        		}

		//call simulations functions for current time step
			simOut = simulation_functions(
						silent,
						ts,
						Ecrit,
						Thresh,
						state,
						thetaSoil,
						tempSoil,
						ti,
						u_ref,
						t_ref,
						D_ref,
						precip,
						Qpar,
						t_canopy,
						D_canopy,
						p_atm,
						CO2_atm,
						t_surface,
						t_soil,
						t_root,
						Zw,
						treesParams,
						EcState,
						EcStateInner,
						GSI_vector,
						EcState[ts-1],
						ModelStatus,
						LeafPsiInner,
                                              	KpInner,
                                              	InnerStep,
						bgc,
						hydraulicModel,
						ecrit_k_psi, k_p_e, ecritOut, pc, klpred, ppredOut, epredOut,
         					nodeFail, nodeTyp, psi, rflux, soilpsiavg, nts, evap, kroot,
         					axr, latr, kshoot, dslat, dsax, drlat, drax, ll, ksat, bsat,
         					ccsat, newb, newc, ip, b, b1, n, n1, r, vu, vl, cf, pe, ws_, ks2,
         					p, dpp, jl, wnu, wu, wnl, wl, cpu, cpl, cpp, ku, kl, f, cc, psinode,
         					psimin, jmatrix, jmatrix2, percent, rssoil, rs, dp, ff, col, row,
         					indxx, subtract, pressure, plc, plcweib, rsquare, soilpsi, al, ar,
         					saturatedKs, ptarg, e, rr, rzw, rd, einc, SoilRhizElements, dt, gmd,
         					gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract, kla, aral,
         					latot, ratot, shootElements, ktotal, soilpsimin, ShootModules, rootElements,
         					RootModules, totmodules, axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b,
         					axRoot_c, latRoot_b, latRoot_c, plco, pdecrement, sumy1, sumy2, sumprod,
         					cincrement, soilvol, rlateral, rLat_base, slateral, sLat_base, raxial,
         					pleafave, tnode, fac, modelledKL, HydraulicModelFailCond);

                        state.set_val_at(xylem_capacity_multiplier, XYLEM);
//Aggregate x and y variables for the data frame

			//if (NEEobs != INVALID && NEEobs > -20 && NEEobs < 10)
			if (NEEobs != INVALID)
			{
				aggNEEobs += NEEobs;
				aggNEEsim += simOut.NEE;
				//aggNEEsim += simOut.A;
				isvalidNEE++;
				validNEEcount += 1.0;
			}
			if (Ecobs > 0 && i > 0)
			{
				aggETobs += Ecobs;
				aggETsim += simOut.Ec;
				isvalidET++;
				validETcount += 1.0;
			}

		} // end if
//clear vector
	   } //end aggregation loop
	   if (isvalidNEE != INVALID)
	   {
		NEEframe.set_y_value(k, aggNEEsim/validNEEcount);
		NEEframe.set_x_value(k, aggNEEobs/validNEEcount);
	   }
	   else
	   {
		NEEframe.set_y_value(k, INVALID);
	   }
	   if (isvalidET != INVALID)
	   {
		Ecframe.set_y_value(k, aggETsim/validETcount);
		Ecframe.set_x_value(k, aggETobs/validETcount);
	   }
	   else
	   {
		Ecframe.set_y_value(k, INVALID);
		Ecframe.set_x_value(k, INVALID);
	   }
	} //end time loop

	GSI_vector.clear();
        EcState.clear();
        EcStateInner.clear();
        LeafPsiInner.clear();
        KpInner.clear();
     	InnerStep.clear();
	state.clear_values();
	delete hydraulicModel;
	delete[] nodeTyp;
/*
        if (treesParams.useHydraulics == true)
        {
                free3DArray(psi, NREC, MD);
                free2DArray(rflux, NREC);
        }
*/
} //end simulator


//---------------------------------------------------------------------------------------------------
//This is called for single simulations with fixed parameters
// Note: The full water balance model can be used here if data availability or problem justifies it
//---------------------------------------------------------------------------------------------------
void simulator(Data_Store& in_data,
		 State_Store& state,
		 Parameter* current_params,
		 trees_params treesParams,
		 Data_Frame& Ecframe,
		 Data_Frame& NEEframe,
		 ofstream& fluxout,
		 ofstream& hydrout, ofstream& lfareaout)
{
	int n_tsteps = in_data.no_of_steps();  //# of time steps in input data
	int i, j, k;
    int l;
	Time ti;
	double u_ref, t_ref, D_ref, precip, Qpar, t_canopy, D_canopy, CO2_atm, Ecobs;
	double p_atm, t_surface, t_soil, t_root, Zw, NEEobs;
	double xylem_capacity_multiplier, Al_init;
	double thetaSoil[ULAT], tempSoil[ULAT];
	int year, yday, hour, min;
	struct sim_out simOut;
        std::vector<double> GSI_vector(22);
//These variables interface with the HydraulicModel model
        std::vector<double> EcState;
        std::vector<double> EcStateInner(50);
        std::vector<double> LeafPsiInner(50);
        std::vector<double> KpInner(50);
        std::vector<int> InnerStep(50);
        int ModelStatus = 0;
        double Ecrit = 0.0;
        double Thresh = 0.0;
        int ts;

	double p_in[MD][MD], ks_in[MD][MD], newb_in[MD][MD], newc_in[MD][MD], psimin_in[MD][MD];

	for (i = 1; i <= 21; ++i)
	{
                GSI_vector[i] = 0.0;
	}

//declare and assign parameters
//Note: position in the .p file is important as name matching is not done
        i = 0;
//Radiation parameters
        treesParams.altitude = current_params[i].get_value(); i++;  //altitude (meters)
        treesParams.lat = current_params[i].get_value(); i++;  //latitude (degrees)
        treesParams.longi = current_params[i].get_value(); i++; //longitude (degrees)
        treesParams.z_ref = current_params[i].get_value(); i++;  //ref height, m
        treesParams.lai = current_params[i].get_value(); i++;  //single sided lai

	treesParams.live_lai = treesParams.lai;
	treesParams.Al = treesParams.lai;
	treesParams.lai_at_sat_kl = treesParams.Al;
	//treesParams.Al = max(1.0, treesParams.lai);
	Al_init = treesParams.Al;
        treesParams.canopy_ht = current_params[i].get_value(); i++; //height of canopy, m
        treesParams.lai_at_full_canopy_height = current_params[i].get_value(); i++; //fraction cloud cover
        treesParams.l_angle = current_params[i].get_value(); i++; //leaf angle distribution
        treesParams.canopy_e = current_params[i].get_value(); i++;  //canopy emissivity
        treesParams.fPAR_beam = current_params[i].get_value(); i++;  //fraction of PAR in beam radiation
        treesParams.fPAR_diff = current_params[i].get_value(); i++;  //fraction of PAR in diffuse radiation
        treesParams.alpha_PAR = current_params[i].get_value(); i++;  //alpha for PAR & NIR typically 0.8 & 0.2
        treesParams.alpha_NIR = current_params[i].get_value(); i++;
        treesParams.omega = current_params[i].get_value(); i++;  //parameter used to adjust ext. coefficient (0 to 1)
        treesParams.p_crown = current_params[i].get_value(); i++;  //parameter used in clumping factor calculation (1 to 3.34)
//Aerodynamic parameters
        treesParams.d_factor = current_params[i].get_value(); i++;  //d = d_factor*canopy_ht
        treesParams.zm_factor = current_params[i].get_value(); i++;  //zm = zm_factor*canopy_ht
        treesParams.zh_factor = current_params[i].get_value(); i++;  //zh = zh_factor*zm
	treesParams.ps_model = (int) current_params[i].get_value(); i++;  //photosynthesis model to select (1,2, or 3)
//Photosynthesis parameters
        treesParams.Rd_mult = current_params[i].get_value(); i++;  //Rd = Rd_mult * Vmax
        treesParams.Jmax_mult = current_params[i].get_value(); i++;  //Jmax = Jmax_mult * Vmax
        treesParams.thetaJ = current_params[i].get_value(); i++;  //J curvature parameter
        treesParams.phiJ_sun = current_params[i].get_value(); i++;  //effective quantum yield of sunlit PSII system, e-/umol
        treesParams.phiJ_shd = current_params[i].get_value(); i++;  //effective quantum yield of shaded PSII system, e-/umol
        treesParams.Nleaf = current_params[i].get_value(); i++;  //leaf N concentration (kg/m2)
        treesParams.N_fixed_proportion = current_params[i].get_value(); i++;  //leaf N concentration (kg/m2)
        treesParams.Nrubisco = current_params[i].get_value(); i++;  //leaf proportion of N in rubisco
        treesParams.Kc25 = current_params[i].get_value(); i++;  //(Pa) MM const carboxylase, 25 deg C
        treesParams.q10Kc = current_params[i].get_value(); i++;  //(DIM) Q_10 for kc
        treesParams.Ko25 = current_params[i].get_value(); i++;  //(Pa) MM const oxygenase, 25 deg C
        treesParams.q10Ko = current_params[i].get_value(); i++;  //(DIM) Q_10 for ko
        treesParams.act25 = current_params[i].get_value(); i++;  //(umol/mgRubisco/min) Rubisco activity
        treesParams.q10act = current_params[i].get_value(); i++;  //(DIM) Q_10 for Rubisco activity
//Added for C4 - DSM August 2019
        treesParams.Vcmax25 = current_params[i].get_value(); i++;  //maximum Rubisco activity at 25 C, umol m-2 s-1
        treesParams.Vpmax25 = current_params[i].get_value(); i++;  //maximum PEP carbolylase activity at 25 C, umol m-2 s-1
        treesParams.Jmax25 = current_params[i].get_value(); i++;  //maximum electron transport rate at 25 C, umol m-2 s-1
	treesParams.gammaStar25 = current_params[i].get_value(); i++;  //compensation point, ubar
        treesParams.Kp25 = current_params[i].get_value(); i++;  //Michaelis constant of PEP carboxylase for CO2 at 25 C, ubar
        treesParams.Vpr = current_params[i].get_value(); i++;  //PEP regeneration rate, umol m-2 s-1
        treesParams.f = current_params[i].get_value(); i++;  //correction for spectral quality of light
        treesParams.x = current_params[i].get_value(); i++;  //partitioning factor of electron transport rate
        treesParams.absorptance = current_params[i].get_value(); i++;  //fraction of irradiance absorbed
        treesParams.E_Vcmax = current_params[i].get_value(); i++;  //activation energy, maximum carboxylation rate, kJ mol-1
        treesParams.E_Vpmax = current_params[i].get_value(); i++;  //activation energy, maximum PEP rate, kJ mol-1
        treesParams.E_Jmax = current_params[i].get_value(); i++;  //activation energy, electron transport, kJ mol-1
        treesParams.E_Kp = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of PEP, kJ mol-1
	treesParams.E_Kc = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of carboxylation, kJ mol-1
        treesParams.E_Ko = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of oxygenation, kJ mol-1
        treesParams.E_Rd = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of mitochondrial respiration, kJ mol-1
        treesParams.E_gammaStar = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of compensation point, kJ mol-1
        treesParams.gm = current_params[i].get_value(); i++;  //mesophyll conductance to CO2, mol m-2 s-1
        treesParams.gbs = current_params[i].get_value(); i++;  //conductance of the bundle sheath, mol m-2 s-1
        treesParams.alphaGmax = current_params[i].get_value(); i++;  //fraction of glycolate carbon diverted to glycine during photorespiration
        treesParams.alphaSmax = current_params[i].get_value(); i++;  //fraction of glycolate carbon diverted to serine during photorespiration
        treesParams.Nmax = current_params[i].get_value(); i++;  //maximum rate of de novo nitrogen supply to the chloroplast, umol N m-2 s-1

//Plant hydraulic Parameters
        treesParams.Gsref0 = current_params[i].get_value(); i++;  //maximum gs, molm-2s-1 (leaf hemisurface area)
        double m = current_params[i].get_value(); i++;  //kPa-1 as D in kPa
	treesParams.isAmphistomatous = (bool) current_params[i].get_value(); i++; //stomata on two sides
        treesParams.Md = current_params[i].get_value(); i++;
        treesParams.midday_at_sat_kl = current_params[i].get_value(); i++;
        treesParams.e_at_saturated_kl= current_params[i].get_value(); i++;
        treesParams.rhizosphere_width= current_params[i].get_value(); i++;
        //treesParams.E_inc= current_params[i].get_value(); i++;
	treesParams.E_inc=0.0001;
        treesParams.soilshells= (int) current_params[i].get_value(); i++;
        treesParams.GMP= current_params[i].get_value(); i++;
        treesParams.GSD= current_params[i].get_value(); i++;
        treesParams.BD= current_params[i].get_value(); i++;
        treesParams.porosity = current_params[i].get_value(); i++;  //
        treesParams.silt_fraction= current_params[i].get_value(); i++;
        treesParams.clay_fraction= current_params[i].get_value(); i++;
        treesParams.residual= current_params[i].get_value(); i++;
        treesParams.frac_absorbing_length= current_params[i].get_value(); i++;
        treesParams.Capacitance= current_params[i].get_value(); i++;
        treesParams.axK_latKl_shoot_modules= current_params[i].get_value(); i++;
        treesParams.axKr_latKr_root_modules= current_params[i].get_value(); i++;
        treesParams.per_total_R_in_root_system= current_params[i].get_value(); i++;
        treesParams.pd_at_sat_kl= current_params[i].get_value(); i++;
//calculate kl from hydraulic parameters, kl=e/(pd-md)
	treesParams.saturated_kl_for_whole_plant = treesParams.e_at_saturated_kl/(treesParams.pd_at_sat_kl-treesParams.midday_at_sat_kl);
        treesParams.ax_Shoot_b_value= current_params[i].get_value(); i++;
        treesParams.ax_Shoot_c_value= current_params[i].get_value(); i++;
        treesParams.lat_Shoot_b_value= current_params[i].get_value(); i++;
        treesParams.lat_Shoot_c_value= current_params[i].get_value(); i++;
        treesParams.ax_Root_b_value= current_params[i].get_value(); i++;
        treesParams.ax_Root_c_value= current_params[i].get_value(); i++;
        treesParams.lat_Root_b_value= current_params[i].get_value(); i++;
        treesParams.lat_Root_c_value= current_params[i].get_value(); i++;
        treesParams.initial_conductivity_root= current_params[i].get_value(); i++;
        treesParams.decrement_root= current_params[i].get_value(); i++;
        treesParams.initial_conductivity_shoot= current_params[i].get_value(); i++;
        treesParams.decrement_shoot= current_params[i].get_value(); i++;
//Biogeochemical cycling parameters
        treesParams.theta_opt = current_params[i].get_value(); i++;  //
        treesParams.optimal_soil_T = current_params[i].get_value(); i++;  //
        treesParams.growth_resp_proportion = current_params[i].get_value(); i++;  //
        treesParams.resp_coef_root = current_params[i].get_value(); i++;  //
        treesParams.resp_coef_stem = current_params[i].get_value(); i++;  //
        treesParams.resp_coef_leaf = current_params[i].get_value(); i++;  //
        treesParams.resp_coefficient = current_params[i].get_value(); i++;  //
        treesParams.EaSx = current_params[i].get_value(); i++;  //
        treesParams.kMsx = current_params[i].get_value(); i++;  //
        treesParams.xASx = current_params[i].get_value(); i++;  //

	treesParams.kd = current_params[i].get_value(); i++;  //
        treesParams.kn = current_params[i].get_value(); i++;  //
        treesParams.kea = current_params[i].get_value(); i++;  //
        treesParams.kes = current_params[i].get_value(); i++;  //
        treesParams.kl = current_params[i].get_value(); i++;  //
        treesParams.kh = current_params[i].get_value(); i++;  //

	treesParams.fr_minCN = current_params[i].get_value(); i++;  //
        treesParams.fr_maxCN = current_params[i].get_value(); i++;  //
        treesParams.leaf_minCN = current_params[i].get_value(); i++;  //
        treesParams.leaf_maxCN = current_params[i].get_value(); i++;  //

        treesParams.Cbelowground = current_params[i].get_value(); i++;  //
        treesParams.Clitter_frac = current_params[i].get_value(); i++;  //
        treesParams.Croot_frac = current_params[i].get_value(); i++;  //
        treesParams.Clitter = treesParams.Clitter_frac * treesParams.Cbelowground;
        treesParams.Croot = treesParams.Croot_frac * treesParams.Cbelowground;
        treesParams.Cstem = current_params[i].get_value(); i++;  //
        treesParams.Csapwood = current_params[i].get_value(); i++;  //
        treesParams.Croot_coarse_frac = current_params[i].get_value(); i++;//
        treesParams.Croot_coarse = treesParams.Croot_coarse_frac * treesParams.Cbelowground;
        treesParams.Csoil = (1.0-treesParams.Clitter_frac - treesParams.Croot_frac - treesParams.Croot_coarse_frac) * treesParams.Cbelowground;
	treesParams.interception_per_leafArea = current_params[i].get_value(); i++;//
        treesParams.litter_capacity = current_params[i].get_value(); i++;//
	treesParams.litter_capacity_init = treesParams.litter_capacity;
        treesParams.theta_deep0 = current_params[i].get_value(); i++;//
        treesParams.theta_mid0 = current_params[i].get_value(); i++;//
        treesParams.theta_shallow0 = current_params[i].get_value(); i++;//
        treesParams.litter_store0 = current_params[i].get_value(); i++;//
        treesParams.SLA = current_params[i].get_value(); i++;  //
	treesParams.SLA_instant = treesParams.SLA;
        treesParams.SRL1 = current_params[i].get_value(); i++;  //
        treesParams.minRootDiam = current_params[i].get_value(); i++;  //
        treesParams.maxRootDiam = current_params[i].get_value(); i++;  //
	treesParams.rootDiamMultiplier = pow(treesParams.maxRootDiam/treesParams.minRootDiam,1.0/9.0);
        treesParams.minRootLifespan = current_params[i].get_value(); i++;  //
        treesParams.LWP_spring_minimum = current_params[i].get_value(); i++;  //
        treesParams.LWP_stomatal_closure = current_params[i].get_value(); i++;  //
        treesParams.is_bryophyte = (int) current_params[i].get_value(); i++;  //
        treesParams.capRiseScalar = current_params[i].get_value(); i++;  //
//How to to multiple precipitation by for experimental drought designs
        double precipReduction = current_params[i].get_value(); i++;  //
        if (precipReduction < 0.0)
	{
                precipReduction = 0.0;
	}
//Threshold for Ec convergence
        treesParams.drainScalar= current_params[i].get_value(); i++;
        treesParams.leafNSCscalar = current_params[i].get_value(); i++;
//What are these for?

        treesParams.usePhenology= (bool) current_params[i].get_value(); i++;
        treesParams.leafLifeSpan= current_params[i].get_value(); i++;
//Maximum number of iterations to achieve convergence of DELTA
        treesParams.max_iterations = (int) current_params[i].get_value(); i++;
//Conductance to use for Darcy's Law, 1=WholePlant,2=AxialComponents, 3=Shoot,4=AxialRoot, 5=LaterialRoot
        treesParams.microbiomeScalar = current_params[i].get_value(); i++;
// MCH 07082020
	treesParams.microbialrainrate = current_params[i].get_value(); i++; //
//MCH 23092020
        treesParams.raininAmmonium = current_params[i].get_value(); i++; //
        treesParams.raininNitrate = current_params[i].get_value(); i++; //
        treesParams.raininMineralN = current_params[i].get_value(); i++; //
        treesParams.raininLabileC = current_params[i].get_value(); i++; //
//Parameters for snowpack accumulation and melt
        treesParams.snowpack_water_equivalent = current_params[i].get_value(); i++; //
        treesParams.snowpack_E_deficit_max = current_params[i].get_value(); i++; //
        treesParams.melt_Rcoef = current_params[i].get_value(); i++; //
//Run with full hydraulics or not (1 = yes, 0 = no)
        treesParams.useHydraulics = (bool) current_params[i].get_value(); i++; //
        treesParams.useInputStress = (bool) current_params[i].get_value(); i++; //
        treesParams.useInputWaterTable = (bool) current_params[i].get_value(); i++; //
        treesParams.dayToStopMaizeRefilling = current_params[i].get_value(); i++; //
	treesParams.allowLeafRefilling = true;
//Leaf growth module parameters
        treesParams.useLeafModule = (bool) current_params[i].get_value(); i++; //
        treesParams.leafAreaMax = current_params[i].get_value(); i++; //
        treesParams.initialLeafSize = current_params[i].get_value(); i++; //
        treesParams.leafArea_Rate = current_params[i].get_value(); i++; //
        treesParams.dur_LeafExpansion= current_params[i].get_value(); i++; //
        treesParams.SLA_max = current_params[i].get_value(); i++; //
    	treesParams.SLA_min = current_params[i].get_value(); i++; //
        treesParams.leaf_insertAngle = current_params[i].get_value(); i++; //
    	treesParams.leaf_len_to_width = current_params[i].get_value(); i++; //
    	treesParams.proportion_CD = current_params[i].get_value(); i++; //
        treesParams.phyllochron = current_params[i].get_value(); i++; //
        treesParams.floweringTime = current_params[i].get_value(); i++; //
        treesParams.Tbase = current_params[i].get_value(); i++; //
        treesParams.therm_plant = current_params[i].get_value(); i++; //
    	treesParams.projectedArea_init = current_params[i].get_value(); i++; //
    	treesParams.pot_size = current_params[i].get_value(); i++; //
    	treesParams.root_to_shoot = current_params[i].get_value(); i++; //
    	treesParams.leaf_to_stem = current_params[i].get_value(); i++; //
        treesParams.useLeafGamma = (bool) current_params[i].get_value(); i++; //
        treesParams.Kalpha = current_params[i].get_value(); i++; //
        treesParams.Kbeta = current_params[i].get_value(); i++; //
        treesParams.Nalpha = current_params[i].get_value(); i++; //
        treesParams.Nbeta = current_params[i].get_value(); i++; //
        treesParams.ralpha = current_params[i].get_value(); i++; //
        treesParams.rbeta = current_params[i].get_value(); i++; //
//MCMC Parameters
        double sd_err1 = current_params[i].get_value(); i++; //
        double sd_err2 = current_params[i].get_value(); i++; //
        double sd_err1_wt = current_params[i].get_value();

        treesParams.delta = m*treesParams.Gsref0; // sensitivity of stomata to VPD, mol m-2 s-1

//Calculate Brooks-Corey ks (cm/hr), bubbling pressure (cm), pore-size distribution index, and residual water content
	double ks, bubbling_pressure, pore_size_index, residual, fieldCapacity;
	double pClay = 100.0*treesParams.clay_fraction;
        double pSand = 100.0*(1.0 - treesParams.silt_fraction - treesParams.clay_fraction);
	soil_texture_parameters(treesParams.porosity, pClay, pSand, ks, bubbling_pressure, pore_size_index, residual);
        treesParams.ks = ks;
        treesParams.bubbling_pressure = bubbling_pressure;
        treesParams.pore_size_index = pore_size_index;
        //treesParams.residual = residual;
	residual = treesParams.residual;

	cout << pSand << '\t' << pClay << '\t' << ks << '\t' << bubbling_pressure << '\t' << pore_size_index << '\t' << residual << endl;
	cout << endl;
	if (treesParams.useHydraulics == true)
	{
		cout << ">>> TREES, v. 3.1.4, deterministic mode with plant water balance <<< \n\n";
	}
	else
	{
		cout << ">>> TREES, v. 3.1.4, deterministic mode without plant water balance <<< \n\n";
	}
	cout << endl;
//Calculate field capacity
	fieldCapacity = field_capacity(bubbling_pressure, pore_size_index, residual, treesParams.porosity);
	treesParams.field_capacity = fieldCapacity;
	if (treesParams.is_bryophyte == 1)
	{
		cout << "Running with Bryophytes\n";
	}
	//cout << "ks = " << ks << endl;
	//cout << "bubbling_pressure = " << bubbling_pressure << endl;
	//cout << "pore_size_index = " << pore_size_index << endl;
	//cout << "residual = " << residual << endl;
	//cout << "fc = " << fieldCapacity << endl;

//Initialize state variables here
	Var_Dictionary vd;
	Var_name deep_name = "deep";
        vd.insert(deep_name);
	Var_name mid_name = "mid";
        vd.insert(mid_name);
        Var_name shallow_name = "shallow";
        vd.insert(shallow_name);
	Var_name litter_name = "litter";
	vd.insert(litter_name);
	Var_name canopy_name = "canopy";
	vd.insert(canopy_name);
        Var_name lwp_name = "lwp";
        Var_name kl_name = "kl";
        Var_name ecrit_name = "ecrit";
        Var_name psicrit_name = "psicrit";
        Var_name snowpack_name = "snowpack";
        Var_name snowEdef_name = "snowEdef";
	Var_name xylemCap_name = "xylem";

//Initialize state variables for carbon pools
	Var_name photosynthate_name = "photosynthate";
	vd.insert(photosynthate_name);
	Var_name nsc_name = "nsc";
	vd.insert(nsc_name);
	Var_name root_name = "root";
	vd.insert(root_name);
        vd.insert(lwp_name);
        vd.insert(kl_name);
        vd.insert(ecrit_name);
        vd.insert(psicrit_name);
        vd.insert(snowpack_name);
        vd.insert(snowEdef_name);

	vd.insert(xylemCap_name);

	Var_name regen_name = "regen";
	vd.insert(regen_name);

	Var_name currentTargetLai_name = "currenttargetlai";
	vd.insert(currentTargetLai_name);
	Var_name forecastTargetLai_name = "forecasttargetlai";
	vd.insert(forecastTargetLai_name);

        state.allocate(vd);   //allocate memory
	state.set_val_at(treesParams.theta_deep0, DEEPSW);
	state.set_val_at(treesParams.theta_mid0, MIDSW);
        state.set_val_at(treesParams.theta_shallow0, SHALLOWSW);
	state.set_val_at(treesParams.litter_store0, LITTER);
	state.set_val_at(0.0, CANOPY);
	state.set_val_at(0.0, PSN);
	state.set_val_at(treesParams.lai/treesParams.SLA*10000.0*((0.14/treesParams.leafLifeSpan)+0.10)+treesParams.Croot*0.10 + treesParams.Croot_coarse*0.10 + 0.04*treesParams.Csapwood, NSC);
	state.set_val_at(treesParams.Cbelowground*treesParams.Croot_frac, ROOTC);
	state.set_val_at(treesParams.snowpack_water_equivalent, SNOWPACK);
        state.set_val_at(treesParams.snowpack_E_deficit_max, SNOWEDEF);
	state.set_val_at(treesParams.saturated_kl_for_whole_plant, KL);

	state.set_val_at(0.0, REGEN);

        state.set_val_at(treesParams.lai, CURRENTTARGETLAI);
        state.set_val_at(0.0, FORECASTTARGETLAI);


//Initializing state variables for HydraulicModel
        for(j=0; j < n_tsteps; j++)
        {
                EcState.push_back(0.0);
        }
        //set to 50 because model should take nowhere near this long to converge.
        for(j=0; j < 50; j++)
        {
                EcStateInner.push_back(0.0);
                LeafPsiInner.push_back(0.0);
                KpInner.push_back(0.0);
        }
        for(j=0; j < 50; j++)
        {
                InnerStep.push_back(0.0);
        }

//Initialize soil theta in each layer (up to ULAT layers)
        for (j = 0; j < treesParams.rmodules; j++)
        {
		if (j == 0 || j == 1)
                {       
                        thetaSoil[j] = treesParams.theta_shallow0;
                }
                else if (j == 2 || j == 3)
                {       
                        thetaSoil[j] = treesParams.theta_mid0;
                }
                else
                {       
                        thetaSoil[j] = treesParams.theta_deep0;
                }      
        }


//define a BiogeochemicalCycles object
	BiogeochemicalCycles bgc(treesParams);
//Define the canopy closure based on crown diameter
	treesParams.canopy_cover = bgc.computeCanopyCover(treesParams);

	double rootArea = bgc.computeRootArea(treesParams);
	bgc.computeLateralRootWeibulls(treesParams);

cout << "rootDiamMultiplier = " << treesParams.rootDiamMultiplier << endl;
cout << "rootArea = " << rootArea << endl;

        treesParams.Ar_Al = rootArea / treesParams.lai / treesParams.canopy_cover;
	treesParams.aral_at_sat_kl = treesParams.Ar_Al;
        treesParams.Ar_Al_init = treesParams.Ar_Al;

//define a HydraulicModel object
        HydraulicModel *hydraulicModel = new HydraulicModel();

//define parameters for setup()
	bool reset;
	HydraulicModel::outputStruct ecrit_k_psi, k_p_e;
	double ecritOut[MD], pc[MD], klpred[MD], ppredOut[MD], epredOut[MD], nodeFail[MD];
	char *nodeTyp = new char[MD];
	double ***psi, **rflux;
	if (treesParams.useHydraulics == true)
        {
    		psi = get3DArray(NREC,MD,ULAT);
    		rflux = get2DArray(NREC,NMAX);
	}
	double soilpsiavg[NREC], soilpsiavg_in[NREC], nts[NREC], evap[NREC], kroot[MD], axr[MD], latr[MD];
	double kshoot[MD], dslat[MD], dsax[MD], drlat[MD],drax[MD], ll[MD];
	double ksat[MD][MD], bsat[MD][MD], ccsat[MD][MD], newb[MD][MD], newc[MD][MD];
	double ip[MD][ULAT], b[MD][ULAT], b1[MD][ULAT], n[MD][ULAT], n1[MD][ULAT];
	double r[MD][ULAT], vu[MD][ULAT], vl[MD][ULAT], cf[MD][ULAT], pe[MD][ULAT], ws_[MD][ULAT];
	double ks2[MD][ULAT], p[MD][ULAT], dpp[MD][ULAT], jl[MD][ULAT], wnu[MD][ULAT], wu[MD][ULAT];
	double wnl[MD][ULAT], wl[MD][ULAT], cpu[MD][ULAT], cpl[MD][ULAT], cpp[MD][ULAT], ku[MD][ULAT];
	double kl[MD][ULAT], f[MD][ULAT], cc[MD][2], psinode[NMAX][NMAX], psimin[NMAX][NMAX];
	double jmatrix[NMAX][NMAX], jmatrix2[NMAX][NMAX];
	double percent[NMAX], rssoil[NMAX], rs[NMAX], dp[NMAX], ff[NMAX];
	int col[NMAX], row[NMAX],  indxx[NMAX];
	double subtract[NMAX], pressure[NKINC], plc[NKINC], plcweib[NKINC];
	double rsquare[NKINC], soilpsi[ULAT], al[MD], ar[MD], saturatedKs[MD][4];
	double ptarg, e=0.0, rr, rzw, rd, einc;
	int SoilRhizElements;
	double dt, gmd, gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract;
	double kla, aral, latot, ratot;
	int shootElements = treesParams.smodules, ktotal;
	double soilpsimin=0.0;
	double soilpsimin_in;
	int ShootModules, rootElements = treesParams.rmodules, RootModules, totmodules;
	double axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b, axRoot_c, latRoot_b, latRoot_c;
	double plco, pdecrement, sumy1, sumy2, sumprod, cincrement, soilvol=0.0;
	double rlateral, rLat_base, slateral, sLat_base, raxial, pleafave;
	double mortality_threshold, fertilization0, fertilization, RLnew, time_to_drop;
	int tnode;
	double fac = 0.0, modelledKL;
	int HydraulicModelFailCond=0;
	double totalRootArea;

//call setup() to initialize HydraulicModel objects
	if (treesParams.useHydraulics == true)
        {
		ratot = rootArea;
		latot = treesParams.lai*treesParams.canopy_cover;
		reset = true;
        	hydraulicModel->setup(reset, ecrit_k_psi, k_p_e, ecritOut, pc, klpred, ppredOut, epredOut,
                	nodeFail, nodeTyp, psi, rflux, soilpsiavg, nts, evap, kroot,
                	axr, latr, kshoot, dslat, dsax, drlat, drax, ll, ksat, bsat,
                	ccsat, newb, newc, ip, b, b1, n, n1, r, vu, vl, cf, pe, ws_, ks2,
                	p, dpp, jl, wnu, wu, wnl, wl, cpu, cpl, cpp, ku, kl, f, cc, psinode,
                	psimin, jmatrix, jmatrix2, percent, rssoil, rs, dp, ff, col, row,
                	indxx, subtract, pressure, plc, plcweib, rsquare, soilpsi, al, ar,
                	saturatedKs, ptarg, e, rr, rzw, rd, einc, SoilRhizElements, dt, gmd,
                	gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract, kla, aral,
                	latot, ratot, shootElements, ktotal, soilpsimin, ShootModules, rootElements,
                	RootModules, totmodules, axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b,
                	axRoot_c, latRoot_b, latRoot_c, plco, pdecrement, sumy1, sumy2, sumprod,
                	cincrement, soilvol, rlateral, rLat_base, slateral, sLat_base, raxial,
                	pleafave, tnode, fac, modelledKL, HydraulicModelFailCond, treesParams);
		latot = treesParams.lai;
		ratot = bgc.computeRootArea(treesParams);
		bgc.updateLateralRootWeibulls(bsat, ccsat, treesParams);
		for (int mm=1; mm<=treesParams.smodules; mm++)
		{
			al[mm] = treesParams.al[mm] * latot;
		}
		for (int mm=1+treesParams.smodules+1; mm <=(treesParams.rmodules+treesParams.smodules+1); mm++)
		{
			ar[mm] = treesParams.ar[mm] * ratot;
			ll[mm] = (ar[mm] / (2.0 * M_PI * rr));
		}
	}

        state.set_val_at(treesParams.pd_at_sat_kl, LPSI);
        state.set_val_at(treesParams.saturated_kl_for_whole_plant, KL);
        state.set_val_at(treesParams.e_at_saturated_kl, ECRIT);

	int mm;
//time loop
	for (i = 0; i < n_tsteps; i++)
	{
		if(in_data.check_step(i) == 0) //valid step w/ valid data
                {

//++++++++++++++++
//Switch parameters from old stand (at death) to regenerating stand with dead needles reducing radiation
// - modification for AGU 2013 DSM
//
// Step 1. Redefine parameters
// Step 2. Re-call setup() for hydraulicModel
// Step 3. Set a flag, treesParams.regen = true
// Step 4. In simulation functions need a modified canopy radiation model

			treesParams.xylemScalar = 1.0;
			time_to_drop = 0.0;
			if (state.get_val_at(REGEN) == 0.0 && treesParams.xylemScalar < 0.06)
			{
				state.set_val_at(0.0, CANOPY);
//record the status of the old canopy dead leaf area index
				treesParams.dead_lai = treesParams.Al;
				treesParams.dead_lai_drop_rate = treesParams.dead_lai / 365.25 / 48.0;
				time_to_drop = treesParams.dead_lai/treesParams.dead_lai_drop_rate;

				treesParams.Gsref0 *= 1.0;
				treesParams.delta = m*treesParams.Gsref0;
				treesParams.e_at_saturated_kl *= 1.0;

//reset appropriate parameters
				treesParams.canopy_ht = 1.0;
				treesParams.Croot_frac *= 1/treesParams.lai_at_sat_kl;
				treesParams.Croot *= 1.0/treesParams.lai_at_sat_kl;
				treesParams.Cstem *= 0.05;
				treesParams.Csapwood = 550.0/treesParams.lai_at_sat_kl;
				treesParams.lai = treesParams.live_lai = treesParams.Al = treesParams.lai_at_sat_kl = 1.0;
				treesParams.Al = 1.01;

				state.set_val_at(treesParams.lai/treesParams.SLA*10000*0.10+treesParams.Cbelowground*treesParams.Croot_frac*0.10 + 0.04*treesParams.Csapwood, NSC);
				RLnew = treesParams.Ar_Al;
				treesParams.Ar_Al = treesParams.Ar_Al_init = treesParams.aral_at_sat_kl = RLnew;
				treesParams.xylemScalar = 1.0;

				for (mm=1; mm<=treesParams.smodules; mm++)
				{
/*
					treesParams.al[mm] = 1.0;
					treesParams.dslat[mm] = 0.4;
*/
					treesParams.dsax[mm] = 1.0;
				}

				e = fac = soilvol = soilpsimin = 0.0;

//Reset K's to saturation
                        	for (mm = 1; mm <= totmodules+1; mm++)
                        	{
                                	if ( mm != ShootModules+1)  // skip over "hidden" module
                                	{
                                        	ks2[mm][1] = ksat[mm][1];
                                        	ks2[mm][0] = ksat[mm][0];
                                	}
                        	}

				//re-call setup() to initialize HydraulicModel  bjects for new canopy
        			if (treesParams.useHydraulics == true)
        			{
					reset = true;
                			hydraulicModel->setup(reset, ecrit_k_psi, k_p_e, ecritOut, pc, klpred, ppredOut, epredOut,
                        			nodeFail, nodeTyp, psi, rflux, soilpsiavg, nts, evap, kroot,
                        			axr, latr, kshoot, dslat, dsax, drlat, drax, ll, ksat, bsat,
                        			ccsat, newb, newc, ip, b, b1, n, n1, r, vu, vl, cf, pe, ws_, ks2,
                        			p, dpp, jl, wnu, wu, wnl, wl, cpu, cpl, cpp, ku, kl, f, cc, psinode,
                        			psimin, jmatrix, jmatrix2, percent, rssoil, rs, dp, ff, col, row,
                        			indxx, subtract, pressure, plc, plcweib, rsquare, soilpsi, al, ar,
                        			saturatedKs, ptarg, e, rr, rzw, rd, einc, SoilRhizElements, dt, gmd,
                        			gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract, kla, aral,
                        			latot, ratot, shootElements, ktotal, soilpsimin, ShootModules, rootElements,
                        			RootModules, totmodules, axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b,
                        			axRoot_c, latRoot_b, latRoot_c, plco, pdecrement, sumy1, sumy2, sumprod,
                        			cincrement, soilvol, rlateral, rLat_base, slateral, sLat_base, raxial,
                        			pleafave, tnode, fac, modelledKL, HydraulicModelFailCond, treesParams);
        			}

				treesParams.forceRefilling = true;

				state.set_val_at(1.0, REGEN);

			} //end REGEN reset
			else if (state.get_val_at(REGEN) == 1.0)
			{
				treesParams.forceRefilling = false;
				treesParams.xylemScalar = 1.0;
			}
			else
			{
				treesParams.forceRefilling = false;
			}
			if (state.get_val_at(REGEN) == 1.0)
			{
				if (time_to_drop <= 1.0)
				{
					treesParams.Nleaf *= (1.0+fertilization);
					fertilization = 0.0;
					treesParams.dead_lai -= treesParams.dead_lai_drop_rate;
					if (treesParams.litter_capacity_init < 0.02)
						treesParams.litter_capacity_init += 8.56E-7;
					if (treesParams.dead_lai < 0.0)
						treesParams.dead_lai = 0.0;
					if (treesParams.dead_lai < 0.01)
					{
						treesParams.Nleaf /= (1.0+fertilization0);
						fertilization0 = 0.0;
					}
				}
				if (time_to_drop > 0)
				{
					time_to_drop -= 0.000057078; //years 30min-1
				}
			}


//Update root area to leaf area ratio from rootCarbon stores
//DSM - July 2015

			rootArea = bgc.computeRootArea(treesParams);
        		//treesParams.Ar_Al = rootArea / max(treesParams.lai, 1.0);
        		treesParams.Ar_Al = rootArea / treesParams.lai / treesParams.canopy_cover;
        		treesParams.aral_at_sat_kl = treesParams.Ar_Al;

			ts = i+1; //HydraulicModel initialization check, ts=0 initializes Gs
//get input data
			ti = in_data.get_time(i);
			j=0;
			u_ref = in_data.get_val_at(j,i); j++;
                	t_ref = in_data.get_val_at(j,i); j++;
                	D_ref = in_data.get_val_at(j,i); j++;
                	precip = in_data.get_val_at(j,i); j++;
			precip *= 0.001; //convert mm to m
			precip *= precipReduction;
                	Qpar = in_data.get_val_at(j,i); j++;
                	t_canopy = in_data.get_val_at(j,i); j++;
                	D_canopy = in_data.get_val_at(j,i); j++;
                	p_atm = in_data.get_val_at(j,i); j++;
                	CO2_atm = in_data.get_val_at(j,i); j++;
                	t_surface = in_data.get_val_at(j,i); j++;
                	t_soil = in_data.get_val_at(j,i); j++;
                	t_root = in_data.get_val_at(j,i); j++;
                	Zw = in_data.get_val_at(j,i); j++;
                	xylem_capacity_multiplier = in_data.get_val_at(j,i); j++;
			if (state.get_val_at(REGEN) == 0.0)
			{
				treesParams.xylemScalar = xylem_capacity_multiplier;
			}
			if (treesParams.useInputStress == true)
			{
                		treesParams.stressScalar = in_data.get_val_at(j,i); j++;
			}

                	//ignore = in_data.get_val_at(j,i); j++;
                	NEEobs = in_data.get_val_at(j,i); j++;
                	Ecobs = in_data.get_val_at(j,i);

//Initialize soil temp in each layer (up to ULAT layers)
                        for (j = 0; j < treesParams.rmodules; j++)
                        {
				if (j == 0)
                                {
                                        tempSoil[j] = 0.5*(t_surface+t_soil);
                                }
                                else if (j == 1)
                                {
                                        tempSoil[j] = t_soil;
                                }
                                else
                                {
                                        tempSoil[j] = t_root;
                                }
                        }

			ti.get_time(year, yday, hour, min);

//This next set of code is used to recall the plant water balance model setup
//  allowing for adjustments in root volume, and xylem refilling
			treesParams.updatedHydraulics = false;
			treesParams.forceRefilling = false;
			if (treesParams.useHydraulics == true)
			{	
//
//Concept: root pressure is sufficient to allow refilling when soil water potential is higher than -0.3 MPa
//         added for maize, refilling occurs at 11:30 pm
//	   Mackay, September 23 2019
//Stages: Early vegetative growth: to day 157; late vegetative growth to 200; 
//  	  early reproductive to day 225
//
				if (treesParams.xylemScalar == 0.99 || (treesParams.useLeafModule == false && treesParams.usePhenology == false && i % 48 == 47 && soilpsiavg[1] > -0.3))
				{
					treesParams.forceRefilling = true;
					treesParams.updatedHydraulics = true;
					//if (yday < 200)
					//if (yday < 252)
					if (yday < treesParams.dayToStopMaizeRefilling)
					{
						treesParams.allowLeafRefilling = true;
					}
					else
					{
						treesParams.allowLeafRefilling = false;
					}
					if (yday < 100 || treesParams.usePhenology == false)
					{
//at refilling reset the potential lai
						treesParams.lai_at_sat_kl = treesParams.live_lai; 
					}
				}
				if ((i > 1 && i % 48 == 0)) //every day
				{
					latot = treesParams.lai*treesParams.canopy_cover;
					ratot = bgc.computeRootArea(treesParams);
					bgc.updateLateralRootWeibulls(bsat, ccsat, treesParams);
					for (mm=1; mm<=treesParams.smodules; mm++)
					{
						al[mm] = treesParams.al[mm] * latot;
					}
					for (mm=1+treesParams.smodules+1; mm <=(treesParams.rmodules+treesParams.smodules+1); mm++)
					{
						ar[mm] = treesParams.ar[mm] * ratot;
						ll[mm] = (ar[mm] / (2.0 * M_PI * rr));
					}
//recompute the maximum saturated hydraulic conductance as the root area and leaf area changes
//for perennials call this once every 10 days; for annuals call every other day
					if ( (treesParams.usePhenology == true && treesParams.forceRefilling == true) || 
						(treesParams.usePhenology == true && treesParams.leafLifeSpan >= 1.0 && i % 480 == 0) ||
						(treesParams.usePhenology == false) )
					{
						for (mm=1; mm <=(treesParams.rmodules+treesParams.smodules+1); mm++)
						{
							psimin_in[mm][0] = psimin[mm][0];
							psimin_in[mm][1] = psimin[mm][1];
						}
						soilpsiavg_in[1] = soilpsiavg[1];
						soilpsimin_in = soilpsimin;
						reset = false;

    						hydraulicModel->setup(reset, ecrit_k_psi, k_p_e, ecritOut, pc, klpred, ppredOut, epredOut,
         						nodeFail, nodeTyp, psi, rflux, soilpsiavg, nts, evap, kroot,
         						axr, latr, kshoot, dslat, dsax, drlat, drax, ll, ksat, bsat,
         						ccsat, newb, newc, ip, b, b1, n, n1, r, vu, vl, cf, pe, ws_, ks2,
         						p, dpp, jl, wnu, wu, wnl, wl, cpu, cpl, cpp, ku, kl, f, cc, psinode,
         						psimin, jmatrix, jmatrix2, percent, rssoil, rs, dp, ff, col, row,
         						indxx, subtract, pressure, plc, plcweib, rsquare, soilpsi, al, ar,
         						saturatedKs, ptarg, e, rr, rzw, rd, einc, SoilRhizElements, dt, gmd,
         						gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract, kla, aral,
         						latot, ratot, shootElements, ktotal, soilpsimin, ShootModules, rootElements,
         						RootModules, totmodules, axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b,
         						axRoot_c, latRoot_b, latRoot_c, plco, pdecrement, sumy1, sumy2, sumprod,
         						cincrement, soilvol, rlateral, rLat_base, slateral, sLat_base, raxial,
         						pleafave, tnode, fac, modelledKL, HydraulicModelFailCond, treesParams);

						soilpsiavg[1] = soilpsiavg_in[1];
						soilpsimin = soilpsimin_in;
						treesParams.updatedHydraulics = true;
					
						for (mm=1; mm <=(treesParams.rmodules+treesParams.smodules+1); mm++)
						{
                                        		psimin[mm][0] = psimin_in[mm][0];
                                        		psimin[mm][1] = psimin_in[mm][1];
						}
						latot = treesParams.lai;
						ratot = bgc.computeRootArea(treesParams);
						bgc.updateLateralRootWeibulls(bsat, ccsat, treesParams);
						for (mm=1; mm<=treesParams.smodules; mm++)
						{
							al[mm] = treesParams.al[mm] * latot;
						}
						for (mm=1+treesParams.smodules+1; mm <=(treesParams.rmodules+treesParams.smodules+1); mm++)
						{
							ar[mm] = treesParams.ar[mm] * ratot;
							ll[mm] = (ar[mm] / (2.0 * M_PI * rr));
						}
					}
				}
			}

			bool silent = false;

//call simulations functions for current time step
                        simOut = simulation_functions(
						silent, ts, Ecrit, Thresh, state, thetaSoil,
						tempSoil, ti, u_ref, t_ref, D_ref,
                                                precip, Qpar, t_canopy, D_canopy, p_atm,
                                                CO2_atm, t_surface, t_soil, t_root, Zw,
                                                treesParams, EcState, EcStateInner, GSI_vector,
                                                EcState[ts-1], ModelStatus, LeafPsiInner, KpInner,
                                                InnerStep, bgc, hydraulicModel,
						ecrit_k_psi, k_p_e, ecritOut, pc, klpred, ppredOut, epredOut,
                                                nodeFail, nodeTyp, psi, rflux, soilpsiavg, nts, evap, kroot,
                                                axr, latr, kshoot, dslat, dsax, drlat, drax, ll, ksat, bsat,
                                                ccsat, newb, newc, ip, b, b1, n, n1, r, vu, vl, cf, pe, ws_, ks2,
                                                p, dpp, jl, wnu, wu, wnl, wl, cpu, cpl, cpp, ku, kl, f, cc, psinode,
                                                psimin, jmatrix, jmatrix2, percent, rssoil, rs, dp, ff, col, row,
                                                indxx, subtract, pressure, plc, plcweib, rsquare, soilpsi, al, ar,
                                                saturatedKs, ptarg, e, rr, rzw, rd, einc, SoilRhizElements, dt, gmd,
                                                gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract, kla, aral,
                                                latot, ratot, shootElements, ktotal, soilpsimin, ShootModules, rootElements,
                                                RootModules, totmodules, axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b,
                                                axRoot_c, latRoot_b, latRoot_c, plco, pdecrement, sumy1, sumy2, sumprod,
                                                cincrement, soilvol, rlateral, rLat_base, slateral, sLat_base, raxial,
                                                pleafave, tnode, fac, modelledKL, HydraulicModelFailCond);

			state.set_val_at(xylem_capacity_multiplier, XYLEM);


                 	fluxout << ti << '\t' << simOut.Et << '\t' << simOut.WPlant_K << '\t';
			fluxout << simOut.Soil_Psi << '\t' << simOut.Leaf_Psi << '\t' << simOut.Psi_Crit << '\t' << simOut.Ecrit << '\t' << simOut.Ec << '\t';
			for (j = 0; j < treesParams.rmodules; j++)
                        {
                                fluxout << simOut.rhizFlux[j] << '\t';
                        }
			fluxout << simOut.Gs << '\t'<< treesParams.lai << '\t' << treesParams.SLA << '\t' << treesParams.live_lai << '\t' << simOut.rmaint << '\t' << simOut.rgrowth << '\t';
			//fluxout  << simOut.nsc << '\t' << simOut.waterStress << '\t' << state.get_val_at(LITTER) << '\t';
			fluxout << bgc.getFruitCarbon() << '\t';
			fluxout  << bgc.getLeafNSC() << '\t' << bgc.getStemNSC() << '\t' << bgc.getRootNSC() << '\t' << bgc.getChloroplastStarch() << '\t' << bgc.getChloroplastSugar() << '\t' << simOut.waterStress << '\t' << state.get_val_at(LITTER) << '\t';
			for (j = 0; j < treesParams.rmodules; j++)
                        {
                                fluxout << thetaSoil[j] << '\t';
                        }
			fluxout  << simOut.thetaRoot << '\t';

			fluxout << simOut.canopy_evaporation << '\t' << simOut.snowpack << '\t' << simOut.snowpack_E_deficit << '\t' << simOut.Vcmax25 << '\t' << simOut.Vcmax_sun << '\t' << simOut.Vcmax_shd << '\t' << simOut.Jmax25 << '\t' << simOut.J_sun << '\t' << simOut.J_shd << '\t';
			fluxout << simOut.A_sun << '\t' << simOut.A_shd << '\t' << simOut.L_sun << '\t' << simOut.L_shd << '\t' << simOut.T_sun << '\t' << simOut.T_shd << '\t' << simOut.D0_sun << '\t' << simOut.D0_shd << '\t' << simOut.Ci_sun << '\t' << simOut.Ci_shd << '\t' << simOut.PAR_sun << '\t' << simOut.PAR_shd << '\t' << simOut.gs_sun << '\t' << simOut.gs_shd << '\t' << simOut.NEE << '\t';
			//fluxout << simOut.NPP << '\t' << simOut.R_total << '\t' << simOut.R_ag << '\t' << simOut.R_bg << '\t' << simOut.Rd_sun << '\t' << simOut.Rd_shd << '\t' << treesParams.Croot << '\t' << treesParams.Csapwood << endl; //output
			fluxout << simOut.NPP << '\t' << simOut.R_total << '\t' << simOut.R_ag << '\t' << simOut.R_bg << '\t' << simOut.Rd_sun << '\t' << simOut.Rd_shd << '\t';
			fluxout << bgc.getLiveStemCarbon() << '\t';

			double rc_sum;
			for (j = 0; j < treesParams.rmodules; j++)
			{
				k = 0;
				fluxout << bgc.getRootCarbon(j,k) << '\t';
			}
			for (j = 0; j < treesParams.rmodules; j++)
			{
				k = 1;
				fluxout << bgc.getRootCarbon(j,k) << '\t';
			}
			for (j = 0; j < treesParams.rmodules; j++)
			{
				rc_sum = 0.0;
				for (k = 0; k < 10; k++)
				{
					rc_sum += bgc.getRootCarbon(j,k);
				}
				fluxout << rc_sum << '\t';
			}
/*
			for (j = 0; j < treesParams.rmodules; j++)
			{
				rc_sum = 0.0;
				for (k = 5; k < 10; k++)
				{
					rc_sum += bgc.getRootCarbon(j,k);
				}
				fluxout << rc_sum << '\t';
			}
*/
			for (j = 0; j < treesParams.rmodules; j++)
			{
				rc_sum = bgc.getFineRootCarbon(j)/bgc.getFineRootBiomassN(j);
				fluxout << rc_sum << '\t';
			}
			fluxout << bgc.getLeafBiomassCarbon()/bgc.getLeafBiomassN() << '\t';
			//fluxout << bgc.getHumus() << '\t';
/*
			for (k = 0; k <10; k++)
			{
				rc_sum = 0.0;
				for (j = 0; j < treesParams.rmodules; j++)
				{
					rc_sum += bgc.getRootCarbon(j,k);;
				}
				fluxout << rc_sum << '\t';
			}
*/
			for (j = 0; j < treesParams.rmodules; j++)
                        {
				fluxout << bgc.getHumus(j) << '\t';
			}
			for (j = 0; j < treesParams.rmodules; j++)
                        {
				fluxout << bgc.getRhizosphereCl(j) << '\t';
			}
			for (j = 0; j < treesParams.rmodules; j++)
                        {
				fluxout << bgc.getRhizosphereNl(j) << '\t';
			}
			for (j = 0; j < treesParams.rmodules; j++)
                        {
				fluxout << bgc.getAminoAcidExudateC(j) << '\t';
			}
			for (j = 0; j < treesParams.rmodules; j++)
                        {
				fluxout << bgc.getSugarExudateC(j) << '\t';
			}
	//		fluxout << bgc.getRhizosphereCl() << '\t' << bgc.getRhizosphereNl() << '\t';
			for (j = 0; j < treesParams.rmodules; j++)
                        {
				fluxout << bgc.getRhizosphereLiveMicrobialCarbon(j) << '\t';
			}
			for (j = 0; j < treesParams.rmodules; j++)
                        {
				fluxout << bgc.getRhizosphereMicrobialNitrogen(j) << '\t';
			}
/*
			for (j = 0; j < treesParams.rmodules; j++)
                        {
				fluxout << bgc.getNitrogenLeaching(j) << '\t';
			}
*/
			//fluxout << bgc.getRhizosphereLiveMicrobialCarbon() << '\t' << bgc.getRhizosphereMicrobialNitrogen() << '\t';
			fluxout << bgc.getRhizosphereNitrateNitrogen() << '\t';
			fluxout << bgc.getRhizosphereAmmoniumNitrogen() << '\t';
			fluxout << bgc.getPlantN() << '\t' << bgc.plantNstatus[0] << '\t';
			fluxout << treesParams.Ar_Al << '\t';
			for (j = 0; j < treesParams.rmodules; j++)
                        {
				fluxout << treesParams.ar[j+3] << '\t';
			}
			fluxout << endl; //output

                 	hydrout << ti << '\t' << k_p_e.latK[1] << '\t';
			for (j = 0; j < treesParams.rmodules; j++)
			{
				hydrout << k_p_e.latK[2+j] << '\t';
			}
                	hydrout << psimin[1][0] << "\t" << psimin[1][1] << "\t";
                	for (j = 0; j < treesParams.rmodules; j++)
                	{
                        	hydrout << psimin[j+3][0] << "\t";
                	}
                	for (j = 0; j < treesParams.rmodules; j++)
                	{
                        	hydrout << psimin[j+3][1] << "\t";
                	}
                	hydrout << ks2[1][0] << "\t" << ks2[1][1] << "\t";
                	for (j = 0; j < treesParams.rmodules; j++)
                	{
                        	hydrout << ks2[j+3][0] << "\t";
                	}
                	for (j = 0; j < treesParams.rmodules; j++)
                	{
                        	hydrout << ks2[j+3][1] << "\t";
                	}
                	hydrout << b[1][0] << "\t" << b[1][1] << "\t";
                	for (j = 0; j < treesParams.rmodules; j++)
                	{
                        	hydrout << b[j+3][0] << "\t";
                	}
                	for (j = 0; j < treesParams.rmodules; j++)
                	{
                        	hydrout << b[j+3][1] << "\t";
                	}
                	hydrout << cc[1][0] << "\t" << cc[1][1] << "\t";
                	for (j = 0; j < treesParams.rmodules; j++)
                	{
                        	hydrout << cc[j+3][0] << "\t";
                	}
                	for (j = 0; j < treesParams.rmodules; j++)
                	{
                        	hydrout << cc[j+3][1] << "\t";
                	}
                	hydrout << endl;

//clear vector
		} //end if
        
        //get leaf area data for this ts
        if(treesParams.useLeafModule== 1)
            lfareaout << ti << '\t';
        {
            for(j=0; j<(bgc.Lf_idx[0]) ; ++j )
            {
                lfareaout << bgc.getSingleLeafArea(j) << "\t" ;
            }
            lfareaout << endl;
            
        }
	} //end time loop
 	EcStateInner.clear();
       	LeafPsiInner.clear();
        KpInner.clear();
        InnerStep.clear();
       	GSI_vector.clear();
       	EcState.clear();
	delete hydraulicModel;
/*
        if (treesParams.useHydraulics == true)
        {
		free3DArray(psi, NREC, MD);	
		free2DArray(rflux, NREC);
        }
*/

} //end simulator

//------------------------------------------------------------------------------------------------------------
//This is called for single simulations with fixed parameters, returning results in frames
//Use this function for simulating replicates, such as those drawn from MCMC posterior parameter distributions
//This function is called by getPIL.cpp
//I don't recommend calling the full plant water balance routines from this call.
//------------------------------------------------------------------------------------------------------------
void simulator(Data_Store& in_data,
		 State_Store& state,
		 Parameter* current_params,
		 trees_params treesParams,
		 Data_Frame& Ecframe,
		 Data_Frame& NEEframe,
		 double& sd1,
		 double& sd2)
{
	int n_tsteps = in_data.no_of_steps();  //# of time steps in input data
	int i, j, ts;
	Time ti;
	double u_ref, t_ref, D_ref, precip, Qpar, t_canopy, D_canopy, CO2_atm, Ecobs;
	double p_atm, t_surface, t_soil, t_root, Zw, NEEobs;
	double xylem_capacity_multiplier;
	double thetaSoil[ULAT], tempSoil[ULAT];
	struct sim_out simOut;
        std::vector<double> GSI_vector(22);
//These variables interface with the HydraulicModel model
        std::vector<double> EcState;
        std::vector<double> EcStateInner(50);
        std::vector<double> LeafPsiInner(50);
        std::vector<double> KpInner(50);
        std::vector<int> InnerStep(50);
        int ModelStatus = 0;
        double Ecrit = 0.0;
        double Thresh = 0.0;

	for (i = 1; i <= 21; ++i)
	{
                GSI_vector[i] = 0.0;
	}

//declare and assign parameters
//Note: position in the .p file is important as name matching is not done
        i = 0;
//Radiation parameters
        treesParams.altitude = current_params[i].get_value(); i++;  //altitude (meters)
        treesParams.lat = current_params[i].get_value(); i++;  //latitude (degrees)
        treesParams.longi = current_params[i].get_value(); i++; //longitude (degrees)
        treesParams.z_ref = current_params[i].get_value(); i++;  //ref height, m
        treesParams.lai = current_params[i].get_value(); i++;  //single sided lai
	treesParams.live_lai = treesParams.lai;
	treesParams.Al = max(1.0, treesParams.lai);
	treesParams.lai_at_sat_kl = treesParams.Al;
        treesParams.canopy_ht = current_params[i].get_value(); i++; //height of canopy, m
        treesParams.lai_at_full_canopy_height = current_params[i].get_value(); i++; //fraction cloud cover
        treesParams.l_angle = current_params[i].get_value(); i++; //leaf angle distribution
        treesParams.canopy_e = current_params[i].get_value(); i++;  //canopy emissivity
        treesParams.fPAR_beam = current_params[i].get_value(); i++;  //fraction of PAR in beam radiation
        treesParams.fPAR_diff = current_params[i].get_value(); i++;  //fraction of PAR in diffuse radiation
        treesParams.alpha_PAR = current_params[i].get_value(); i++;  //alpha for PAR & NIR typically 0.8 & 0.2
        treesParams.alpha_NIR = current_params[i].get_value(); i++;
        treesParams.omega = current_params[i].get_value(); i++;  //parameter used to adjust ext. coefficient (0 to 1)
        treesParams.p_crown = current_params[i].get_value(); i++;  //parameter used in clumping factor calculation (1 to 3.34)
//Aerodynamic parameters
        treesParams.d_factor = current_params[i].get_value(); i++;  //d = d_factor*canopy_ht
        treesParams.zm_factor = current_params[i].get_value(); i++;  //zm = zm_factor*canopy_ht
        treesParams.zh_factor = current_params[i].get_value(); i++;  //zh = zh_factor*zm
	treesParams.ps_model = (int) current_params[i].get_value(); i++;  //photosynthesis model to select (1,2, or 3)
        treesParams.Rd_mult = current_params[i].get_value(); i++;  //Rd = Rd_mult * Vmax
//Photosynthesis parameters
        treesParams.Jmax_mult = current_params[i].get_value(); i++;  //Jmax = Jmax_mult * Vmax
        treesParams.thetaJ = current_params[i].get_value(); i++;  //J curvature parameter
        treesParams.phiJ_sun = current_params[i].get_value(); i++;  //effective quantum yield of sunlit PSII system, e-/umol
        treesParams.phiJ_shd = current_params[i].get_value(); i++;  //effective quantum yield of shaded PSII system, e-/umol
        treesParams.Nleaf = current_params[i].get_value(); i++;  //leaf N concentration (kg/m2)
        treesParams.N_fixed_proportion = current_params[i].get_value(); i++;  //leaf N concentration (kg/m2)
        treesParams.Nrubisco = current_params[i].get_value(); i++;  //leaf proportion of N in rubisco
        treesParams.Kc25 = current_params[i].get_value(); i++;  //(Pa) MM const carboxylase, 25 deg C
        treesParams.q10Kc = current_params[i].get_value(); i++;  //(DIM) Q_10 for kc
        treesParams.Ko25 = current_params[i].get_value(); i++;  //(Pa) MM const oxygenase, 25 deg C
        treesParams.q10Ko = current_params[i].get_value(); i++;  //(DIM) Q_10 for ko
        treesParams.act25 = current_params[i].get_value(); i++;  //(umol/mgRubisco/min) Rubisco activity
        treesParams.q10act = current_params[i].get_value(); i++;  //(DIM) Q_10 for Rubisco activity
//Added for C4 - DSM August 2019
        treesParams.Vcmax25 = current_params[i].get_value(); i++;  //maximum Rubisco activity at 25 C, umol m-2 s-1
        treesParams.Vpmax25 = current_params[i].get_value(); i++;  //maximum PEP carbolylase activity at 25 C, umol m-2 s-1
        treesParams.Jmax25 = current_params[i].get_value(); i++;  //maximum electron transport rate at 25 C, umol m-2 s-1
	treesParams.gammaStar25 = current_params[i].get_value(); i++;  //compensation point, ubar
        treesParams.Kp25 = current_params[i].get_value(); i++;  //Michaelis constant of PEP carboxylase for CO2 at 25 C, ubar
        treesParams.Vpr = current_params[i].get_value(); i++;  //PEP regeneration rate, umol m-2 s-1
        treesParams.f = current_params[i].get_value(); i++;  //correction for spectral quality of light
        treesParams.x = current_params[i].get_value(); i++;  //partitioning factor of electron transport rate
        treesParams.absorptance = current_params[i].get_value(); i++;  //fraction of irradiance absorbed
        treesParams.E_Vcmax = current_params[i].get_value(); i++;  //activation energy, maximum carboxylation rate, kJ mol-1
        treesParams.E_Vpmax = current_params[i].get_value(); i++;  //activation energy, maximum PEP rate, kJ mol-1
        treesParams.E_Jmax = current_params[i].get_value(); i++;  //activation energy, electron transport, kJ mol-1
        treesParams.E_Kp = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of PEP, kJ mol-1
	treesParams.E_Kc = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of carboxylation, kJ mol-1
        treesParams.E_Ko = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of oxygenation, kJ mol-1
        treesParams.E_Rd = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of mitochondrial respiration, kJ mol-1
        treesParams.E_gammaStar = current_params[i].get_value(); i++;  //activation energy, Michaelis reaction of compensation point, kJ mol-1
        treesParams.gm = current_params[i].get_value(); i++;  //mesophyll conductance to CO2, mol m-2 s-1
        treesParams.gbs = current_params[i].get_value(); i++;  //conductance of the bundle sheath, mol m-2 s-1
        treesParams.alphaGmax = current_params[i].get_value(); i++;  //fraction of glycolate carbon diverted to glycine during photorespiration
        treesParams.alphaSmax = current_params[i].get_value(); i++;  //fraction of glycolate carbon diverted to serine during photorespiration
        treesParams.Nmax = current_params[i].get_value(); i++;  //maximum rate of de novo nitrogen supply to the chloroplast, umol N m-2 s-1

//Plant hydraulic Parameters
        treesParams.Gsref0 = current_params[i].get_value(); i++;  //maximum gs, molm-2s-1 (leaf hemisurface area)
        double m = current_params[i].get_value(); i++;  //kPa-1 as D in kPa
	treesParams.isAmphistomatous = (bool) current_params[i].get_value(); i++; //stomata on two sides
        treesParams.Md = current_params[i].get_value(); i++;
        //treesParams.Al = current_params[i].get_value(); i++;
        treesParams.midday_at_sat_kl = current_params[i].get_value(); i++;
        treesParams.e_at_saturated_kl= current_params[i].get_value(); i++;
        treesParams.rhizosphere_width= current_params[i].get_value(); i++;
        //treesParams.E_inc= current_params[i].get_value(); i++;
	treesParams.E_inc=0.0001;
        treesParams.soilshells= (int) current_params[i].get_value(); i++;
        treesParams.GMP= current_params[i].get_value(); i++;
        treesParams.GSD= current_params[i].get_value(); i++;
        treesParams.BD= current_params[i].get_value(); i++;
        treesParams.porosity = current_params[i].get_value(); i++;  //
        treesParams.silt_fraction= current_params[i].get_value(); i++;
        treesParams.clay_fraction= current_params[i].get_value(); i++;
	treesParams.residual= current_params[i].get_value(); i++;
        treesParams.frac_absorbing_length= current_params[i].get_value(); i++;
        treesParams.Capacitance= current_params[i].get_value(); i++;
        treesParams.axK_latKl_shoot_modules= current_params[i].get_value(); i++;
        treesParams.axKr_latKr_root_modules= current_params[i].get_value(); i++;
        treesParams.per_total_R_in_root_system= current_params[i].get_value(); i++;
        //treesParams.lai_at_sat_kl= current_params[i].get_value(); i++;
        treesParams.pd_at_sat_kl= current_params[i].get_value(); i++;
//calculate kl from hydraulic parameters, kl=e/(pd-md)
	treesParams.saturated_kl_for_whole_plant = treesParams.e_at_saturated_kl/(treesParams.pd_at_sat_kl-treesParams.midday_at_sat_kl);

        treesParams.ax_Shoot_b_value= current_params[i].get_value(); i++;
        treesParams.ax_Shoot_c_value= current_params[i].get_value(); i++;
        treesParams.lat_Shoot_b_value= current_params[i].get_value(); i++;
        treesParams.lat_Shoot_c_value= current_params[i].get_value(); i++;
        treesParams.ax_Root_b_value= current_params[i].get_value(); i++;
        treesParams.ax_Root_c_value= current_params[i].get_value(); i++;
        treesParams.lat_Root_b_value= current_params[i].get_value(); i++;
        treesParams.lat_Root_c_value= current_params[i].get_value(); i++;
        treesParams.initial_conductivity_root= current_params[i].get_value(); i++;
        treesParams.decrement_root= current_params[i].get_value(); i++;
        treesParams.initial_conductivity_shoot= current_params[i].get_value(); i++;
        treesParams.decrement_shoot= current_params[i].get_value(); i++;
//Biogeochemical cycling parameters
        treesParams.theta_opt = current_params[i].get_value(); i++;  //
        treesParams.optimal_soil_T = current_params[i].get_value(); i++;  //
        treesParams.growth_resp_proportion = current_params[i].get_value(); i++;  //
        treesParams.resp_coef_root = current_params[i].get_value(); i++;  //
        treesParams.resp_coef_stem = current_params[i].get_value(); i++;  //
        treesParams.resp_coef_leaf = current_params[i].get_value(); i++;  //
        treesParams.resp_coefficient = current_params[i].get_value(); i++;  //
        treesParams.EaSx = current_params[i].get_value(); i++;  //
        treesParams.kMsx = current_params[i].get_value(); i++;  //
        treesParams.xASx = current_params[i].get_value(); i++;  //

	treesParams.kd = current_params[i].get_value(); i++;  //
        treesParams.kn = current_params[i].get_value(); i++;  //
        treesParams.kea = current_params[i].get_value(); i++;  //
        treesParams.kes = current_params[i].get_value(); i++;  //
        treesParams.kl = current_params[i].get_value(); i++;  //
        treesParams.kh = current_params[i].get_value(); i++;  //

	treesParams.fr_minCN = current_params[i].get_value(); i++;  //
        treesParams.fr_maxCN = current_params[i].get_value(); i++;  //
        treesParams.leaf_minCN = current_params[i].get_value(); i++;  //
        treesParams.leaf_maxCN = current_params[i].get_value(); i++;  //

        treesParams.Cbelowground = current_params[i].get_value(); i++;  //
        treesParams.Clitter_frac = current_params[i].get_value(); i++;  //
        treesParams.Croot_frac = current_params[i].get_value(); i++;  //
        treesParams.Clitter = treesParams.Clitter_frac * treesParams.Cbelowground;
        treesParams.Croot = treesParams.Croot_frac * treesParams.Cbelowground;
        treesParams.Cstem = current_params[i].get_value(); i++;  //
        treesParams.Csapwood = current_params[i].get_value(); i++;  //
        treesParams.Croot_coarse_frac = current_params[i].get_value(); i++;//
        treesParams.Croot_coarse = treesParams.Croot_coarse_frac * treesParams.Cbelowground;
        treesParams.Csoil = (1.0-treesParams.Clitter_frac - treesParams.Croot_frac - treesParams.Croot_coarse_frac) * treesParams.Cbelowground;
	treesParams.interception_per_leafArea = current_params[i].get_value(); i++;//
        treesParams.litter_capacity = current_params[i].get_value(); i++;//
	treesParams.litter_capacity_init = treesParams.litter_capacity;
        treesParams.theta_deep0 = current_params[i].get_value(); i++;//
        treesParams.theta_mid0 = current_params[i].get_value(); i++;//
        treesParams.theta_shallow0 = current_params[i].get_value(); i++;//
        treesParams.litter_store0 = current_params[i].get_value(); i++;//
        treesParams.SLA = current_params[i].get_value(); i++;  //
	treesParams.SLA_instant = treesParams.SLA;
        treesParams.SRL1 = current_params[i].get_value(); i++;  //
        treesParams.minRootDiam = current_params[i].get_value(); i++;  //
        treesParams.maxRootDiam = current_params[i].get_value(); i++;  //
	treesParams.rootDiamMultiplier = pow(treesParams.maxRootDiam/treesParams.minRootDiam,1.0/9.0);
        treesParams.minRootLifespan = current_params[i].get_value(); i++;  //
        treesParams.LWP_spring_minimum = current_params[i].get_value(); i++;  //
        treesParams.LWP_stomatal_closure = current_params[i].get_value(); i++;  //
        treesParams.is_bryophyte = (int) current_params[i].get_value(); i++;  //
        treesParams.capRiseScalar = current_params[i].get_value(); i++;  //
//How to to multiple precipitation by for experimental drought designs
        double precipReduction = current_params[i].get_value(); i++;  //
        if (precipReduction < 0.0)
                precipReduction = 0.0;
//Threshold for Ec convergence
        treesParams.drainScalar= current_params[i].get_value(); i++;
        treesParams.leafNSCscalar= current_params[i].get_value(); i++;
//What are these for?
        treesParams.usePhenology= (bool) current_params[i].get_value(); i++;
        treesParams.leafLifeSpan= current_params[i].get_value(); i++;
//Maximum number of iterations to achieve convergence of DELTA
        treesParams.max_iterations = (int) current_params[i].get_value(); i++;
//Conductance to use for Darcy's Law, 1=WholePlant,2=AxialComponents, 3=Shoot,4=AxialRoot, 5=LaterialRoot
        treesParams.microbiomeScalar = current_params[i].get_value(); i++;
// MCH 07082020
	treesParams.microbialrainrate = current_params[i].get_value(); i++; //
//MCH 23092020
        treesParams.raininAmmonium = current_params[i].get_value(); i++; //
        treesParams.raininNitrate = current_params[i].get_value(); i++; //
        treesParams.raininMineralN = current_params[i].get_value(); i++; //
        treesParams.raininLabileC = current_params[i].get_value(); i++; //
//Parameters for snowpack accumulation and melt
        treesParams.snowpack_water_equivalent = current_params[i].get_value(); i++; //
        treesParams.snowpack_E_deficit_max = current_params[i].get_value(); i++; //
        treesParams.melt_Rcoef = current_params[i].get_value(); i++; //
//Run with full hydraulics or not (1 = yes, 0 = no)
        treesParams.useHydraulics = (bool) current_params[i].get_value(); i++; //
        treesParams.useInputStress = (bool) current_params[i].get_value(); i++; //
        treesParams.useInputWaterTable = (bool) current_params[i].get_value(); i++; //
        treesParams.dayToStopMaizeRefilling = current_params[i].get_value(); i++; //
//Leaf growth module parameters
        treesParams.useLeafModule = (bool) current_params[i].get_value(); i++; //
        treesParams.leafAreaMax = current_params[i].get_value(); i++; //
        treesParams.initialLeafSize = current_params[i].get_value(); i++; //
        treesParams.leafArea_Rate = current_params[i].get_value(); i++; //
        treesParams.dur_LeafExpansion= current_params[i].get_value(); i++; //
        treesParams.SLA_max = current_params[i].get_value(); i++; //
    	treesParams.SLA_min = current_params[i].get_value(); i++; //
        treesParams.leaf_insertAngle = current_params[i].get_value(); i++; //
    	treesParams.leaf_len_to_width = current_params[i].get_value(); i++; //
    	treesParams.proportion_CD = current_params[i].get_value(); i++; //
        treesParams.phyllochron = current_params[i].get_value(); i++; //
        treesParams.floweringTime = current_params[i].get_value(); i++; //
        treesParams.Tbase = current_params[i].get_value(); i++; //
        treesParams.therm_plant = current_params[i].get_value(); i++; //
    	treesParams.projectedArea_init = current_params[i].get_value(); i++; //
    	treesParams.pot_size = current_params[i].get_value(); i++; //
    	treesParams.root_to_shoot = current_params[i].get_value(); i++; //
    	treesParams.leaf_to_stem = current_params[i].get_value(); i++; //
        treesParams.useLeafGamma = (bool) current_params[i].get_value(); i++; //
        treesParams.Kalpha = current_params[i].get_value(); i++; //
        treesParams.Kbeta = current_params[i].get_value(); i++; //
        treesParams.Nalpha = current_params[i].get_value(); i++; //
        treesParams.Nbeta = current_params[i].get_value(); i++; //
        treesParams.ralpha = current_params[i].get_value(); i++; //
        treesParams.rbeta = current_params[i].get_value(); i++; //
//MCMC Parameters
        sd1 = current_params[i].get_value(); i++; //
        sd2 = current_params[i].get_value(); i++; //
        double sd_err1_wt = current_params[i].get_value();

        treesParams.delta = m*treesParams.Gsref0; // sensitivity of stomata to VPD, mol m-2 s-1

	//treesParams.live_lai = treesParams.live_lai * (1.0 - 0.15/treesParams.leafLifeSpan);

//Calculate Brooks-Corey ks (cm/hr), bubbling pressure (cm), pore-size distribution index, and residual water content
	double ks, bubbling_pressure, pore_size_index, residual, fieldCapacity;
	double pClay = 100.0*treesParams.clay_fraction;
        double pSand = 100.0*(1.0 - treesParams.silt_fraction - treesParams.clay_fraction);
	soil_texture_parameters(treesParams.porosity, pClay, pSand, ks, bubbling_pressure, pore_size_index, residual);
        treesParams.ks = ks;
        treesParams.bubbling_pressure = bubbling_pressure;
        treesParams.pore_size_index = pore_size_index;
        //treesParams.residual = residual;
	residual = treesParams.residual;
//Calculate field capacity
	fieldCapacity = field_capacity(bubbling_pressure, pore_size_index, residual, treesParams.porosity);
	treesParams.field_capacity = fieldCapacity;

//Initialize state variables here
	Var_Dictionary vd;
	Var_name deep_name = "deep";
	Var_name mid_name = "mid";
        Var_name shallow_name = "shallow";
	Var_name litter_name = "litter";
	Var_name canopy_name = "canopy";
        Var_name lwp_name = "lwp";
        Var_name kl_name = "kl";
        Var_name ecrit_name = "ecrit";
        Var_name psicrit_name = "psicrit";
        Var_name snowpack_name = "snowpack";
        Var_name snowEdef_name = "snowEdef";
	Var_name xylemCap_name = "xylem";

        vd.insert(deep_name);
        vd.insert(mid_name);
        vd.insert(shallow_name);
	vd.insert(litter_name);
	vd.insert(canopy_name);
        vd.insert(lwp_name);
        vd.insert(kl_name);
        vd.insert(ecrit_name);
        vd.insert(psicrit_name);
	vd.insert(snowpack_name);
        vd.insert(snowEdef_name);
	vd.insert(xylemCap_name);

        Var_name regen_name = "regen";
        vd.insert(regen_name);

        Var_name currentTargetLai_name = "currenttargetlai";
        vd.insert(currentTargetLai_name);
        Var_name forecastTargetLai_name = "forecasttargetlai";
        vd.insert(forecastTargetLai_name);

        state.allocate(vd);   //allocate memory
	state.set_val_at(treesParams.theta_deep0, DEEPSW);
	state.set_val_at(treesParams.theta_mid0, MIDSW);
        state.set_val_at(treesParams.theta_shallow0, SHALLOWSW);
	state.set_val_at(treesParams.litter_store0, LITTER);
	state.set_val_at(0.0, CANOPY);
        state.set_val_at(0.0, PSN);
	state.set_val_at(treesParams.lai/treesParams.SLA*10000.0*((0.14/treesParams.leafLifeSpan)+0.10)+treesParams.Croot*0.10 + treesParams.Croot_coarse*0.10 + 0.04*treesParams.Csapwood, NSC);
        state.set_val_at(treesParams.Cbelowground*treesParams.Croot_frac, ROOTC);
        state.set_val_at(treesParams.snowpack_water_equivalent, SNOWPACK);
        state.set_val_at(treesParams.snowpack_E_deficit_max, SNOWEDEF);
	state.set_val_at(treesParams.saturated_kl_for_whole_plant, KL);
        state.set_val_at(treesParams.lai, CURRENTTARGETLAI);
        state.set_val_at(0.0, FORECASTTARGETLAI);

//Initializing state variables for HydraulicModel
        for(j=0; j < n_tsteps; j++)
        {
                EcState.push_back(0.0);
        }
        //set to 50 because model should take nowhere near this long to converge.
        for(j=0; j < 50; j++)
        {
                EcStateInner.push_back(0.0);
                LeafPsiInner.push_back(0.0);
                KpInner.push_back(0.0);
                //RInner.push_back(0.0);

        }
        for(j=0; j < 50; j++)
        {
                InnerStep.push_back(0.0);
        }

//Initialize soil theta in each layer (up to ULAT layers)
        for (j = 0; j < treesParams.rmodules; j++)
        {
		if (j == 0 || j == 1)
                {       
                        thetaSoil[j] = treesParams.theta_shallow0;
                }
                else if (j == 2 || j == 3)
                {       
                        thetaSoil[j] = treesParams.theta_mid0;
                }
                else
                {       
                        thetaSoil[j] = treesParams.theta_deep0;
                }      
        }

//define a BiogeochemicalCycles object
	BiogeochemicalCycles bgc(treesParams);
//Define the canopy closure based on crown diameter
	treesParams.canopy_cover = bgc.computeCanopyCover(treesParams);
	double rootArea = bgc.computeRootArea(treesParams);
	bgc.computeLateralRootWeibulls(treesParams);

//set root to leaf area ratio
        treesParams.Ar_Al = rootArea / treesParams.lai / treesParams.canopy_cover;
	treesParams.aral_at_sat_kl = treesParams.Ar_Al;
        treesParams.Ar_Al_init = treesParams.Ar_Al;

//define a HydraulicModel object
        HydraulicModel *hydraulicModel = new HydraulicModel();

//define parameters for setup()
	bool reset = true;
        HydraulicModel::outputStruct ecrit_k_psi, k_p_e;
	double ecritOut[MD], pc[MD], klpred[MD], ppredOut[MD], epredOut[MD], nodeFail[MD];
        char *nodeTyp = new char[MD];
        double ***psi, **rflux;
        if (treesParams.useHydraulics == true)
        {
                psi = get3DArray(NREC,MD,ULAT);
                rflux = get2DArray(NREC,NMAX);
        }
        double soilpsiavg[NREC],  nts[NREC], evap[NREC], kroot[MD], axr[MD], latr[MD];
        double kshoot[MD], dslat[MD], dsax[MD], drlat[MD],drax[MD], ll[MD];
        double ksat[MD][MD], bsat[MD][MD], ccsat[MD][MD], newb[MD][MD], newc[MD][MD];
        double ip[MD][ULAT], b[MD][ULAT], b1[MD][ULAT], n[MD][ULAT], n1[MD][ULAT];
        double r[MD][ULAT], vu[MD][ULAT], vl[MD][ULAT], cf[MD][ULAT], pe[MD][ULAT], ws_[MD][ULAT];
        double ks2[MD][ULAT], p[MD][ULAT], dpp[MD][ULAT], jl[MD][ULAT], wnu[MD][ULAT], wu[MD][ULAT];
        double wnl[MD][ULAT], wl[MD][ULAT], cpu[MD][ULAT], cpl[MD][ULAT], cpp[MD][ULAT], ku[MD][ULAT];
        double kl[MD][ULAT], f[MD][ULAT], cc[MD][2], psinode[NMAX][NMAX], psimin[NMAX][NMAX];
        double jmatrix[NMAX][NMAX], jmatrix2[NMAX][NMAX];
        double percent[NMAX], rssoil[NMAX], rs[NMAX], dp[NMAX], ff[NMAX];
        int col[NMAX], row[NMAX],  indxx[NMAX];
        double subtract[NMAX], pressure[NKINC], plc[NKINC], plcweib[NKINC];
        double rsquare[NKINC], soilpsi[ULAT], al[MD], ar[MD], saturatedKs[MD][4];
        double ptarg, e=0.0, rr, rzw, rd, einc;
        int SoilRhizElements;
        double dt, gmd, gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract;
        double kla, aral, latot, ratot;
        int shootElements = treesParams.smodules, ktotal;
        double soilpsimin=0.0;
        int ShootModules, rootElements = treesParams.rmodules, RootModules, totmodules;
        double axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b, axRoot_c, latRoot_b, latRoot_c;
        double plco, pdecrement, sumy1, sumy2, sumprod, cincrement, soilvol=0.0;
        double rlateral, rLat_base, slateral, sLat_base, raxial, pleafave;
        int tnode;
        double fac = 0.0, modelledKL;
        int HydraulicModelFailCond=0;



//call setup() to initialize HydraulicModel objects
        if (treesParams.useHydraulics == true)
        {
        	hydraulicModel->setup(reset, ecrit_k_psi, k_p_e, ecritOut, pc, klpred, ppredOut, epredOut,
                nodeFail, nodeTyp, psi, rflux, soilpsiavg, nts, evap, kroot,
                axr, latr, kshoot, dslat, dsax, drlat, drax, ll, ksat, bsat,
                ccsat, newb, newc, ip, b, b1, n, n1, r, vu, vl, cf, pe, ws_, ks2,
                p, dpp, jl, wnu, wu, wnl, wl, cpu, cpl, cpp, ku, kl, f, cc, psinode,
                psimin, jmatrix, jmatrix2, percent, rssoil, rs, dp, ff, col, row,
                indxx, subtract, pressure, plc, plcweib, rsquare, soilpsi, al, ar,
                saturatedKs, ptarg, e, rr, rzw, rd, einc, SoilRhizElements, dt, gmd,
                gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract, kla, aral,
                latot, ratot, shootElements, ktotal, soilpsimin, ShootModules, rootElements,
                RootModules, totmodules, axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b,
                axRoot_c, latRoot_b, latRoot_c, plco, pdecrement, sumy1, sumy2, sumprod,
                cincrement, soilvol, rlateral, rLat_base, slateral, sLat_base, raxial,
                pleafave, tnode, fac, modelledKL, HydraulicModelFailCond, treesParams);
	}

        state.set_val_at(treesParams.pd_at_sat_kl, LPSI);
        state.set_val_at(treesParams.saturated_kl_for_whole_plant, KL);
        state.set_val_at(treesParams.e_at_saturated_kl, ECRIT);

//time loop
	for (i = 0; i < n_tsteps; i++)
	{
		if(in_data.check_step(i) == 0) //valid step w/ valid data
                {
			ts = i+1; //HydraulicModel initialization check, ts=0 initializes Gs
//get input data
			ti = in_data.get_time(i);
			j=0;
			u_ref = in_data.get_val_at(j,i); j++;
                	t_ref = in_data.get_val_at(j,i); j++;
                	D_ref = in_data.get_val_at(j,i); j++;
                	precip = in_data.get_val_at(j,i); j++;
			precip *= 0.001; //convert mm to m
			precip *= precipReduction;
                	Qpar = in_data.get_val_at(j,i); j++;
                	t_canopy = in_data.get_val_at(j,i); j++;
                	D_canopy = in_data.get_val_at(j,i); j++;
                	p_atm = in_data.get_val_at(j,i); j++;
                	CO2_atm = in_data.get_val_at(j,i); j++;
                	t_surface = in_data.get_val_at(j,i); j++;
                	t_soil = in_data.get_val_at(j,i); j++;
                	t_root = in_data.get_val_at(j,i); j++;
                	Zw = in_data.get_val_at(j,i); j++;
                	xylem_capacity_multiplier = in_data.get_val_at(j,i); j++;
			treesParams.xylemScalar = xylem_capacity_multiplier;
			if (treesParams.useInputStress == true)
			{
				treesParams.stressScalar = in_data.get_val_at(j,i); j++;
			}

			//treesParams.e_at_saturated_kl *= xylem_capacity_multiplier;
                	//ignore = in_data.get_val_at(j,i); j++;
                	NEEobs = in_data.get_val_at(j,i); j++;
                	Ecobs = in_data.get_val_at(j,i);

			bool silent = true;

//Initialize soil temp in each layer (up to ULAT layers)
                        for (j = 0; j < treesParams.rmodules; j++)
                        {
				if (j == 0)
                                {
                                        tempSoil[j] = 0.5*(t_surface+t_soil);
                                }
                                else if (j == 1)
                                {
                                        tempSoil[j] = t_soil;
                                }
                                else
                                {
                                        tempSoil[j] = t_root;
                                }

                        }

//call simulations functions for current time step
                        simOut = simulation_functions(
						silent,
                                                ts,
                                                Ecrit,
                                                Thresh,
                                                state,
						thetaSoil,
						tempSoil,
                                                ti,
                                                u_ref,
                                                t_ref,
                                                D_ref,
                                                precip,
                                                Qpar,
                                                t_canopy,
                                                D_canopy,
						p_atm,
                                                CO2_atm,
                                                t_surface,
                                                t_soil,
                                                t_root,
                                                Zw,
                                                treesParams,
                                                EcState,
                                                EcStateInner,
                                                GSI_vector,
                                                EcState[ts-1],
                                                ModelStatus,
                                                LeafPsiInner,
                                                KpInner,
                                                InnerStep,
						bgc,
						hydraulicModel,
						ecrit_k_psi, k_p_e, ecritOut, pc, klpred, ppredOut, epredOut,
                                                nodeFail, nodeTyp, psi, rflux, soilpsiavg, nts, evap, kroot,
                                                axr, latr, kshoot, dslat, dsax, drlat, drax, ll, ksat, bsat,
                                                ccsat, newb, newc, ip, b, b1, n, n1, r, vu, vl, cf, pe, ws_, ks2,
                                                p, dpp, jl, wnu, wu, wnl, wl, cpu, cpl, cpp, ku, kl, f, cc, psinode,
                                                psimin, jmatrix, jmatrix2, percent, rssoil, rs, dp, ff, col, row,
                                                indxx, subtract, pressure, plc, plcweib, rsquare, soilpsi, al, ar,
                                                saturatedKs, ptarg, e, rr, rzw, rd, einc, SoilRhizElements, dt, gmd,
                                                gsd, bkd, fs, fc, Fabs, cpplant, saxlat, raxlat, rfract, kla, aral,
                                                latot, ratot, shootElements, ktotal, soilpsimin, ShootModules, rootElements,
                                                RootModules, totmodules, axShoot_b, axShoot_c, latShoot_b, latShoot_c, axRoot_b,
                                                axRoot_c, latRoot_b, latRoot_c, plco, pdecrement, sumy1, sumy2, sumprod,
                                                cincrement, soilvol, rlateral, rLat_base, slateral, sLat_base, raxial,
                                                pleafave, tnode, fac, modelledKL, HydraulicModelFailCond);

                        state.set_val_at(xylem_capacity_multiplier, XYLEM);

//assign y variables in the data frames
			Ecframe.set_y_value(i, simOut.Ec);
			NEEframe.set_y_value(i, simOut.NEE);
			//NEEframe.set_y_value(i, simOut.A);
//clear vector
		} //end if
	} //end time loop
        GSI_vector.clear();
        EcState.clear();
        EcStateInner.clear();
        LeafPsiInner.clear();
        KpInner.clear();
        InnerStep.clear();
	state.clear_values();
        delete hydraulicModel;
        delete[] nodeTyp;
        if (treesParams.useHydraulics == true)
        {
                free3DArray(psi, NREC, MD);
                free2DArray(rflux, NREC);
        }
} //end simulator


