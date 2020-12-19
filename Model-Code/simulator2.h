//
//function definitions for calculating
//incoming global radiation,
//transmission within canopy a layer,
//and canopy conductances  - SS July 2002
//refer data_store.h for object definitions for intended use
//this is a consolidated version containing Ec simulation
//with various parameter configs
//Sim_Ec_t - constraint on gc : D and gsmax
//Sim_Ec_q - constraint on gc : light, D and gsmax
//Sim_Ec_hc - my heat cap version, shh... :)
//***********************************************************
#ifndef _SIMULATOR_H
#define _SIMULATOR_H

#include <vector>
#include <fstream>
#include <math.h>
//#include <random>
#include <algorithm>
#include "parameter.h"
#include "state_store.h"
#include "data.h"
#include "constants.h"
#include "data_store.h"

#define DEEPSW 0
#define MIDSW 1
#define SHALLOWSW 2
#define LITTER 3
#define CANOPY 4
#define PSN 5
#define NSC 6
#define ROOTC 7
#define LPSI 8
#define KL 9
#define ECRIT 10
#define PSICRIT 11
#define SNOWPACK 12
#define SNOWEDEF 13
#define XYLEM 14
#define REGEN 15
#define CURRENTTARGETLAI 16
#define FORECASTTARGETLAI 17

#define MD 21   //maximum number of root+canopy layers
#define SLAT 2  //maximum number of lateral shoot elements
#define ULAT 20     //maximum number of below ground lateral elements
#define IM 1.0   //set high enough to saturate limiting flux and pressure result
#define PI 3.141592653589793
#define MAXITS 50   //maximum number of iterations allowed in solvewater
#define NMAX 121    //maximum number of nodes
#define NREC 100000   //maximum number of output lines
#define NKINC 15000   //maximum number of increments of K for solving Weibulls
#define STPP 2.50662827465
#define ITMAX 100
#define EPS 0.0000003
#define MAXH 400   //maximum number of steady-state loops
#define LFUSION 334000.0    // Latent heat of fusion (KJ mm-1 m-2) 
#define LVAPOUR 2495000.0   // Latent heat of vapourization (KJ mm-1 m-2) 

using namespace std;

const double MAXDOUBLE = 1.7 * pow(10.0,300.0);
const double MINDOUBLE = -1.0 * 1.0/MAXDOUBLE;


double*** get3DArray(int x, int y, int z);
double** get2DArray(int x, int y);

void free3DArray(double ***array, int x, int y);
void free2DArray(double **array, int x);

struct trees_params {
	double altitude;
	double lat;
	double longi;
	double z_ref;
	double lai;
	double live_lai;
	double dead_lai;
	double dead_lai_drop_rate;
	double canopy_cover;
	double canopy_ht;
	double lai_at_full_canopy_height;
	double l_angle;
	double canopy_e;
	double fPAR_beam;
	double fPAR_diff;
	double alpha_PAR;
	double alpha_NIR;
	double omega;
	double p_crown;
	double d_factor;
	double zm_factor;
	double zh_factor;
	double Gsref0;
	double delta;
	bool isAmphistomatous;
	int    ps_model;
	double Rd_mult;
	double Jmax_mult;
	double thetaJ;
	double phiJ_sun;
	double phiJ_shd;
	double Nleaf;
	double N_fixed_proportion;
	double Nrubisco;
	double Kc25;
	double q10Kc;
	double Ko25;
	double q10Ko;
	double act25;
	double q10act;
	double Vcmax25;
	double Vpmax25;
	double Jmax25;
	double gammaStar25;
	double Kp25;
	double Vpr;
	double f;
	double x;
	double absorptance;
	double E_Vcmax;
	double E_Vpmax;
	double E_Jmax;
	double E_Kp;
	double E_Kc;
	double E_Ko;
	double E_Rd;
	double E_gammaStar;
	double gm;
	double gbs;
	double alphaGmax;
	double alphaSmax;
	double Nmax;
	double Md;
        double Al;
        double Ar_Al;
        double Ar_Al_init;
        double midday_at_sat_kl;
        double e_at_saturated_kl;
        double rhizosphere_width;
        double E_inc;
        int soilshells;
        double GMP;
        double GSD;
        double BD;
        double silt_fraction;
        double clay_fraction;
        double frac_absorbing_length;
        double Capacitance;
        double axK_latKl_shoot_modules;
        double axKr_latKr_root_modules;
        double per_total_R_in_root_system;
        double saturated_kl_for_whole_plant;
        double aral_at_sat_kl;
        double lai_at_sat_kl;
        double pd_at_sat_kl;
        int smodules;
        int rmodules;
        double al[MD];
        double dslat[MD];
        double dsax[MD];
        double ar[MD];
        double drlat[MD];
        double drax[MD];
        double ax_Shoot_b_value;
        double ax_Shoot_c_value;
        double lat_Shoot_b_value;
        double lat_Shoot_c_value;
        double ax_Root_b_value;
        double ax_Root_c_value;
        double lat_Root_b_value;
        double lat_Root_c_value;
        double initial_conductivity_root;
        double decrement_root;
        double initial_conductivity_shoot;
        double decrement_shoot;
	double drainScalar;
	double xylemScalar;
	double leafNSCscalar;
	double stressScalar;
	bool usePhenology;
	double leafLifeSpan;
	int max_iterations;
	double microbiomeScalar;
//Microbioal rain in rate MCH 05082020
	double microbialrainrate;
//..
//MCH 23092020
        double raininAmmonium;
        double raininNitrate;
        double raininMineralN;
        double raininLabileC;
//..
	double theta_opt;
	double optimal_soil_T;
	double growth_resp_proportion;
	double resp_coef_root;
	double resp_coef_stem;
	double resp_coef_leaf;
	double resp_coefficient;
	double EaSx; //for DAMM model
	double kMsx; //for DAMM model
	double xASx; //for DAMM model
	double kd; //for rhizosphere BGC model
	double kn; //for rhizosphere BGC model
	double kea; //for rhizosphere BGC model
	double kes;
	double kl; //for rhizosphere BGC model
	double kh; //for rhizosphere BGC model
	double fr_minCN; //minimum fine root C:N ratio
	double fr_maxCN; //maximum fine root C:N ratio
	double leaf_minCN; //minimum leaf C:N ratio
	double leaf_maxCN; //maximum leaf C:N ratio
	double Cbelowground;
	double Clitter_frac;
	double Croot_frac;
	double Clitter;
	double Croot;
	double Croot_coarse;
	double Csoil;
	double Cstem;
	double Csapwood;
	double Croot_coarse_frac;
	double interception_per_leafArea;
	double litter_capacity;
	double litter_capacity_init;
	double theta_deep0;
	double theta_mid0;
	double theta_shallow0;
	double litter_store0;
	double SLA;
	double SLA_instant;
	double SRL1;
	double minRootDiam;
	double maxRootDiam;
	double minRootLifespan;
	double rootDiamMultiplier;
	double porosity;
	double ks;
	double bubbling_pressure;
	double pore_size_index;
	double residual;
	double field_capacity;
	double LWP_spring_minimum;
	double LWP_stomatal_closure;
	int is_bryophyte;
	double capRiseScalar;
	double snowpack_water_equivalent;
	double snowpack_E_deficit_max;
	double melt_Rcoef;
	bool useHydraulics;
	bool useInputStress;
	bool useRefilling;
	bool useInputWaterTable;
	bool forceRefilling;
	int dayToStopMaizeRefilling;
	bool updatedHydraulics;
	bool allowLeafRefilling;
//NEW leaf module
        bool useLeafModule;
        double leafAreaMax; // K
        double initialLeafSize; //A_pot_in
        double leafArea_Rate; //r
        double dur_LeafExpansion;//d_exp
        double SLA_max; //SLA_max
        double SLA_min; //SLA_min
        double leaf_insertAngle; // leaf insertion angle
        double leaf_len_to_width; // leaf length to width ratio
        double proportion_CD; //a
        double phyllochron; //phyllochron
        double floweringTime; //TTF
        double Tbase; //Tb
        double therm_plant; //therm_plant
        double projectedArea_init; // projected shoot area at initiation
        double pot_size; //max projected area
        double root_to_shoot;
        double leaf_to_stem;
        bool useLeafGamma;
        double Kalpha;
        double Kbeta;
        double Nalpha;
        double Nbeta;
        double ralpha;
        double rbeta;

};

//TREES - proposed new biogeochemical cycling routines
//define a class of objects to keep track of soil and plant nutrients
//we will start with carbon and nitrogen, both structural and non-structural
//we will assume soil pH is vertically uniform
//starting with some simplifications, e.g. one leaf and multiple root NSC pools,
//  phloem transport depends on xylem K and soil-sink strength,
//  microbial activity is concentrated in the rhizosphere,
//  labile carbon exits the roots as a function root K and source-sink strength,
//  N uptake by roots as a function of root K and soil-sink strength
//conceptual elements based on Porporato et al 2003 Advances in Water Resources


class BiogeochemicalCycles
{
	public:
		BiogeochemicalCycles();
		BiogeochemicalCycles(trees_params& treesParams);
		~BiogeochemicalCycles();

	private:
	//methods to allocate and initialize pools
		void allocatePools();
		void initializePools(trees_params treesParams);
		void clearPools();

	public:
	//methods to compute photosynthesis
		double C4photosynthesis(trees_params treesParams,
                                        double Ca,
                                        double Iin,
					double phiJ,
                                        double Tl,
                                        double gc,
                                        double& Ci,
                                        double& Cm);
		void coupledA4_gc(trees_params treesParams,
                                  double Ca,
                                  double Iin,
				  double phiJ,
                                  double Tl,
                                  double gc0,
                                  double& gc,
                                  double& An,
                                  double& Ci,
                                  double& Cm);
		double C3photosynthesis(trees_params treesParams,
                                        double Ca,
					double Pa,
                                        double Iin,
					double phiJ,
                                        double Tl,
                                        double gc,
					double Nmax,
                                        double& Ci,
                                        double& Cc,
					double& NG,
					double& NS);
		void coupledA3_gc(trees_params treesParams,
                                  double Ca,
				  double Pa,
                                  double Iin,
				  double phiJ,
                                  double Tl,
                                  double gc0,
				  double Nmax,
                                  double& gc,
                                  double& An,
                                  double& Ci,
                                  double& Cc,
				  double& NG,
				  double& NS);
		void storeGlycineAndSerine(double NGsun,
					     double NGshd,
					     double NSsun,
					     double NSshd,
					     double Lsun,
					     double Lshd);
		double photosynthesis(trees_params treesParams,
                                      double Jmax_mult,
                                      double thetaJ,
                                      double phiJ,
                                      struct farqin in,
                                      struct farqout& out);

	//methods to compute absorbing root area
		double computeRootArea(trees_params treesParams);
		double computeRootArea(trees_params treesParams,
				       int j);
		double computeRootArea(trees_params treesParams,
				       double SRA,
				       int j,
				       int k);

	//methods to compute absorbing fine root area
		double computeFineRootArea(trees_params treesParams);
		double computeFineRootArea(trees_params treesParams,
					   int j);
		double computeFineRootArea(trees_params treesParams,
					   double SRA,
					   int j,
					   int k);

	//method to compute maximum allowable carbon allocation to each root
		void computeMaximumRootBiomassCarbon(trees_params treesParams);

	//method to compute canopy cover as a proportion of the maximum ground area based on lateral root extension
		double computeCanopyCover(trees_params treesParams);

	//method to compute lateral root Weibull parameters based on root area
		void computeLateralRootWeibulls(trees_params& treesParams);
		void computeLateralRootWeibulls(int j, 
                                                        trees_params treesParams,
                                                        double& bsum,
                                                        double& csum);
		void updateLateralRootWeibulls(double bsat[][MD], 
                                                        double ccsat[][MD],
                                                        trees_params treesParams);
	//method to query leaf carbon
		double getLeafBiomassCarbon();
	//method to query fruit carbon
		double getFruitCarbon();

	//methods to query root carbon
		double getRootCarbon();
		double getRootCarbon(int j);
		double getRootCarbon(int j,
				     int k);
		double getFineRootCarbon();
		double getFineRootCarbon(int j);
		double getFineRootCarbon(int j,
				         int k);

	//method to query live stem carbon
		double getLiveStemCarbon();

	//method to query dead stem carbon
		double getDeadStemCarbon();

	//method to query leaf NSC
		double getLeafNSC();

	//method to query chloroplast starch
		double getChloroplastStarch();

	//method to query chloroplast sugar
		double getChloroplastSugar();

	//method to query stem NSC
		double getStemNSC();

	//methods to query root NSC
		double getRootNSC();
		double getRootNSC(int j);
		double getRootNSC(int j,
				  int k);

	//methods to query fine root NSC
		double getFineRootNSC();
		double getFineRootNSC(int j);
		double getFineRootNSC(int j,
				  int k);

	//method to assign a value to leaf NSC
		void putLeafNSC(double nsc);

	//method to modify the starch content of the chloroplast
		void putChloroplastStarch(double nsc);

	//method to modify the sugar content of the chloroplast
		void putChloroplastSugar(double nsc);

	//method to assign a value to stem NSC
		void putStemNSC(double nsc);

	//method to update values in root NSC
		void updateRootNSC(double delta_nsc);
		void updateRootNSC(double delta_nsc,
				   int j,
				   int k);

	//method to query humus carbon
		double getHumus();
		double getHumus(int j);

	//method to query rhizosphere litter carbon
		double getRhizosphereCl();
		double getRhizosphereCl(int j);
		double getRhizosphereCl(int j,
					int k);

	//method to query rhizosphere litter nitrogen
		double getRhizosphereNl();
		double getRhizosphereNl(int j);
		double getRhizosphereNl(int j,
					int k);

	//method to query amino acid exudate carbon pool
		double getAminoAcidExudateC();
		double getAminoAcidExudateC(int j);
		double getAminoAcidExudateC(int j,
					    int k);

	//method to query amino acid exudate nitrogen pool
		double getAminoAcidExudateN();
		double getAminoAcidExudateN(int j);
		double getAminoAcidExudateN(int j,
					    int k);

	//method to query sugar exudate carbon pool
		double getSugarExudateC();
		double getSugarExudateC(int j);
		double getSugarExudateC(int j,
					int k);

        //method to query rhizosphere live microbial carbon
                double getRhizosphereLiveMicrobialCarbon();
                double getRhizosphereLiveMicrobialCarbon(int j);
                double getRhizosphereLiveMicrobialCarbon(int j,
                                        		 int k);

        //method to query rhizosphere microbial nitrogen
                double getRhizosphereMicrobialNitrogen();
                double getRhizosphereMicrobialNitrogen(int j);
                double getRhizosphereMicrobialNitrogen(int j,
                                        		int k);

	//method to query rhizosphere dead microbial carbon
                double getRhizosphereDeadMicrobialCarbon();
		double getRhizosphereDeadMicrobialCarbon(int j);
                double getRhizosphereDeadMicrobialCarbon(int j,
                                        		 int k);

        //method to query rhizosphere mineral nitrogen
                double getRhizosphereMineralNitrogen();
                double getRhizosphereMineralNitrogen(int j);
                double getRhizosphereMineralNitrogen(int j,
                                        	     int k);

        //method to query rhizosphere ammonium nitrogen
                double getRhizosphereAmmoniumNitrogen();
                double getRhizosphereAmmoniumNitrogen(int j);
                double getRhizosphereAmmoniumNitrogen(int j,
                                        	     int k);

        //method to query rhizosphere nitrate nitrogen
                double getRhizosphereNitrateNitrogen();
                double getRhizosphereNitrateNitrogen(int j);
                double getRhizosphereNitrateNitrogen(int j,
                                        	     int k);

	//method to query nitrogen leaching flux variable
		double getNitrogenLeaching(int j);

	//method to compute NSC fluxes
		void computeNSCfluxes(trees_params treesParams,
				      double kratio,
				      double* kratio_vector);

	//method to compute N fluxes
		void computeNfluxes(trees_params treesParams,
				    double kratio,
			            double* kratio_vector);

	//method to compute root exudates
		void computeRootExudates(trees_params treesParams,
					 double* thetaSoil,
					 double* kratio_vector);
		void computeRootExudates(trees_params treesParams,
					 double* thetaSoil,
					 double* kratio_vector,
                                         int j);
		void computeRootExudates(trees_params treesParams,
					 double* thetaSoil,
					 double* kratio_vector,
                                         int j,
                                         int k);

	//method to update leaf carbon pools
		void computeLeafAllocation(trees_params& treesParams,
                                           double newc[][MD],
                                           double N_avail_rate_plant,
                                           double RL,
                                           double nsc,
                                           double psn,
                                           double nscRatio,
                                           double& rgrowth,
                                           double& leafCfraction,
                                           double lai,
					   int yday,
                                           double& stressedLeafLifeSpan,
                                           double SLA);

	//method to update lef carbon pools
		void updateLeafCarbonNitrogenPools(int k,
						   trees_params& treesParams,
					   	   double delta_lai,
						   double RL,
						   double N_neg_fract,
                                                   double& N_neg_demand,
                                                   double& N_pos_demand);
		void updateLeafCarbonNitrogenPools(trees_params& treesParams,
					   	   double delta_lai,
						   double RL,
						   double N_neg_fract,
                                                   double& N_neg_demand,
                                                   double& N_pos_demand);

	//method to update stem carbon pools
		void updateStemCarbonNitrogenPools(trees_params treesParams,
						   double N_avail_rate_plant,
						   double kratio,
					   	   double nscRatio,
                                           	   double rgrowth,
						   double& stemAllocation,
						   double N_neg_fract,
						   double& N_neg_demand,
						   double& N_pos_demand,
                                           	   double leafCfraction,
						   double t_canopy,
                                           	   double lai,
                                           	   double lifeSpan,
                                           	   double SLA);

	//method to update root carbon and nitrogen pools
		void updateRootCarbonNitrogenPools(trees_params treesParams,
					   	   double* tempSoil,
                                           	   double rgrowth,
						   double N_neg_fract,
						   double& N_neg_demand,
						   double& N_pos_demand,
						   double* kratio_vector,
                                           	   double stressedLeafLifeSpan,
                                           	   double unstressedLeafLifeSpan,
                                           	   double N_avail_rate_plant,
                                           	   double leafCfraction,
                                           	   double stemAllocation);

	//method to compute root maintenance respiration
		double root_respiration_rate(double resp_coef_root,
                               		     double resp_coefficient,
                                             double t_soil,
                                             double Croot,
                                             double transport_rate);

	//method to compute stem maintenance respiration
		double stem_respiration_rate(double resp_coef_stem,
                                             double resp_coefficient,
                                             double t_canopy,
                                             double Cstem);
	//method to compute leaf maintenance respiration
		double leaf_respiration_rate(double resp_coef_leaf,
                                             double resp_coefficient,
                                             double tLeaf,
                                             double lai,
                                             double SLA);

	//methods to get rhizosphere N totals
		void getRhizosphereN(double& N_neg,
				     double& N_pos);
		void getRhizosphereN(double& N_neg,
				     double& N_pos,
				     int j);


		double getRhizosphereN_neg(int j,
					   int k);

		double getRhizosphereN_pos(int j,
					   int k);

	//methods to get plant N totals
		
		double getPlantN();

		double getLeafBiomassN();
		double getLeafMineralN();
		double getRootN();
		double getRootN(int j);
		double getRootN(int j,
				int k);
		double getFineRootN();
		double getFineRootN(int j);
		double getFineRootN(int j,
				    int k);

		double getRootBiomassN();
		double getRootBiomassN(int j);
		double getRootBiomassN(int j,
				       int k);
		double getFineRootBiomassN();
		double getFineRootBiomassN(int j);
		double getFineRootBiomassN(int j,
					   int k);

		void computeLeafNdemand(trees_params treesParams,
					double& delta_lai,
					double N_neg_fract,
					double& N_neg_demand,
					double& N_pos_demand);

	//methods to compute nitrogen transformations, NH4+, NO3-, N2O
		void computeRootNitrogenUptake(double UP_neg[][10],
                               		     double UP_pos[][10],
                               		     trees_params treesParams,
					     double* thetaSoil,
                               		     double* Rflux,
                               		     double ratot,
					     double field_capacity,
                               		     double DEM_neg,
                               		     double DEM_pos);
		void computeLeaching(double LE_neg[][10],
				     double LE_pos[][10],
				     trees_params treesParams,
				     double* thetaSoil,
				     double* depthSoil,
				     double* bypassFlow);

	//methods to compute the volume of rhizosphere
		double computeRhizosphereVolume(int rootNum,
						double porosity,
						double rootRadius,
						double rhizosphereWidth);
		double computeRhizosphereVolume(trees_params treesParams,
					        int rootNum);
		double computeRhizosphereVolume(trees_params treesParams,
					        int rootNum,
					        int rootOrder);
		double computeRhizosphereVolume(int rootNum,
						int rootOrder,
						double porosity,
						double rootRadius,
						double rhizosphereWidth);

		void updateRhizospherePools(trees_params treesParams,
					    double* thetaSoil,
					    double* tempSoil,
					    double* depthSoil,
					    double UP_neg[][10],
					    double UP_pos[][10],
					    double LE_neg[][10],
					    double LE_pos[][10]);
		double computeMineralizationImmobilization(double& MIN,
							   double& PHI,
							   double& IMM_pos,
							   double& IMM_neg,
							   double& DECl,
							   double& DECh,
							   double& DECea,
							   double& DECes,
							   double& rh, 
                                                   	   double& rr,
                                                   	   double& kd,
                                                   	   double& kn, 
                                                   	   double& fns,
							   int root_num,
							   int root_order,
							   double t_soil,
							   double theta,
							   double AtoV,
							   double AtoVbulk,
							   double depth,
							   trees_params);
/*
		double computeAmmonification();
		double computeNitrification();
*/
		double nitrification_water_function(double,
					    	    double);
		double slow_temp_decomp_rate(double t_soil,
					     double optimal_soil_T,
					     double slow_mineral_rate_parm);
		double fast_temp_decomp_rate(double t_soil,
					     double optimal_soil_T,
					     double fast_mineral_rate_parm);
		double slow_decomposition_rate( double theta_opt,
                                		double theta_soil,
                                		double t_soil,
                                		double optimal_soil_T,
                                		double slow_mineral_rate_parm,
                                		double slow_decomp_proportion);
		double fast_decomposition_rate(double theta_soil,
                                	       double t_soil,
					       trees_params);

		double DAMM_Cflux(trees_params treesParams,
                                double t_soil,
                                double theta,
                                double porosity,
                                double depth,
                                double Clitter);

		double DAMM_reactionVelocity(trees_params treesParams,
                       			   double t_soil,
                       			   double theta,
                       			   double porosity,
                       			   double depth,
                       			   double Clitter);

	//methods to compute nitrogen fluxes and update pools
		double computeDenitrification();
		void computeRootNitrogenUptake();
		double computeNitrogenFlowToCanopy();
		double computerNitrogenFlowToRubisco();

	//methods to compute photosynthesis, NSC pools
		double computePhotosynthesis();
		double updateLeafNSC();
		double computePhloemNSC();
		double updateRootNSC();
		double computeLabileCarbonFlux();

	//methods to compute growth allocation, biomass turnover
		double computeLeafGrowth();
		double computeRootGrowth();
		double computeStemGrowth();
		double updateLitterPools();

	//methods to compute plant maintenance respiration
		double computeLeafMaintenance();
		double computeRootMaintenance();
		double computeStemMaintenance();

	//methods to compute microbial pools and litter decomposition
		double computeLitterDecomposition();
		double computeHumusDecomposition();


	//methods for nitrogen cycle

	public:
	//Proportion of N demand is currently being met by supply
		double* plantNstatus;
	//Heterotrophic respiration, CO2 flux
		double* heterotrophicRespiration;

	private:
	//Number of soil-root layers
		int nRoots;
		int nFineRootOrders;
		int nRootOrders;
	//State variables for plant residue inputs
		double*  leafResidueCarbon;
		double*  leafResidueNitrogen;
		double*  stemResidueCarbon;
		double*  stemResidueNitrogen;
		double** rootResidueCarbon;
		double** rootResidueNitrogen;
	//Flux variables for root residue inputs
		double** dCrootResidue;
		double** dNrootResidue;
	//State variables for stabilized soil carbon
		double* humusCarbon;
		double* humusNitrogen;
	//State variables for mineral nitrogen in the bulk soil
		double* soilAmmoniumNitrogen;
		double* soilNitrateNitrogen;
	//State variables for litter C and litter N
		double** rhizosphereCl;
		double** rhizosphereNl;
	//State variables for rhizosphere microbes, labile carbon, and immobile N, and mineral N
		double** rhizosphereLiveMicrobialCarbon;
		double** rhizosphereDeadMicrobialCarbon;
		double** rhizosphereLabileCarbon;
		double** rhizosphereMicrobialNitrogen;
		double** rhizosphereMineralNitrogen;
		double** rhizosphereAmmoniumNitrogen;
		double** rhizosphereNitrateNitrogen;
	//State variables for root exudates
		double** rootExudateSugarCarbon;
		double** rootExudateAminoAcidCarbon;
		double** rootExudateAminoAcidNitrogen;
	//State variables for roots
		double** maximumRootBiomassCarbon;
		double* rootAreaScalar;
		double** fineRootBiomassCarbon;
		double** fineRootBiomassNitrogen;
		double** coarseRootBiomassCarbon;
		double** coarseRootBiomassNitrogen;
		double** rootNSC;
		double** rootMineralNitrogen;
		double** rootArea;
	//State variables for stem
		double*  liveStemCarbon;
		double*  stemNSC;
		double*  liveStemNitrogen;
		double*  deadStemCarbon;
		double*  deadStemNitrogen;
	//State variables for leaf
		double*  leafBiomassCarbon;
		double*  leafBiomassNitrogen;
		double*  leafNSC;
		double*  chloroplastStarch;
		double*  chloroplastSugar;
		double*  glycineNitrogen;
		double*  glycineCarbon;
		double*  serineNitrogen;
		double*  serineCarbon;
		double*  leafStoredNitrogen;
		double*  leafRubiscoNitrogen;
	//State variables for fruit production
		double*  fruitCarbon;
		double*  fruitNitrogen;
	//Flux variable for nitrogen leaching
		double* nitrogenLeaching;
	//Initial input lateral root Weibull parameters
		double* lat_Root_b_value_init;
		double* lat_Root_c_value_init;

    // methods for leaf growth module
	public:
    	//leaf number trackers
    		int* Lf_idx;
    	//plant phyllochron tracker
    		double* phyllo_tracker;
    	// plant average SLA
    		double* SLA_avg;
    	// projected leaf area
    		double* Projected_LeafArea_total;
        // K leaf parameter array
            double* Karray;
        // N leaf parameter array
            double* Narray;
        // r leaf parameter array
            double* rarray;

    //methods to get leaf-level state variables
    		double getSingleLeafBiomassCarbon(int k);
    		double getSingleLeafBiomassNitrogen(int k);
    		double getSingleLeafNSC(int k);
    		double getSingleLeafchloroplastStarch(int k);
    		double getSingleLeafchloroplastSugar(int k);
    		double getSingleLeafStoredNitrogen(int k);
    		double getSingleLeafRubiscoNitrogen(int k);
    		double getSingleLeafArea(int k);
    		double getSingleLeafAreaPotential(int k);
    		double getSingleLeafThermTm(int k);
    //methods to put leaf-level state variables
    		void putSingleLeafAreaPotential(double area_pot,
						int idx);
    //methods to convert between leaf-level units to plant-level units
    		double cm2_to_m2(double value);
    		double cm2_to_ha(double value);
    		double g_to_kg(double value);
    //methods to compute relative growth rates
    		double calcRAR(double thermCD,
			       trees_params& treesParams, double* Karray, double* Narray, double* rarray, int idx);
    		double calcRER(double thermEXP,
			       trees_params& treesParams, double* Karray, double* Narray, double* rarray, int idx);
    //methods to update single leaf state variables
		void updateLeaf(int k,
				double delta_thermTime,
				trees_params& treesParams,
				double ProjGrndArea,
				double RL,
				double N_neg_fract,
				double& N_neg_demand,
				double& N_pos_demand, double* Karray, double* Narray, double* rarray, int idx, double thetaRoot);
    		double calcdeltaC(double SLA, 
				  double deltaAreaL);
    		double calcdeltaN(double deltaC, 
				  trees_params treesParams);
    //methods to compute biophysical constraints on relative growth rates (PLACEHOLDERS)
    		double calcHydRAR(double hydraulicStatus);
    		double calcNRAR(double NStatus);
    		double calcHydRER(double waterStat);
    		double calcNSCRER(double NSCStatus, 
				  double deltaC);
    //methods to compute SLA from N status
    		double calcSLA(double plantNitrogenStatus, 
			       double SLA_max, 
			       double SLA_min);
    //methods for computing growth respiration associated with division and expansion
    		double calcRespCD(double deltaAPot, 
				  double pseudoC_coef);
    		double calcRespEXP(double deltaC);
    //method for computing plant-level SLA (sum(single leaf areas)/sum(single leaf carbon))
    		double calcPlantSLA(int Lf_idx, 
				    trees_params treesParams);
    //methods to calculate plant LAI
    		double calcLAI(int Lf_idx, 
			       trees_params treesParams, double* Karray);
    //methods to calculate plant maximum projected ground area
    		double calc_projGndArea(trees_params treesParams, double* Karray);
    //methods to calculate plant total leaf area
    		double calc_projArea(int Lf_idx, 
			       trees_params treesParams);
    //methods for computing plant-level state variables from leaf-level state variables
    		double getSumLeafBiomassCarbon(int Lf_idx, 
					       double ProjGrndArea_instant);
    		double getSumLeafBiomassNitrogen(int Lf_idx, 
						 double ProjGrndArea_instant);
    		double getSumLeafNSC(int Lf_idx);
    		double getSumLeafchloroplastStarch(int Lf_idx);
    		double getSumLeafchloroplastSugar(int Lf_idx);
    		double getSumLeafStoredNitrogen(int Lf_idx);
    		double getSumLeafRubiscoNitrogen(int Lf_idx);
    //method to compute total leaf growth respiration
    		double getTotLeafGrowthRespiration(int Lf_idx, 
						   double ProjGrndArea_instant);
    //method to get upper and lower limits for truncated gamma distributions
            //double* getQuantileVals(double alpha, double beta, double min_quant, double max_quant);
            std::vector<double> getQuantileVals(double alpha, double beta, double min_quant, double max_quant);
    //method to sample from truncated gamma distribution
            double sampleTruncatedGamma(double alpha, double beta, double minlim, double maxlim);
    
	private:
    	//leaf-level state variables that directly link to plant level leaf state variables
    		double* SingleLeafBiomassCarbon;
    		double* SingleLeafBiomassNitrogen;
    		double* SingleLeafNSC;
    		double* SingleLeafchloroplastStarch;
    		double* SingleLeafchloroplastSugar;
    		double* SingleLeafStoredNitrogen;
    		double* SingleLeafRubiscoNitrogen;
    		double* SingleLeafGrowthRespiration;
    	//leaf-level state variables unique to single leaves
    		double* SingleLeafArea;
    		double* SingleLeafAreaPotential;
    		double* SingleLeafThermTm;
};

class HydraulicModel {

public:

        //Data_Store *data_store;

        struct outputStruct {
                int run[209];
                int module[209];
                double axialK[209];
                double latK[209];
                double rhizK[209];
                double rhizFlux[209];
        };
        

public:

        void solveWater(bool silent,
		bool callRewet,
		double Ypd[],
                double&  ev,
		double&  Kl,
		double&  leafpsi,
		double&  Ecrit,
		int &HydraulicModelFailCond,
		double &PsiCrit,
		int hydraulicModel_flag,
		int ts,
		trees_params& treesParams,
		outputStruct &ecrit_k_psi, outputStruct &k_p_e, double ecritOut[], 
		double pc[], double klpred[],
                double ppredOut[], double epredOut[], double nodeFail[], 
		char nodeTyp[], double ***psi,
                double **rflux, double soilpsiavg[], double nts[], double evap[], double kroot[],
                double axr[], double latr[], double kshoot[], double dslat[], double dsax[], double drlat[],
                double drax[], double l[], double ksat[][MD], double bsat[][MD], double ccsat[][MD],
                double newb[][MD], double newc[][MD], double ip[][ULAT], double b[][ULAT], double b1[][ULAT],
                double n[][ULAT], double n1[][ULAT], double r[][ULAT], double vu[][ULAT], double vl[][ULAT],
                double cf[][ULAT], double pe[][ULAT], double ws_[][ULAT], 
		double ks[][ULAT], double p[][ULAT],
                double dpp[][ULAT], double jl[][ULAT], double wnu[][ULAT], 
		double wu[][ULAT], double wnl[][ULAT],
                double wl[][ULAT], double cpu[][ULAT], double cpl[][ULAT], 
		double cpp[][ULAT], double ku[][ULAT],
                double kl[][ULAT], double f[][ULAT], double cc[][2], 
		double psinode[][NMAX], double psimin[][NMAX],
                double jmatrix[][NMAX], double jmatrix2[][NMAX], 
		double percent[], double rssoil[], double rs[],
                double dp[], double ff[], int col[], int row[], 
		int indxx[], double subtract[], double pressure[],
                double plc[], double plcweib[], double rsquare[], 
		double soilpsi[], double al[], double ar[],
                double saturatedKs[][4], double &ptarg, double &e, 
		double &rr, double &rzw, double &rd, double &einc,
                int &SoilRhizElements, double &dt, double &gmd, 
		double &gsd, double &bkd, double &fs, double &fc,
                double &Fabs, double &cpplant, double &saxlat, 
		double &raxlat, double &rfract, double &kla,
                double &aral, double &latot, double &ratot, int &shootElements, 
		int &ktotal, double &soilpsimin,
                int &ShootModules, int &rootElements, int &RootModules, 
		int &totmodules, double &axShoot_b,
                double &axShoot_c, double &latShoot_b, double &latShoot_c, 
		double &axRoot_b, double &axRoot_c,
                double &latRoot_b, double &latRoot_c, double &plco, double &pdecrement, double &sumy1,
                double &sumy2, double &sumprod, double &cincrement, double &soilvol, double &rlateral,
                double &rLat_base, double &slateral, double &sLat_base, double &raxial, double &pleafave,
                int &tnode, double fac, double &modelledKL);

public:  
      
        double getData(std::ifstream &in);

        void ksolve( double r[][ULAT], double pe[][ULAT], double b[][ULAT], double b1[][ULAT],
                    double n[][ULAT], double n1[][ULAT], double ks[][ULAT], double ws_[][ULAT],
                    double rssoil[], double cf[][ULAT], double vu[][ULAT], double vl[][ULAT],
                    double saturatedKs[][4], double cpp[][ULAT],
                    double wnu[][ULAT], double wu[][ULAT], double wnl[][ULAT], double wl[][ULAT],
                    double ip[][ULAT], double p[][ULAT], double al[],
                    double l[], double &rlateral, double &rLat_base,
                    double slateral, double sLat_base, double ratot,
                    double &raxial, double raxlat, const int smodules,
                    const int totmodules, double ar[], double axr[],
                    double drax[], double latr[], double drlat[],
                    double ksat[][MD], double soilpsiavg[], double rfract,
                    double &ptarg, double saxlat, double latot,
                    double dsax[], double dslat[], double &pleafave,
                    double psimin[][NMAX], double jl[][ULAT], double soilpsi[], double cpu[][ULAT],
                    double cpl[][ULAT], double f[][ULAT], double ff[NMAX], double jmatrix[][NMAX],
                    double dpp[][ULAT], double dp[NMAX], int col[], int indxx[], double subtract[],
                    const int rmodules, double rr, double rd,
                    double Fabs, double cpplant, double gmd,
                    double gsd, double fc, double fs, double bkd,
                    int &tnode, double e, double dt, double fac, const int SoilRhizElements,
                    double cc[][2], double kl[][ULAT], double ku[][ULAT], double **rflux,
		    double &modelledKL, double &soilvol, double initial_conductivity_root, 
		    double decrement_root, double initial_conducitivty_shoot, 
		    double decrement_shoot, int &HydraulicModelFailCond );

        //ifstream paramFile;

        void preSetup(double md_raw[], double Md, double e_raw[], double ev,  double la_raw[], double Al,
                      double arAl[], double Ar_Al,  double **l_raw, double Ypd[]);

   
public:

   	void setup( bool reset,
		    outputStruct &ecrit_k_psi,
                    outputStruct &k_p_e,
                    double ecritOut[],
                    double pc[],
                    double klpred[],  
                    double ppredOut[], 
                    double epredOut[], 
                    double nodeFail[], 
                    char nodeTyp[],   
                    double ***psi,
                    double **rflux,    
                    double soilpsiavg[],
                    double nts[],
                    double evap[],
                    double kroot[],  
                    double axr[],
                    double latr[],
                    double kshoot[],  
                    double dslat[], 
                    double dsax[],
                    double drlat[],
                    double drax[],
                    double l[], 
                    double ksat[][MD],
                    double bsat[][MD], 
                    double ccsat[][MD],
                    double newb[][MD], 
                    double newc[][MD],
                    double ip[][ULAT],
                    double b[][ULAT],  
                    double b1[][ULAT],  
                    double n[][ULAT],
                    double n1[][ULAT],
                    double r[][ULAT],
                    double vu[][ULAT],
                    double vl[][ULAT],
                    double cf[][ULAT],
                    double pe[][ULAT],
                    double ws_[][ULAT],
                    double ks[][ULAT],
                    double p[][ULAT],
                    double dpp[][ULAT],
                    double jl[][ULAT],
                    double wnu[][ULAT],
                    double wu[][ULAT], 
                    double wnl[][ULAT],
                    double wl[][ULAT],
                    double cpu[][ULAT],
                    double cpl[][ULAT],
                    double cpp[][ULAT], 
                    double ku[][ULAT],
                    double kl[][ULAT],
                    double f[][ULAT],
                    double cc[][2],   
                    double psinode[][NMAX],
                    double psimin[][NMAX],
                    double jmatrix[][NMAX],
	      	    double jmatrix2[][NMAX],
                    double percent[], 
                    double rssoil[], 
                    double rs[],
                    double dp[],
                    double ff[],
                    int col[],
                    int row[],
                    int indxx[],
                    double subtract[], 
                    double pressure[], 
                    double plc[],
                    double plcweib[], 
                    double rsquare[], 
                    double soilpsi[],
                    double al[],      
                    double ar[],
                    double saturatedKs[][4],
                    double &ptarg,
           	    double &e,
                    double &rr,
                    double &rzw,
                    double &rd, 
                    double &einc,
                    int &SoilRhizElements,
                    double &dt,
                    double &gmd,
                    double &gsd,
                    double &bkd,
                    double &fs,
                    double &fc,  
                    double &Fabs,
                    double &cpplant,  
                    double &saxlat,  
                    double &raxlat,   
                    double &rfract,
                    double &kla,
                    double &aral, 
            	    double &latot,
                    double &ratot,
                    int &shootElements,
                    int &ktotal,
                    double &soilpsimin,
                    int &smodules,
                    int &rootElements,
                    int &rmodules,
                    int &totmodules,
                    double &axShoot_b,
                    double &axShoot_c,
                    double &latShoot_b,
                    double &latShoot_c,
                    double &axRoot_b, 
                    double &axRoot_c,
                    double &latRoot_b,
                    double &latRoot_c,
                    double &plco,
                    double &pdecrement,
                    double &sumy1,
                    double &sumy2,
                    double &sumprod,   
                    double &cincrement,
                    double &soilvol,   
                    double &rlateral,
                    double &rLat_base,
                    double &slateral,
                    double &sLat_base,
                    double &raxial,   
                    double &pleafave, 
                    int &tnode,
                    double fac,  
                    double &modelledKL,
		    int &HydraulicModelFailCond,
		    trees_params treesParams);

    void solvek( double r[][ULAT], double pe[][ULAT], double b[][ULAT], double b1[][ULAT], double n[][ULAT],
                  double n1[][ULAT], double ks[][ULAT], double ws_[][ULAT], double rssoil[], double cf[][ULAT],
                  double vu[][ULAT], double vl[][ULAT], double saturatedKs[][4],
                  double cpp[][ULAT], double wnu[][ULAT], double wu[][ULAT], double wnl[][ULAT],
                  double wl[][ULAT], double ip[][ULAT], double p[][ULAT], double al[], double l[],
                  double jl[][ULAT], double soilpsi[], double cpu[][ULAT], double cpl[][ULAT],
                  double cc[][2], double kl[][ULAT], double ku[][ULAT], double **rflux,
                  double ar[], double f[][ULAT], double ff[NMAX], double jmatrix[][NMAX],
                  double dpp[][ULAT], double dp[NMAX], int col[], int indxx[], double subtract[],
                  double latot, const int rmodules, double rr, double rd,
                  double Fabs, double cpplant, double gmd,
                  double gsd, double fc, double fs, double bkd,
                  int &tnode, double &pleafave,
                  double e, double dt, double fac, const int SoilRhizElements,
                  const int smodules, const int totmodules, double &soilvol, int &HydraulicModelFailCond );

    void initialize( double r[MD][ULAT], double pe[MD][ULAT], double b[MD][ULAT], double b1[MD][ULAT],
                     double n[MD][ULAT], double n1[MD][ULAT], double ks[MD][ULAT], double ws_[MD][ULAT],
                     double rssoil[], double cf[MD][ULAT], double vu[MD][ULAT], double vl[MD][ULAT],
                     double saturatedKs[][4], double cpp[MD][ULAT], double wnu[MD][ULAT], double wu[MD][ULAT],
                     double wnl[MD][ULAT], double wl[MD][ULAT], double ip[MD][ULAT], double p[MD][ULAT],
                     double al[], double l[], int SoilRhizElements, int smodules, int
totmodules,
                     double latot, int rmodules, double rr, double rd, double Fabs,
                     double cpplant, double gsd, double fc, double fs,
                     double bkd, double gmd, double &soilvol );

     bool pgamma( double aa, double nn, double fac, double &x,
                              double &gamma, double &gammap);


     void gauss( int smodules, int SoilRhizElements, int tnode, double jmatrix[][NMAX],
                 int col[], double ff[], int indxx[], double subtract[], double dp[]);

     bool verifyPowCondition(double base, double exponent);

     bool verifyExpCondition( double exponent);

     void gcf( double aa, double fac, double nn, double &gammcf, double &gln, double &x );

     void gammap( double aa, double &gammln, double &x );

     bool gser( double aa, double nn, double &gamser, double &x, double &gln );

     void readin(bool silent, bool callRewet, double soilpsi[], double soilpsiavg[2000], double ks[][ULAT],
                 double pressure[NKINC], double plc[NKINC], double **l_raw,
                 double ksat[MD][MD], double psimin[NMAX][NMAX], double bsat[MD][MD],
                 double ccsat[MD][MD], double newc[MD][MD], double plcweib[NKINC],
                 double rsquare[NKINC], double b[MD][ULAT], double newb[MD][MD],
                 double p[MD][ULAT], double md_raw[], double e_raw[],
                 double la_raw[], double arAl[], double cc[MD][2],
                 const int smodules, const int totmodules, const int rmodules,
                 const int SoilRhizElements, const int rootElements, double &soilpsimin, const int t,
                 double &e, double &latot, double &aral,
                 double &ratot, double &ptarg,
		 trees_params& treesParams);

    void print( double final_pe, double final_ks, double final_b,
                double final_ws, double final_soilvol, double final_rd,
                outputStruct k_p_e, outputStruct ecrit_k_psi,
                double epredOut[], double ecritOut[],
                double ppredOut[], double saturatedKs[][4], double klpred[],
                double md_raw[], double la_raw[], double arAl[], double e_raw[],
                double pc[] );

 public:
        HydraulicModel();



};


//Coupled CO2 and H2O, based on Farquhar and Katul

struct sim_out
        {
	double thetaRoot;
	double waterStress;
        double Gs_est;
        double Gs;
	double WPlant_K;
	double Ecrit;
	double Psi_Crit;
	double Leaf_Psi;
	double Soil_Psi;
	double OutFlag;
        double Ec;
	double rhizFlux[ULAT];
	double rhizFlux1;
	double rhizFlux2;
	double rhizFlux3;
        double Ec_t;
        double Et;
	double snowpack;
	double snowpack_E_deficit;
	double canopy_evaporation;
	double Vcmax25;
	double Vcmax_sun;
	double Vcmax_shd;
	double Jmax25;
	double J_sun;
	double J_shd;
        double A;
        double A_sun;
        double A_shd;
	double Rd_sun;
	double Rd_shd;
        double L_sun;
        double L_shd;
        double T_sun;
        double T_shd;
        double G_sun;
        double G_shd;
        double D0_sun;
        double D0_shd;
        double Ci_sun;
        double Ci_shd;
        double PAR_sun;
        double PAR_shd;
        double gs_sun;
        double gs_shd;
	double NPP;
        double NEE;
        double R_total;
	double R_ag;
        double R_bg;
	double R_leaf;
	double psn;
	double nsc;
	double rmaint;
	double rgrowth;
	double min_temp;
	double pd_lwp; // predawn lwp
	double sum_vpd;
	double GSI[21];
	double GSImean;
	double leaf_growth;
	double leaf_growth_respiration;
        };

const int MAX_PARAMS = 45;
typedef Parameter Param_Array[MAX_PARAMS];


void simulator(Data_Store& in_data,
                 State_Store& state,
                 Parameter* current_params,
                 trees_params treesParams,
                 Data_Frame& Ecframe,
                 Data_Frame& NEEframe,
                 Data_Frame& RchamberFrame,
                 Data_Frame& RstemFrame,
                 Data_Frame& RleafFrame,
                 double& sd_err1,
                 double& sd_err2,
                 double& sd_err3,
                 double& sd_err4,
                 double& sd_err5,
                 double& sd_err1_wt,
                 int n_agg);

void simulator(Data_Store& in_data,
                 State_Store& state,
                 Parameter* current_params,
                 trees_params treesParams,
                 Data_Frame& Ecframe,
                 Data_Frame& NEEframe,
                 Data_Frame& RchamberFrame,
                 Data_Frame& RstemFrame,
                 Data_Frame& RleafFrame,
                 double& sd1,
                 double& sd2,
                 double& sd3,
                 double& sd4,
                 double& sd5);

void simulator(Data_Store& in_data,
                 State_Store& state,
                 Parameter* current_params,
                 trees_params treesParams,
                 Data_Frame& Ecframe,
                 Data_Frame& NEEframe,
                 double& sd_err1,
                 double& sd_err2,
                 double& sd_err1_wt,
                 int n_agg);

void simulator(Data_Store& in_data,
                 State_Store& state,
                 Parameter* current_params,
                 trees_params treesParams,
                 Data_Frame& Ecframe,
                 Data_Frame& NEEframe,
                 std::ofstream& fluxout,
		 std::ofstream& hydrout, std::ofstream& lfareaout);

void simulator(Data_Store& in_data,
                 State_Store& state,
                 Parameter* current_params,
                 trees_params treesParams,
                 Data_Frame& Ecframe,
                 Data_Frame& NEEframe,
                 double& sd1,
                 double& sd2);


struct sim_out simulation_functions(
				bool silent,
				int ts,  //time step
				double Ecrit, //Rate of transpiration above which hydraulic failure occurs
				double &Thresh,     //
				State_Store& state, //current state
				double thetaSoil[],
				double tempSoil[],
				const Time ti,     //current time
                                double u_ref,       //wind speed @ ref height, m/s
                                double t_ref,       //deg C @ ref ht
                                double D_ref,       //vapor pressure deficit @ ref ht, kPa
                                double precip,
                                double Qpar,        //umol m-2 s-1 above canopy, horizontal plane
                                double t_canopy,    //deg C inside canopy
                                double D_canopy,
				double p_air,	    //atm pressure, kPa
                                double co2_atm,     //PPM atmospheric CO2
				double t_surface,	//Soil temperature at air-soil interface, deg. C
				double t_soil,	//Soil temperature, deg. C
				double t_root,	//Rooting depth average temperature, deg. C
				double Zw,		//Water table depth, m
				trees_params &treesParms,
				std::vector<double> EcState,
                 		std::vector<double> EcStateInner,
                 		std::vector<double> &GSI_vector,
				double Ec_,
				int ModelStatus,
				vector<double> LeafPsiInner,
                 		vector<double> KpInner,
                 		vector<int> InnerStep,
				BiogeochemicalCycles& bgc,
				HydraulicModel *hydraulicModel,
				HydraulicModel::outputStruct &ecrit_k_psi, 
				HydraulicModel::outputStruct &k_p_e, double ecritOut[], 
				double pc[], double klpred[],
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
				double ks[][ULAT], double p[][ULAT],
                    		double dpp[][ULAT], double jl[][ULAT], double wnu[][ULAT], 
				double wu[][ULAT], double wnl[][ULAT],
                    		double wl[][ULAT], double cpu[][ULAT], double cpl[][ULAT], 
				double cpp[][ULAT], double ku[][ULAT],
                    		double kl[][ULAT], double f[][ULAT], double cc[][2], 
				double psinode[][NMAX], double psimin[][NMAX],
                    		double jmatrix[][NMAX], double jmatrix2[][NMAX], 
				double percent[], double rssoil[], double rs[],
                    		double dp[], double ff[], int col[], int row[], 
				int indxx[], double subtract[], double pressure[],
                    		double plc[], double plcweib[], double rsquare[], 
				double soilpsi[], double al[], double ar[],
                    		double saturatedKs[][4], double &ptarg, double &e, 
				double &rr, double &rzw, double &rd, double &einc,
                    		int &SoilRhizElements, double &dt, double &gmd, 
				double &gsd, double &bkd, double &fs, double &fc,
                    		double &Fabs, double &cpplant, double &saxlat, 
				double &raxlat, double &rfract, double &kla,
                    		double &aral, double &latot, double &ratot, 
				int &shootElements, int &ktotal, double &soilpsimin,
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
                    		int &tnode, double fac, double &modelledKL, int &HydraulicModelFailCond);

/*
 * conversion from (mmol s-1 m-2) to (mm/sec)
 */
double EcConversion_StoT(double Ec, double LAI);

/*
 * conversion from (mm/sec) to (mmol s-1 m-2)
 */
double EcConversion_TtoS(double Ec, double LAI);

//break in atmos stuff, distribution in canopy, and layer element level stuff

//******* calculate physical properties @ current temp ********
//these functions should use measurements in canopy for calcs.

//atmospheric pressure based on altitude, kPa
double air_pr(double altitude); //m above msl

//molar density of air @ pressure pressure (kPa), temp tcelsius (C)
//in molm-3
double mol_density(double pressure, double tcelsius);

//calculates density of air in kgm-3 using above
double air_density(double pressure, double tcelsius);

//caculates latent heat of vaporization at tcelsius in Jkg-1
double lheat_vap(double tcelsius);

//calculates s (delta/pa)n C-1, slope of saturation mole fraction function
//used in Penman-Monteith, C&N 3.8, 3.9, 3.10
double calc_s(	double temp, 	//air temp, C
		double pa);	//atmospheric pressure, kPa

//calculate saturation vapore pressure in kPa
//Source: Eq. 7-4 in Dingman, 2nd Ed.
double sat_vapor_pressure(double Tair);   //air temp, C

//*************** global radiation and its components ***********
//returns zenith angle(radians)
double zenith(const Time& t,    //current time
              double lat,       //latitude
              double longi);    //longitude

//returns daylength (hours)
double daylength(const Time& t,
                 double lat,
                 double longi);

//calculates leaf area index (LAI) using the GSI algorithm
void compute_GSI_LAI(int ts,
                     int hour,
                     int min,
                     int yday,
                     double daylen,
                     double RL,
		     double D_ref,
                     double t_canopy,
                     double N_neg_fract,
                     double N_neg_demand,
                     double N_post_demand,
                     BiogeochemicalCycles& bgc,
                     State_Store& state,
                     struct sim_out& simOut,
                     trees_params& treesParams,
                     vector<double> &GSI_vector);
//returns GSI, index of leaf phenological potential based on climate indicators
/*
double calc_GSI(double photoPeriod,
                double temperature,
                double vpd);
*/
double calc_GSI(double photoPeriod,
                double temperature,
                double gsi_lwp);

//sky emissivity - Campbell and Norman eq 10.11, 10.12
//Swinbank (1963) & Monteith & Unsworth (1990)
//NOTE: should fcloud remain parameter?, should we use eq 10.10? Brutsaert (1984)?
double sky_emissivity(double tK,        //observed temp @ ref ht, K
                      double fcloud);   //fraction cloud cover - parameter ??

//************* flux partition above canopy *********************
//from Spitters et al., 1986
//returns Qdt - diffuse flux W m-2 (both Qt and Qdt are above canopy values)
//NOTE: Qbt - beam flux at the top of the canopy can be obtained by Qt - Qbt
double diffuse_flux(const Time& t,      //current time
                    double Qt,          //measured radiation flux W m-2 (ground)
                    double zenith);     //zenith angle (radians)

//from Spitters et al., 1986
//partition total flux above canopy in diffuse and beam components (Wm-2)
//Qob - beam above canopy on horizontal surface, returned by ref
//Qod - diffuse above canopy on horizontal surface, returned by ref
void rad_partition(const Time& t,      //current time
                   double Qt,          //Qt, measured radiation flux W m-2 (ground)
                   double zenith,       //zenith angle (radians)        
                   double& Qob,		//ref args
                   double& Qod);

//*******************************************************
//*							*
//*		TWO COMPARTMENT BIG LEAF MODEL		*
//*							*
//*******************************************************
//---- 	Extinction coefficients of canopies 	-----
//Kbe (black leaves) - using eq 15.4 in Campbell & Norman
double cnpy_beam_ext(double zenith,	//radians
                     double l_angle);	//leaf angle distribution, parameter

double cnpy_beam_ext(double zenith,     //radians
                     double l_angle,  //leaf angle distribution, parameter
		     double omega);   //canopy openness, parameter

double cnpy_beam_ext(double zenith,     //radians
                     double l_angle,  //leaf angle distribution, parameter
                     double omega,      //canopy clumpiness, zenith=0, parameter
                     double p_crown);   //canopy depth&diameter, parameter

//diffuse radiation, Kd
//numerical integration 15.5, C&N
double cnpy_diff_ext(double l_angle,
		     double lai_total,
		     double omega,
		     double p_crown);

//returns LAI in shade/sun, Rabs in shade/sun by ref
//alpha differs in PAR & NIR
void absorb_rad(const Time& t,      //current time
		double Qt,          //Qt, measured radiation flux W m-2 (ground)
		double lai_total,   //measured lai
		double zenith,      //zenith angle (radians)
		double Kd,          //need not calc inside
	//parameters
		double l_angle,	    //leaf angle distribution
		double fPAR_beam,   //fraction of PAR in beam and diffuse radiation
		double fPAR_diff,
		double alpha_PAR,   //alpha for PAR & NIR typically 0.8 & 0.2
		double alpha_NIR,
		double omega,	    //parameter for canopy openness
	//ref args
		double& Rabs_sun,
		double& Rabs_shd,
		double& lai_sun,
		double& lai_shd);

//042004
//adding overloaded function to return incident PAR in sun and shd compartments
void absorb_rad(const Time& t,      //current time
                double Qt,          //Qt, measured radiation flux W m-2 (ground)
                double lai_total,   //measured lai
                double zenith,      //zenith angle (radians)
                double Kd,          //need not calc inside
        //parameters
                double l_angle,     //leaf angle distribution
                double fPAR_beam,   //fraction of PAR in beam and diffuse radiation
                double fPAR_diff,
                double alpha_PAR,   //alpha for PAR & NIR typically 0.8 & 0.2
                double alpha_NIR,
		double omega,	    //parameter for canopy openness
		double p_crown,
        //ref args
                double& PAR_sun,        //incident PAR in sun, added 032004 to get gvc
                double& PAR_shd,        //incident PAR in shade
                double& Rabs_sun,
                double& Rabs_shd,
                double& lai_sun,
                double& lai_shd);

//042004
//adding overloaded function to return incident PAR in sun and shd compartments
//adding above and below canopies
void absorb_rad(const Time& t,      //current time
                double Qt,          //Qt, measured radiation flux W m-2 (ground)
                double lai_above,   //measured lai
                double lai_below,   //measured lai
                double zenith,      //zenith angle (radians)
                double Kd,          //need not calc inside
        //parameters
                double l_angle,     //leaf angle distribution
                double fPAR_beam,   //fraction of PAR in beam and diffuse radiation
                double fPAR_diff,
                double alpha_PAR,   //alpha for PAR & NIR typically 0.8 & 0.2
                double alpha_NIR,
                double omega,       //parameter for canopy openness
                double p_crown,
        //ref args
                double& PAR_sun,        //incident PAR in sun, added 032004 to get gvc
                double& PAR_shd,        //incident PAR in shade
                double& Rabs_sun,
                double& Rabs_shd,
                double& lai_sun,
                double& lai_shd);


//*******************************************************
//*                                                     *
//*		CALCULATING CONDUCTANCES		*
//*                                                     *
//*******************************************************
//*******************************************************
//	BOUNDARY LAYER & RADIATIVE CONDUCTANCES
//*******************************************************
//calculate gr 12.7 in C&N
//constants used cp_air, sigma
double calc_gr( double ta,      //air temp, C (ref ht used, may try canopy F
                double es);     //canopy emissivity - param (F - leave for now 6/6)

//calculating gHa to get gHr = gHa + gr, 7.28 in C&N
//uses stability correction obained by successive substitution
//through diabatic correction factors
//as the vpd used in PM is @ ref height, gHa (molm-2s-1) is @ same
//F - can use more correct profile for multiple layers later
//uses von karman const from constants.h
double calc_gHa(double uz, 	//wind speed @ ref height, m/s
		double zee,	//ref height, m
		double canopy_ht,	//height of canopy, m
//the followings are parameters
		double d_factor,	//d = d_factor*canopy_ht
		double zm_factor,	//zm = zm_factor*canopy_ht
		double zh_factor,	//zh = zh_factor*zm
//the followings are calculated values
		double molar_dens,	//molar density of air, molm-3
		double psi_m,		//diabatic correction factors
		double psi_h);

//calculating psi_h and psi_m at a known zeta (stability)
//these are functions for straightforward calculations
//that do not need numerical or iterative solutions
//eq. 7.26 & 7.27 -  Campbell & Norman
double calc_psim(double zeta);
double calc_psih(double zeta);

//returns stability by successive subst 
//refer 7.21, 7.24, 7.26 & 7.27, C&N
//conditions for termination hardcoded in terms of change in u_star value
//in successive substituted value
double stability_sucs(  double tr,      //temp @ ref ht, C
                        double tc,      //temp in canopy, C
                        double z_ref,   //ref ht, m
//parameters
                        double zm_factor,
                        double zh_factor,
                        double d_factor,
//observations
                        double ur,      	//wind speed @ ref ht, ms-1
                        double h_canopy, 	//canopy ht, m
                        double pressure);       //atmospheric pressure, kPa

//*******************************************************
//                      VAPOR CONDUCTANCES
//*******************************************************
//assumption: gva = gHa

//canopy conductance calculation
//Jarvis model, gs = gsmax*(1-delta*D),   gvc = gs*LAI
//returns gvc molm-2s-1

double calc_gvc(double lai,	//observed, single sided
		double D,	//vapor pressure deficit @ ref ht (try in canopy?) F
		double gsmax,   //maximum gs, molm-2s-1 (leaf hemisurface area?)
		double delta);  //kPa-1 as D in kPa

//added 042004
//Jarvis model, gs = gsmax*(1-delta*D)*(PAR/(PAR+A))
double calc_gvc(double lai,     //observed, single sided
                double D,       //vapor pressure deficit @ ref ht (try in canopy?) F
                double PAR,     //constraining with PAR, W/m2
                double gsmax,   //maximum gs, molm-2s-1 (leaf hemisurface area?)
                double delta,   //kPa-1 as D in kPa
                double abse);   //absorbtion efficiency, A, Nishida et al JGR 2003

//Canopy average stomatal conductance from Darcy's Law
//returns in mol m-2 s-1 (Ewers, Mackay, Samanta 2000) Equation 1.
double calc_gvc(double Kl,
                double soilpsi,
                double leafpsi,
                double lai,
                double canopy_ht,
                trees_params &treesParams,
		double E,
                double D,
                double t_canopy,
                double pair);


//**************** GET gHr and gv *********************************
//this one function returns both conductances using above functions
//conductances are molm-2s-1 ground area basis returned in ref args
//a flag is returned
int get_conductances(double& gHr,	//molm-2s-1
		     double& gv,	//molm-2s-1
//observations
                     double h_canopy,    //canopy ht, m
		     double t_ref,	//temp @ ref ht, C
		     double t_canopy,	//temp in canopy, C
		     double z_ref,	//ref ht, m
                     double u_ref,	//wind speed @ ref ht, ms-1
                     double D,		//vapor pressure deficit @ ref ht
		     double lai,	//single sided lai
                     double pressure,	//atmospheric pressure, kPa
//parameters
                     double es,		//canopy emissivity
                     double zm_factor,
                     double zh_factor,
                     double d_factor,
                     double gsmax,   //maximum gs, molm-2s-1
                     double delta,   //kPa-1 as D in kPa
		     double zeta);   //as input, do not want to calculate over

//*************** calculate leaf boundary layer conductance, given wind speed and leaf dimension
//Reference: Magnani et al, 1999, Plant, Cell and Environment
double calc_gvb(
                double zeta,         //Monin-Obukhov stability parameter, already calculated
                double d_factor,     //
                double h_canopy,     //Canopy height
                double zm_factor,
                double z_ref,
                double ur,           // Reference height wind speed
                double alpha,        // Wind speed attenuation coefficient (parameter)
                double lai,
                double d_leaf);      // Characteristic dimension of the leaf (parameter)


//*********SOIL WATER***********
//Water potential
double water_potential(	double porosity,
			double bubbling_pressure,
			double pore_size_index,
			double residual,
			double theta);


//Calculate the proportion of maximum hydraulic conductance of water through the plant
//as a function of soil water content
//return values from 0.001 (hydraulic failure) to 1.0 (no soil water limitation)
double theta_gsref_reduction(   double porosity,
                                double bubbling_pressure,
                                double pore_size_index,
				double residual,
                                double LWP_spring_minimum,
                                double LWP_stomatal_closure,
                                double theta,
				double field_capacity);

//Calculate the proportion of maximum hydraulic conductance of water through the plant
//as a function of soil water content
//return values from 0.001 (hydraulic failure) to 1.0 (no soil water limitation)
//includes a saturated or near-saturated reduction for handling wetlands
double theta_gsref_reduction_wet(   double porosity,
                                double bubbling_pressure,
                                double pore_size_index,
                                double residual,
                                double LWP_spring_minimum,
                                double LWP_stomatal_closure,
                                double theta,
                                double field_capacity);

//Calculate the proportion of potential evaporation from the soil surface
//as a function of soil water content
//returns values from 0.001 (dry soil) to 1.0 (saturated soil)
double theta_soil_evap_reduction(double field_capacity,
                                 double residual,
                                 double theta);

//calculate soil moisture limited conductance rate for evaporation
//DSM - 02/22/2010
double soil_cond_rate(double dZ,
			  double VPD,
			  double Tsoil,
                          double residual,
                          double theta,
                          double field_capacity,
                          double ko,
                          double pore_size_index,
                          double air_entry_pressure,
                          double porosity,
                          double mineral_fraction);

//Calculate Brooks-Corey bubbling pressure (cm), pore-size distribution index, and residual water content
//Source: Rawls, W.J., L.R. Ahuja, and D.L. Brakensiek. 1992. Estimating soil hydraulic properties from
//            soils data. In van Genuchten, M.T. and F.J. Leij (Eds.). Indirection Methods for Estimating
//            the Hydraulic Properties of Unsaturated Soils, Proceedings of the International Workshop,
//            Riverside, California, October 11-13, 1989.

void soil_texture_parameters(   double por,
                                double pClay,
                                double pSand,
				double& ks,
                                double& bubbling_pressure,
                                double& pore_size_index,
                                double& residual);

//Calculate root zone soil moisture
//Added October 24, 2014 DSM
void root_zone_moisture(double Zw,
                        double capFringe,
                        double Droot,
                        int rmodules,
                        double rootDepth[],
                        double thetaSoil[],
                        double porosity,
                        double mineral_fraction);

//Added July 14, 2013 DSM
void root_zone_moisture(double Zw,
                        double capFringe,
                        double Dshallow,
                        double Dmid,
                        double Ddeep,
                        double theta_shallow,
                        double theta_mid,
                        double theta_deep,
                        double porosity,
                        double mineral_fraction,
                        double& theta_root_shallow,
                        double& theta_root_mid,
                        double& theta_root_deep);


//Capillary rise function
//Uses a steady-state solution of Richards equation
//References Gardner (1958), Eagleson (1978), Famiglietti and Wood (1994)
//Warning: becomes undefined as water table depth + suction head approach zero
//         and upper bound has been added to avoid this problem
//Units:
//
double capillary_rise(double water_table_depth,
                      double suction_head,
                      double pore_size_index,
                      double sat_hydraulic_cond);


/* Soil water retention functions using the van Genuchten (1980) equation */
/*       Added 6/16/99 DSM */
double field_capacity(double air_entry_pressure,
                      double pore_size_index,
                      double residual,
                      double porosity);


//
//compute_effective_precipitation()
//
double compute_effective_precipitation(State_Store& state,
				       trees_params treesParams,
                                       double canopy_store_max,
                                       double precip,
                                       double t_ref);

//Infiltration
//
void infiltration(double input,
                  State_Store& state,
                  double* rootDepth,
                  double* thetaSoil,
                  double* bypassFlow,
		  int rmodules,
                  double litter_capacity,
                  double porosity,
                  double pore_size_index,
                  double field_capacity,
                  double residual,
                  double capFringe,
                  double mineral_fraction,
		  double drainScalar,
                  double ks);

// Unsaturated drainage depends on unsaturated hydraulic conductivity
// Use Darcy's Law
double unsat_zone_drainage(double dZ1,
			   double dZ2,
                           double pore_size_index,
                           double ko,
                           double theta1,
                           double theta2,
                           double field_capacity,
                           double residual,
                           double air_entry_pressure,
                           double porosity,
                           double mineral_fraction);

// Allow unsaturated drainage as long as the soil moisture content is above field capacity 
double unsat_zone_drainage(double Zw,
                           double pore_size_index,
                           double ko,
                           double theta,
                           double field_capacity,
                           double residual,
                           double air_entry_pressure,
                           double porosity,
			   double mineral_fraction);


void DoHydraulicModel(bool silent,
		HydraulicModel *hydraulicModel, 
		int ts,
		bool callRewet,
		double Ysoil[],
                double  ev,
                double& Kl,
                double& leafpsi,
                double& Ecrit,
                int &HydraulicModelFailCond,
                double &PsiCrit,
                int hydraulicModel_flag,
                trees_params& treesParams,
		HydraulicModel::outputStruct &ecrit_k_psi, HydraulicModel::outputStruct &k_p_e, 
		double ecritOut[], double pc[], double klpred[],
                double ppredOut[], double epredOut[], double nodeFail[], char nodeTyp[], double ***psi,
                double **rflux, double soilpsiavg[], double nts[], double evap[], double kroot[],
                double axr[], double latr[], double kshoot[], double dslat[], double dsax[], double drlat[],
                double drax[], double l[], double ksat[][MD], double bsat[][MD], double ccsat[][MD],
                double newb[][MD], double newc[][MD], double ip[][ULAT], double b[][ULAT], double b1[][ULAT],
                double n[][ULAT], double n1[][ULAT], double r[][ULAT], double vu[][ULAT], double vl[][ULAT],
                double cf[][ULAT], double pe[][ULAT], double ws_[][ULAT], 
		double ks[][ULAT], double p[][ULAT],
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
                double plc[], double plcweib[], double rsquare[], double soilpsi[], double al[], double ar[],
                double saturatedKs[][4], double &ptarg, double &e, 
		double &rr, double &rzw, double &rd, double &einc,
                int &SoilRhizElements, double &dt, double &gmd, double &gsd, 
		double &bkd, double &fs, double &fc,
                double &Fabs, double &cpplant, double &saxlat, double &raxlat, double &rfract, double &kla,
                double &aral, double &latot, double &ratot, int &shootElements, 
		int &ktotal, double &soilpsimin,
                int &ShootModules, int &rootElements, int &RootModules, int &totmodules, double &axShoot_b,
                double &axShoot_c, double &latShoot_b, double &latShoot_c, 
		double &axRoot_b, double &axRoot_c,
                double &latRoot_b, double &latRoot_c, double &plco, double &pdecrement, double &sumy1,
                double &sumy2, double &sumprod, double &cincrement, double &soilvol, double &rlateral,
                double &rLat_base, double &slateral, double &sLat_base, double &raxial, double &pleafave,
                int &tnode, double fac, double &modelledKL);


//
//compute_canopy_evaporation()
//
double compute_canopy_evaporation(State_Store& state,
                                   double canopy_store_max,
                                   double precip,
                                   double& canopy_wetness,
                                   double Rnet_sun,
                                   double Rnet_shd,
                                   double T_sun,
                                   double T_shd,
                                   double D,
                                   double p_air,
                                   double gHr,
                                   double gva,
				   double canopy_cover);

//*************** calculate ET using PM, given conductance and Rabs values
//uses Penman Monteith to get Ecanopy (mm s-1 m-2 gr area) - terms as 14.12
double do_pm(
                double Rnet,       //absorbed radiation W/m2
                double t_ref,      //deg C @ ref ht (measured)
                double Dee,        //vapor deficit, D (kPa) (measured)
                double pa,         //atmospheric pressure, (kPa)
                double gHr,        //combined boundary layer & radiative cond.
                double gv);        //vapor cond. (both cond. mol m-2 s-1)


//***********************************************************
//function definitions for calculating
//photosynthesis
//refer data_store.h for object definitions for intended use
//***********************************************************


/* structure for input into Farquhar routine  */

struct farqin
	{
        double pa;
	double D;
	double co2;
	double t;
	double irad;
	double g;
	double Rd;
	double lnc;
	double flnr;
	};

/* structure for output from Farquhar routine  */

struct farqout
	{
	 double E;
	 double g;
	 double O2;
	 double Ca;
	 double Ci;
	 double gammaStar;
	 double Kc;
	 double Ko;
	 double act;
	 double Vcmax25;
	 double Vcmax;
	 double Jmax25;
	 double Jmax;
	 double J;
	 double Av;
	 double Aj;
	 double A;
	 double Rd;
        };

double photosynthesis(double Rd_mult, 
			double Jmax_mult, 
			double thetaJ, 
			double psiJ, 
			struct farqin in,
			struct farqout& out, 
			int verbose);


double snow_E_acc(double old_acc,
                  double tavg,
                  double maxdeficit,
                  double snow);

double latentMelt(double snow_E_deficit,
            double Tair,
            double snowpack,
            double vpd,
            double wind_speed);

double sensibleMelt(double snow_E_deficit,
            double Tair,
            double Pair,
            double snowpack,
            double vpd,
            double wind_speed);

double radMelt(double snow_E_deficit,
            double Rad_floor,
            double Tair,
            double melt_Rcoef,
            double snowpack);

double rainMelt(double dewpoint,
                double rainfall,
                double snowpack,
                double &snow_E_deficit);

double sublimate(double snow_E_deficit,
                 double Rad_floor,
                 double Tair,
                 double melt_Rcoef,
                 double snowpack);

double tempMelt(double snow_E_deficit,
                double Tair,
                double melt_Tcoef);



#endif
