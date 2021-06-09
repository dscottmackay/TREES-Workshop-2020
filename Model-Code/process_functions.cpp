//
//process_functions.cpp
//
//This file contains many of the functions called by simulator.cpp and simulation_functions.cpp
//    as part of the TREES model
//
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
//Addition of root phenology and rhizosphere volume dynamics - 2015 DSM
//
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "simulator2.h"


/*
 * Notes: E is per unit leaf area in sperry
 *        E is per unit ground area in trees
 *        E Conversion (1 mol of H20 = 18.2g ->
 *                      1 mm = 1000g ->
 *                      1000g / 18.02g mol-1 = 55.49 mol ->
 *                      55490 mmol /m2 s
 */


//
// conversion from (mmol s-1 m-2) to (mm/sec)
//
double EcConversion_StoT(double Ec, double LAI) 
{
	double Ec_mm;

	Ec_mm = Ec*LAI/55490.0;
	return (Ec_mm);
}

//
// conversion from (mm/sec) to (mmol s-1 m-2)
//
double EcConversion_TtoS(double Ec, double LAI) 
{
	double Ec_molar;

	Ec_molar = Ec/LAI*55490.0;
	return (Ec_molar);
}

//
//-------------- calculate physical properties @ current temp -----------------------
//C&N 3.7
double air_pr(double altitude)
{
	double P;

	P = 101.3*exp(-(altitude/8200.0));
        return (P);
}
 
//
//molar density of air @ pressure pressure (kPa), temp tcelsius (C)
//
double mol_density(double pressure, double tcelsius)
{
	double density;

	density = 44.6*pressure*273.15/(101.3*C2K(tcelsius));
        return (density);
}
 
//mm_air is molecular mass gmol-1 from constants.h
double air_density(double pressure, double tcelsius)
{
	double density;

	density = mol_density(pressure, tcelsius)*mm_air*pow(10.0,-3.0);
        return (density);
}

double lheat_vap(double tcelsius) //Jkg-1
{
	double lv;

	lv = (2.501 - (2.361*pow(10.0,-3.0)*tcelsius))*pow(10.0,6.0);
        return (lv);
}

double calc_s(double temp, double pa) //C-1
{
        double a = 0.611; //kPa pg 41 C&N
        double b = 17.502;
        double c = 240.97;
        double s, esT, delta;

        esT = a*(exp(b*temp/(temp+c)));
        delta = b*c*esT/((c+temp)*(c+temp));
        s = delta/pa;
        return (s);
}

double sat_vapor_pressure(double Tair)
{
	double svp;

	svp = 0.611*exp((17.3*Tair)/(Tair+237.3));
	return (svp);
}


//--------------- global radiation and its components -----------
double zenith(  const Time& t,
                double lat,
                double longi)
{
//calculations are done in radians
        int year, jday, hour, min;
        double declin; //solar declination
        double solnoon;//corrected solar noon
        double corr; //total correction
        double temp, ztemp, zenith;
 
//get the data from the time object
        t.get_time(year, jday, hour, min);

//longitude correction
        corr = long_corr(lat, longi);
//for calculations of corrections due to eq of time
        temp = deg2rad(279.575 + (0.9856*(double)jday));
//add the eqt in hours - eq. 11.4 in C&N
        corr += (-104.7*sin(temp) +596.2*sin(2*temp) +4.3*sin(3*temp) -12.7*sin(4*temp)
                -429.3*cos(temp) -2.0*cos(2.0*temp) +19.3*cos(3*temp))/3600.0;
        solnoon = 12.0 - corr; //eq 11.3 C&N
        lat = deg2rad(lat); //convert to radian
//eq 11.2 - C&N
        temp = deg2rad(278.97 + 0.9856*(double)jday
                        + 1.9165*sin(deg2rad(356.6+0.9856*(double)jday)));
        declin = asin(0.39785*sin(temp));
//eq 11.1 - C&N
        ztemp = sin(lat)*sin(declin)
                + cos(lat)*cos(declin)*
                cos(deg2rad(15.0*((double)hour + (double)min/(double)MINS - solnoon)));
	zenith = acos(ztemp); //the zenith angle
 
        return (zenith);
}


//daylength	Campbell and Norman eq. 11.6, 11.1
double daylength(const Time& t,
		 double lat,
		 double longi)
{
//calculations are done in radians
        int year, jday, hour, min;
        double declin; //solar declination
        double corr; //total correction
        double temp, ztemp, half_day_deg, daylen;

//get the data from the time object
        t.get_time(year, jday, hour, min);

//longitude correction
        corr = long_corr(lat, longi);
        lat = deg2rad(lat); //convert to radian
//eq 11.2 - C&N
        temp = deg2rad(278.97 + 0.9856*(double)jday
                        + 1.9165*sin(deg2rad(356.6+0.9856*(double)jday)));
        declin = asin(0.39785*sin(temp));
//eq 11.1 - C&N
/*
        double solnoon;//corrected solar noon
        ztemp = sin(lat)*sin(declin)
                + cos(lat)*cos(declin)*
                cos(deg2rad(15*((double)hour + (double)min/(double)MINS - solnoon)));
*/
//biologically significant twilight hours assumed, so use 96% zenith
	ztemp = deg2rad(96.0);
//eq 11.6 - C&N - time in radian
	half_day_deg = acos((cos(ztemp) - sin(lat)*sin(declin))/(cos(lat)*cos(declin)));
//eq 11.7 - C&N - time in hours, divide by 15 degrees
	daylen = 2.0*rad2deg(half_day_deg)/15.0;

	return (daylen); //the daylength
}

/*
state, ts, simOut, treesParams, daylen, hour, min, GSI_vector, yday
State_Store &state
int ts
struct sim_out& simOut
trees_params &treesParams
double daylen
int hour
int min
vector<double> &GSI_vector
int yday
double RL
double N_neg_fract
double N_neg_demand
double N_post_demand

*/

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
		     double N_pos_demand,
		     BiogeochemicalCycles& bgc,
		     State_Store& state,
		     struct sim_out& simOut,
		     trees_params& treesParams,
		     vector<double> &GSI_vector)
{
        int i, daystep, dayindex;
        double GSI_sum, GSI_value, lai_scalar, old_lai, delta_lai, delta_nsc, daily_temp, daily_vpd, gsi;
        double daily_pd_lwp, nsc;

        double lai_pot = state.get_val_at(CURRENTTARGETLAI);

	if (ts == 1)
        {
        	simOut.sum_vpd = 0.0;
                simOut.pd_lwp = 0.0; // predawn leaf water potential
                simOut.min_temp = 999.0; //set to an initially high value
                simOut.leaf_growth = 0.0;
                simOut.leaf_growth_respiration = 0.0;
                simOut.GSImean = 0.0;
        }
	if (yday == 364 && (hour+min) == 0.0)
        {
        	lai_pot = (1.0-1.0/treesParams.leafLifeSpan)*state.get_val_at(CURRENTTARGETLAI) +
                                        state.get_val_at(FORECASTTARGETLAI)/treesParams.leafLifeSpan;
                treesParams.lai_at_sat_kl = treesParams.live_lai = lai_pot;
                state.set_val_at(lai_pot, CURRENTTARGETLAI);
                state.set_val_at(0.0, FORECASTTARGETLAI);
        }
        if (ts < 48)
        {
                treesParams.lai = lai_pot * (1.0 - 1.0/treesParams.leafLifeSpan);
        }
        if (treesParams.lai < 0.01)
        {
                treesParams.lai = 0.01;
        }
//added to handled leaf water potential
        if ((hour+min) == 3.0 && treesParams.useHydraulics == true)
        {
                //simOut.pd_lwp = leafpsi;
                simOut.pd_lwp = state.get_val_at(LPSI);
        }

        simOut.sum_vpd += D_ref;
        if (t_canopy < simOut.min_temp)
        {
                simOut.min_temp = t_canopy;
        }
        delta_nsc = 0.0;
	if (hour+(double)min/100.0 == 23.5 && ts >= 48)
        {
                daily_temp = simOut.min_temp;
                daily_pd_lwp = simOut.pd_lwp; // predawn LWP
                simOut.pd_lwp = 0.0;
                daily_vpd = simOut.sum_vpd/48.0;
                simOut.sum_vpd = 0.0;
                simOut.min_temp = 999.0; //set to an initially high value for the next day
                gsi = calc_GSI(daylen*3600.0, daily_temp, daily_pd_lwp);
                daystep = (int) ts/48-1;
                dayindex = (int) daystep % 21 + 1;
                GSI_vector[dayindex] = gsi;
                GSI_sum = 0.0;
                for (i = 1; i <= 21; ++i)
                {
                        GSI_sum += GSI_vector[i];
                }
                if (daystep < 21)
                {
                        GSI_value = GSI_sum/(double)dayindex;
                }
                else
                {
                        GSI_value = GSI_sum / 21.0;
                }
//incipient leaf expansion at GSI = 0.5
//next set of lines ensures a monotonic increase and decrease in LAI during spring and fall, respectively - DSM Oct 2014
                if (yday < 210 && GSI_value < simOut.GSImean)
                {
                        GSI_value = simOut.GSImean;
                }
                else if (yday >= 210 && GSI_value > simOut.GSImean)
                {
                        GSI_value = simOut.GSImean;
                }
                else
                {
                        simOut.GSImean = GSI_value;
                }
                if (GSI_value <= 0.5)
                {
                        lai_scalar = 0.0;
                }
                else if (GSI_value >= 1.0)
                {
                        lai_scalar = 1.0;
                }
                else
                {
                        lai_scalar = (GSI_value-0.5)/0.5;
                }
		nsc = bgc.getLeafNSC();
                delta_nsc = 0.0;
                delta_lai = 0.0;
//when leaf expansion is occurring
                if (lai_scalar > 0.0)
                {
                        lai_scalar = 1.0-1.0/treesParams.leafLifeSpan+lai_scalar/treesParams.leafLifeSpan;
                        old_lai = treesParams.lai;
//first determine the change in LAI without water stress
                        delta_lai = lai_pot * lai_scalar - old_lai;
                }
                        bgc.updateLeafCarbonNitrogenPools(treesParams, delta_lai, RL,
                                                N_neg_fract, N_neg_demand, N_pos_demand);

                        treesParams.lai = bgc.getLeafBiomassCarbon()*treesParams.SLA/10000.0;
                        delta_nsc = nsc - bgc.getLeafNSC();
                        simOut.leaf_growth = 0.86*delta_nsc/48.0;
                        simOut.leaf_growth_respiration = 0.14*delta_nsc /48.0;
	}
}


/*
double calc_GSI(double photoPeriod,
		double temperature,
		double vpd)
*/
double calc_GSI(double photoPeriod,
		double temperature,
		double gsi_lwp)
{
	double gsi, tmin, tmax, temp_indicator;
	double lwpmin, lwpmax, lwp_indicator; //predawn LWP
	double photomin, photomax, photo_indicator;

/*
	tmin = -2.0;
	tmax = 5.0;
	tmin = 2.0;
	tmax = 9.0;
*/
	tmin = -5.0;
	tmax = 10.0;
	if (temperature <= tmin)
	{
		temp_indicator = 0.0;
	}
	else if (temperature >= tmax)
	{
		temp_indicator = 1.0;
	}
	else
	{
		temp_indicator = (temperature-tmin)/(tmax-tmin);
	}

/*
	double vpd_indicator;
        double vpdmin = 3.9;
        double vpdmax = 7.1;
        if (vpd <= vpdmin)
	{
                vpd_indicator = 1.0;
	}
        else if (vpd >= vpdmax)
	{
                vpd_indicator = 0.0;
	}
        else
	{
                vpd_indicator = 1.0-(vpd-vpdmin)/(vpdmax-vpdmin);
	}
*/

	lwpmin = -6.0;
	lwpmax = -1.0;
	if (gsi_lwp <= lwpmin)
	{
		lwp_indicator = 0.0;
	}
	else if (gsi_lwp >= lwpmax)
	{
		lwp_indicator = 1.0;
	}
	else
	{
		lwp_indicator = (gsi_lwp-lwpmin)/(lwpmax-lwpmin);
	}

        photomin = 36000.0;
        photomax = 39600.0;
        if (photoPeriod <= photomin)
	{
                photo_indicator = 0.0;
	}
        else if (photoPeriod >= photomax)
	{
                photo_indicator = 1.0;
	}
        else
	{
                photo_indicator = (photoPeriod-photomin)/(photomax-photomin);
	}

	//gsi = temp_indicator*vpd_indicator*photo_indicator;
	gsi = temp_indicator*lwp_indicator*photo_indicator;

	return (gsi);
}
 
//sky emissivity    Campbell and Norman eq 10.11, 10.12
//depricated function replaced with diffuse_flux() function
double sky_emissivity(double tK,        //ref temp in K - observed
                      double fcloud)    //fraction cloud cover
{
	double emm;

	emm = (1.0 - 0.84*fcloud)*(0.0000092*tK*tK) + (0.84 * fcloud);
        return (emm);
}

//------------- flux partition above canopy ------------------------------------
double diffuse_flux(const Time& t,
                    double Qt,
                    double zenith)
{
        int year, jday, hour, min;
        double elev; //elevation angle in radians
        double R, K, Qo, atm_trns, fdiffuse, Qdiffuse; //variables used here
        t.get_time(year, jday, hour, min); //get values from t
        elev = M_PI_2 - zenith; //in radians
 
//Eq.1 - Spitters et al 1986
        Qo = Ssc * sin(elev) * (1.0 + 0.033*cos(deg2rad(360.0*(double)jday/365.0)));
 
//the follwing is a set of eq. 20a, b, c, d from Spitters et al 1986 appendix
        atm_trns = Qt/Qo;
        R = 0.847 - (1.61*sin(elev)) + (1.04 * sin(elev)*sin(elev));
        K = (1.47 - R)/1.66;
 
//use the conditions to get fdiffuse
        if (atm_trns <= 0.22) 
	{
		fdiffuse = 1.0;
	}
        if (atm_trns > 0.22 && atm_trns <= 0.35) 
	{
		fdiffuse = 1.0-(6.4*(atm_trns-0.22)*(atm_trns-0.22));
	}
        if (atm_trns > 0.35 && atm_trns <= K) 
	{
		fdiffuse = 1.47 - (1.66*atm_trns);
	}
        if (atm_trns > K) 
	{
		fdiffuse = R;
	}
	Qdiffuse = fdiffuse*Qt;

        return (Qdiffuse);
}


void rad_partition(const Time& t,      //current time
                   double Qt,          //Qt, measured radiation flux W m-2 (ground)
                   double zenith,       //zenith angle (radians)        
                   double& Qob,         //ref args
                   double& Qod)
{
	Qod = diffuse_flux(t, Qt, zenith);
	Qob = Qt - Qod;
}

//*******************************************************
//*                                                     *
//*             TWO COMPARTMENT BIG LEAF MODEL          *
//*                                                     *
//*******************************************************

//----  Extinction coefficients of canopies     -----
//Kbe (black leaves) - using eq 15.4 in Campbell & Norman
double cnpy_beam_ext(double zenith,
                     double l_angle)
{
        double numerator, denominator, Kbe;

        numerator = sqrt((l_angle*l_angle)+tan(zenith)*tan(zenith));
        denominator = l_angle + 1.774*pow((l_angle+1.182),-0.733);
	Kbe = numerator/denominator;

        return (Kbe);
}

//Kbe (black leaves) with adjustment for effective canopy openness
double cnpy_beam_ext(double zenith,     //radians
                     double l_angle,  //leaf angle distribution, parameter
                     double omega)   //canopy openness, parameter
{
 	double numerator, denominator, Kbe;

	numerator = sqrt((l_angle*l_angle)+tan(zenith)*tan(zenith));
	denominator = l_angle + 1.774*pow((l_angle+1.182),-0.733);
	Kbe = numerator/denominator*omega;

	return (Kbe);
}

//Kbe (black leaves) with adjustment for effective canopy openness
//adjust omega by zenith angle - DSM 8/22/2008
double cnpy_beam_ext(double zenith,     //radians
                     double l_angle,  //leaf angle distribution, parameter
                     double omega,	//canopy clumpiness, zenith=0, parameter
		     double p_crown)   //canopy depth&diameter, parameter
{
        double numerator, denominator, omega_zen, Kbe;

        numerator = sqrt((l_angle*l_angle)+tan(zenith)*tan(zenith));
        denominator = l_angle + 1.774*pow((l_angle+1.182),-0.733);
	omega_zen = omega/(omega+(1.0-omega)*exp(-2.2*pow(zenith,p_crown))); //Eq. 15.35 C&N 1998
	Kbe = (numerator/denominator)*omega_zen;

        return (Kbe);
}

 
//diffuse radiation, Kd
//numerical integration 15.5, C&N
double cnpy_diff_ext(double l_angle,
                     double lai_total,
		     double omega,
		     double p_crown)
{
        double Kd, tau_d;
//variables for numerical integration
        double integral; //value of sum
        double steps;
        double psi, dpsi, max_psi;

//initialize
        steps = 90.0; //about 1 degree steps
        max_psi = M_PI_2; //pi/2
        dpsi = max_psi/steps;
        psi = 0.0;
        integral = 0.0;
        while(psi < max_psi) //left sum
        {
                integral += exp(-cnpy_beam_ext(psi, l_angle, omega, p_crown)*lai_total)
							*sin(psi)*cos(psi)*dpsi;
                psi += dpsi;
        }
 
        tau_d = 2.0*integral;
        Kd = -log(tau_d)/lai_total;

        return (Kd);
}
 

//absorbed radiation in sun and shade parts of canopy, big leaf, alpha diff in PAR & NIR
//returns LAI in shade/sun, Rabs in shade/sun by ref
//ch 15.8 C&N
void absorb_rad(const Time& t,      //current time
                double Qt,          //Qt, measured radiation flux W m-2 (ground)
                double lai_total,   //measured lai
                double zenith,      //zenith angle (radians)
		double Kd,	    //need not calc inside
	//parameters
                double l_angle,     //leaf angle distribution
                double fPAR_beam,   //fraction of PAR in beam and diffuse radiation
                double fPAR_diff,
                double alpha_PAR,   //alpha for PAR & NIR typically 0.8 & 0.2
                double alpha_NIR,
		double omega,       //parameter for canopy openness
        //ref args
                double& Rabs_sun,
                double& Rabs_shd,
                double& lai_sun,
                double& lai_shd)
{
//for above canopy
	double Qod, Qob, Qod_NIR, Qod_PAR, Qob_NIR, Qob_PAR;
//in canopy averages: diffuse, unintercepted beam, total beam, scattered radiation
	double Qd_NIR, Qd_PAR, Qb_NIR, Qb_PAR, Qbt_NIR, Qbt_PAR, Qsc_NIR, Qsc_PAR;
//average absorbed flux density
	double Qsun, Qshd;
//extinction coefficient
	double Kb;

	Kb = cnpy_beam_ext(zenith, l_angle, omega);

//break down in diffuse and beam components
	rad_partition(t, Qt, zenith, Qob, Qod);
	Qob_PAR = fPAR_beam*Qob;
	Qob_NIR = Qob - Qob_PAR;
	Qod_PAR = fPAR_diff*Qod;
	Qod_NIR = Qod - Qod_PAR;

//calculate averages
	//for PAR
	Qd_PAR = Qod_PAR*(1.0-exp(-sqrt(alpha_PAR)*Kd*lai_total))/(sqrt(alpha_PAR)*Kd*lai_total);
	Qbt_PAR = Qob_PAR*exp(-sqrt(alpha_PAR)*Kb*lai_total);
	Qb_PAR = Qob_PAR*exp(-Kb*lai_total);
	Qsc_PAR = 0.5*(Qbt_PAR - Qb_PAR);    //average (scattered 0 @ top)

        //for NIR
        Qd_NIR = Qod_NIR*(1.0-exp(-sqrt(alpha_NIR)*Kd*lai_total))/(sqrt(alpha_NIR)*Kd*lai_total);
        Qbt_NIR = Qob_NIR*exp(-sqrt(alpha_NIR)*Kb*lai_total);
        Qb_NIR = Qob_NIR*exp(-Kb*lai_total);
        Qsc_NIR = 0.5*(Qbt_NIR - Qb_NIR);    //average (scattered 0 @ top)

//absorbed flux density
	Qshd = alpha_PAR*(Qd_PAR+Qsc_PAR) + alpha_NIR*(Qd_NIR+Qsc_NIR);
	Qsun = alpha_PAR*(Kb*Qob_PAR+Qd_PAR+Qsc_PAR) + alpha_NIR*(Kb*Qob_NIR+Qd_NIR+Qsc_NIR);

//lai partition
	lai_sun = (1.0-exp(-Kb*lai_total))/Kb;
	lai_shd = lai_total-lai_sun;

//absorbed flux
	Rabs_sun = Qsun*lai_sun;
	Rabs_shd = Qshd*lai_shd;
}

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
		double omega,       //parameter for canopy openness
		double p_crown,
        //ref args
                double& PAR_sun,        //incident PAR in sun, added 032004 to get gvc
                double& PAR_shd,        //incident PAR in shade
                double& Rabs_sun,
                double& Rabs_shd,
                double& lai_sun,
                double& lai_shd)
{
//for above canopy
        double Qod, Qob, Qod_NIR, Qod_PAR, Qob_NIR, Qob_PAR;
//in canopy averages: diffuse, unintercepted beam, total beam, scattered radiation
        double Qd_NIR, Qd_PAR, Qb_NIR, Qb_PAR, Qbt_NIR, Qbt_PAR, Qsc_NIR, Qsc_PAR;
//average absorbed flux density
        double Qsun, Qshd;
//extinction coefficient
        double Kb;

        Kb = cnpy_beam_ext(zenith, l_angle, omega, p_crown);

//break down in diffuse and beam components
        rad_partition(t, Qt, zenith, Qob, Qod);
        Qob_PAR = fPAR_beam*Qob;
        Qob_NIR = Qob - Qob_PAR;
        Qod_PAR = fPAR_diff*Qod;
        Qod_NIR = Qod - Qod_PAR;
 
//calculate averages
        //for PAR
        Qd_PAR = Qod_PAR*(1.0-exp(-sqrt(alpha_PAR)*Kd*lai_total))/(sqrt(alpha_PAR)*Kd*lai_total);
        Qbt_PAR = Qob_PAR*exp(-sqrt(alpha_PAR)*Kb*lai_total);
        Qb_PAR = Qob_PAR*exp(-Kb*lai_total);
        Qsc_PAR = 0.5*(Qbt_PAR - Qb_PAR);    //average (scattered 0 @ top)
 
        //for NIR
        Qd_NIR = Qod_NIR*(1.0-exp(-sqrt(alpha_NIR)*Kd*lai_total))/(sqrt(alpha_NIR)*Kd*lai_total);
        Qbt_NIR = Qob_NIR*exp(-sqrt(alpha_NIR)*Kb*lai_total);
        Qb_NIR = Qob_NIR*exp(-Kb*lai_total);
        Qsc_NIR = 0.5*(Qbt_NIR - Qb_NIR);    //average (scattered 0 @ top)
 
//absorbed flux density
        Qshd = alpha_PAR*(Qd_PAR+Qsc_PAR) + alpha_NIR*(Qd_NIR+Qsc_NIR);
        Qsun = alpha_PAR*(Kb*Qob_PAR+Qd_PAR+Qsc_PAR) + alpha_NIR*(Kb*Qob_NIR+Qd_NIR+Qsc_NIR);
 
//incident flux (this is the mod 042004)
        PAR_sun = Kb*Qob_PAR+Qd_PAR+Qsc_PAR; //average
        PAR_shd = Qd_PAR+Qsc_PAR;
	if (PAR_sun > 2200.0)
	{
		PAR_sun = 2200.0;
	}

//lai partition
        lai_sun = (1.0-exp(-Kb*lai_total))/Kb;
        lai_shd = lai_total-lai_sun;

//absorbed flux
        Rabs_sun = Qsun*lai_sun;
        Rabs_shd = Qshd*lai_shd;
}
//042004
//adding overloaded function to return incident PAR in sun and shd compartments
void absorb_rad(const Time& t,      //current time
                double Qt,          //Qt, measured radiation flux W m-2 (ground)
                double lai_above,   //measured lai emergent canopy
                double lai_below,   //measured lai supressed canopy
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
                double& lai_shd)
{
//for above canopy
        double Qod, Qob, Qod_NIR, Qod_PAR, Qob_NIR, Qob_PAR;
//in canopy averages: diffuse, unintercepted beam, total beam, scattered radiation
        double Qd_NIR, Qd_PAR, Qb_NIR, Qb_PAR, Qbt_NIR, Qbt_PAR, Qsc_NIR, Qsc_PAR;
        double tQd_NIR, tQd_PAR, tQb_NIR, tQb_PAR, tQbt_NIR, tQbt_PAR, tQsc_NIR, tQsc_PAR;
        double dQd_NIR, dQd_PAR, dQb_NIR, dQb_PAR, dQbt_NIR, dQbt_PAR, dQsc_NIR, dQsc_PAR;
//average absorbed flux density
        double Qsun, Qshd, lai_total;
//extinction coefficient
        double Kb;

        Kb = cnpy_beam_ext(zenith, l_angle, omega, p_crown);

//break down in diffuse and beam components
        rad_partition(t, Qt, zenith, Qob, Qod);
        Qob_PAR = fPAR_beam*Qob;
        Qob_NIR = Qob - Qob_PAR;
        Qod_PAR = fPAR_diff*Qod;
        Qod_NIR = Qod - Qod_PAR;

	lai_total = lai_above + lai_below;

//calculate averages for overstory
        //for PAR
        Qd_PAR = Qod_PAR*(1.0-exp(-sqrt(alpha_PAR)*Kd*lai_above))/(sqrt(alpha_PAR)*Kd*lai_above);
        Qbt_PAR = Qob_PAR*exp(-sqrt(alpha_PAR)*Kb*lai_above);
        Qb_PAR = Qob_PAR*exp(-Kb*lai_above);
        Qsc_PAR = 0.5*(Qbt_PAR - Qb_PAR);    //average (scattered 0 @ top)
 
        //for NIR
        Qd_NIR = Qod_NIR*(1.0-exp(-sqrt(alpha_NIR)*Kd*lai_above))/(sqrt(alpha_NIR)*Kd*lai_above);
        Qbt_NIR = Qob_NIR*exp(-sqrt(alpha_NIR)*Kb*lai_above);
        Qb_NIR = Qob_NIR*exp(-Kb*lai_above);
        Qsc_NIR = 0.5*(Qbt_NIR - Qb_NIR);    //average (scattered 0 @ top)
 
//calculate remaining diffuse and direct components not absorbed by overstory
	tQd_PAR = Qod_PAR*(1.0-exp(-sqrt(alpha_PAR)*Kd*lai_total))/(sqrt(alpha_PAR)*Kd*lai_total);
	tQbt_PAR = Qob_PAR*exp(-sqrt(alpha_PAR)*Kb*lai_total);
	tQb_PAR = Qob_PAR*exp(-Kb*lai_total);
        tQsc_PAR = 0.5*(Qbt_PAR - Qb_PAR);    //average (scattered 0 @ top)

	tQd_NIR = Qod_NIR*(1.0-exp(-sqrt(alpha_NIR)*Kd*lai_total))/(sqrt(alpha_NIR)*Kd*lai_total);
	tQbt_NIR = Qob_NIR*exp(-sqrt(alpha_NIR)*Kb*lai_total);
	tQb_NIR = Qob_NIR*exp(-Kb*lai_total);
        tQsc_NIR = 0.5*(Qbt_NIR - Qb_NIR);    //average (scattered 0 @ top)

	dQd_PAR = tQd_PAR - Qd_PAR;
	dQbt_PAR = tQbt_PAR - Qbt_PAR;

	dQb_PAR = Qb_PAR;

	dQsc_PAR = tQsc_PAR - Qsc_PAR;

	dQd_NIR = tQd_NIR - Qd_NIR;
	dQbt_NIR = tQbt_NIR - Qbt_NIR;

	dQb_NIR = Qb_NIR;

	dQsc_NIR = tQsc_NIR - Qsc_NIR;


//absorbed flux density
        Qshd = alpha_PAR*(tQd_PAR+tQsc_PAR) + alpha_NIR*(tQd_NIR+tQsc_NIR);
        Qsun = alpha_PAR*(Kb*dQb_PAR+tQd_PAR+tQsc_PAR) + alpha_NIR*(Kb*dQb_NIR+tQd_NIR+tQsc_NIR);
 
//incident flux (this is the mod 042004)
        PAR_sun = Kb*dQb_PAR+tQd_PAR+tQsc_PAR; //average
        PAR_shd = tQd_PAR+tQsc_PAR;
	if (PAR_sun > 2200.0)
	{
		PAR_sun = 2200.0;
	}

//lai partition
        lai_sun = (1.0-exp(-Kb*(lai_total)))/Kb - (1.0-exp(-Kb*(lai_above)))/Kb;
	if (lai_sun < 0.01)
	{
		lai_sun = 0.01;
	}
        lai_shd = lai_below-lai_sun;

//absorbed flux
        Rabs_sun = Qsun*lai_sun;
        Rabs_shd = Qshd*lai_shd;
}

//*******************************************************
//*                                                     *
//*             CALCULATING CONDUCTANCES                *
//*                                                     *
//*******************************************************
//calculate gr 12.7 in C&N
//constants used cp_air, sigma
//changed cp_air to molm-2s-1 in constants.h 042004
double calc_gr( double ta,      //air temp, C
                double es)      //canopy emissivity - param
{
        double gr;

	gr = 4.0*es*sigma*pow(C2K(ta), 3.0)/cp_air; //corrected unit to molm-2s-1

        return (gr);
}

//calculating gHa to get gHr = gHa + gr
//using 7.28 in C&N
double calc_gHa(double uz,      //wind speed @ ref height, m/s
                double zee,     //ref height, m
                double canopy_ht,       //height of canopy, m
        //the followings are parameters
                double d_factor,        //d = d_factor*canopy_ht
                double zm_factor,       //zm = zm_factor*canopy_ht
                double zh_factor,       //zh = zh_factor*zm
        //the followings are calculated values
                double molar_dens,      //molar density of air, molm-3
                double psi_m,           //diabatic correction factors
                double psi_h)
{
	double gHa, lnz_d;

	lnz_d = log(zee - d_factor*canopy_ht);
	gHa = (vkk*vkk*molar_dens*uz)/((lnz_d - log(zm_factor*canopy_ht) + psi_m)*
		(lnz_d - log(zh_factor*zm_factor*canopy_ht) + psi_h)); 
//enable for check
/*
cout << vkk << '\t' << molar_dens << '\t' << uz  << '\t' 
<< (lnz_d - log(zm_factor*canopy_ht) + psi_m) << '\t' 
<< (lnz_d - log(zh_factor*zm_factor*canopy_ht) + psi_h) << '\t' 
<< (vkk*vkk*molar_dens*uz) << '\t'  << gHa << '\t' 
<< (log((zee - d_factor*canopy_ht)/(zm_factor*canopy_ht))+psi_m) << '\t' 
<< (log((zee - d_factor*canopy_ht)/(zh_factor*zm_factor*canopy_ht))+psi_h) << endl;
*/
	return (gHa);
}

//calculating psi_h and psi_m at a known zeta (stability)
//eq. 7.26 & 7.27 -  Campbell & Norman
double calc_psim(double zeta)
{
        double psi_m;
 
        if(zeta < 0.0)
	{
                psi_m = -1.2*log((1.0 + sqrt(1.0 - 16.0*zeta))/2.0);
	}
        else
	{
                psi_m = 6.0*log1p(zeta);
	}
 
        return (psi_m);
}
 
double calc_psih(double zeta)
{
	double psi_h;

        if(zeta < 0.0)
	{
                psi_h = calc_psim(zeta)/0.6;
	}
        else
	{
                psi_h = calc_psim(zeta);
	}

	return (psi_h);
}

//trying to find by successive subst
//refer 7.21, 7.24, 7.26 & 7.27, C&N
double stability_sucs(  double tr,      //temp @ ref ht, C
                        double tc,      //temp in canopy, C
                        double z_ref,   //ref ht, m
                //parameters
                        double zm_factor,
                        double zh_factor,
                        double d_factor,
                //observations
                        double ur,              //wind speed @ ref ht, ms-1
                        double h_canopy,        //canopy ht, m
                        double pressure)       //atmospheric pressure, kPa
{
//true constants used are vkk and cp_air & gr_acc
//the fixed terms
        double dee = d_factor*h_canopy;
        double zee_m = zm_factor*h_canopy;
        double zee_h = zh_factor*zee_m;
        double rho = mol_density(pressure, tr); //molm-3
        double lnm = log((z_ref - dee)/zee_m);
        double lnh = log((z_ref - dee)/zee_h);
//temps in K
        double trK = C2K(tr);
        double tcK = C2K(tc);
 
//the vars
        double u_star, h_flux, zeta, psi_m, psi_h;
        double u_star_last, u_star_diff;        //u* is what we need, so check this

//********** check effects F *************
        double u_star_accuracy = 0.001; //accuracy of ur is <= 0.01 m/s 042004 - a change @ 4th pl dec
        int i;
        int n_times = 50;  //give up after this
	int min_n = 10; //do at least 10 iterations 042004 - a change @ 4th pl dec
 
//estimates using naive values
        zeta = psi_m = psi_h = 0.0;
        u_star = u_star_last = ur*vkk/(lnm + psi_m);
        h_flux = (tcK - trK)*vkk*rho*cp_air*u_star/(lnh + psi_h); //left as is 042004 mistake before?
/*
cout << "stability : " << zeta
<< "psi : " << psi_m << " " << psi_h
<< " U*: " << u_star
<< " H: " << h_flux << endl;
*/
//iterate
        u_star_diff = 1.0; //force loop
        i = 0;
        while((u_star_diff > u_star_accuracy || i <= min_n) && i < n_times)
        {
                zeta = -vkk*gr_acc*z_ref*h_flux/(rho*cp_air*trK*u_star*u_star*u_star); //7.21 
									//left as is 042004 mistake before?
									//obtains same values
                psi_m = calc_psim(zeta); //7.26 or 7.27
                psi_h = calc_psih(zeta); //7.26 or 7.27
                u_star = ur*vkk/(lnm + psi_m); //7.24
                h_flux = (tcK - trK)*vkk*rho*cp_air*u_star/(lnh + psi_h); //unnumbered in C&N
									  //left as is 042004 mistake before?
                u_star_diff = u_star_last - u_star;
                u_star_diff = absdbl(u_star_diff);
                u_star_last = u_star;
//check
/*
cout << " " << i ;
cout << "\nstability : " << zeta
<< "psi : " << psi_m << " " << psi_h
<< " U*: " << u_star
<< " H: " << h_flux << endl;
*/
                i++;
        }//end for
/*
if(i >= n_times)
        cout << "WARNING! " << n_times << " iterations exceeded.\n";
*/
        return (zeta);
}

//canopy conductance calculation
//returns gvc molm-2s-1
double calc_gvc(double lai,     //observed, single sided
                double D,       //vapor pressure deficit @ ref ht (try in canopy?) F
                double gsmax,   //maximum gs, molm-2s-1 (leaf hemisurface area?)
                double delta)	//kPa-1 as D in kPa
{
	double gvc;

	if (D > 0.1)
	{
		gvc = (gsmax-delta*log(D))*lai;
	}
	else
	{
		gvc = (gsmax-delta*log(0.1))*lai;
	}

	return (gvc);
}
//042004
double calc_gvc(double lai,     //observed, single sided
                double D,       //vapor pressure deficit @ ref ht (try in canopy?) F
                double PAR,     //constraining with PAR, W/m2
                double gsmax,   //maximum gs, molm-2s-1 (leaf hemisurface area?)
                double delta,   //kPa-1 as D in kPa
                double abse)    //absorbtion efficiency, Nishida et al JGR 2003
{
        double gvc;

        gvc = gsmax*(1.0-delta*D)*(PAR/(PAR+(abse*0.235)))*lai; //converting mumolm-2s-1 to Wm-2
	if (gvc < 0.00000001) 
	{
		gvc = 0.00000001;
	}

        return (gvc);
}

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
                double pair)
{
       	double gvc;
       	double gs;
       	double D_molar; //converted to molar fraction mol/mol
	double KG; //conductance coefficient kPa m3 kg-1
//change lai to leaf area. total leaf surface area one sided.

	if (soilpsi > 0.5*(treesParams.pd_at_sat_kl+soilpsi))
        {
                soilpsi = 0.5*(treesParams.pd_at_sat_kl+soilpsi);
        }

//comment out if you don't want to include canopy height
	//canopy_ht = 0.0;

	if (soilpsi < 0.0 && Kl > 0.0 && leafpsi < soilpsi && !isnan(leafpsi))
	{
		D_molar = D/pair;
//Assume mean height to pull water against gravity is at 2/3 canopy height
		canopy_ht=0.67*canopy_ht;
		gs = (Kl/D_molar) * (soilpsi - leafpsi-(canopy_ht*(9.81)*0.001));
	}
	else
	{
		KG = 115.8 + 0.4236*t_canopy;
		gs = KG*(E/lai)/D; //m s-1 per unit leaf area basis
		gs *= 42000.0; //to mmols
	}

//must return gvc in mols s-1 instead of mmols m-2 s-1 as calculated above.
        gvc = gs * 0.001; //to mols
	if (gvc < 0.000001)
	{
		gvc = 0.000001;
	}

        return (gvc);
}


//this one function computes both conductances using above functions
//conductances are molm-2s-1 ground area basis returned in ref args
//a flag is returned
//NOT USED
int get_conductances(double& gHr,       //molm-2s-1
                     double& gv,        //molm-2s-1
                //observations
                     double h_canopy,    //canopy ht, m
                     double t_ref,      //temp @ ref ht, C
                     double t_canopy,   //temp in canopy, C
                     double z_ref,      //ref ht, m
                     double u_ref,      //wind speed @ ref ht, ms-1
                     double D,          //vapor pressure deficit @ ref ht
                     double lai,        //single sided lai
                     double pressure,   //atmospheric pressure, kPa
                //parameters
                     double es,         //canopy emissivity
                     double zm_factor,
                     double zh_factor,
                     double d_factor,
                     double gsmax,   //maximum gs, molm-2s-1
                     double delta,   //kPa-1 as D in kPa
		     double zeta)    //as input, do not want to calculate over
{
	int flag = 1;
//intermediate variables
	double psi_m, psi_h, rho_mol;
//conductance components
	double gr, gHa, gva, gvc;

//calculate required values
	rho_mol = mol_density(pressure, t_ref);
        psi_m = calc_psim(zeta);
        psi_h = calc_psih(zeta);

//calculate the components
	gr = calc_gr(t_ref, es);
	gHa = calc_gHa(u_ref, z_ref, h_canopy, d_factor, zm_factor, zh_factor, rho_mol, psi_m, psi_h);
	gva = gHa; //14.9 C&N
	gvc = calc_gvc(lai, D, gsmax, delta);

//calculate gHr
	gHr = gHa + gr; //14.5 C&N
	gv = 1.0/(1.0/gva + 1.0/gvc); //pg 229 C&N

	return (flag);
}


//calculate leaf boundary layer conductance
double calc_gvb(
		double zeta, 	     //Monin-Obukhov stability parameter, already calculated
		double d_factor,     //
		double h_canopy,     //Canopy height
		double zm_factor,
		double z_ref,
		double ur,           // Reference height wind speed
		double alpha,        // Wind speed attenuation coefficient (parameter)
		double lai,
		double d_leaf)       // Characteristic dimension of the leaf (parameter)
{
	double dee, zee_m, lnm, lnm_canopy, psi_m, u_star, u_H, u_avg, gvb;

        dee = d_factor*h_canopy;
        zee_m = zm_factor*h_canopy;
        lnm = log((z_ref - dee)/zee_m);
        lnm_canopy = log((h_canopy - dee)/zee_m);
        psi_m = calc_psim(zeta); //7.26 or 7.27
        u_star = ur*vkk/(lnm + psi_m); //7.24
	u_H = u_star/vkk*(lnm_canopy - psi_m);	     // wind speed at top of canopy
	u_avg = u_H*(1.0-exp(0.5*alpha)/alpha);  // average wind speed attenduated through the canopy
	if (u_avg < 0.0001)
	{
		u_avg = 0.0001;
	}
	gvb = lai*diff_h2o*sqrt(u_avg/d_leaf);

	return(gvb);
}

//Calculate water potential, MPa
//Assues bubbling pressure in cm
double water_potential( double porosity,
                        double bubbling_pressure,
                        double pore_size_index,
                        double residual,
                        double theta)
{
	double psi_soil, n, m, S;

	S = (theta - residual) / (porosity - residual);
        if (S < 0.001)
	{
                S = 0.001;
	}
        else if (S > 1.0)
	{
                S = 1.0;
	}
        n = pore_size_index + 1.0;
        m = pore_size_index / n;

//Use van Ganuchten model of soil water potential
        psi_soil = -0.0001019977334*bubbling_pressure * pow(pow(S,-1.0/m)-1.0, 1.0/n);

	if (psi_soil < -10.0)
	{
		psi_soil = -10.0;
	}

	return(psi_soil);
}


//Calculate the proportion of maximum hydraulic conductance of water through the plant
//as a function of soil water content
//return values from 0.001 (hydraulic failure) to 1.0 (no soil water limitation)
double theta_gsref_reduction(   double por,
                                double bubbling_pressure,
                                double pore_size_index,
				double residual,
                                double LWP_spring_minimum,
                                double LWP_stomatal_closure,
                                double theta,
				double field_capacity)
{

        double psi_soil, n, m, S, LWP, water_stress;
	double Sfc, PSI_crit;

        S = (theta - residual) / (por - residual);
        if (S < 0.0000001)
	{
                S = 0.0000001;
	}
        else if (S > 1.0)
	{
                S = 1.0;
	}
        n = pore_size_index + 1.0;
        m = pore_size_index / n;

//Use van Ganuchten model of soil water potential
        psi_soil = 0.0001*bubbling_pressure * pow(pow(S,-1.0/m)-1.0, 1.0/n);

//Convert water potential to MPa (positive values refer to negative water potentials)
        LWP = psi_soil;

//Determine water potential equivalent of field capacity
        Sfc = (field_capacity - residual) / (por - residual);
        if (Sfc < 0.0000001)
	{
                Sfc = 0.0000001;
	}
        else if (Sfc > 1.0)
	{
                Sfc = 1.0;
	}
	PSI_crit = 0.0001*bubbling_pressure;

//Water stress for both low and high soil moisture (experimental)
	if (LWP > PSI_crit)
	{
		LWP = max(LWP, LWP_spring_minimum);
        	water_stress = min(1.0, (LWP_stomatal_closure - LWP)/(LWP_stomatal_closure - LWP_spring_minimum));
	}
	else
	{
//Based on Chen et al 2005, Journal of Hydrology
//		water_stress = max(0.7, 1.0-0.3*(PSI_crit-LWP)/PSI_crit);
		water_stress = 1.0; //turn off this reduction effect
	}

        if (water_stress > 0.0000001)
	{
                return(water_stress);
	}
        else
	{
                return(0.0000001);
	}
}

//Calculate Brooks-Corey bubbling pressure (cm), pore-size distribution index, and residual water content
//Sources: 
//	Rawls, W.J. and D.L. Brakensiek. 1985. Prediction of soil water properties for hydrologic modeling. 
//		In Jones, E.B. and T.J. Ward (Eds.). Watershed Management in the Eighties, Proceedings of the 
//		Symposium sponsored by the Committee on Watershed Management of the Irrigatio and Drainage 
//		Division of the American Society of Civil Engineers in conjunction with the ASCE Convention in 
//		Denver, Colorado, April 30-May 1, 1985, American Society of Civil Engineers, New York, 293-299.
//
//	Rawls, W.J., L.R. Ahuja, and D.L. Brakensiek. 1992. Estimating soil hydraulic properties from
//            soils data. In van Genuchten, M.T. and F.J. Leij (Eds.). Indirection Methods for Estimating
//            the Hydraulic Properties of Unsaturated Soils, Proceedings of the International Workshop,
//            Riverside, California, October 11-13, 1989.
//
// Added 11/2008 and modified 04/2009 by DSM to minimize on the number of parameters needed to characterize
// soil texture. Parameter inputs are porosity, percent clay, and percent sand. 

void soil_texture_parameters(   double por,
                                double pClay,
                                double pSand,
				double& ks,
                                double& bubbling_pressure,
                                double& pore_size_index,
                                double& residual)
{

//Define boundary values for percent Clay and percent Sand
//Currently assumes that pClay + pSand don't exceed 100 percent
        if (pClay < 5.000001)
	{
                pClay = 5.000001;
	}
        if (pClay > 59.999999)
	{
                pClay = 59.999999;
	}
        if (pSand < 5.000001)
	{
                pSand = 5.000001;
	}
        if (pSand > 69.999999)
	{
                pSand = 69.999999;
	}

        double por2 = por*por;
        double pClay2 = pClay*pClay;
        double pSand2 = pSand*pSand;

	ks = exp(19.52348*por - 8.96847 - 0.028212*pClay + 0.00018107*pSand2 - 0.0094125*pClay2 - 8.395215*por2 + 0.077718*pSand*por - 0.00298*pSand2*por2 - 0.019492*pClay2*por2 + 0.0000173*pSand2*pClay + 0.02733*pClay2*por + 0.001434*pSand2*por - 0.0000035*pClay2*pSand);

        bubbling_pressure = exp(5.33967 + 0.1845 * pClay -2.483945*por - 0.00213853*pClay2 - 0.04356*pSand*por - 0.61745*pClay*por + 0.00143598*pSand2*por2 - 0.00855375*pClay2*por2 - 0.00001282*pSand2*pClay + 0.00895359*pClay2*por - 0.00072472*pSand2*por + 0.0000054*pClay2*pSand + 0.50028*por2*pClay);

        pore_size_index = exp(-0.7842831 + 0.0177544*pSand - 1.062498*por - 0.00005304*pSand2 - 0.00273493*pClay2 + 1.111349*por2 - 0.03088295*pSand*por + 0.00026587*pSand2*por2 - 0.00610522*pClay2*por2 - 0.00000235*pSand2*pClay + 0.00798746*pClay2*por - 0.00674491*por2*pClay);

        residual = -0.0182482 + 0.00087269*pSand + 0.00513488*pClay + 0.02939286*por - 0.00015395*pClay2 - 0.0010827*pSand*por - 0.00018233*pClay2*por2 + 0.00030703*pClay2*por - 0.0023584*por2*pClay;

}

//Calculate root zone soil moisture
//Added October 24, 2014 DSM
void root_zone_moisture(double Zw,
			double capFringe,
			double Droot,
			int rmodules,
			double rootDepth[],
			double thetaSoil[],
			double porosity,
			double mineral_fraction)
{
	double cumDepth, suzRoot[ULAT];

	if ((Zw-capFringe) < Droot)
	{
		cumDepth = 0.0;
		for (int i = 0; i < rmodules; i++)
		{
			cumDepth += rootDepth[i];
			if ((Zw-capFringe) >= cumDepth)
			{
				suzRoot[i] = thetaSoil[i]*rootDepth[i]*mineral_fraction;
			}
			else if ((Zw-capFringe) >= ((cumDepth-rootDepth[i]) && (Zw-capFringe) < cumDepth))
			{
				suzRoot[i] = (cumDepth-(Zw-capFringe))*porosity + 
					(Zw-capFringe-(cumDepth-rootDepth[i]))*thetaSoil[i]*mineral_fraction;
				//thetaSoil[i] = suzRoot[i]/rootDepth[i]*(1/mineral_fraction);
			}
			else
			{
				suzRoot[i] = porosity*rootDepth[i]*mineral_fraction;
				thetaSoil[i] = suzRoot[i]/rootDepth[i]*(1/mineral_fraction);
			}
		}
	}
}

//Added June 18, 2009 DSM
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
                        double& theta_root_deep)
{
        double Droot, suz_root_shallow, suz_root_mid, suz_root_deep;
	
	Droot = Dshallow + Dmid + Ddeep;
        if ((Zw-capFringe) >= Droot)
        {
                suz_root_shallow = theta_shallow*Dshallow*mineral_fraction;
                suz_root_mid = theta_mid*Dmid*mineral_fraction;
                suz_root_deep = Ddeep*theta_deep*mineral_fraction;
        }
	else if ((Zw-capFringe) >= (Dshallow+Dmid) && (Zw-capFringe) < Droot)
        {
                suz_root_shallow = theta_shallow*Dshallow*mineral_fraction;
                suz_root_mid = theta_mid*Dmid*mineral_fraction;
                suz_root_deep = ((Droot-(Zw-capFringe))*porosity + 
					(Zw-capFringe-Dshallow-Dmid)*theta_deep)*mineral_fraction;
	}
        else if ((Zw-capFringe) >= Dshallow && (Zw-capFringe) < (Dshallow+Dmid))
        {
                suz_root_shallow = theta_shallow*Dshallow*mineral_fraction;
                suz_root_mid = ((Dshallow+Dmid-(Zw-capFringe))*porosity + 
					(Zw-capFringe-Dshallow-Dmid)*theta_mid)*mineral_fraction;
		suz_root_deep = Ddeep*porosity*mineral_fraction;
        }
        else if ((Zw-capFringe) > 0.0 && (Zw-capFringe) < Dshallow)
        {
                suz_root_shallow = ((Dshallow-(Zw-capFringe))*porosity + 
					(Zw-capFringe)*theta_shallow)*mineral_fraction;
                suz_root_mid = Dmid*porosity*mineral_fraction;
                suz_root_deep = Ddeep*porosity*mineral_fraction;
        }
        else
        {
                suz_root_shallow = Dshallow*porosity*mineral_fraction;
                suz_root_mid = Dmid*porosity*mineral_fraction;
                suz_root_deep = Ddeep*porosity*mineral_fraction;
        }

        theta_root_shallow = suz_root_shallow/Dshallow*(1.0/mineral_fraction);
        theta_root_mid = suz_root_mid/Dmid*(1.0/mineral_fraction);
        theta_root_deep = suz_root_deep/Ddeep*(1.0/mineral_fraction);
}


//Calculate resistance to evaporation from the soil surface
//as a function of soil water content
//returns resistance in m2 s mol-1
double theta_soil_evap_reduction(double reference_moisture,
                                 double residual,
                                 double theta)
{
        double Sr;

        Sr = 4104.0*(reference_moisture - theta) - 805.0; //units are s m-1

	if (Sr < 42.0)
	{
		Sr = 42.0;
	}

	Sr /= 42.0;

	return(Sr); //convert to molar units
}

//calculate soil moisture limited conductance rate for evaporation
//uses Darcy's Law + Fick's Law
//depends on:
//    1) vapor pressure deficit at the air-soil interface
//    2) soil water potential
//    3) water flow rate to evaporating site - unsaturated hydraulic conductivity
//rate returned in mol m-2 s-1
//DSM - 02/22/2010
//assumes air_entry_pressure in meters
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
			double mineral_fraction)
{
	double cond, ku, S, n, m, psi1, psi2, VPD_relative;
	double store, min_store;

	if (theta <= residual)
	{
		return(1.0E-10);
	}

	VPD_relative = VPD/sat_vapor_pressure(Tsoil); //unitless, relative measure of VPD
	if (VPD_relative < 0.000001)
	{
		VPD_relative = 0.000001;
	}

	store = dZ*theta*mineral_fraction;
	min_store = dZ*residual*mineral_fraction;

//calculate hydraulic conductivity for evaporating soil layer
        n = pore_size_index + 1.0;
        m = pore_size_index / n;
        S = (theta - residual) / (porosity - residual);
        if (S < 0.000001)
	{
                S = 0.000001;
	}
        else if (S > 1.0)
	{
                S = 1.0;
	}
        ku = ko * pow(S, n) * pow(1.0 - pow(1.0 - pow(S, 1.0/m), m), 2.0);
//calculate soil water potential for evaporating soil layer
        psi1 = -air_entry_pressure*pow(pow(S,-1.0/m)-1.0,1.0/n);
//calculate water potential at air-soil surface, assuming atmospheric VPD just above soil surface
	psi2 = -2243.95;
//calculate flux rate using Darcy's Law + Fick's Law analogy
//evap = cond * D && evap = ku * dpsi/dz => cond = ku/D * dpsi/dz
//                                          m/s = m/s/(kPa/kPa) * m/m
	cond = ku/VPD_relative*((psi1-psi2)/dZ);
	cond *= 42.0; //convert from m s-1 to mol m-2 s-1
        if (cond < 1.0E-10)
	{
               cond = 1.0E-10;
	}

    	return(cond);
}


//Capillary rise function
//Uses a steady-state solution of Richards equation
//References Gardner (1958), Eagleson (1978), Famiglietti and Wood (1994)
//Warning: becomes undefined as water table depth + suction head approach zero
//         an upper bound has been added to avoid this problem
//Units: 
//
double capillary_rise(double water_table_depth,
                      double suction_head,
                      double pore_size_index,
                      double sat_hydraulic_cond)
{

    	double C, m, a, b, w, zb;

// conversion of units,  m to cm 

    	water_table_depth  = water_table_depth  * 100.0;

	if (water_table_depth <= suction_head)
	{
		return(0.0);
	}

// m/day to cm/s

    	sat_hydraulic_cond = sat_hydraulic_cond * 100.0 / 86400.0;

    	m = pore_size_index;
    	b = 2.0 + 3.0*m;
    	C = 1.0 + (1.5/(1.0+3.0*b));
    	a = sat_hydraulic_cond*pow(suction_head, b);


// guard against potential explosion
// zb = upper bound of water movement (ie solve equation for z.
// when zb > water_table_depth then the amount of water that can
// enter the unsaturated zone is restricted  to that amount

    	zb = pow(C*a/water_table_depth, 1.0/b) + suction_head;

    	if (zb > water_table_depth)
	{
        	w = water_table_depth/86400.0;  // cm/s
	}
    	else
	{
        	w = (C * a)/pow((water_table_depth-suction_head), b);
	}

// convert units back */

	if (w < 0.0)
	{
		w = 0.0;
	}
	else
	{
    		w = w * 86400.0 / 100.0;      // cm/s to m/day
		w = w / 48.0;		    // m/day to m/30-min
	}

    	return(w);
}


/* Soil water retention functions using the van Genuchten (1980) equation */
/*       Added 6/16/99 DSM */
double field_capacity(double air_entry_pressure,
                      double pore_size_index,
                      double residual,
                      double porosity)
{
	double fc, m, n, fc_potential;

    	n = pore_size_index + 1.0;
    	m = pore_size_index / n;
    	fc_potential = 340.0;

    	fc = (porosity - residual) * pow(1.0 + pow(fc_potential/air_entry_pressure, n), -m) + residual;

    	return(fc);
}


//
//compute_effective_precipitation()
//
double compute_effective_precipitation(State_Store& state,
			      		trees_params treesParams,
                                        double canopy_store_max,
                                        double precip,
				      	double t_ref)
{
	double canopy_store, interception_capacity, interception, eff_precip, snowpack, snowmelt, snowpack_E_deficit;

	canopy_store = state.get_val_at(CANOPY);
	interception_capacity = canopy_store_max - canopy_store;
        if(treesParams.useLeafModule == true){
            interception_capacity = 0.0;
        }
        if (precip > 0.0 && interception_capacity > 0.0)
        {
                interception = min(precip, interception_capacity);
                eff_precip = precip - interception;
                state.update_val_at(interception, CANOPY);
        }
        else
        {
                eff_precip = precip;
        }

//SNOWPACK - add to snowpack if precipitation falls as snow, set eff_precip for drainage to zero
        snowpack = state.get_val_at(SNOWPACK);
        snowpack_E_deficit = state.get_val_at(SNOWEDEF);
        snowpack_E_deficit = snow_E_acc(snowpack_E_deficit, t_ref,
                                        treesParams.snowpack_E_deficit_max, snowpack);
        state.set_val_at(snowpack_E_deficit, SNOWEDEF);
//if precipitation is snow
        if (eff_precip > 0.0 && t_ref <= 0.0)
        {
                snowpack += eff_precip;
                eff_precip = 0.0;
                state.set_val_at(snowpack, SNOWPACK);
        }
//if rain falls on snowpack
        if (eff_precip > 0.0 && snowpack > 0.0 && t_ref > 0.0)
        {
                snowmelt = rainMelt(t_ref, eff_precip, snowpack, snowpack_E_deficit);
                eff_precip += snowmelt;
                snowpack -= snowmelt;
                if (snowpack < 0.0)
                {
                        eff_precip -= snowpack;
                        snowpack = 0.0;
                }
                state.set_val_at(snowpack, SNOWPACK);
        }
	return(eff_precip);
}
				   


//
//DRAINAGE: precipitation --> litter layer --> shallow soil layer --> deep soil layer
//
//Note: drainage at this point occurs for one-half the time step, or at present 15 min
//      so that other processes are influenced by the mean soil moisture at this time step;
//Future: Add flow to water table, currently using measured water table as a flux boundary
//        update soil moisture with precipitation and drainage
//        bypassFlow elements are accumulators, so need to initialize them to zero first
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
		  double ks)
{
	double bypass, bpp, drainMax, retentionMax, updatedTheta;
	double litter_store, drain, delta_Sr;
	double kinter;
	int kiter;

        for (int i = 0; i < rmodules; i++)
        {
                bypassFlow[i] = 0.0;
        }

//1) drainage fron litter layer
//allow a small amount of precipitation to bypass the litter zone
        if (input > 0.0)
        {
                drain = min(0.00005, input);
                input -= drain;
        }
        else
        {
                drain = 0.0;
        }
        litter_store = state.get_val_at(LITTER);
        litter_store += input;
        if (litter_store > litter_capacity)
        {
                drain += litter_store - litter_capacity;
                litter_store = litter_capacity;
        }
        state.set_val_at(litter_store, LITTER);

//2) drainage through shallow layer, add drainage from litter layer and then drain the shallow layer
        ks = ks * 0.01 * 0.5 * 0.5;   //convert saturated hydraulic conductivity from cm/hr to m/15-min
//define a by-pass flow to the next soil layer - assumes a film flow in a macropore
        bpp = 1.2; //exponent to control rate of change of bypass flow with input "drain"

        delta_Sr = drain;
        drain = 0.0;

        if (delta_Sr > 0.0) //Set small time steps for infiltration of new water
        {
                kiter = 30; //intervals of 30 seconds
                kinter = 30.0;
        }
        else //Set 30 second time steps if soils are wetter than field capacity
        {
                kiter = 0; //if soils are dry then no need to call drainage
                kinter = 1.0;
                for (int i = 0; i < rmodules-1; i++)
                {
                        if (thetaSoil[i] > field_capacity)
                        {
                                kiter = 30; //intervals of 30 seconds
                                kinter = 30.0;
                        }
                }
        }
        while (kiter > 0)
        {
                kiter -= 1;
                for (int i = 0; i < rmodules; i++)
                {
                        bypass = 0.0;
                        if (delta_Sr > 0.0)
                        {
                                drainMax = 0.001*pow(80.0/(bpp+1.0),1.0/bpp);
                                retentionMax = drainMax - (pow(1000.0*drainMax,bpp)*0.01+0.2)*drainMax;
                                if (delta_Sr > retentionMax)
                                {
                                        bypass = delta_Sr - retentionMax;
                                }
                                else
                                {
                                        bypass = (pow(1000.0*delta_Sr,bpp)*0.01+0.2)*delta_Sr;
                                }
                        }
//remove bypass flux from ith layer
                        delta_Sr -= bypass;
//update soil theta for the ith layer
                        updatedTheta = thetaSoil[i]+delta_Sr/(rootDepth[i]*mineral_fraction);
                        drain = 0.0;
                        if (updatedTheta > porosity)
                        {
                                bypass += (updatedTheta - porosity)*(rootDepth[i]*mineral_fraction);
                                updatedTheta = porosity;
                        }
                        if (i < rmodules-1)
                        {
                                drain += unsat_zone_drainage(rootDepth[i], rootDepth[i+1], pore_size_index,
                                                                ks/kinter, updatedTheta, thetaSoil[i+1],
                                                                field_capacity, residual, capFringe,
                                                                porosity, mineral_fraction);
                        }
                        else
                        {
                                drain += unsat_zone_drainage(rootDepth[i], rootDepth[i], pore_size_index,
                                                                ks/kinter, updatedTheta, field_capacity,
                                                                field_capacity, residual, capFringe,
                                                                porosity, mineral_fraction);
                                drain *= drainScalar;
                        }
                        delta_Sr -= drain;
                        thetaSoil[i] += delta_Sr/(rootDepth[i]*mineral_fraction);
			if (thetaSoil[i] > porosity)
			{
				bypass += (thetaSoil[i] - porosity)*(rootDepth[i]*mineral_fraction);
				thetaSoil[i] = porosity;
			}
                        delta_Sr = drain + bypass;
                        bypassFlow[i] += bypass;
                }
                drain = delta_Sr = 0.0;
        }
        if (thetaSoil[rmodules-1] > porosity)
        {
                thetaSoil[rmodules-1] = porosity;
        }
}



// Unsaturated drainage 
//    Uses:
//         (1) unsaturated hydraulic conductivity
//	   (2) soil water potential gradient
// Use Darcy's Law + Conservation equation
// A quasi-explicit numerical approximation of Richards Eqn, 
//    probably best if layers are thin and short time steps are used
// DSM 2009/2010
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
                           double mineral_fraction)
{
    	double drain, totalDrain, ku, S, n, m, fc, suz1, suz2, psi1, psi2;
	double Z1, Z2;
	int iter;

	if (dZ2 <= 0.0)
	{
		return(0.0);
	}
	if (dZ1 < 0.0001)
	{
		return(0.0001);
	}
	if (theta1 <= field_capacity)
	{
		return(0.0);
	}

    	drain = totalDrain = 0.0;
	if (theta1 > 0.9*porosity)
	{
		iter=30; //step time down as low as 100 msec
		ko = ko/30.0;
	}
	else if (theta1 > 0.67*porosity)
	{
		iter=6;
		ko = ko/6.0;
	}
	else if (theta1 > field_capacity)
	{
		iter = 3;
		ko = ko/3.0;
	}
	else
	{
		iter = 1;
		ko = ko/1.0;
	}

	Z1 = 0.5*dZ1;
	Z2 = 0.5*(dZ1 + dZ2);

	fc = field_capacity*dZ1*mineral_fraction;
	for (int i = 0; i < iter; ++i)
	{
//no drainage if draining layer is already at field capacity;
		suz1 = theta1*dZ1*mineral_fraction;
		suz2 = theta2*dZ2*mineral_fraction;

//calculate hydraulic conductivity for draining soil layer
		n = pore_size_index + 1.0;
		m = pore_size_index / n;
		S = (theta1 - residual) / (porosity - residual);
		if (S <= 0.001)
		{
			S = 0.001;
		}
		else if (S > 1.0)
		{
			S = 1.0;
		}
        	//ku = ko * pow(S, n) * pow(1.0 - pow(1.0 - pow(S, 1.0/m), m), 2.0);
        	ku = ko * pow(S, 0.5) * pow(1.0 - pow(1.0 - pow(S, 1.0/m), m), 2.0);
//calculate soil water potential for draining soil layer
		psi1 = air_entry_pressure*pow(pow(S,-1.0/m)-1,1.0/n);
//calculate soil water potential for receiving soil layer
		S = (theta2 - residual) / (porosity - residual);
		if (S <= 0.001)
		{
			S = 0.001;
		}
		else if (S > 1.0)
		{
			S = 1.0;
		}
		psi2 = air_entry_pressure*pow(pow(S,-1.0/m)-1.0,1.0/n);
//calculate drainage using Darcy's Law
		drain = ku*((psi2-psi1)/(0.5*(dZ1+dZ2))+1.0);
		//drain = ku*((psi2-psi1)/(Z2-Z1)+1.0);
//update total drainage and moisture contents
		totalDrain += drain;
		suz1 -= drain;
		suz2 += drain;
		theta1 = suz1/(dZ1*mineral_fraction);
		theta2 = suz2/(dZ2*mineral_fraction);
	}
	if (totalDrain > (suz1-fc))
	{
		totalDrain = suz1-fc;
	}
	if (totalDrain < 0.0)
	{
		totalDrain = 0.0;
	}

	return(totalDrain);
}

// Unsaturated drainage depends on unsaturated hydraulic conductivity
double unsat_zone_drainage(double Zw,
                           double pore_size_index,
                           double ko,
                           double theta,
                           double field_capacity,
                           double residual,
                           double air_entry_pressure,
                           double porosity,
			   double mineral_fraction)
{
    	double drain, ku, S, n, m, fc, suz;

    	drain = 0.0;

    	S = (theta - residual) / (porosity - residual);
    	if (S > 1.0)
	{
		S = 1.0;
	}
    	fc = field_capacity * Zw * mineral_fraction;
    	suz = theta * Zw * mineral_fraction;

    	if (S > 0.0 && suz > fc)
    	{
        	n = pore_size_index + 1.0;
        	m = pore_size_index / n;
        	ku = ko * pow(S, n) * pow(1.0 - pow(1.0 - pow(S, 1.0/m), m), 2.0);
        	drain = min(ku, suz - fc);
    	}

    	return(drain);
}

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
				   double canopy_cover)
{
	double canopy_evaporation, canopy_store, can_evap_sun, can_evap_shd;

//get canopy storage and wetness
        canopy_store = state.get_val_at(CANOPY);
        canopy_wetness = canopy_store/canopy_store_max;
//calculate canopy evaporation
        canopy_evaporation = 0.0;
        if (precip < 0.0005)
        {
                if (canopy_wetness > 1.0E-6)
                {
                        can_evap_sun = do_pm(Rnet_sun, T_sun, D, p_air, gHr, gva);
                        can_evap_shd = do_pm(Rnet_shd, T_shd, D, p_air, gHr, gva);
                        canopy_evaporation = canopy_wetness * canopy_cover * (can_evap_sun + can_evap_shd);
                        if (canopy_evaporation > canopy_store)
                        {
                                canopy_evaporation = canopy_store;
                        }
                        if (canopy_evaporation < 0.0)
                        {
                                canopy_evaporation = 0.0;
                        }
                        state.update_val_at(-canopy_evaporation, CANOPY);
                }
                else
                {
                        canopy_evaporation = canopy_store;
                        state.set_val_at(0.0, CANOPY);
                }
        }
	return(canopy_evaporation);
}

//uses Penman Monteith to get Ecanopy, mm s-1 (m-2 gr area) - terms as 14.12
//calculates s, lambda inside
//assume pa = 96 kPa @ chequamegon 450m el
double do_pm(
                double Rnet,       //net absorbed radiation W/m2
                double t_ref,      //deg C @ ref ht (measured)
                double Dee,        //vapor deficit, D (kPa)
                double pa,         //atmospheric pressure, (kPa)
                double gHr,        //combined boundary layer & radiative cond
                double gv)         //vapor cond. (both cond. mol m-2 s-1)
{
//canopy transpiration, mm s-1 m-2 gr area
//latent heat of vaporization, Jkg-1, Jmol-1
//apparent psychrometer const, C-1
        double  Ecanopy, lambda, l_mol, gamma_star; 
        double ess;         //delta/pa, C-1
        double num_R, num_D; //tem vars

//from constants.h cp_air (J mol-1 C-1) -042004, sigma(W m-2 K-4)
//calculations
        ess = calc_s(t_ref, pa);
        lambda = lheat_vap(t_ref);   //Jkg-1
        l_mol  = lambda*18.02*pow(10.0,-3.0);       //convert to Jmol-1
        gamma_star = cp_air*gHr/(l_mol*gv); //C-1, C&N, pg232 042004

//enable for check
/*
        double temp1 = emissivity*sigma*pow(t_ref, 4.0); //emission from canopy (W m-2)
	cout << temp1 << '\t';
*/
	//053104
	if(Rnet > 0.0)
	{
        	num_R = ess*Rnet;         //radiation driven term (W m-2 K-1)
	}
	else
	{
		num_R = 0.0;
	}
        num_D = gamma_star*l_mol*gv*Dee/pa;       //drying power of air driven term
                                                  //W m-2 C-1
        Ecanopy = (num_R + num_D)/(ess + gamma_star); //W m-2
        Ecanopy = Ecanopy/lambda;  //kg m-2 s-1 or mm s-1 (m-2 gr area)

//enable for check
/*
	double Rec = ess*Rnet/(ess + gamma_star); //Wm-2
	Rec = Rec/lambda; //kg m-2 s-1 or mm s-1

	double Dec = (gamma_star*gv*Dee*l_mol)/(pa*(ess + gamma_star)); //Wm-2
	Dec = Dec/lambda; //kg m-2 s-1 or mm s-1

cout << lambda << '\t' << l_mol << '\t' << gamma_star  << '\t' 
	<< ess << '\t' << Rec << '\t' << Dec << '\t' << Ecanopy << '\t'; 
*/

        return (Ecanopy);
}


//accumulate snowpack degree-day deficit (deg. C days)
//tavg/48 to convert from 30-min to days
double snow_E_acc(double old_acc,
                  double tavg,
                  double maxdeficit,
                  double snow)
{
    	double new_acc;

    	if (snow >  0.000000001)
	{
        	new_acc =max((old_acc+tavg/48.0),maxdeficit);
	}
    	else
	{
        	new_acc  = 0.0001;
	}
    	return(new_acc);
}


//calculate latent heat melting of snowpack water equivalent in m 30-min-1
double latentMelt(double snow_E_deficit,
            double Tair,
            double snowpack,
            double vpd,
            double wind_speed)
{
    	double latent, turb_melt, esd, ea;

    	esd = 6.1078 * exp((17.269*Tair)/(237.3 + Tair));

    	ea = esd - vpd;

    	wind_speed = 0.1*wind_speed;

/* Source: Berris and Harr, 1987 */
    	latent = 0.0011 * pow(30.0*30.0, -1.0/6.0) * (1000.0*ea - 611.0) * wind_speed;
//convert from mm hour-1 to m 30-min-1
    	latent = latent * 0.5 * 0.001;

    	if (latent < 0.0) 
	{
		turb_melt = max(-snowpack, latent);
	}
    	else 
	{
		turb_melt = min(snowpack, latent);
	}

    	return(turb_melt);
}


//melt snowpack as a function of radiant energy
double radMelt(double snow_E_deficit,
            double Rad_floor,
            double Tair,
            double melt_Rcoef,
            double snowpack)
{
    	double melt, Rad_net;

    	snow_E_deficit = min(0.0, snow_E_deficit);
    	snow_E_deficit = 0.0;

    /*  RADIATION  MELT */

    	Rad_net = Rad_floor * 1800.0 / 1000.0; //W m-2 converted to KJ m-2 over 30 minutes
    	if (Rad_net < 0.0) 
	{
		Rad_net = 0.0;
	}
    	if ((Tair > 0.0) && (snow_E_deficit >= 0.0))
	{
        	melt = max((Rad_net*melt_Rcoef/LFUSION), 0.0);
	}
    	else
	{
        	melt = 0.0;
	}
    	melt = min(melt, snowpack);

    	return(melt);
}

/* Heat transfer from rain falling on snowpack */
/* Q = density_water * specific_heat_water * wet_bulb_temperature * rainfall_input */
/*        kg m-3             kJ kg-1 C-1             C                  m3 m-2 */

double rainMelt(double dewpoint,
                double rainfall,
                double snowpack,
                double &snow_E_deficit)
{
    	double melt, temp, Tsnow, Train;

/* need rain melt - refreezing, so need snowpack temperature here as well */

    	dewpoint = max(dewpoint, 1.0);
    	Tsnow = snow_E_deficit + 273.1;
    	Train = dewpoint + 273.1;
    	temp = (Tsnow * snowpack + Train*rainfall)/ (rainfall + snowpack) - 273.1;
    	melt = 1000.0*4.218*dewpoint*rainfall/LFUSION;

/* Negative melt is returned if refreezing occurs */

    	return(melt);
}


double sensibleMelt(double snow_E_deficit,
            double Tair,
	    double Pair,
            double snowpack,
            double vpd,
            double wind_speed)
{
    	double sensible, Po, turb_melt;

    	Po = 101.0;
    	snow_E_deficit = min(0.0, snow_E_deficit);
    	snow_E_deficit = 0.0;
    	wind_speed = 0.1*wind_speed;

/* Source: Berris and Harr, 1987 */
    	sensible = 0.061 * (Pair/Po)*pow(30.0*30.0,-1.0/6.0) * (Tair - snow_E_deficit)*wind_speed;
//convert from mm hour-1 to m 30-min-1
    	sensible = sensible * 0.5 * 0.001;

    	turb_melt = min(snowpack, max(0.0,sensible));

    	return(turb_melt);
}


double sublimate(double snow_E_deficit,
                 double Rad_floor,
                 double Tair,
                 double melt_Rcoef,
                 double snowpack)
{
    	double sublimation, Rad_net;

    	Rad_net = Rad_floor*1800.0/1000.0; //convert from W m-2 to KJ m-2 over 30-min period
    	if (Rad_net < 0.0) 
	{
		Rad_net = 0.0;
	}
    	sublimation = max((Rad_net*melt_Rcoef/(LFUSION+LVAPOUR)),0.0);

    /* CORRECT FOR NEGATIVE SNOWMELT */

    	if (sublimation > snowpack)
	{
        	return(snowpack);
	}
    	else
	{
        	return(sublimation);
	}
}


double tempMelt(double snow_E_deficit,
                double Tair,
                double melt_Tcoef)
{
    	double Tmelt;

    	if ((Tair > 0.0) && (snow_E_deficit >= 0.0))
	{
        	Tmelt = melt_Tcoef * Tair;
	}
    	else 
	{
		Tmelt = 0.0;
	}
    	return(Tmelt);
}


