//*******************************************
//values of physical constants used in TREES
//SS - June, 2002
//*******************************************
#ifndef CONSTANTS_H
#define CONSTANTS_H

using namespace std;

const int DAYS = 365;
const int HOURS = 24;
const int MINS = 60;
const int SECONDS = 60;

const double sigma   = 5.67e-8; //stefan-boltzman constant W m-2 K-4
const double rho_w = 1000.0;  // density of water  kg m-3
const double mm_h2o = 18.02;	//molecular mass of H2O, g/mol
const double mm_air = 28.97;	//molecular mass of air, g/mol
const double vkk     = 0.4;     // von karmen constant
const double b_rb     = 0.0067;  // proportionality constant used to calc rb, m s^-1/2 (magnani)
const double gr_acc     = 9.8;     // gravitational acceleration  m s-2
const double cp_air    = 29.3;  // specific heat of air  J mol-1 C-1
const double Ssc = 1370.0;  // solar constant  W m-2

const double zr    = 30.0;    // reference height m, data specific

const double R	= 8.314472;	//Gas Law constant, m3Pa/mol K (or J/mol/K)

const double diff_h2o = 0.147;  //Leaf boundary layer diffisivity for water (mol m-2 s-1)
const double diff_co2 = 0.110;  //Leaf boundary layer diffisivity for carbon dioxide (mol m-2 s-1)

//used in the photosynthesis routine
/* the weight proportion of Rubisco to its nitrogen content, fnr, is
        calculated from the relative proportions of the basic amino acids
        that make up the enzyme, as listed in the Handbook of Biochemistry,
        Proteins, Vol III, p. 510, which references:
        Kuehn and McFadden, Biochemistry, 8:2403, 1969 */

const double fnr = 7.16;   /* kg Rub/kg NRub */

/* the following constants are from:
        Woodrow, I.E., and J.A. Berry, 1980. Enzymatic regulation of photosynthetic
                CO2 fixation in C3 plants. Ann. Rev. Plant Physiol. Plant Mol. Biol.,
                39:533-594.
        Note that these values are given in the units used in the paper, and that
        they are converted to units appropriate to the rest of this function before
        they are used.
        */
        /* I've changed the values for Kc and Ko from the Woodrow and Berry
        reference, and am now using the values from De Pury and Farquhar,
        1997. Simple scaling of photosynthesis from leaves to canopies
        without the errors of big-leaf models. Plant, Cell and Env. 20: 537-557.
         All other parameters, including the q10's for Kc and Ko are the same
         as in Woodrow and Berry. */

//const double Kc25 = 40.4;    /* (Pa) MM const carboxylase, 25 deg C */
//const double q10Kc = 2.1;    /* (DIM) Q_10 for kc */
//const double Ko25 = 24800.0; /* (Pa) MM const oxygenase, 25 deg C */
//const double q10Ko = 1.2;    /* (DIM) Q_10 for ko */
//const double act25 = 3.6;    /* (umol/mgRubisco/min) Rubisco activity */
//const double q10act = 2.4;   /* (DIM) Q_10 for Rubisco activity */
const double pabs = 0.85;    /* (DIM) photons absorbed by PSII per e- transported (0.15 is spectral correction) */
const double ppe = 2.0;      /* (DIM) fPAR effectively absorbed by PSII */

#endif
