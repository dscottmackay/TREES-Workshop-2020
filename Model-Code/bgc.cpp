//***********************************************************
//Implementation of carbon and nitrogen cycline (BGC) class
//
//This code is part of the TREES model
//  but a small effort needed to adapt to other models
// 2015-2018 DSM
//
// This object manages states and fluxes associated with:
//   - photosynthetic carbon sequestration
//   - allocation of carbon to structural pools
//   - allocation of non-structural carbon pools
//   - carbon and nitrogen cycles
//   - rhizosphere C and N dynamics
//
//******************** DSM 2015 *****************************

#include <random>
#include "simulator2.h"

BiogeochemicalCycles::BiogeochemicalCycles()
{
	leafResidueCarbon = NULL;
	leafResidueNitrogen = NULL;
	stemResidueCarbon = NULL;
	stemResidueNitrogen = NULL;
	rootResidueCarbon = NULL;
	rootResidueNitrogen = NULL;
	dCrootResidue = NULL;
	dNrootResidue = NULL;
	humusCarbon = NULL;
	humusNitrogen = NULL;
	soilAmmoniumNitrogen = NULL;
	soilNitrateNitrogen = NULL;
	rhizosphereCl = NULL;
	rhizosphereNl = NULL;
	rhizosphereLiveMicrobialCarbon = NULL;
	rhizosphereDeadMicrobialCarbon = NULL;
	rhizosphereLabileCarbon = NULL;
	rhizosphereMicrobialNitrogen = NULL;
	rhizosphereMineralNitrogen = NULL;
	rhizosphereAmmoniumNitrogen = NULL;
	rhizosphereNitrateNitrogen = NULL;
	rootExudateSugarCarbon = NULL;
	rootExudateAminoAcidCarbon = NULL;
	rootExudateAminoAcidNitrogen = NULL;
	maximumRootBiomassCarbon = NULL;
	rootAreaScalar = NULL;
	fineRootBiomassCarbon = NULL;
	fineRootBiomassNitrogen = NULL;
	coarseRootBiomassCarbon = NULL;
	coarseRootBiomassNitrogen = NULL;
	rootNSC = NULL;
	rootMineralNitrogen = NULL;
	rootArea = NULL;
	liveStemCarbon = NULL;
	stemNSC = NULL;
	liveStemNitrogen = NULL;
	deadStemCarbon = NULL;
	deadStemNitrogen = NULL;
	leafBiomassCarbon = NULL;
	leafBiomassNitrogen = NULL;
	leafNSC = NULL;
	chloroplastStarch = NULL;
	chloroplastSugar = NULL;
	glycineNitrogen = NULL;
	glycineCarbon = NULL;
	serineNitrogen = NULL;
	serineCarbon = NULL;
	leafStoredNitrogen = NULL;
	leafRubiscoNitrogen = NULL;
	fruitCarbon = NULL;
	fruitNitrogen = NULL;
	plantNstatus = NULL;
	nitrogenLeaching = NULL;
	lat_Root_b_value_init = NULL;
	lat_Root_c_value_init = NULL;
	heterotrophicRespiration = NULL;
	Lf_idx = NULL;
	phyllo_tracker = NULL;
	SLA_avg = NULL;
	Projected_LeafArea_total = NULL;
    	Karray = NULL;
    	Narray = NULL;
    	rarray = NULL;
	SingleLeafBiomassCarbon = NULL;
	SingleLeafBiomassNitrogen = NULL;
	SingleLeafNSC = NULL;
	SingleLeafchloroplastStarch = NULL;
	SingleLeafchloroplastSugar = NULL;
	SingleLeafStoredNitrogen = NULL;
	SingleLeafRubiscoNitrogen = NULL;
	SingleLeafGrowthRespiration = NULL;
	SingleLeafArea = NULL;
	SingleLeafAreaPotential = NULL;
	SingleLeafThermTm = NULL;
}

BiogeochemicalCycles::BiogeochemicalCycles(trees_params& treesParams)
{
//Number of soil-root layers not to exceed ULAT
	nRoots = treesParams.rmodules;
//Number of root sizes or orders
	nFineRootOrders = 5;
	nRootOrders = 10;

	allocatePools();
	initializePools(treesParams);
}

BiogeochemicalCycles::~BiogeochemicalCycles()
{
	clearPools();
}

//
//allocate memory to store BGC pools
//
void BiogeochemicalCycles::allocatePools()
{
	assert(nRoots > 0);
	assert(nRootOrders > 0);

//allocate pool for leaf residue C and N
	leafResidueCarbon = new double;
	assert(leafResidueCarbon != NULL);
	leafResidueNitrogen = new double;
	assert(leafResidueNitrogen != NULL);

//allocate pool for stem residue C and N
	stemResidueCarbon = new double;
	assert(stemResidueCarbon != NULL);
	stemResidueNitrogen = new double;
	assert(stemResidueNitrogen != NULL);

//allocate pools for root residue C and N, assuming we will lump fine and coarse root necromass
	rootResidueCarbon = new double*[nRoots];
	assert(rootResidueCarbon != NULL);
	rootResidueNitrogen = new double*[nRoots];
	assert(rootResidueNitrogen != NULL);
	dCrootResidue = new double*[nRoots];
	assert(dCrootResidue != NULL);
	dNrootResidue = new double*[nRoots];
	assert(dNrootResidue != NULL);

	for (int j = 0; j < nRoots; j++)
	{
		rootResidueCarbon[j] = new double[nRootOrders];
		assert(rootResidueCarbon[j] != NULL);
		rootResidueNitrogen[j] = new double[nRootOrders];
		assert(rootResidueNitrogen[j] != NULL);
		dCrootResidue[j] = new double[nRootOrders];
		assert(dCrootResidue[j] != NULL);
		dNrootResidue[j] = new double[nRootOrders];
		assert(dNrootResidue[j] != NULL);
	}

//allocate pool for recalcitrant (humus) C
	humusCarbon = new double[nRoots];
	assert(humusCarbon != NULL);
	humusNitrogen = new double[nRoots];
	assert(humusNitrogen != NULL);
	soilAmmoniumNitrogen = new double[nRoots];
	assert(soilAmmoniumNitrogen != NULL);
	soilNitrateNitrogen = new double[nRoots];
	assert(soilNitrateNitrogen != NULL);

//State variables for litter C and litter N
	rhizosphereCl = new double*[nRoots];
	assert(rhizosphereCl != NULL);
	rhizosphereNl = new double*[nRoots];
	assert(rhizosphereNl != NULL);

//allocate microbial pools
	rhizosphereLiveMicrobialCarbon = new double*[nRoots];
	assert(rhizosphereLiveMicrobialCarbon != NULL);
	rhizosphereDeadMicrobialCarbon = new double*[nRoots];
	assert(rhizosphereDeadMicrobialCarbon != NULL);
	rhizosphereLabileCarbon = new double*[nRoots];
	assert(rhizosphereLabileCarbon != NULL);
	rhizosphereMicrobialNitrogen = new double*[nRoots];
	assert(rhizosphereMicrobialNitrogen != NULL);
	rhizosphereMineralNitrogen = new double*[nRoots];
	assert(rhizosphereMineralNitrogen != NULL);
	rhizosphereAmmoniumNitrogen = new double*[nRoots];
	assert(rhizosphereMineralNitrogen != NULL);
	rhizosphereNitrateNitrogen = new double*[nRoots];
	assert(rhizosphereMineralNitrogen != NULL);

	for (int j = 0; j < nRoots; j++)
	{
		rhizosphereCl[j] = new double[nRootOrders];
		rhizosphereNl[j] = new double[nRootOrders];
		rhizosphereLiveMicrobialCarbon[j] = new double[nRootOrders];
		assert(rhizosphereLiveMicrobialCarbon[j] != NULL);
		rhizosphereDeadMicrobialCarbon[j] = new double[nRootOrders];
		assert(rhizosphereDeadMicrobialCarbon[j] != NULL);
		rhizosphereLabileCarbon[j] = new double[nRootOrders];
		assert(rhizosphereLabileCarbon[j] != NULL);
		rhizosphereMicrobialNitrogen[j] = new double[nRootOrders];
		assert(rhizosphereMicrobialNitrogen[j] != NULL);
		rhizosphereMineralNitrogen[j] = new double[nRootOrders];
		assert(rhizosphereMineralNitrogen[j] != NULL);
		rhizosphereAmmoniumNitrogen[j] = new double[nRootOrders];
		assert(rhizosphereAmmoniumNitrogen[j] != NULL);
		rhizosphereNitrateNitrogen[j] = new double[nRootOrders];
		assert(rhizosphereNitrateNitrogen[j] != NULL);
	}

//allocate heterotrophic respiration
	heterotrophicRespiration = new double[nRoots];
	assert(heterotrophicRespiration != NULL);

//allocate root pools
	rootExudateSugarCarbon = new double*[nRoots];
	assert(rootExudateSugarCarbon != NULL);
	rootExudateAminoAcidCarbon = new double*[nRoots];
	assert(rootExudateAminoAcidCarbon != NULL);
	rootExudateAminoAcidNitrogen = new double*[nRoots];
	assert(rootExudateAminoAcidNitrogen != NULL);
	maximumRootBiomassCarbon = new double*[nRoots];
	assert(maximumRootBiomassCarbon != NULL);
	rootAreaScalar = new double[nRoots];
	assert(rootAreaScalar != NULL);
	fineRootBiomassCarbon = new double*[nRoots];
	assert(fineRootBiomassCarbon != NULL);
	fineRootBiomassNitrogen = new double*[nRoots];
	assert(fineRootBiomassNitrogen != NULL);
	coarseRootBiomassCarbon = new double*[nRoots];
	assert(coarseRootBiomassCarbon != NULL);
	coarseRootBiomassNitrogen = new double*[nRoots];
	assert(coarseRootBiomassNitrogen != NULL);
	rootNSC = new double*[nRoots];
	assert(rootNSC != NULL);
	rootMineralNitrogen = new double*[nRoots];
	assert(rootMineralNitrogen != NULL);
	rootArea = new double*[nRoots];
	assert(rootArea != NULL);

	for (int j = 0; j < nRoots; j++)
	{
		rootExudateSugarCarbon[j] = new double[nRootOrders];
		assert(rootExudateSugarCarbon[j] != NULL);
		rootExudateAminoAcidCarbon[j] = new double[nRootOrders];
		assert(rootExudateAminoAcidCarbon[j] != NULL);
		rootExudateAminoAcidNitrogen[j] = new double[nRootOrders];
		assert(rootExudateAminoAcidNitrogen[j] != NULL);
		maximumRootBiomassCarbon[j] = new double[nRootOrders];
		assert(maximumRootBiomassCarbon[j] != NULL);
		fineRootBiomassCarbon[j] = new double[nRootOrders];
		assert(fineRootBiomassCarbon[j] != NULL);
		fineRootBiomassNitrogen[j] = new double[nRootOrders];
		assert(fineRootBiomassNitrogen[j] != NULL);
		coarseRootBiomassCarbon[j] = new double[nRootOrders];
		assert(coarseRootBiomassCarbon[j] != NULL);
		coarseRootBiomassNitrogen[j] = new double[nRootOrders];
		assert(coarseRootBiomassNitrogen[j] != NULL);
		rootNSC[j] = new double[nRootOrders];
		assert(rootNSC[j] != NULL);
		rootMineralNitrogen[j] = new double[nRootOrders];
		assert(rootMineralNitrogen[j] != NULL);
		rootArea[j] = new double[nRootOrders];
		assert(rootArea[j] != NULL);
	}

//allocate stem pools
	liveStemCarbon = new double;
	assert(liveStemCarbon != NULL);
	stemNSC = new double;
	assert(stemNSC != NULL);
	liveStemNitrogen = new double;
	assert(liveStemNitrogen != NULL);
	deadStemCarbon = new double;
	assert(deadStemCarbon != NULL);
	deadStemNitrogen = new double;
	assert(deadStemNitrogen != NULL);

//allocate leaf pools
	leafBiomassCarbon = new double;
	assert(leafBiomassCarbon != NULL);
	leafBiomassNitrogen = new double;
	assert(leafBiomassNitrogen != NULL);
	leafNSC = new double;
	assert(leafNSC != NULL);
	chloroplastStarch = new double;
	assert(chloroplastStarch != NULL);
	chloroplastSugar = new double;
	assert(chloroplastSugar != NULL);

	glycineNitrogen = new double;
	assert(glycineNitrogen != NULL);
	glycineCarbon = new double;
	assert(glycineCarbon != NULL);
	serineNitrogen = new double;
	assert(serineNitrogen != NULL);
	serineCarbon = new double;
	assert(serineCarbon != NULL);

	leafStoredNitrogen = new double;
	assert(leafStoredNitrogen != NULL);
	leafRubiscoNitrogen = new double;
	assert(leafRubiscoNitrogen != NULL);

//allocation reproductive pools
	fruitCarbon = new double;
	assert(fruitCarbon != NULL);
	fruitNitrogen = new double;
	assert(fruitNitrogen != NULL);

//allocation for plant nitrogen status
	plantNstatus = new double;
	assert(plantNstatus != NULL);

//allocation for the nitrogen leaching flux variable
	nitrogenLeaching = new double[nRoots];
	assert(nitrogenLeaching != NULL);

//allocation for lateral root Weibulls
	lat_Root_b_value_init = new double;
	assert(lat_Root_b_value_init != NULL);
	lat_Root_c_value_init = new double;
	assert(lat_Root_c_value_init != NULL);

//allocation for leaf growth module variables
    	Lf_idx = new int;
    	assert(Lf_idx != NULL);
    	phyllo_tracker = new double;
    	assert(phyllo_tracker != NULL);
    	SLA_avg = new double;
    	assert(SLA_avg != NULL);
    	Projected_LeafArea_total = new double;
    	assert(Projected_LeafArea_total != NULL);
    
        Karray = new double[500];
        assert(Karray != NULL);
        Narray = new double[500];
        assert(Narray != NULL);
        rarray = new double[500];
        assert(rarray != NULL);
    
    	SingleLeafBiomassCarbon = new double[500];
    	assert(SingleLeafBiomassCarbon != NULL);
    	SingleLeafBiomassNitrogen = new double[500];
    	assert(SingleLeafBiomassNitrogen != NULL);
    	SingleLeafNSC = new double[500];
    	assert(SingleLeafNSC != NULL);
    	SingleLeafchloroplastStarch = new double[500];
    	assert(SingleLeafchloroplastStarch != NULL);
    	SingleLeafchloroplastSugar = new double[500];
    	assert(SingleLeafchloroplastSugar != NULL);
    	SingleLeafStoredNitrogen = new double[500];
    	assert(SingleLeafStoredNitrogen != NULL);
    	SingleLeafRubiscoNitrogen = new double[500];
    	assert(SingleLeafRubiscoNitrogen != NULL);
    	SingleLeafGrowthRespiration = new double[500];
    	assert(SingleLeafGrowthRespiration != NULL);
    	SingleLeafArea = new double[500];
    	assert(SingleLeafArea != NULL);
    	SingleLeafAreaPotential = new double[500];
    	assert(SingleLeafAreaPotential != NULL);
    	SingleLeafThermTm = new double[500];
    	assert(SingleLeafThermTm != NULL);
}

//
//set BGC pools to initial values
//the assumption with this code is that many state variables are correlated
//   and can be estimates from a few input paramaters such as SLA and leaf lifespan
//Common C:N for meterials:
//    > deciduous leaves 32
//        > litter       64
//    > needle leaf     110
//        > litter      220
//    > microbes          8
//    > soil humus       10
//    > amino acids      24/14
//
void BiogeochemicalCycles::initializePools(trees_params treesParams)
{
	double totalRootArea, scalar, rootScalar, sumCl, leafN;
	double SLAscalar, SLAreference, totDepth, cumDepth;
	double CN, sumMaxC, numLayers;

	SLAreference = 22.0;
	SLAscalar = SLAreference / treesParams.SLA;
	if (SLAscalar > 3.0)
	{
		SLAscalar = 3.0;
	}
	if (SLAscalar < 1.0)
	{
		SLAscalar = 1.0;
	}

	totDepth = cumDepth = numLayers = 0.0;
	for (int j = 0; j < treesParams.rmodules; j++)
	{
		if (treesParams.ar[j+3] > 0.0)
		{
			totDepth += treesParams.drax[j+3];
		}
		numLayers += 1.0;
	}
//
//leaf-level state variables to initiate only if leafmodule is turned on
//
    	if (treesParams.useLeafModule == 1)
    	{
// initialize leaf module state variables at 0's
// the Single Leaf array holds leaf number 3 and higher
// the two emerged leaves at start of simulation are accounted in plant-level state variables

        	Lf_idx[0] = 0;
        	phyllo_tracker[0] = 0.0;
            
        	for (int j = 0; j < 500; j++)
        	{
                    //arrays that hold leaf growth parameters
                    Karray[j] = treesParams.leafAreaMax;
                    Narray[j] = treesParams.initialLeafSize;
                    rarray[j] = treesParams.leafArea_Rate;
                
                    //single leaf arrays
            		SingleLeafBiomassCarbon[j] = 0.0;
            		SingleLeafBiomassNitrogen[j] =0.0;
            		SingleLeafNSC[j] =0.0;
            		SingleLeafchloroplastStarch[j] =0.0;
            		SingleLeafchloroplastSugar[j] =0.0;
            		SingleLeafStoredNitrogen[j] =0.0;
            		SingleLeafRubiscoNitrogen[j] =0.0;
            		SingleLeafGrowthRespiration[j] =0.0;
            		SingleLeafArea[j] =0.0;
            		SingleLeafAreaPotential[j] =0.0;
            		SingleLeafThermTm[j] =0.0;
        	}

// add code here to modify initial values for plant-level state variables for leaf and stem
// based on initial projected shoot area parameter.
// leaf initialization
        	double shoot_init, leaf_init, stem_init;

        	shoot_init = 1.0/(treesParams.SLA_max*10.0)*treesParams.projectedArea_init; // in grams

        	leaf_init = 0.8*g_to_kg(shoot_init)/cm2_to_ha(treesParams.pot_size) ; // kg/ha
        	stem_init = 0.2*g_to_kg(shoot_init)/cm2_to_ha(treesParams.pot_size) ; // kg/ha

// leaf carbon
        	leafBiomassCarbon[0] = leaf_init/(1.0+treesParams.leafNSCscalar);
cout << "leafBiomassCarbon = " << leafBiomassCarbon[0] << endl;
        	leafNSC[0] = leafBiomassCarbon[0]*treesParams.leafNSCscalar; //retained

// leaf nitrogen
 //       	leafN = treesParams.Nleaf*treesParams.lai*10000.0; //retained
//		leafN *= treesParams.projectedArea_init / treesParams.pot_size;
		//double leafNscalar = min(3.0,log10(treesParams.microbiomeScalar+1))/3.0;
		//leafN = leafBiomassCarbon[0] / (treesParams.leaf_maxCN - leafNscalar*(treesParams.leaf_maxCN-treesParams.leaf_minCN));

		leafN = leafBiomassCarbon[0] / (0.5 * (treesParams.leaf_minCN + treesParams.leaf_maxCN));

		leafBiomassNitrogen[0] = leafN;
cout << "leafBiomassNitrogen = " << leafBiomassNitrogen[0] << endl;
		CN = leafBiomassCarbon[0]/leafBiomassNitrogen[0];
// rubisco N, chloroplast storage, and stored N eqns are retained
        	leafRubiscoNitrogen[0] = leafBiomassNitrogen[0] * treesParams.Nrubisco;

        	chloroplastStarch[0] = 0.01*leafRubiscoNitrogen[0]*leafBiomassCarbon[0]/leafBiomassNitrogen[0];
        	chloroplastSugar[0] = 0.0;

		//leafStoredNitrogen[0] = leafNSC[0]/CN;
		leafStoredNitrogen[0] = 0.0;

// stem initialization; very young stem tissue --> keep same ratios as leaf tissue
        	liveStemCarbon[0] = stem_init/(1+ treesParams.leafNSCscalar*0.4) ;
        	stemNSC[0] = treesParams.leafNSCscalar*0.4 * liveStemCarbon[0];
        	liveStemNitrogen[0] = 0.003 * liveStemCarbon[0];

// negligible dead stem tissue
        	deadStemCarbon[0] = 0.000001;
        	deadStemNitrogen[0] = 0.000001;

// leaf and stem residues set to close to 0
        	leafResidueCarbon[0] = 0.000001;
        	leafResidueNitrogen[0] = 0.000001;
        	stemResidueCarbon[0] = 0.000001;
        	stemResidueNitrogen[0] = 0.000001;
    	}

//Set initial leaf states for "big leaf" type simulation
	else
	{
//assume initial leaf residue is proportion to leaf biomass and increasing with declining SLA
//units are kgC ha-1 and kgN ha-1
//set leaf residue N to 50% of live leaf N, assuming retranslocation
		if (treesParams.usePhenology == true) //best for perennial plants
		{
        		leafResidueCarbon[0] = (1.0/treesParams.leafLifeSpan)*treesParams.lai / treesParams.SLA * 10000.0 * SLAscalar;
        		leafResidueNitrogen[0] = (0.5 * treesParams.Nleaf * treesParams.lai * 10000.0)*SLAscalar;
//set stem residue to 5% of standing stem biomass
//stem [N] is assumed to be low (0.3% of carbon)
        		stemResidueCarbon[0] = 0.05 * treesParams.Cstem;
        		stemResidueNitrogen[0] = 0.003 * stemResidueCarbon[0];
			if (treesParams.leafLifeSpan > 1.0)
			{
        			leafBiomassCarbon[0] = (1.0 - 1.0/treesParams.leafLifeSpan)*treesParams.lai/treesParams.SLA*10000.0+0.01;
			}
			else
			{
        			leafBiomassCarbon[0] = (1.0 - 0.9/treesParams.leafLifeSpan)*treesParams.lai/treesParams.SLA*10000.0+0.01;
			}
		}
		else //best for annual plants
		{
			leafBiomassCarbon[0] = treesParams.lai/treesParams.SLA*10000.0+0.01;
			liveStemCarbon[0] = treesParams.Csapwood;
			stemNSC[0] = treesParams.leafNSCscalar*0.4 * liveStemCarbon[0];
			liveStemNitrogen[0] = 0.003 * liveStemCarbon[0];
			deadStemCarbon[0] = treesParams.Cstem;
			deadStemNitrogen[0] = 0.003 * deadStemCarbon[0];
			leafResidueCarbon[0] = 0.00001;
			leafResidueNitrogen[0] = 0.000001;
        		stemResidueCarbon[0] = 0.00001;
        		stemResidueNitrogen[0] = 0.000001;
		}
        	leafNSC[0] = leafBiomassCarbon[0]*treesParams.leafNSCscalar;
		//leafN = min(leafBiomassCarbon[0]/treesParams.leaf_minCN, treesParams.Nleaf*treesParams.lai*10000.0);
		//leafN = treesParams.Nleaf*treesParams.lai*10000.0;
		leafN = leafBiomassCarbon[0] / (0.5 * (treesParams.leaf_minCN + treesParams.leaf_maxCN));
		leafBiomassNitrogen[0] = leafN;
        	//leafBiomassNitrogen[0] = max(0.1,(treesParams.leafLifeSpan-1.0)/treesParams.leafLifeSpan)*leafN;
        	leafBiomassNitrogen[0] = leafN;
		leafRubiscoNitrogen[0] = leafBiomassNitrogen[0] * treesParams.Nrubisco;
		chloroplastStarch[0] = 0.01*leafRubiscoNitrogen[0]*leafBiomassCarbon[0]/leafBiomassNitrogen[0];
		chloroplastSugar[0] = 0.0;
		CN = leafBiomassCarbon[0]/leafBiomassNitrogen[0];
cout << "leafBiomassCarbon[0] = " << leafBiomassCarbon[0] << "; " << leafBiomassNitrogen[0] << "; CN = " << CN << endl;
		leafStoredNitrogen[0] = leafNSC[0]/treesParams.leaf_maxCN;
	}

	glycineNitrogen[0] = 0.0;
	glycineCarbon[0] = 0.0;
	serineNitrogen[0] = 0.0;
	serineCarbon[0] = 0.0;

//set the upper limit for carbon allocation to each root
//  as a function of the soil volume defined in the param modules
	computeMaximumRootBiomassCarbon(treesParams);


//set root and rhizosphere state variables
//   there are nRoots = number of soil-root layers in TREES, defined in "param_mod" input file
//   there are 10 root orders, assuming the first 5 are fine root (from 1/4 mm up to 4 mm diameter)
        for (int j = 0; j < nRoots; j++)
        {
		scalar = 0.25 / max(0.03125, treesParams.minRootLifespan);
		rootScalar = 1.0;
		sumCl = 0.0;

		heterotrophicRespiration[j] = 0.0;

//maximumRootBiomassCarbon[j][k]
		sumMaxC = 0.0;
		for (int k = 0; k < nFineRootOrders; k++)
		{
			sumMaxC += maximumRootBiomassCarbon[j][k];
		}

		for (int k = 0; k < nFineRootOrders; k++)
		{
//State variables for litter C and litter N
                	fineRootBiomassCarbon[j][k] = treesParams.Croot * treesParams.ar[j+3] * maximumRootBiomassCarbon[j][k]/sumMaxC;
			if (fineRootBiomassCarbon[j][k] > maximumRootBiomassCarbon[j][k])
			{
				fineRootBiomassCarbon[j][k] = maximumRootBiomassCarbon[j][k];
			}
                	//fineRootBiomassNitrogen[j][k] = fineRootBiomassCarbon[j][k] / (20.0*rootScalar*SLAscalar);
                	fineRootBiomassNitrogen[j][k] = fineRootBiomassCarbon[j][k] / 
						(leafBiomassCarbon[0]/leafBiomassNitrogen[0]*rootScalar);
                	coarseRootBiomassCarbon[j][k] = 0.0;
                	coarseRootBiomassNitrogen[j][k] = 0.0;

                	rhizosphereCl[j][k] = (scalar * SLAscalar)*
				(fineRootBiomassCarbon[j][k] + coarseRootBiomassCarbon[j][k]);
			sumCl += rhizosphereCl[j][k];
			if (j < 3)
			{
				rhizosphereCl[j][k] += 0.2/3.0*treesParams.Clitter_frac * treesParams.Cbelowground;
			}
			CN = fineRootBiomassCarbon[j][k]/fineRootBiomassNitrogen[j][k];
                	rhizosphereNl[j][k] = rhizosphereCl[j][k] / (CN+20.0*(1.0-1.0/SLAscalar));
			rootResidueCarbon[j][k] = scalar*fineRootBiomassCarbon[j][k];
			rootResidueNitrogen[j][k] = scalar*fineRootBiomassNitrogen[j][k];
                	rhizosphereLiveMicrobialCarbon[j][k] = 0.12*rhizosphereCl[j][k];
//cout << rhizosphereLiveMicrobialCarbon[j][k] << '\t';
                	rhizosphereDeadMicrobialCarbon[j][k] = 0.2*rhizosphereLiveMicrobialCarbon[j][k];
			sumCl += rhizosphereLiveMicrobialCarbon[j][k] + rhizosphereDeadMicrobialCarbon[j][k];
                	rhizosphereMicrobialNitrogen[j][k] = rhizosphereLiveMicrobialCarbon[j][k] / (8.0+2.0*(1.0-1.0/SLAscalar));
                	rhizosphereMineralNitrogen[j][k] = 0.1*rhizosphereDeadMicrobialCarbon[j][k] / (8.0+2.0*(1.0-1.0/SLAscalar));
                	rhizosphereAmmoniumNitrogen[j][k] = 0.01 * rhizosphereMineralNitrogen[j][k];
                	rhizosphereNitrateNitrogen[j][k] = 0.99 * rhizosphereMineralNitrogen[j][k];

                	rootNSC[j][k] = treesParams.leafNSCscalar * fineRootBiomassCarbon[j][k];
                	rhizosphereLabileCarbon[j][k] = 0.1*rootNSC[j][k];
			rootMineralNitrogen[j][k] = 0.1*fineRootBiomassNitrogen[j][k];
//allow for a scaling of the initial nutrient levels in the microbiome
                        if (treesParams.microbiomeScalar < 0.0)
                        {
				treesParams.microbiomeScalar = 0.0;
			}
                        rhizosphereMineralNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereAmmoniumNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereNitrateNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereLiveMicrobialCarbon[j][k] *= treesParams.microbiomeScalar;
			rhizosphereDeadMicrobialCarbon[j][k] *= treesParams.microbiomeScalar;
			rhizosphereMicrobialNitrogen[j][k] *= treesParams.microbiomeScalar;

			rootExudateSugarCarbon[j][k] = 0.0;
			rootExudateAminoAcidCarbon[j][k] = 0.0;
			rootExudateAminoAcidNitrogen[j][k] = 0.0;

			scalar *= 0.80;
			rootScalar *= 1.25;
		}
		sumMaxC = 0.0;
		for (int k = nFineRootOrders; k < nRootOrders; k++)
		{
			sumMaxC += maximumRootBiomassCarbon[j][k];
		}
//cout << endl;
		for (int k = nFineRootOrders; k < nRootOrders; k++)
		{

                	fineRootBiomassCarbon[j][k] = 0.0;
                	fineRootBiomassNitrogen[j][k] = 0.0;
                	coarseRootBiomassCarbon[j][k] = treesParams.Croot_coarse * treesParams.ar[j+3] * maximumRootBiomassCarbon[j][k]/sumMaxC * 0.99;
                	coarseRootBiomassCarbon[j][k] += treesParams.Croot_coarse * treesParams.drax[j+3]/totDepth * 0.01;
			if (coarseRootBiomassCarbon[j][k] > maximumRootBiomassCarbon[j][k])
			{
				coarseRootBiomassCarbon[j][k] = maximumRootBiomassCarbon[j][k];
			}
                	coarseRootBiomassNitrogen[j][k] = coarseRootBiomassCarbon[j][k] / 
						(leafBiomassCarbon[0]/leafBiomassNitrogen[0]*rootScalar);

                	rhizosphereCl[j][k] = (scalar * SLAscalar)*
					(fineRootBiomassCarbon[j][k] + coarseRootBiomassCarbon[j][k]);
			sumCl += rhizosphereCl[j][k];
			if (j < 3)
			{
				rhizosphereCl[j][k] += 0.2/3.0*treesParams.Clitter_frac * treesParams.Cbelowground;
			}
			//CN = coarseRootBiomassCarbon[j][k]/coarseRootBiomassNitrogen[j][k];
                	rhizosphereNl[j][k] = rhizosphereCl[j][k] / (CN+20.0*(1.0-1.0/SLAscalar));
			rootResidueCarbon[j][k] = scalar*coarseRootBiomassCarbon[j][k];
			rootResidueNitrogen[j][k] = scalar*coarseRootBiomassNitrogen[j][k];
                	rhizosphereLiveMicrobialCarbon[j][k] = 0.12 * rhizosphereCl[j][k];
//cout << rhizosphereLiveMicrobialCarbon[j][k] << '\t';
                	rhizosphereDeadMicrobialCarbon[j][k] = 0.2*rhizosphereLiveMicrobialCarbon[j][k];
			sumCl += rhizosphereLiveMicrobialCarbon[j][k] + rhizosphereDeadMicrobialCarbon[j][k];
                	rhizosphereMicrobialNitrogen[j][k] = rhizosphereLiveMicrobialCarbon[j][k] / (8.0+2.0*(1.0-1.0/SLAscalar));
                	rhizosphereMineralNitrogen[j][k] = 0.1*rhizosphereDeadMicrobialCarbon[j][k] / (8.0+2.0*(1.0-1.0/SLAscalar));
                	rhizosphereAmmoniumNitrogen[j][k] = 0.01 * rhizosphereMineralNitrogen[j][k];
                	rhizosphereNitrateNitrogen[j][k] = 0.99 * rhizosphereMineralNitrogen[j][k];

                	rootNSC[j][k] = treesParams.leafNSCscalar * coarseRootBiomassCarbon[j][k];
                	rhizosphereLabileCarbon[j][k] = 0.1*rootNSC[j][k];
			rootMineralNitrogen[j][k] = 0.1*fineRootBiomassNitrogen[j][k];
//allow for a scaling of the initial nutrient levels in the microbiome
                        if (treesParams.microbiomeScalar < 0.0001)
                        {
				treesParams.microbiomeScalar = 0.0001;
			}
                        rhizosphereMineralNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereAmmoniumNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereNitrateNitrogen[j][k] *= treesParams.microbiomeScalar;
			rhizosphereLiveMicrobialCarbon[j][k] *= treesParams.microbiomeScalar;
			rhizosphereDeadMicrobialCarbon[j][k] *= treesParams.microbiomeScalar;
			rhizosphereMicrobialNitrogen[j][k] *= treesParams.microbiomeScalar;

			rootExudateSugarCarbon[j][k] = 0.0;
			rootExudateAminoAcidCarbon[j][k] = 0.0;
			rootExudateAminoAcidNitrogen[j][k] = 0.0;

			scalar *= 0.80;
		}
//cout << endl;
//State variables for stabilized soil carbon
		cumDepth += treesParams.drax[j+3];
		if (treesParams.ar[j+3] > 0.0)
		{
                	//humusCarbon[j] = treesParams.Csoil * 0.5*(treesParams.ar[j+3]+treesParams.drax[j+3]/totDepth) - sumCl;
                	humusCarbon[j] = treesParams.Csoil * treesParams.drax[j+3]/totDepth - sumCl;
		}
		else
		{
			humusCarbon[j] = 0.1;
		}
                humusNitrogen[j] = 0.09 * humusCarbon[j];

//State variables for mineral nitrogen in the bulk soil
		soilAmmoniumNitrogen[j] = 0.1*humusNitrogen[j];
		soilNitrateNitrogen[j] = 0.9*humusNitrogen[j];
//Flux variable for nitrogen leaching set to zero
		nitrogenLeaching[j] = 0.0;

        }

	totalRootArea = computeRootArea(treesParams);

        liveStemCarbon[0] = treesParams.Csapwood;
        stemNSC[0] = treesParams.leafNSCscalar*0.4 * liveStemCarbon[0];
        liveStemNitrogen[0] = 0.003 * liveStemCarbon[0];
        deadStemCarbon[0] = treesParams.Cstem;
        deadStemNitrogen[0] = 0.003 * deadStemCarbon[0];

        //leafStoredNitrogen[0] *= treesParams.microbiomeScalar;

        fruitCarbon[0] = 0.0;
        fruitNitrogen[0] = 0.0;

	plantNstatus[0] = 1.0; //values are between 0 and 1
//keep a record of the input lateral root Weibull curve parameters
	lat_Root_b_value_init[0] = treesParams.lat_Root_b_value;
	lat_Root_c_value_init[0] = treesParams.lat_Root_c_value;


}


//clear and delete BGC pool variables
void BiogeochemicalCycles::clearPools()
{
	if (leafResidueCarbon != NULL)
	{
		delete leafResidueCarbon;
		leafResidueCarbon = NULL;
	}
	if (leafResidueNitrogen != NULL)
	{
		delete leafResidueNitrogen;
		leafResidueNitrogen = NULL;
	}
	if (stemResidueCarbon != NULL)
	{
		delete stemResidueCarbon;
		stemResidueCarbon = NULL;
	}
	if (stemResidueNitrogen != NULL)
	{
		delete stemResidueNitrogen;
		stemResidueNitrogen = NULL;
	}
	if (rootResidueCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootResidueCarbon[j] != NULL)
			{
				delete[] rootResidueCarbon[j];
				rootResidueCarbon[j] = NULL;
			}
		}
		delete[] rootResidueCarbon;
		rootResidueCarbon = NULL;
	}
	if (rootResidueNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootResidueNitrogen[j] != NULL)
			{
				delete[] rootResidueNitrogen[j];
				rootResidueNitrogen[j] = NULL;
			}
		}
		delete[] rootResidueNitrogen;
		rootResidueNitrogen = NULL;
	}
	if (humusCarbon != NULL)
	{
		delete[] humusCarbon;
		humusCarbon = NULL;
	}
	if (humusNitrogen != NULL)
	{
		delete[] humusNitrogen;
		humusNitrogen = NULL;
	}
	if (heterotrophicRespiration != NULL)
	{
		delete[] heterotrophicRespiration;
		heterotrophicRespiration = NULL;
	}
	if (rhizosphereCl != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereCl[j] != NULL)
			{
				delete[] rhizosphereCl[j];
				rhizosphereCl[j] = NULL;
			}
		}
		delete[] rhizosphereCl;
		rhizosphereCl = NULL;
	}
	if (rhizosphereNl != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereNl[j] != NULL)
			{
				delete[] rhizosphereNl[j];
				rhizosphereNl[j] = NULL;
			}
		}
		delete[] rhizosphereNl;
		rhizosphereNl = NULL;
	}
	if (rhizosphereLiveMicrobialCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereLiveMicrobialCarbon[j] != NULL)
			{
				delete[] rhizosphereLiveMicrobialCarbon[j];
				rhizosphereLiveMicrobialCarbon[j] = NULL;
			}
		}
		delete[] rhizosphereLiveMicrobialCarbon;
		rhizosphereLiveMicrobialCarbon = NULL;
	}
	if (rhizosphereDeadMicrobialCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereDeadMicrobialCarbon[j] != NULL)
			{
				delete[] rhizosphereDeadMicrobialCarbon[j];
				rhizosphereDeadMicrobialCarbon[j] = NULL;
			}
		}
		delete[] rhizosphereDeadMicrobialCarbon;
		rhizosphereDeadMicrobialCarbon = NULL;
	}
	if (rhizosphereLabileCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereLabileCarbon[j] != NULL)
			{
				delete[] rhizosphereLabileCarbon[j];
				rhizosphereLabileCarbon[j] = NULL;
			}
		}
		delete[] rhizosphereLabileCarbon;
		rhizosphereLabileCarbon = NULL;
	}
	if (rhizosphereMicrobialNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereMicrobialNitrogen[j] != NULL)
			{
				delete[] rhizosphereMicrobialNitrogen[j];
				rhizosphereMicrobialNitrogen[j] = NULL;
			}
		}
		delete[] rhizosphereMicrobialNitrogen;
		rhizosphereMicrobialNitrogen = NULL;
	}
	if (rhizosphereMineralNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereMineralNitrogen[j] != NULL)
			{
				delete[] rhizosphereMineralNitrogen[j];
				rhizosphereMineralNitrogen[j] = NULL;
			}
		}
		delete[] rhizosphereMineralNitrogen;
		rhizosphereMineralNitrogen = NULL;
	}
	if (rhizosphereAmmoniumNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereAmmoniumNitrogen[j] != NULL)
			{
				delete[] rhizosphereAmmoniumNitrogen[j];
				rhizosphereAmmoniumNitrogen[j] = NULL;
			}
		}
		delete[] rhizosphereAmmoniumNitrogen;
		rhizosphereAmmoniumNitrogen = NULL;
	}
	if (rhizosphereNitrateNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rhizosphereNitrateNitrogen[j] != NULL)
			{
				delete[] rhizosphereNitrateNitrogen[j];
				rhizosphereNitrateNitrogen[j] = NULL;
			}
		}
		delete[] rhizosphereNitrateNitrogen;
		rhizosphereNitrateNitrogen = NULL;
	}
	if (rootExudateSugarCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootExudateSugarCarbon[j] != NULL)
			{
				delete[] rootExudateSugarCarbon[j];
				rootExudateSugarCarbon[j] = NULL;
			}
		}
		delete[] rootExudateSugarCarbon;
		rootExudateSugarCarbon = NULL;
	}
	if (rootExudateAminoAcidCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootExudateAminoAcidCarbon[j] != NULL)
			{
				delete[] rootExudateAminoAcidCarbon[j];
				rootExudateAminoAcidCarbon[j] = NULL;
			}
		}
		delete[] rootExudateAminoAcidCarbon;
		rootExudateAminoAcidCarbon = NULL;
	}
	if (rootExudateAminoAcidNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootExudateAminoAcidNitrogen[j] != NULL)
			{
				delete[] rootExudateAminoAcidNitrogen[j];
				rootExudateAminoAcidNitrogen[j] = NULL;
			}
		}
		delete[] rootExudateAminoAcidNitrogen;
		rootExudateAminoAcidNitrogen = NULL;
	}
	if (maximumRootBiomassCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (maximumRootBiomassCarbon[j] != NULL)
			{
				delete[] maximumRootBiomassCarbon[j];
				maximumRootBiomassCarbon[j] = NULL;
			}
		}
		delete[] maximumRootBiomassCarbon;
		maximumRootBiomassCarbon = NULL;
	}
	if (rootAreaScalar != NULL)
	{
		delete[] rootAreaScalar;
		rootAreaScalar = NULL;
	}
	if (fineRootBiomassCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (fineRootBiomassCarbon[j] != NULL)
			{
				delete[] fineRootBiomassCarbon[j];
				fineRootBiomassCarbon[j] = NULL;
			}
		}
		delete[] fineRootBiomassCarbon;
		fineRootBiomassCarbon = NULL;
	}
	if (fineRootBiomassNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (fineRootBiomassNitrogen[j] != NULL)
			{
				delete[] fineRootBiomassNitrogen[j];
				fineRootBiomassNitrogen[j] = NULL;
			}
		}
		delete[] fineRootBiomassNitrogen;
		fineRootBiomassNitrogen = NULL;
	}
	if (coarseRootBiomassCarbon != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (coarseRootBiomassCarbon[j] != NULL)
			{
				delete[] coarseRootBiomassCarbon[j];
				coarseRootBiomassCarbon[j] = NULL;
			}
		}
		delete[] coarseRootBiomassCarbon;
		coarseRootBiomassCarbon = NULL;
	}
	if (coarseRootBiomassNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (coarseRootBiomassNitrogen[j] != NULL)
			{
				delete[] coarseRootBiomassNitrogen[j];
				coarseRootBiomassNitrogen[j] = NULL;
			}
		}
		delete[] coarseRootBiomassNitrogen;
		coarseRootBiomassNitrogen = NULL;
	}
	if (rootNSC != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootNSC[j] != NULL)
			{
				delete[] rootNSC[j];
				rootNSC[j] = NULL;
			}
		}
		delete[] rootNSC;
		rootNSC = NULL;
	}
	if (rootMineralNitrogen != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootMineralNitrogen[j] != NULL)
			{
				delete[] rootMineralNitrogen[j];
				rootMineralNitrogen[j] = NULL;
			}
		}
		delete[] rootMineralNitrogen;
		rootMineralNitrogen = NULL;
	}
	if (rootArea != NULL)
	{
		for (int j = 0; j < nRoots; j++)
		{
			if (rootArea[j] != NULL)
			{
				delete[] rootArea[j];
				rootArea[j] = NULL;
			}
		}
		delete[] rootArea;
		rootArea = NULL;
	}
	if (liveStemCarbon != NULL)
	{
		delete liveStemCarbon;
		liveStemCarbon = NULL;
	}
	if (stemNSC != NULL)
	{
		delete stemNSC;
		stemNSC = NULL;
	}
	if (liveStemNitrogen != NULL)
	{
		delete liveStemNitrogen;
		liveStemNitrogen = NULL;
	}
	if (deadStemCarbon != NULL)
	{
		delete deadStemCarbon;
		deadStemCarbon = NULL;
	}
	if (deadStemNitrogen != NULL)
	{
		delete deadStemNitrogen;
		deadStemNitrogen = NULL;
	}
	if (leafBiomassCarbon != NULL)
	{
		delete leafBiomassCarbon;
		leafBiomassCarbon = NULL;
	}
	if (leafBiomassNitrogen != NULL)
	{
		delete leafBiomassNitrogen;
		leafBiomassNitrogen = NULL;
	}
	if (leafNSC != NULL)
	{
		delete leafNSC;
		leafNSC = NULL;
	}
	if (chloroplastStarch != NULL)
	{
		delete chloroplastStarch;
		chloroplastStarch = NULL;
	}
	if (chloroplastSugar != NULL)
	{
		delete chloroplastSugar;
		chloroplastSugar = NULL;
	}
	if (glycineNitrogen != NULL)
	{
		delete glycineNitrogen;
		glycineNitrogen = NULL;
	}
	if (glycineCarbon != NULL)
	{
		delete glycineCarbon;
		glycineCarbon = NULL;
	}
	if (serineNitrogen != NULL)
	{
		delete serineNitrogen;
		serineNitrogen = NULL;
	}
	if (serineCarbon != NULL)
	{
		delete serineCarbon;
		serineCarbon = NULL;
	}
	if (leafStoredNitrogen != NULL)
	{
		delete leafStoredNitrogen;
		leafStoredNitrogen = NULL;
	}
	if (leafRubiscoNitrogen != NULL)
	{
		delete leafRubiscoNitrogen;
		leafRubiscoNitrogen = NULL;
	}
	if (fruitCarbon != NULL)
	{
		delete fruitCarbon;
		fruitCarbon = NULL;
	}
	if (fruitNitrogen != NULL)
	{
		delete fruitNitrogen;
		fruitNitrogen = NULL;
	}
	if (plantNstatus != NULL)
	{
		delete plantNstatus;
		plantNstatus = NULL;
	}
	if (nitrogenLeaching != NULL)
	{
		delete nitrogenLeaching;
		nitrogenLeaching = NULL;
	}
	if (lat_Root_b_value_init != NULL)
	{
		delete lat_Root_b_value_init;
		lat_Root_b_value_init = NULL;
	}
	if (lat_Root_c_value_init != NULL)
	{
		delete lat_Root_c_value_init;
		lat_Root_c_value_init = NULL;
	}
}

//
//C4photosynthesis()
//
//August 2019 D.S. Mackay
//
//This functions leaf photosynthesis for C4 plants
//	Inputs are: 
//		Leaf surface CO2 partial pressure (Ca), ubar 
//		Incoming photosynthetically active radiation (Iin), umol m-2 s-1
//		Leaf temperature (Tl), degrees C
//		Stomatal conductance to CO2 (gc), mol m-2 s-1
//	Outputs are:
//		Intercellular CO2 partial pressure (Ci), ubar
//		Mesophyll CO2 partial pressure (Cm), ubar
//		C4 net photosynthesis (An), umol m-2 s-1
//
//	Reference: von Caemmerer (2013) Plant, Cell & Environment
//
double BiogeochemicalCycles::C4photosynthesis(trees_params treesParams,
						double Ca,
						double Iin,
						double phiJ,
						double Tl,
						double gc,
						double& Ci,
						double& Cm)
{
	double Vcmax25, Vpmax25, Jmax25, thetaJ, Kp25, Vpr, f, x, absorptance, I2, Rd, Rm, Rs;
	double E_Vcmax, E_Vpmax, E_Jmax, E_Kp;
	double gm, gt, gbs, Tlk, Vcmax, Vpmax, Kp, Jmax;
	double Aref, An, Ac, Aj, Jt, Vp;
	double Cm_max, Cm_min;
	bool stop;

	Vcmax25 = treesParams.Vcmax25; //maximum Rubisco activity at 25 C, umol m-2 s-1
	Vpmax25 = treesParams.Vpmax25; //maximum PEP carbolylase activity at 25 C, umol m-2 s-1
	Jmax25 = treesParams.Jmax25;   //maximum electron transport rate at 25 C, umol m-2 s-1
	thetaJ = treesParams.thetaJ;   //unitless curvature parameter
	Kp25 = treesParams.Kp25;       //Michaelis constant of PEP carboxylase for CO2 at 25 C, ubar
	Vpr = treesParams.Vpr;         //PEP regeneration rate, umol m-2 s-1
	f = treesParams.f;              //correction for spectral quality of light	
	x = treesParams.x;             //partitioning factor of electron transport rate
	absorptance = treesParams.absorptance; //fraction of irradiance absorbed
	//I2 = Iin * absorptance * (1.0 - f) / 2.0; //PAR irradiance that is useful for photosynthesis
	I2 = Iin * absorptance * (1.0 - f) * phiJ; //PAR irradiance that is useful for photosynthesis
	E_Vcmax = treesParams.E_Vcmax; //activation energy, maximum carboxylation rate, kJ mol-1
	E_Vpmax = treesParams.E_Vpmax; //activation energy, maximum PEP rate, kJ mol-1
	E_Jmax = treesParams.E_Jmax;   //activation energy, electron transport, kJ mol-1
	E_Kp = treesParams.E_Kp;       //activation energy, Michaelis reaction of PEP, kJ mol-1
	gm = treesParams.gm;	       //mesophyll conductance to CO2, mol m-2 s-1
	gt = 1.0 / (1.0/gc + 1.0/gm);  //diffusive conductance of CO2 from mesophyll to leaf surface, mol m-2 s-1
	gbs = treesParams.gbs;         //conductance of the bundle sheath, mol m-2 s-1

	Tlk = Tl + 273.1;              //absolute temperature
	Vcmax = Vcmax25 * exp((Tl-25.0)*E_Vcmax*1000.0/(298.1*R*Tlk));
	Vpmax = Vpmax25 * exp((Tl-25.0)*E_Vpmax*1000.0/(298.1*R*Tlk));
	Kp = Kp25 * exp((Tl-25.0)*E_Kp*1000.0/(298.1*R*Tlk));
	Jmax = Jmax25 * exp((Tl-25.0)*E_Jmax*1000.0/(298.1*R*Tlk));

	Rd = treesParams.Rd_mult * Vcmax; //leaf mitochondrial respiration, umol m-2 s-1
	Rm = 0.5 * Rd;                 //mesophyll mitochondrial respiration, umol m-2 s-1
	Rs = Rd - Rm;                  //bundle sheath mitochondrial respiration, umol m-2 s-1

	stop = false;
	Cm_max = Ca;
	Cm_min = 0.0;

//Iterative solution of coupled Cm, Ci, and A
//Correction solution is found when gas exchange Aref equals the C4 biochemical A
//Aref declines as Cm increaes, while C4 A increases with Cm, and so there is an
//   equilibrium point, which is found efficiently with a binary search

	while (stop == false)
	{
		Cm = 0.5 * (Cm_max + Cm_min); //using the mid-point between the current limits

//Note: from leaf surface to intercellular, A = gc(Ca - Ci)
//      and from intercellular to mesophyll, A = gm(Ci - Cm), and so
//      from leaf surface to mesophyll, A = gt(Ca - Cm),
//      where gt = 1/(1/gc + 1/gm)

//(1) gas exchange solution for photosynthesis
		Aref = gt * (Ca - Cm);

//(2) biochemical solution for photosynthesis, after von Caemmerer (2013)
//(2a) Enzyme-limited photosynthesis
		Vp = min(Cm*Vpmax/(Cm+Kp), Vpr);
		Ac = min((Vp - gbs*Cm - Rm), (Vcmax - Rd));

//(2b) Transport-limited photosynthesis
		Jt = (I2 + Jmax - sqrt(pow(I2+Jmax,2.0) - 4.0 * thetaJ * I2 * Jmax))/(2.0*thetaJ);
		Aj = min((x*Jt/2.0 - Rm + gbs*Cm), ((1.0-x)*Jt/3.0 - Rd));

//Photosynthesis is the minimum of Ac and Aj
		An = min(Ac, Aj);

//Back out Ci as a function of A and Cm
		Ci = Cm + An/gm;

		if (Aref > An)
		{
			Cm_min = Cm; //solution must not be lower than current Cm
		}
		else
		{
			Cm_max = Cm; //solution must not be higher than current Cm
		}

		if (((An/Aref) > 0.9999 && (An/Aref) < 1.0001) || An < 0.0)
		{
			stop = true; //solution has been found
		}
	}
	return An;
}

//
//coupledA4_gc()
//
//August 2019 D.S. Mackay
//
//Solves the coupled C4 photosynthesis & stomatal conductance
//
//This algorithm finds the lowest gc that maximizes An. It starts with an initial gc, which can obtained
//  from the soil-xylem hydraulics, and finds the smallest gc for the reference An, which is either
//  enzyme- or transport-limited. This function calls C4photosynthesis() to compute An, Ci, and Cm.
//  The algorithm employs a binary search for efficiency.
//
//Why does this matter? At low light, An is transport-limited and so gc can be reduced from its maximum
//  without lowering An. At high atmospheric CO2 concentration, An may be enzyme-limited, and
//  so gc can be reduced without lowering An.
//
void BiogeochemicalCycles::coupledA4_gc(trees_params treesParams,
                                          double Ca,
                                          double Iin,
					  double phiJ,
                                          double Tl,
					  double gc0,
                                          double& gc,
					  double& An,
                                          double& Ci,
                                          double& Cm)
{
	bool stop;
	double gc_min, gc_max, A_ref;

	stop = false;
	gc_min = 0.0001;
	gc_max = gc0;

//compute reference photosynthesis at the initial, maximum, stomatal conductance
	A_ref = C4photosynthesis(treesParams, Ca, Iin, phiJ, Tl, gc_max, Ci, Cm);
	if (A_ref <= 0.0)
	{
		An = A_ref;
		gc = gc_min;
		stop = true;
	}

//find the lowest gc that gives the highest An
	while (stop == false)
	{
		gc = 0.5 * (gc_min + gc_max);
		An = C4photosynthesis(treesParams, Ca, Iin, phiJ, Tl, gc, Ci, Cm);
		if (An < A_ref)
		{
			gc_min = gc;  //solution cannot be lower than current gc
		}
		else
		{
			gc_max = gc;  //solution cannot be higher than current gc
		}
		if ((gc_min/gc_max) > 0.99999 && (gc_min/gc_max) < 1.00001)
		{
			stop = true; // solution found
		}
	}
}
	
//
//C3photosynthesis()
//
//(c) August 2019, July 2020 D.S. Mackay
//
//This functions leaf photosynthesis for C3 plants
//	Inputs are: 
//		Leaf surface CO2 partial pressure (Ca), ubar 
//		Incoming photosynthetically active radiation (Iin), umol m-2 s-1
//		Leaf temperature (Tl), degrees C
//		Stomatal conductance to CO2 (gc), mol m-2 s-1
//		de novo nitrogen NO3 available (must be determined at host model end)
//	Outputs are:
//		Intercellular CO2 partial pressure (Ci), ubar
//		Chloroplast CO2 partial pressure (Cc), ubar
//		C3 net photosynthesis (An), umol m-2 s-1
//		Serine and glycine nitrogen assimilated (host model would need to do something with these)
//
//	Reference: von Caemmerer (2013) Plant, Cell & Environment
//		   Busch et al (2017) Nature Plants
//
//This modified C3 photosynthesis model can assimilate transiently available nitrogen into organic
//  compounds (serine and glycine). It utilizes the photorespiratory pathway, and can take up
//  additional CO2 under both TPU and Rubisco limitation, as well as assimilate large amounts of N.
//
double BiogeochemicalCycles::C3photosynthesis(trees_params treesParams,
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
						double& NS)
{
	double Vcmax25, Jmax25, thetaJ, Kc25, Ko25, f, x, absorptance, I2, Rd25, Rd, Rm;
	double E_Vcmax, E_Jmax, E_Kc, E_Ko, E_Rd, E_gammaStar;
	double gm, gt, Tlk, Vcmax, O, Jmax, gammaStar, gammaStar25;
	double Aref, An, Ac, Aj, Ap, Jt, Tp, Kc, Ko;
	double Cc_max, Cc_min;
	double alphaGmax, alphaSmax, beta, PHI, alphaGc, alphaSc, alphaGj, alphaSj, alphaGp, alphaSp, gammaStarG;
	double Voc, Voj, Vop;
	bool stop;

	O = 0.21 * Pa;                 //O2 partial pressure
	Vcmax25 = treesParams.Vcmax25; //maximum Rubisco activity at 25 C, umol m-2 s-1
	Jmax25 = treesParams.Jmax25;   //maximum electron transport rate at 25 C, umol m-2 s-1
	thetaJ = treesParams.thetaJ;   //unitless curvature parameter
	Kc25 = 10.0*treesParams.Kc25;  //Michaelis-Menten constant of Rubisco carboxylation at 25 C, ubar
	Ko25 = 0.01*treesParams.Ko25;  //Michaelis-Menten constant of Ribusco oxygenation at 25 C, mbar
	gammaStar25 = treesParams.gammaStar25; // , ubar
	f = treesParams.f;              //correction for spectral quality of light	
	x = treesParams.x;             //partitioning factor of electron transport rate
	absorptance = treesParams.absorptance; //fraction of irradiance absorbed
	//I2 = Iin * absorptance * (1.0 - f) / 2.0; //PAR irradiance that is useful for photosynthesis
	I2 = Iin * absorptance * (1.0 - f) * phiJ; //PAR irradiance that is useful for photosynthesis
	Rd25 = treesParams.Rd_mult * Vcmax25; //leaf mitochondrial respiration, umol m-2 s-1
	E_Vcmax = treesParams.E_Vcmax; //activation energy, maximum carboxylation rate, kJ mol-1
	E_Jmax = treesParams.E_Jmax;   //activation energy, electron transport, kJ mol-1
	E_Kc = treesParams.E_Kc;       //activation energy, M-M carboxylation, kJ mol-1
	E_Ko = treesParams.E_Ko;       //activation energy, M-M oxygenation, kJ mol-1
	E_Rd = treesParams.E_Rd;       //activation energy for leaf mitochondrial respiration, kJ mol-1
	E_gammaStar = treesParams.E_gammaStar; //activation energy, CO2 compensation point, kJ mol-1
	gm = treesParams.gm;	       //mesophyll conductance to CO2, mol m-2 s-1
	gt = 1.0 / (1.0/gc + 1.0/gm);  //diffusive conductance of CO2 from mesophyll to leaf surface, mol m-2 s-1
	alphaGmax = treesParams.alphaGmax; //fraction of glycolate carbon diverted to glycine during photorespiration
	alphaSmax = treesParams.alphaSmax; //fraction of glycolate carbon diverted to serine during photorespiration

	Tlk = Tl + 273.1;              //absolute temperature
	Vcmax = Vcmax25 * exp((Tlk-298.0)*E_Vcmax*1000.0/(298.0*R*Tlk));
	Kc = Kc25 * exp((Tlk-298.0)*E_Kc*1000.0/(298.0*R*Tlk));
	Ko = Ko25 * exp((Tlk-298.0)*E_Ko*1000.0/(298.0*R*Tlk));
	Rd = Rd25 * exp((Tlk-298.0)*E_Rd*1000.0/(298.0*R*Tlk));
	Rm = 0.5 * Rd;                 //mesophyll mitochondrial respiration, umol m-2 s-1
	gammaStar = gammaStar25 * exp((Tlk-298.0)*E_gammaStar*1000.0/(298.0*R*Tlk));
	Jmax = Jmax25 * exp((Tlk-298.0)*E_Jmax*1000.0/(298.0*R*Tlk));
	Tp = 0.167 * Vcmax; //Triosphosphate limitation constant, umol m-2 s-1

	stop = false;
	Cc_max = Ca;
	Cc_min = 0.0;

//Iterative solution of coupled Cc, Ci, and An
//Correction solution is found when gas exchange Aref equals the C3 biochemical An
//Aref declines as Cc increaes, while C3 An increases with Cc, and so there is an
//   equilibrium point, which is found efficiently with a binary search

	while (stop == false)
	{
		Cc = 0.5 * (Cc_max + Cc_min); //using the mid-point between the current limits

//Note: from leaf surface to intercellular, A = gc(Ca - Ci)
//      and from intercellular to mesophyll, A = gm(Ci - Cc), and so
//      from leaf surface to mesophyll, A = gt(Ca - Cc),
//      where gt = 1/(1/gc + 1/gm)

//(1) gas exchange solution for photosynthesis
		Aref = gt * (Ca - Cc);

//(2) Biochemical solution for C & N assimilation, after von Caemmerer (2013) PC&E and Busch et al (2017) Nature Plants
//
		beta = 0.0;
		if (alphaGmax > 0.0)
		{
			beta = 3.0*alphaGmax / (3.0*alphaGmax + 2.0*alphaSmax);
		}
		PHI = 2.0*gammaStar/Cc;
		Voc = PHI * Cc * Vcmax / (Cc + Kc * (1.0 + O/Ko));
//(2a) Enzyme-limited photosynthesis with N assimilation during photorespiration
		alphaGc = min(alphaGmax, Nmax*beta/Voc);
		alphaSc = min(alphaSmax, 3.0*Nmax*(1.0-beta)/(2.0*Voc));
//Note: enzyme-limited photosynthesis does not change with addition of N assimilation
		Ac = ((Cc - gammaStar)*Vcmax)/(Cc + Kc*(1.0+O/Ko)) - Rd;

//(2b) Transport-limited photosynthesis with N assimilation during photorespiration
		Jt = (I2 + Jmax - sqrt(pow(I2+Jmax,2.0) - 4.0 * thetaJ * I2 * Jmax))/(2.0*thetaJ);
		if (Jt > Nmax * (2.0*beta+6.0))
		{
			alphaGj = min(alphaGmax, 4.0*Nmax*beta*(1.0/PHI + 1)/(Jt - Nmax*(2.0*beta+6.0)));
			alphaSj = min(alphaSmax, 6.0*Nmax*(1.0-beta)*(1.0/PHI + 1)/(Jt - Nmax*(2.0*beta+6.0)));
		}
		else
		{
			alphaGj = alphaGmax;
			alphaSj = alphaSmax;
		}
		gammaStarG = gammaStar * (1.0 - alphaGj);
		Aj = (Cc - gammaStarG) * Jt / (4.0*Cc + 8.0*gammaStar + (1.0 + 2.0*alphaGj + alphaSj)) - Rd;

//(2c) Triosephosphate-limited photosynthesis
		alphaGp = min(alphaGmax, Nmax*beta*(2.0/PHI - 1.0)/(6.0*Tp + 3.0*Nmax*(2.0-beta)));
		alphaSp = min(alphaSmax, (3.0/2.0)*Nmax*(1.0-beta)*(2.0/PHI-1.0)/(6.0*Tp + 3.0*Nmax*(2.0-beta)));
		gammaStarG = gammaStar * (1.0 - alphaGp);
		Ap = 3.0*Tp*(Cc - gammaStarG)/(Cc - gammaStarG * (1.0 + 3.0*alphaGp + 4.0*alphaSp)) - Rd;

//Photosynthesis is the minimum of Ac, Aj, and Ap
		An = min(Ac, Aj);
		An = min(An, Ap);

//Compute export of glycine and serine, umol m-2 s-1
		if (An == Ac)
		{
			NG = alphaGc * Voc;
			NS = 2.0/3.0 * alphaSc * Voc;
		}
		else if (An == Aj)
		{
			Voj = Jt / (4.0/PHI + (4.0 + 8.0*alphaGj + 4.0*alphaSj));
			NG = alphaGj * Voj;
			NS = 2.0/3.0 * alphaSj * Voj;
		}
		else
		{
			Vop = 3.0 * Tp / (1.0/PHI - 0.5 * (1.0 + 3.0 * alphaGp + 4.0 * alphaSp));
			NG = alphaGp * Vop;
			NS = 2.0/3.0 * alphaSp * Vop;
		}

//Back out Ci as a function of A and Cm
		Ci = Cc + An/gm;

		if (Aref > An)
		{
			Cc_min = Cc; //solution must not be lower than current Cm
		}
		else
		{
			Cc_max = Cc; //solution must not be higher than current Cm
		}

		if (((An/Aref) > 0.9999 && (An/Aref) < 1.0001) || An < 0.0)
		{
			stop = true; //solution has been found
		}
	}
	return An;
}

//
//coupledA3_gc()
//
//(c) August 2019, July 202 D.S. Mackay
//
//Solves the coupled C3 photosynthesis & stomatal conductance
//
//This algorithm finds the lowest gc that maximizes An. It starts with an initial gc, which can obtained
//  from the soil-xylem hydraulics, and finds the smallest gc for the reference An, which is either
//  enzyme- or transport-limited. This function calls C3photosynthesis() to compute An, Ci, and Cm.
//  The algorithm employs a binary search for efficiency.
//
//Why does this matter? At low light, An is transport-limited and so gc can be reduced from its maximum
//  without lowering An. At high atmospheric CO2 concentration, An may be enzyme-limited, and
//  so gc can be reduced without lowering An.
//
//Limitations: C3 An is asymptic and so an arbitrary threshold of 99% of An is used to find the 
//             stopping point for computing gc
//
void BiogeochemicalCycles::coupledA3_gc(trees_params treesParams,
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
					  double& NS)
{
	bool stop;
	double gc_min, gc_max, gc_nolimit, A_ref, dAdg;

	stop = false;
//gc = stomatal conductance, mol m-2 s-1
	gc_min = 0.00125;
	gc_max = gc0;
	gc_nolimit = 1.01*gc0;

//compute reference photosynthesis, A_ref, for no stomatal conductance limitation
	A_ref = C3photosynthesis(treesParams, Ca, Pa, Iin, phiJ, Tl, gc_nolimit, Nmax, Ci, Cc, NG, NS);
	if (A_ref <= 0.0)
	{
		An = A_ref;
		gc = gc_min;
		stop = true;
	}
//find the lowest gc that gives the highest An
	while (stop == false)
	{
		gc = 0.5 * (gc_min + gc_max);
		An = C3photosynthesis(treesParams, Ca, Pa, Iin, phiJ, Tl, gc, Nmax, Ci, Cc, NG, NS);
		if (An < 0.99*A_ref)
		{
			gc_min = gc;  //solution cannot be lower than current gc
		}
		else
		{
			gc_max = gc;  //solution cannot be higher than current gc
		}
		if ((gc_min/gc_max) > 0.9999 && (gc_min/gc_max) < 1.0001)
		{
			stop = true; // solution found
		}
	}
}

//
//storeGlycineAndSerine()
//
//Updates to glycine and serine exported from C3 photosynthesis
//
// N exported in umol N m-2 leaf s-1 --> kgN ha-1
// kgN ha-1 30min-1 = (umol N m-2 leaf s-1 ) * (m2 ha-1) * (m2 leaf m-2 ground) * (kg N g-1 N) * (gN mol-1) * (mol umol-1) * (s 30min-1)
//                  = N exported * 10000 * L * 0.001 * 14.0067 * 1e-6 * 1800
//                  = N exported * L * 0.2521206
//
// Glycine: C2H5NO2
// 	kgC ha-1 30min-1 = N exported * 2 * 10000 * L * 0.001 * 12.0107 * 1e-6 * 1800
//                       = N exported * L * 0.4323852
// Serine: C3H7NO3
//	kgC ha-1 30-min01 = N exported * 3 * 10000 * L * 0.001 * 12.0107 * 1e-6 * 1800
//			  = N exported * L * 0.6485778
//

void BiogeochemicalCycles::storeGlycineAndSerine(double NGsun, 
						 double NGshd, 
						 double NSsun, 
						 double NSshd, 
						 double Lsun, 
						 double Lshd)
{
	glycineNitrogen[0] += (NGsun * Lsun + NGshd * Lshd)* 0.2521206;
	glycineCarbon[0] += (NGsun * Lsun + NGshd * Lshd) * 0.4323852;
	serineNitrogen[0] += (NSsun * Lsun + NSshd * Lshd) * 0.2521206;
	serineCarbon[0] += (NSsun * Lsun + NSshd * Lshd) * 0.6485778;
//cout << glycineNitrogen[0] << '\t' << glycineCarbon[0] << '\t' << serineNitrogen[0] << '\t' << serineCarbon[0] << endl;
}
	

//
// photosynthesis()
//     This function computes leaf photosynthesis for C3 plants,
//	returning also stomatal conductance
//
//     August 2017 DSM - Modification to allow noctural stomatal conductance
//
double BiogeochemicalCycles::photosynthesis(trees_params treesParams,
						double Jmax_mult, 
						double thetaJ, 
						double phiJ,
                                                struct farqin in,
						struct farqout& out)
{
	double t;      			// (deg C) temperature
	double tk;     			// (K) absolute temperature
	double g, gin;      		// (umol/m2/s/Pa) conductance to CO2
	double O2;     			// (Pa) atmospheric partial pressure O2
	double Ca;     			// (Pa) atmospheric partial pressure CO2
	double gammaStar;  		// (Pa) co2 compensation point, no dark respiration
	double Kc25, Kc;     		// (Pa) MM constant for carboxylase reaction
	double q10Kc;
	double Ko25, Ko;     		// (Pa) MM constant for oxygenase reaction
	double q10Ko;
	double act25, act;    		// (umol/kgRubisco/s) Rubisco activity
	double q10act;
	double Rd;     			// (umol/m2/s) dark respiration rate
	double Vcmax, Vcmax25;   	// (umol/m2/s) maximum carboxylation velocity
	double Jmax, Jmax25;   		// (umol/m2/s) maximum rate of electron transport
	double J;      			// (umol/m2/s) maximum rate of Rubisco regeneration
	double Av;     			// (umol/m2/s) Rubisco limited assimilation rate
	double Aj;     			// (umol/m2/s) RuBP regeneration limited assim rate
	double As;      		// (umol/m2/s) Sink-limited assimilation rate
	double A;      			// (umol/m2/s) net assimilation rate
	double A_kg;			// (kgC ha-1)
	//double Ile;   		// (umol/m2/s) PAR effectively absorbed by PSII per unit leaf area
	double Ji;      		// ((umol/m2/s) rate of electron transport per unit absorbed PAR
	double absorptance, spectralQuality;
	double aa,bb,cc,det;
	double kappa1, kappa2, kappa3;	//constants used in Katul et al, 2003, Plant, Cell and Env
	double alpha1, alpha2;   	//constants used in Katul et al, 2003, Plant, Cell and Env
	double Ci, D;
	double E;       		//transpiration used to check validity of Ci
	double g_mult;
	double Rd_mult = treesParams.Rd_mult;

/* local variables  */

        t = in.t;               	//temperature should be leaf temperature
        tk = t + 273.15;        	//leaf temperature (deg. K)
        g = gin = in.g*1.0e6/(R*tk);    //convert conductance from m/s -. umol/m2/s/Pa
        if ( g < 0.00000001 )
        {
                g = 0.00000001;
        }
        Ca = in.co2 * 1.0e-6 * in.pa;   //convert atmospheric CO2 from ppm to Pa
        O2 = 0.2095 * in.pa; 		//atmospheric O2 in Pa, assumes 20.95% O2 by volume

// correct kinetic constants for temperature, and do unit conversions
//Future work - replace with Arrhenius kinetics
        Kc25 = treesParams.Kc25;
        q10Kc = treesParams.q10Kc;
        Ko25 = treesParams.Ko25;
        q10Ko = treesParams.q10Ko;
        act25 = treesParams.act25;
        q10act = treesParams.q10act;
        Ko = Ko25 * pow(q10Ko, (t-25.0)/10.0);
        if (t > 15.0)
        {
                Kc = Kc25 * pow(q10Kc, (t-25.0)/10.0);
                act = act25 * pow(q10act, (t-25.0)/10.0);
        }
        else
        {
                Kc = Kc25 * pow(1.8*q10Kc, (t-15.0)/10.0) / q10Kc;
                act = act25 * pow(1.8*q10act, (t-15.0)/10.0) / q10act;
        }

        act = act * 1.0e6 / 60.0;     // umol/mg/min to umol/kg/s

// calculate gammaStar (Pa) - DePury and Farquhar, 1997
        gammaStar = 3.69+0.188*(t-25.0)+0.0036*(t-25.0)*(t-25.0);

// calculate Vcmax from leaf nitrogen data and Rubisco activity

/*          kg Nleaf   kg NRub    kg Rub      umol            umol
           -------- X -------  X ------- X ---------   =   --------
              m2      kg Nleaf   kg NRub   kg Rub * s       m2 * s

             (lnc)  X  (flnr)  X  (fnr)  X   (act)     =    (Vcmax)
*/

        Vcmax = in.lnc * in.flnr * fnr * act;

/* Leaf respiration not including photorespiration calc. as
           Rd = 0.0089Vcmax
           from Leuning et al.  1995     */

        Rd = Rd_mult * Vcmax;
        out.Rd = Rd;

/* calculate Jmax = f(Vcmax), reference:
        Wullschleger, S.D., 1993.  Biochemical limitations to carbon assimilation
                in C3 plants - A retrospective analysis of the A/Ci curves from
                109 species. Journal of Experimental Botany, 44:907:920.
*/

//calculate J = f(Jmax, ppfd), reference: de Pury and Farquhar 1997 Plant Cell and Env.
//Jmax proportion to Vcmax at the reference temperature, 25 deg. C (JUNE 2008 DSM)
        Vcmax25 = in.lnc * in.flnr * fnr * act25 * 1.0e6/60.0;
        Jmax25 = Jmax_mult*Vcmax25;

        double Ea = 37000.0;   //Activation energy, kJ mol-1
        double S = 710.0;      //Electrong transport temperature response parameter, J K-1 mol-1
        double H = 220000.0;   //Electron transport temperature curvature parameter, J mol-1
        Jmax = Jmax25*exp((tk-298.0)*Ea/(R*tk*298.0))*
				(1.0+exp((S*298.0-H)/(R*298.0)))/(1.0+exp((S*tk-H)/(R*tk)));

//phiJ is effective quantum yield, mol electrons mol-1 photons
//assumes phiJ = quantum yield (mol C mol-1 photons) * 4 e- per C * leaf absorptance
//phiJ can be varied between sun and shade leaves to take into consideration direct versus diffuse rad'n
//therefore, Ile no longer used, but embedded in phiJ
        if (in.irad > 2200.0)
        {
                in.irad = 2200.0;
        }
        if (in.irad < 0.0)
        {
                in.irad = 0.0;
        }
	absorptance = 0.85;
	spectralQuality = 0.85;
        //Ji = in.irad * phiJ;
	Ji = in.irad * absorptance * spectralQuality * phiJ;

        aa = thetaJ;
        bb = -Ji -Jmax;
        cc = Ji*Jmax;
        J = (-bb - sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa);

/* solve for Av and Aj using the quadratic equation, substitution for Ci
        Farquhar, G.D., and S. von Caemmerer, 1982.  Modelling of photosynthetic
                response to environmental conditions.  In Encyclopedia of Plant
                Physiology, New Series, Vol. 12B, Physiological Plant Ecology II,
                O.L. Lange, P.S. Nobel, C.B. Osmond, and H. Ziegler, eds, Springer-
                Verlag, Berlin, Germany, pp 549-587.

        from A = g(Ca-Ci) into the equations from Farquhar and von Caemmerer:

               Vcmax (Ci - gammaStar)
        Av =  -------------------   -   Rd
              Ci + Kc (1 + O2/Ko)

        Use Aj equation from dePury and Farquhar, 1997

                 J (Ci - gammaStaR)
        Aj  =  -------------------  -   Rd
              4(Ci +  2*gammaStar)

   	and As as sink-limited photosynthesis from accumulation of starch at the chloroplast
                (Bonan 2008, page 246)

        As = Vxmax / 2

        */

// quadratic solution for Av
        aa = -1.0/g;
        bb = Ca + (Vcmax - Rd)/g + Kc*(1.0 + O2/Ko);
        cc = Vcmax*(gammaStar - Ca) + Rd*(Ca + Kc*(1.0 + O2/Ko));

        if ((det = bb*bb - 4.0*aa*cc) < 0.0)
        {
                Av = 0.0;
        }
        else
        {
                Av = (-bb + sqrt(det)) / (2.0*aa);
        }

// quadratic solution for Aj
        aa = -4.0/g;
        bb = 4.0*Ca + 8.0*gammaStar + J/g - 4.0*Rd/g;
        cc = J*(gammaStar - Ca) + Rd*(4.0*Ca + 8.0*gammaStar);

        if ((det = bb*bb - 4.0*aa*cc) < 0.0)
        {
                Aj = 0.0;
        }
        else
        {
                Aj = (-bb + sqrt(det)) / (2.0*aa);
        }

// sink-limited As
        As = Vcmax / 2.0;

//determine Rubisco activity or electron transport limited
        if (Av < Aj)
        {
                alpha1 = Vcmax;
                alpha2 = Kc*(1.0+O2/Ko);
                g_mult = g;
        }
        else
        {
                alpha1 = J;
                alpha2 = 2.0*gammaStar;
                g_mult = 4.0*g;
        }

//implementation of quadratic solution of Ci
//from Katul et al, 2003, Plant, Cell & Environment
//DSM June 7, 2005 - Added Rd
        kappa1 = -g_mult;
        kappa2 = g_mult * (Ca - alpha2) - alpha1;
        kappa3 = g_mult*alpha2*Ca + alpha1*gammaStar;

        Ci = (-kappa2 - sqrt(kappa2*kappa2 - 4.0*kappa1*kappa3))/(2.0*kappa1);

        if (Ci > (Ca+0.000001))
        {
                Ci = (-kappa2 + sqrt(kappa2*kappa2 - 4.0*kappa1*kappa3))/(2.0*kappa1);
        }
        if (Ci < gammaStar)
        {
                Ci = gammaStar;
        }

//Calculate conductance using most limiting of Rubisco-limited or
//electron-transport limited assimilation, and Ci from Katul

        if (Av > Aj)
        {
                A = Aj;
        }
        else
        {
                A = Av;
        }

        if (A > As)
        {
                A = As;
        }

	A_kg = A / 4.6296; //convert from umol m-2 s-1 to kg ha-1


//Positive photosynthesis
        if (A >= 0.0)
        {
                g = (A)/(Ca-Ci+0.01);
                out.g = g * (R * tk) / 1.0e6;
        }
//At compensation point or lower, hold stomata open
        else
        {
		if (getChloroplastStarch() > (-A_kg))
		{
			putChloroplastStarch(getChloroplastStarch()-(-A_kg));
		}
		else
		{
			Rd = Rd * getChloroplastStarch()/(-A_kg);
			A = -getChloroplastStarch()*4.6296;
			putChloroplastStarch(0.0);
		}
                g = (-A)/(2.0*gammaStar);
                out.g = g * (R * tk) / 1.0e6;
        }
//Cuticular conductance
        if (g < 0.0001/1.6/42.0)
        {
                g = 0.0001/1.6/42.0;
                out.g = g;
        }
        D = in.D;
        E = 1.6*g*D;

        out.E = E;
        out.Ca = Ca; 			//(Pa) atmospheric [CO2]
        out.Ci = Ci*1.0e6/in.pa;	//(Pa) intercellular [CO2]
        out.gammaStar = gammaStar;	//(Pa) CO2 compensation point, no Rd
        out.O2 = O2;			//(Pa) atmospheric [O2]
        out.Kc = Kc;			//(Pa) MM constant carboxylation
        out.Ko = Ko;			//(Pa) MM constant oxygenation
        out.act = act;			//(umol/kg/s) Rubisco activity
        out.Vcmax25 = Vcmax25;		//(umol/m2/s) max rate carboxylation at 25C
        out.Vcmax = Vcmax;		//(umol/m2/s) max rate carboxylation
        out.Jmax25 = Jmax25;		//(umol/m2/s) max rate electron transport at 25C
        out.Jmax = Jmax;		//(umol/m2/s) max rate electron transport
        out.J = J;			//(umol/m2/s) rate of RuBP regeneration
        out.Av = Av;			//(umol/m2/s) carboxylation limited assimilation
        out.Aj = Aj;			//(umol/m2/s) RuBP regen limited assimilation

        out.A = A;			//(umol/m2/s) final assimilation rate
        return out.A;
} //end photosynthesis function


//
//computeRootArea()
//calculate absorbing root area based on the fine root carbon
//Limited data suggests a SRL average in the first two root classes of about
//    30 m gC-1 in conifer and about twice as much in deciduous
//    e.g., Withington et al 2006 Ecological Monographs
//		Meinen et al 2009 Oecologia
//
double BiogeochemicalCycles::computeRootArea(trees_params treesParams)
{
        double totalRootArea;

        totalRootArea = 0.0;

        for (int j = 0; j < nRoots; j++)
        {
		totalRootArea += computeRootArea(treesParams, j);
        }
	return(totalRootArea);
}
double BiogeochemicalCycles::computeRootArea(trees_params treesParams,
					     int j)
{
        double totalRootArea, refRootDiam, rootDiam, SRL, SRA;

	assert(j >= 0);
	assert(j < nRoots);

        totalRootArea = 0.0;
	refRootDiam = 0.00025; //m
	rootDiam = treesParams.minRootDiam; //m

//root segments are approximated as cylinders
//specific root length (m kgC-1 root)
//SRL1 is defined at a root diameter of 0.00025 m
        SRL = 1000.0*treesParams.SRL1*pow(refRootDiam, 2.0)/pow(rootDiam,2.0);
        for (int k = 0; k < nRootOrders; k++)
        {
//specific root area = m2 root m-2 ground area kgC-1 root
        	SRA = SRL * rootDiam * M_PI;
		totalRootArea += computeRootArea(treesParams, SRA, j, k);
        	rootDiam *= treesParams.rootDiamMultiplier;
        	SRL /= pow(treesParams.rootDiamMultiplier, 2.0);
        }
	return(totalRootArea);
}
double BiogeochemicalCycles::computeRootArea(trees_params treesParams,
					     double SRA,
					     int j,
					     int k)
{
        double rootCarbon, totalRootArea;

	assert(j >= 0);
	assert(j < nRootOrders);

	rootCarbon = fineRootBiomassCarbon[j][k] + coarseRootBiomassCarbon[j][k];

//convert kgC ha-1 to kgC m-2
//m2 root area m-2 ground area = 
//      kgC root m-2 ground area * m-2 ground area * m2 root m-2 ground area kgC-1 root
//
	rootArea[j][k] = 0.0000001;
	if (treesParams.drlat[j+3] > 0.0001)
	{
		rootArea[j][k] = 0.0001 * rootCarbon * SRA;
	}
	
        totalRootArea = rootArea[j][k];

	return(totalRootArea);
}

//
//computeFineRootArea()
//calculate absorbing root area based on the fine root carbon only
//Limited data suggests a SRL average in the first two root classes of about
//    30 m gC-1 in conifer and about twice as much in deciduous
//    e.g., Withington et al 2006 Ecological Monographs
//		Meinen et al 2009 Oecologia
//
double BiogeochemicalCycles::computeFineRootArea(trees_params treesParams)
{
        double totalRootArea;

        totalRootArea = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
		totalRootArea += computeFineRootArea(treesParams, j);
        }
	return(totalRootArea);
}
double BiogeochemicalCycles::computeFineRootArea(trees_params treesParams,
						 int j)
{
        double totalRootArea, SRL, SRA, rootDiam, refRootDiam;

	assert(j >= 0);
	assert(j < nRoots);

        totalRootArea = 0.0;
	refRootDiam = 0.00025; //m
	rootDiam = treesParams.minRootDiam; //m

//root segments are approximated as cylinders
//specific root length (m kgC-1 root)
//SRL1 is defined at a root diameter of 0.00025 m
        SRL = 1000.0*treesParams.SRL1*pow(refRootDiam, 2.0)/pow(rootDiam,2.0);

        for (int k = 0; k < nFineRootOrders; k++)
        {
        	SRA = SRL * rootDiam * M_PI;
		totalRootArea += computeFineRootArea(treesParams, SRA, j, k);
        	rootDiam *= treesParams.rootDiamMultiplier;
        	SRL /= pow(treesParams.rootDiamMultiplier, 2.0);
        }
	return(totalRootArea);
}
double BiogeochemicalCycles::computeFineRootArea(trees_params treesParams,
						 double SRA,
					         int j,
						 int k)
{
        double rootCarbon, totalRootArea;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nFineRootOrders);

        totalRootArea = 0.0;

	rootCarbon = fineRootBiomassCarbon[j][k] + coarseRootBiomassCarbon[j][k];

//convert kgC ha-1 to kgC m-2
//m2 root area m-2 ground area = 
//      kgC root m-2 ground area * m-2 ground area * m2 root m-2 ground area kgC-1 root
//
	rootArea[j][k] = 0.0001 * rootCarbon * SRA; //m2 root m2 ground area
        totalRootArea += rootArea[j][k];

	return(totalRootArea);
}

//
//computeMaximumRootBiomassCarbon()
//
//This function sets upper limits on carbon allocation to each root order in each layer,
//    by assuming the maximum available volume is given by the cyclinder defined by
//    root diameter times soil porosity.
//
//Note: root and leaf areas are per unit ground area; mass is per hectare
//
//root segments are approximated as cylinders
//specific root length (m kgC-1 root) is input relationship between carbon and root length
//SRL1 is defined at a root diameter of 0.00025 m
//finest root diameter is an input m; maximum root diameter is an input
//
void BiogeochemicalCycles::computeMaximumRootBiomassCarbon(trees_params treesParams)
{
	double maxCarbon, totalRootAxialLength, rootAxialLength, lateralLength, porosity;
	double refRootDiam, rootDiam, SRL, areaScalar, allocationScalar;
	double poreVolume, rhizWidth, rhizCrossArea, maxRootLength;

	porosity = treesParams.porosity;
	rhizWidth = 0.001*treesParams.rhizosphere_width;

	totalRootAxialLength = 0.0;
	for (int j = 0; j < nRoots; ++j)
	{
		totalRootAxialLength += treesParams.drax[j+3];
	}

	refRootDiam = 0.00025; //m

	for (int j = 0; j < nRoots; ++j)
	{
		rootDiam = treesParams.minRootDiam; //m
        	SRL = 1000.0*treesParams.SRL1*pow(refRootDiam, 2.0)/pow(rootDiam,2.0);
		lateralLength = treesParams.drlat[j+3]*totalRootAxialLength+0.0000001;
		rootAxialLength = treesParams.drax[j+3];
//how many defined cylinders fit into a hectare, m2 m-2
		areaScalar = 10000.0/(M_PI*pow(lateralLength, 2.0));
		rootAreaScalar[j] = 1.0/areaScalar;
		poreVolume = M_PI * pow(lateralLength, 2.0) * rootAxialLength * porosity;
		allocationScalar = 0.5;
		for (int k = 0; k < nRootOrders; ++k)
		{
			rhizCrossArea = M_PI * pow(0.5*rootDiam+0.25*rhizWidth, 2.0)+0.0000001;
			maxRootLength = poreVolume/rhizCrossArea;
			maximumRootBiomassCarbon[j][k] = areaScalar*allocationScalar*maxRootLength/SRL;
        		SRL /= pow(treesParams.rootDiamMultiplier, 2.0);
        		rootDiam *= treesParams.rootDiamMultiplier;
			allocationScalar /= 2.0;
		}
	}
}


//
//computeCanopyCover()
//
//This functions computes the projected crown area divided by ground area
//    defined by root extent. When this number exceed unity it means the
//    roots are confined relative to the canopy extent, and when it is
//    less than unity then the canopy is confined relative to the roots.
//
double BiogeochemicalCycles::computeCanopyCover(trees_params treesParams)
{
	double totalRootAxialLength, rootAxialLength, maximumLateralLength;
	double crownRadius, crownArea, groundArea, canopy_cover;

//Define ground area using longest lateral root length (defined in param_mod)	
	totalRootAxialLength = 0.0;
	for (int j = 0; j < nRoots; ++j)
	{
		totalRootAxialLength += treesParams.drax[j+3];
	}
	maximumLateralLength = 0.0;
	for (int j = 0; j < nRoots; ++j)
	{
		maximumLateralLength = max(maximumLateralLength, treesParams.drlat[j+3]*totalRootAxialLength);;
	}
	groundArea = M_PI * pow(maximumLateralLength, 2.0)+0.0000001;

//Define the canopy closure based on crown diameter
        if (treesParams.usePhenology == false)
        {
                crownRadius = treesParams.dslat[1]*treesParams.dsax[1];
                crownArea = M_PI*pow(crownRadius,2.0);
                canopy_cover = crownArea/groundArea;
cout << "canopy_cover = " << canopy_cover << endl;
        }
        else
        {
                canopy_cover = 1.0;
        }
	return(canopy_cover);
}


//
//computeLateralRootWeibulls(), updateLateralRootWeibulls()
//
//  This code adjusts the lateral root saturated vulnerability curve parameters.
//  Adjustment of b is downward for increases in fine root area and increased
//    fine root CN ratio. The adjustment of c is upward for the same conditions.
//
void BiogeochemicalCycles::computeLateralRootWeibulls(trees_params& treesParams)
{
        double bval, cval, bsum, csum;
        double totalRootArea, totalRootCarbon, rootCarbonRatio;

        totalRootArea = computeRootArea(treesParams);
        for (int j = 0; j < nRoots; j++)
        {
        	bsum = csum = 0.0;
                bval = lat_Root_b_value_init[0];
                cval = lat_Root_c_value_init[0];
		totalRootCarbon = getFineRootCarbon(j)+0.0000001;
		rootCarbonRatio = fineRootBiomassCarbon[j][2]/totalRootCarbon;
                bsum += bval*rootCarbonRatio;
                csum += cval*rootCarbonRatio;
                for (int k = 1; k >= 0; k--)
                {
                        bval *= 0.95;
                        cval /= 0.95;
			rootCarbonRatio = fineRootBiomassCarbon[j][k]/totalRootCarbon;
                        bsum += bval*rootCarbonRatio;
                        csum += cval*rootCarbonRatio;
                }
                bval = lat_Root_b_value_init[0];
                cval = lat_Root_c_value_init[0];
                for (int k = 3; k <= 4; k++)
                {
                        bval /= 0.95;
                        cval *= 0.95;
			rootCarbonRatio = fineRootBiomassCarbon[j][k]/totalRootCarbon;
                        bsum += bval*rootCarbonRatio;
                        csum += cval*rootCarbonRatio;
                }
        }
        treesParams.lat_Root_b_value = bsum;
        treesParams.lat_Root_c_value = csum;
}
void BiogeochemicalCycles::computeLateralRootWeibulls(int j, 
							trees_params treesParams,
							double& bsum,
							double& csum)
{
        double bval, cval;
        double totalRootArea, CN, CNratio, SLAscalar, rootScalar;
	double totalRootCarbon, rootCarbonRatio;
	double bfactor = 0.95;
	double cfactor = 1.0/0.95;
	SLAscalar = 22.0/treesParams.SLA;
	if (SLAscalar > 3.0)
	{
		SLAscalar = 3.0;
	}
	else if (SLAscalar < 1.0)
	{
		SLAscalar = 1.0;
	}
	rootScalar = 1.5625;
        bsum = csum = 0.0;
	totalRootCarbon = getFineRootCarbon(j)+0.0000001;
	CN = fineRootBiomassCarbon[j][2]/fineRootBiomassNitrogen[j][2];
	CN = max(CN, 20.0);
	CNratio = 0.5 * (1.0 + 20.0*rootScalar*SLAscalar/CN);
	CNratio = max(0.8, CNratio);
	CNratio = min(1.2, CNratio);
        bval = lat_Root_b_value_init[0]*CNratio;
        cval = lat_Root_c_value_init[0]/CNratio;
	rootCarbonRatio = fineRootBiomassCarbon[j][2]/totalRootCarbon;
        bsum += bval*rootCarbonRatio;
        csum += cval*rootCarbonRatio;
        for (int k = 1; k >= 0; k--)
        {
		rootScalar /= 1.25;
		CN = fineRootBiomassCarbon[j][k]/fineRootBiomassNitrogen[j][k];
		CN = max(CN, 20.0);
		CNratio = 0.5 * (1.0 + 20.0*rootScalar*SLAscalar/CN);
		CNratio = max(0.8, CNratio);
		CNratio = min(1.2, CNratio);
        	bval = lat_Root_b_value_init[0]*(bfactor*CNratio);
                cval = lat_Root_c_value_init[0]*(cfactor/CNratio);
		rootCarbonRatio = fineRootBiomassCarbon[j][k]/totalRootCarbon;
                bsum += bval*rootCarbonRatio;
                csum += cval*rootCarbonRatio;
		bfactor *= 0.95;
		cfactor /= 0.95;
        }
	rootScalar = 1.5625;
	bfactor = 1.0/0.95;
	cfactor = 0.95;
        bval = lat_Root_b_value_init[0]*bfactor;
        cval = lat_Root_c_value_init[0]*cfactor;
        for (int k = 3; k <= 4; k++)
        {
		rootScalar *= 1.25;
		CN = fineRootBiomassCarbon[j][k]/fineRootBiomassNitrogen[j][k];
		CN = max(CN, 20.0);
		CNratio = 20.0*rootScalar*SLAscalar/CN;
		CNratio = max(0.8, CNratio);
		CNratio = min(1.2, CNratio);
        	bval = lat_Root_b_value_init[0]*(bfactor*CNratio);
                cval = lat_Root_c_value_init[0]*(cfactor/CNratio);
		rootCarbonRatio = fineRootBiomassCarbon[j][k]/totalRootCarbon;
                bsum += bval*rootCarbonRatio;
                csum += cval*rootCarbonRatio;
		bfactor /= 0.95;
		cfactor *= 0.95;
        }
}
void BiogeochemicalCycles::updateLateralRootWeibulls(double bsat[][MD], 
							double ccsat[][MD],
							trees_params treesParams)
{
	double bsum, csum;
        for (int j = treesParams.smodules+2; j < nRoots+treesParams.smodules+2; j++)
        {
        	bsum = csum = 0.0;
		computeLateralRootWeibulls(j-treesParams.smodules-2, treesParams, bsum, csum);
		bsat[j][1] = bsum;
		ccsat[j][1] = csum;
        }
}

//
//getLeafBiomassCarbon()
//
double BiogeochemicalCycles::getLeafBiomassCarbon()
{
	return(leafBiomassCarbon[0]);
}

//
//getFruitCarbon()
//
double BiogeochemicalCycles::getFruitCarbon()
{
	return(fruitCarbon[0]);
}

//
//getRootCarbon()
//This code sums up all the carbon within the root system
//
double BiogeochemicalCycles::getRootCarbon()
{
	double rootCarbon = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		rootCarbon += getRootCarbon(j);
	}
	return (rootCarbon);
}
double BiogeochemicalCycles::getRootCarbon(int j)
{
	double rootCarbon = 0.0;

	for (int k = 0; k < nRootOrders; k++)
	{
		rootCarbon += getRootCarbon(j, k);
	}
	return (rootCarbon);
}
double BiogeochemicalCycles::getRootCarbon(int j, 
					   int k)
{
	double rootCarbon;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	rootCarbon = fineRootBiomassCarbon[j][k] + coarseRootBiomassCarbon[j][k];
	return(rootCarbon);
}

//
//getFineRootCarbon()
//This code sums up all the fine root carbon within the root system
// --Note: Assumption is that the first 5 root orders are fine (0.25 mm to 4.0 mm diameter)
//
double BiogeochemicalCycles::getFineRootCarbon()
{
	double rootCarbon = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		rootCarbon += getFineRootCarbon(j);
	}
	return (rootCarbon);
}
double BiogeochemicalCycles::getFineRootCarbon(int j)
{
	double rootCarbon = 0.0;

	for (int k = 0; k < nFineRootOrders; k++)
	{
		rootCarbon += getFineRootCarbon(j, k);
	}
	return (rootCarbon);
}
double BiogeochemicalCycles::getFineRootCarbon(int j,
					       int k)
{
	double rootCarbon;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	rootCarbon = fineRootBiomassCarbon[j][k];
	return (rootCarbon);
}

//
//getLiveStemCarbon()
//
double BiogeochemicalCycles::getLiveStemCarbon()
{
	double stemCarbon;
	stemCarbon = liveStemCarbon[0];
	return(stemCarbon);
}

//
//getDeadStemCarbon()
//
double BiogeochemicalCycles::getDeadStemCarbon()
{
	double stemCarbon;
	stemCarbon = deadStemCarbon[0];
	return(stemCarbon);
}

//
//getLeafNSC()
//
double BiogeochemicalCycles::getLeafNSC()
{
	double nsc;
	nsc = leafNSC[0];
	return(nsc);
}

//
//getChloroplastStarch()
//
double BiogeochemicalCycles::getChloroplastStarch()
{
	double nsc;
	nsc = chloroplastStarch[0];
	return(nsc);
}

//
//getChloroplastSugar()
//
double BiogeochemicalCycles::getChloroplastSugar()
{
	double nsc;
	nsc = chloroplastSugar[0];
	return(nsc);
}

//
//getStemNSC()
//
double BiogeochemicalCycles::getStemNSC()
{
	double nsc;
	nsc = stemNSC[0];
	return(nsc);
}

//
//getLeafBiomassN()
//
double BiogeochemicalCycles::getLeafBiomassN()
{
        return(leafBiomassNitrogen[0]);
}

//
//getLeafMineralN()
//
double BiogeochemicalCycles::getLeafMineralN()
{
        return(leafStoredNitrogen[0]);
}

//
//getRootNSC()
//
double BiogeochemicalCycles::getRootNSC()
{
	double nsc = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		nsc += getRootNSC(j);
	}
	return(nsc);
}
double BiogeochemicalCycles::getRootNSC(int j)
{
	double nsc = 0.0;

	for (int k = 0; k < nRootOrders; k++)
	{
		nsc += getRootNSC(j, k);
	}
	return(nsc);
}
double BiogeochemicalCycles::getRootNSC(int j,
					int k)
{
	double nsc;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	nsc = rootNSC[j][k];
	return(nsc);
}

//
//getFineRootNSC()
//
double BiogeochemicalCycles::getFineRootNSC()
{       
        double nsc = 0.0;
        
        for (int j = 0; j < nRoots; j++)
        {
                nsc += getFineRootNSC(j);
        }
        return(nsc);
}
double BiogeochemicalCycles::getFineRootNSC(int j)
{       
        double nsc = 0.0;
        
        for (int k = 0; k < nFineRootOrders; k++)
        {
                nsc += getFineRootNSC(j, k);
        }
        return(nsc);
}
double BiogeochemicalCycles::getFineRootNSC(int j,
                                        int k)
{
        double nsc;

        assert(j >= 0);
        assert(j < nRoots);
        assert(k >= 0);
        assert(k < nRootOrders);

        nsc = rootNSC[j][k];
        return(nsc);
}


//
//putLeafNSC()
//
void BiogeochemicalCycles::putLeafNSC(double nsc)
{
	leafNSC[0] = nsc;
}

//
//putChloroplastStarch()
//
void BiogeochemicalCycles::putChloroplastStarch(double nsc)
{
	chloroplastStarch[0] = nsc;
}

//
//putChloroplastSugar()
//
void BiogeochemicalCycles::putChloroplastSugar(double nsc)
{
	chloroplastSugar[0] = nsc;
}

//
//putStemNSC()
//
void BiogeochemicalCycles::putStemNSC(double nsc)
{
	stemNSC[0] = nsc;
}

//
//updateRootNSC()
//This functions serves primiarly as a means of deducting NSC from roots
//   as a cost of root respiration
//Use the first function if you don't have root-specific NSC updates, as this
//   function will distribute the NSC consumption in proportion to the live
//   root carbon.
//The better approach is to have root-specific NSC updates, in which case
//   you would use the second function.
//
void BiogeochemicalCycles::updateRootNSC(double delta_nsc)
{
	double rootC, rootCproportion;

	rootC = getRootCarbon();
	if (rootC > 0.00001)
	{
		for (int j = 0; j < nRoots; j++)
		{
			for (int k = 0; k < nRootOrders; k++)
			{
				rootCproportion = getRootCarbon(j,k)/(rootC+0.00000001);
				rootNSC[j][k] += delta_nsc*rootCproportion;
				if (rootNSC[j][k] < 0.000001)
				{
					rootNSC[j][k] = 0.000001;
				}
			}
		}
	}
}
void BiogeochemicalCycles::updateRootNSC(double delta_nsc,
					 int j,
					 int k)
{

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	rootNSC[j][k] += delta_nsc;

}

//
//getHumus()
//
//Note: humus is currently assumed to be well-mixed by soil-root layer, and
//      so it is undifferentiated at the level of root orders
//
double BiogeochemicalCycles::getHumus()
{
	double humus = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		humus += getHumus(j);
	}
	return(humus);
}
double BiogeochemicalCycles::getHumus(int j)
{
	double humus;

	assert(j >= 0);
	assert(j < nRoots);

	humus = humusCarbon[j];
	return humus;
}

//
//getRhizosphereCl()
//
double BiogeochemicalCycles::getRhizosphereCl()
{
	double Cl = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		Cl += getRhizosphereCl(j);
	}
	return(Cl);
}
double BiogeochemicalCycles::getRhizosphereCl(int j)
{
	double Cl = 0.0;

	for (int k = 0; k < nRootOrders; k++)
	{
		Cl += getRhizosphereCl(j, k);
	}
	return(Cl);
}
double BiogeochemicalCycles::getRhizosphereCl(int j,
 					      int k)
{
	double Cl;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	Cl = rhizosphereCl[j][k];
	return(Cl);
}

//
//getRhizosphereNl()
//
double BiogeochemicalCycles::getRhizosphereNl()
{
	double Nl = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		Nl += getRhizosphereNl(j);
	}
	return(Nl);
}
double BiogeochemicalCycles::getRhizosphereNl(int j)
{
	double Nl = 0.0;

	for (int k = 0; k < nRootOrders; k++)
	{
		Nl += getRhizosphereNl(j, k);
	}
	return(Nl);
}
double BiogeochemicalCycles::getRhizosphereNl(int j,
					      int k)
{
	double Nl;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	Nl = rhizosphereNl[j][k];
	return(Nl);
}

//
//getAminoAcidExudateC()
//
double BiogeochemicalCycles::getAminoAcidExudateC()
{
	double Ce = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		Ce += getAminoAcidExudateC(j);
	}
	return(Ce);
}
double BiogeochemicalCycles::getAminoAcidExudateC(int j)
{
	double Ce = 0.0;

	for (int k = 0; k < nRootOrders; k++)
	{
		Ce += getAminoAcidExudateC(j, k);
	}
	return(Ce);
}
double BiogeochemicalCycles::getAminoAcidExudateC(int j,
						  int k)
{
	double Ce;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	Ce = rootExudateAminoAcidCarbon[j][k];
	return(Ce);
}

//
//getAminoAcidExudateN()
//
double BiogeochemicalCycles::getAminoAcidExudateN()
{
	double Ne = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		Ne += getAminoAcidExudateN(j);
	}
	return(Ne);
}
double BiogeochemicalCycles::getAminoAcidExudateN(int j)
{
	double Ne = 0.0;

	for (int k = 0; k < nRootOrders; k++)
	{
		Ne += getAminoAcidExudateN(j, k);
	}
	return(Ne);
}
double BiogeochemicalCycles::getAminoAcidExudateN(int j,
						  int k)
{
	double Ne;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	Ne = rootExudateAminoAcidNitrogen[j][k];
	return(Ne);
}

//
//getSugarExudateC()
//
double BiogeochemicalCycles::getSugarExudateC()
{
	double Ce = 0.0;

	for (int j = 0; j < nRoots; j++)
	{
		Ce += getSugarExudateC(j);
	}
	return(Ce);
}
double BiogeochemicalCycles::getSugarExudateC(int j)
{
	double Ce = 0.0;

	for (int k = 0; k < nRootOrders; k++)
	{
		Ce += getSugarExudateC(j, k);
	}
	return(Ce);
}
double BiogeochemicalCycles::getSugarExudateC(int j,
					      int k)
{
	double Ce;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

	Ce = rootExudateSugarCarbon[j][k];
	return(Ce);
}

//
//getRhizosphereLiveMicrobialCarbon()
//
double BiogeochemicalCycles::getRhizosphereLiveMicrobialCarbon()
{
        double Cm = 0.0;

        for (int j = 0; j < nRoots; j++)
        {
                Cm += getRhizosphereLiveMicrobialCarbon(j);
        }
        return(Cm);
}
double BiogeochemicalCycles::getRhizosphereLiveMicrobialCarbon(int j)
{
	double Cm = 0.0;

	for (int k = 0; k < nRootOrders; k++)
	{
		Cm += getRhizosphereLiveMicrobialCarbon(j, k);
	}
	return(Cm);
}
double BiogeochemicalCycles::getRhizosphereLiveMicrobialCarbon(int j,
                                            		       int k)
{
        double Cm;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

        Cm = rhizosphereLiveMicrobialCarbon[j][k];
        return(Cm);
}

//
//getRhizosphereMicrobialNitrogen()
//
double BiogeochemicalCycles::getRhizosphereMicrobialNitrogen()
{
        double Nm = 0.0;

        for (int j = 0; j < nRoots; j++)
        {
                Nm += getRhizosphereMicrobialNitrogen(j);
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereMicrobialNitrogen(int j)
{
	double Nm = 0.0;

	for (int k = 0; k < nRootOrders; k++)
	{
		Nm += getRhizosphereMicrobialNitrogen(j, k);
	}
	return(Nm);
}
double BiogeochemicalCycles::getRhizosphereMicrobialNitrogen(int j,
                                            		     int k)
{
        double Nm;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

        Nm = rhizosphereMicrobialNitrogen[j][k];
        return(Nm);
}

//
//getRhizosphereDeadMicrobialCarbon()
//
double BiogeochemicalCycles::getRhizosphereDeadMicrobialCarbon()
{
        double Cm = 0.0;

        for (int j = 0; j < nRoots; j++)
        {                                   
                Cm += getRhizosphereDeadMicrobialCarbon(j);
        }
        return(Cm);
}
double BiogeochemicalCycles::getRhizosphereDeadMicrobialCarbon(int j)
{
        double Cm = 0.0;

        for (int k = 0; k < nRootOrders; k++)
        {
        	Cm += getRhizosphereDeadMicrobialCarbon(j, k);
        }
        return(Cm);
}
double BiogeochemicalCycles::getRhizosphereDeadMicrobialCarbon(int j,
                                            		       int k)
{
        double Cm;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

        Cm = rhizosphereDeadMicrobialCarbon[j][k];
        return(Cm);
}

//
//getRhizosphereMineralNitrogen()
//
double BiogeochemicalCycles::getRhizosphereMineralNitrogen()
{
        double Nm = 0.0;

        for (int j = 0; j < nRoots; j++)
        {
                Nm += getRhizosphereMineralNitrogen(j);
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereMineralNitrogen(int j)
{
        double Nm = 0.0;

        for (int k = 0; k < nRootOrders; k++)
        {
                Nm += getRhizosphereMineralNitrogen(j, k);
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereMineralNitrogen(int j,
                                            		   int k)
{
        double Nm;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

        Nm = rhizosphereMineralNitrogen[j][k];
        return(Nm);
}

//
//getNitrogenLeaching()
//
double BiogeochemicalCycles::getNitrogenLeaching(int j)
{
	double Nl;

	Nl = nitrogenLeaching[j];
	return(Nl);
}

//
//getRhizosphereAmmoniumNitrogen()
//
double BiogeochemicalCycles::getRhizosphereAmmoniumNitrogen()
{
        double Nm = 0.0;

        for (int j = 0; j < nRoots; j++)
        {
                Nm += getRhizosphereAmmoniumNitrogen(j);
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereAmmoniumNitrogen(int j)
{
        double Nm = 0.0;

        for (int k = 0; k < nRootOrders; k++)
        {
                Nm += getRhizosphereAmmoniumNitrogen(j, k);
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereAmmoniumNitrogen(int j,
                                            		    int k)
{
        double Nm;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

        Nm = rhizosphereAmmoniumNitrogen[j][k];
        return(Nm);
}

//
//getRhizosphereNitrateNitrogen()
//
double BiogeochemicalCycles::getRhizosphereNitrateNitrogen()
{
        double Nm = 0.0;

        for (int j = 0; j < nRoots; j++)
        {
                Nm += getRhizosphereNitrateNitrogen(j);
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereNitrateNitrogen(int j)
{
        double Nm = 0.0;

        for (int k = 0; k < nRootOrders; k++)
        {
                Nm += getRhizosphereNitrateNitrogen(j, k);
        }
        return(Nm);
}
double BiogeochemicalCycles::getRhizosphereNitrateNitrogen(int j,
                                            		   int k)
{
        double Nm;

	assert(j >= 0);
	assert(j < nRoots);
	assert(k >= 0);
	assert(k < nRootOrders);

        Nm = rhizosphereNitrateNitrogen[j][k];
        return(Nm);
}

//
//computeNSCfluxes(kratio, kratio_vector)
//Move NSC between leaf-stem and stem-root pools based on concentration gradients
//  and relative hydraulic conductance (kratio)
//
//Assumes phloem loading has already taken place and energy costs consumed
//  -- these costs are taken where net photosynthesis is added to the PSN state variable in 
//     simulation_function()
//
//Other assumptions: 0.16 diffusion term assumes that it would take 6 simulation time
//                   steps or 3 hours to move all NSC from leaf to stem, and another
//		     3 hours to move from stem to roots
//
void BiogeochemicalCycles::computeNSCfluxes(trees_params treesParams,
					    double kratio,
				       	    double* kratio_vector)
{
	double C, flux, leafCfract, stemCfract, rootCfract; 
	double rootCarbon, singleRootC, kratioRoot;
	double cfract;
//flux between canopy and stem
	leafCfract = getLeafNSC()/(getLeafBiomassCarbon()+0.0000001);
	stemCfract = getStemNSC()/(getLiveStemCarbon()+0.0000001);
	//cfract = getFineRootCarbon()/(getLeafBiomassCarbon()+0.0000001);;
	cfract = 2.5;
	if (leafCfract > cfract*stemCfract)
	{
		C = 0.16*getLeafNSC();
	}
	else if (leafCfract < cfract*stemCfract)
	{
		C = 0.16*getStemNSC();
	}
	else
	{
		C = 0.0;
	}
	//cfract = 1.0;
	flux = C*kratio*(leafCfract-cfract*stemCfract);
	if (flux > 0.9*getLeafNSC() || -flux > 0.9*getStemNSC())
	{
		flux = 0.0;
	}
	putLeafNSC(getLeafNSC()-flux);
	putStemNSC(getStemNSC()+flux);
//flux between stem and roots
	rootCarbon = getRootCarbon();
	for (int j = 0; j < nRoots; j++)
	{
		kratioRoot = kratio_vector[j];

//Note: there needs to be a graceful way to handle soil-root zones with near zero fine root
//      to avoid underflow problems that lead to hydraulic failure

/*
		if (treesParams.drlat[j+3] < 0.0001)
		{
			continue;
		}
*/

//nFineRootOrders
		cfract = 0.4;
		for (int k = 0; k < nRootOrders; k++)
		{
			singleRootC = getRootCarbon(j, k);
			rootCfract = getRootNSC(j, k)/(singleRootC+0.0000001);
			stemCfract = getStemNSC()/(getLiveStemCarbon()+0.0000001);
			if (stemCfract > cfract*rootCfract)
			{
				C = 0.16*getStemNSC()*(singleRootC/(rootCarbon+0.0000001));
			}
			else if (stemCfract < cfract*rootCfract)
			{
				C = 0.16*getRootNSC(j, k)*(singleRootC/(rootCarbon+0.0000001));
			}
			else
			{
				C = 0.0;
			}
			flux = C*kratioRoot*(stemCfract - cfract*rootCfract);
			if (flux > 0.9*getStemNSC() || -flux > 0.9*getRootNSC(j, k))
			{
				flux = 0.0;
			}
			updateRootNSC(flux, j, k);
			putStemNSC(getStemNSC()-flux);
			//cfract += 0.03;
		}
	}
}

//
//computeNfluxes(kratio, kratio_vector)
//Move N between leaf-root pools based on concentration gradients
//    and relative hydraulic conductance
//
void BiogeochemicalCycles::computeNfluxes(trees_params treesParams,
					  double kratio,
					  double* kratio_vector)
{
	double C, nfract, nscalar, flux, leafNfract, rootNfract; 
	double rootNitrogen, leafNstore, singleRootN, kratioRoot;

	leafNstore = leafStoredNitrogen[0];
	rootNitrogen = getRootBiomassN();
	//nscalar = rootNitrogen/(getLeafBiomassN()+0.0000001);
	nscalar = 1.0;
	for (int j = 0; j < nRoots; j++)
	{
//Note: there needs to be a graceful way to handle soil-root zones with near zero fine root
//      to avoid underflow problems that lead to hydraulic failure

		if (treesParams.drlat[j+3] < 0.0001)
		{
			continue;
		}

		//nfract = 0.1*nscalar;
		nfract = 1.0*nscalar;
		kratioRoot = kratio_vector[j];
		for (int k = 0; k < nRootOrders; k++)
		{
			singleRootN = getRootBiomassN(j, k);
			leafNfract = leafNstore/(getLeafBiomassN()+0.0000001);
			rootNfract = getRootN(j, k)/(singleRootN+0.0000001);
			if (leafNfract > nfract*rootNfract)
			{
				C = 0.16*leafNstore*(singleRootN/(rootNitrogen+0.00000001));
			}
			else if (leafNfract < nfract*rootNfract)
			{
				C = 0.16*getRootN(j, k)*(singleRootN/(rootNitrogen+0.00000001));
			}
			else
			{
				C = 0.0;
			}
			flux = C*kratioRoot*(leafNfract - nfract*rootNfract);
			if (rootMineralNitrogen[j][k] + flux <= 0.0)
			{
				flux = 0.0;
			}
			rootMineralNitrogen[j][k] += flux;
			leafStoredNitrogen[0] -= flux;
			//nfract += 0.03;
		}
	}
}

//
//computeRootExudates()
//For each rhizosphere define two types of root exudates that can be used to manipulate microbes:
//1) simple sugars
//2) amino acids
//For sugars only carbon units are moved from the root to the rhizosphere, while for amino acids
//  both carbon and nitrogen are output. The assumption is a 2C to 1N promportion or 24 moles of C
//  to 14 moles of N for a C:N ratio of 24/14.
//Export of exudates is assumed to be of excess NSC and N at the end of a diurnal growth cycle
//
void BiogeochemicalCycles::computeRootExudates(trees_params treesParams,
					       double* thetaSoil,
					       double* kratio_vector)
{
	for (int j = 0; j < nRoots; j++)
        {
		computeRootExudates(treesParams, thetaSoil, kratio_vector, j);
	}
}
void BiogeochemicalCycles::computeRootExudates(trees_params treesParams,
					       double* thetaSoil,
					       double* kratio_vector,
					       int j)
{
	for (int k = 0; k < nFineRootOrders; k++)
        {
		computeRootExudates(treesParams, thetaSoil, kratio_vector, j, k);
	}
}
void BiogeochemicalCycles::computeRootExudates(trees_params treesParams,
					       double* thetaSoil,
					       double* kratio_vector,
					       int j,
					       int k)
{
	double residualNSC, residualN, aminoAcidC, aminoAcidN, CN;
	double rhizAAconc, rhizSconc, rootNSCconc, degreeSaturation;
	double dDist, kratioRoot, flux, rootC, rhizC, nsc, porosity;
	double relativeRootArea[5] = { 1.0, 0.25, 0.0625, 0.01625, 0.00390625 };

	assert(k < nFineRootOrders);

//carbon:nitrogen ratio for amino acid
	CN = 24.0/14.0;
//distance from root to the mid-point of the rhizosphere
	dDist = 0.5*0.001*treesParams.rhizosphere_width;
	kratioRoot = kratio_vector[j];

//compute rhizosphere degree of saturation
	porosity = treesParams.porosity;
	degreeSaturation = thetaSoil[j]/porosity;
	rhizC = getRhizosphereCl(j, k) + (getAminoAcidExudateC(j, k) + getSugarExudateC(j, k))/degreeSaturation;
	rhizAAconc = getAminoAcidExudateC(j, k)/degreeSaturation/(rhizC+0.0000001);
	rhizSconc = getSugarExudateC(j, k)/degreeSaturation/(rhizC+0.0000001);

//move NSC and N from the higher order roots to prime the pump for exudates
	rootNSC[j][k] += 0.5*relativeRootArea[k]*rootNSC[j][k+5];
	rootNSC[j][k+5] -= 0.5*relativeRootArea[k]*rootNSC[j][k+5];
	rootMineralNitrogen[j][k] += 0.5*relativeRootArea[k]*rootMineralNitrogen[j][k+5];
	rootMineralNitrogen[j][k+5] -= 0.5*relativeRootArea[k]*rootMineralNitrogen[j][k+5];
	nsc = getRootNSC(j, k);

//compute residual NSC and N that is potentially available for exudates
	rootC = getRootCarbon(j, k);
	residualNSC = nsc - treesParams.leafNSCscalar*rootC;
	if (residualNSC < 0.0)
	{
		residualNSC = 0.0;
	}
	residualN = getFineRootN(j, k) - getFineRootBiomassN(j, k);
	if (residualN < 0.0)
	{
		residualN = 0.0;
	}

//amino acid production is limited first by residualN, then by residual NSC
//assumes that 20% C is a cost to produce the amino acid
//flux is dimensionless, accounts for normalized concentration differences 
//  between root NSC and exudates, and relative hydraulic conductance
	rootNSCconc = nsc/(nsc+rootC+0.0000001);
	flux = ((rootNSCconc-rhizAAconc)/(rootNSCconc+rhizAAconc))*kratioRoot/12.0;
	flux *= relativeRootArea[k];
	if (flux > 1.0)
	{
		flux = 1.0;
	}
	if (residualN > 0.0 && flux > 0.0)
	{
		aminoAcidC = residualN * CN;
		aminoAcidC *= flux;
		if (aminoAcidC*1.2 > residualNSC)
		{
			aminoAcidC = residualNSC/1.2;
			residualNSC = 0.0;
		}
		else
		{
			residualNSC -= 1.2*aminoAcidC;
		}
		aminoAcidN = aminoAcidC / CN;
		heterotrophicRespiration[j] += aminoAcidC * 0.2;
	} 
	else
	{
		aminoAcidC = 0.0;
		aminoAcidN = 0.0;
	}
	rootExudateAminoAcidCarbon[j][k] += aminoAcidC;
	rootExudateAminoAcidNitrogen[j][k] += aminoAcidN;
	updateRootNSC(-aminoAcidC, j, k);
	rootMineralNitrogen[j][k] -= aminoAcidN;

//re-compute concentration of NSC in the root
	nsc = getRootNSC(j, k);
	rootNSCconc = nsc/(nsc+rootC+0.0000001);

//sugar exudates are whatever residualNSC is left over
//sugar export is reduced as its concentration in the rhizosphere exceeds amino acid concentration
	flux = ((rootNSCconc-rhizSconc)/(rootNSCconc+rhizSconc))*kratioRoot/12.0;
	flux *= min(1.0, getAminoAcidExudateC(j, k)/(getSugarExudateC(j, k) + 0.0000001));
	flux *= relativeRootArea[k];
	if (flux > 1.0)
	{
		flux = 1.0;
	}
	if (residualNSC > 0.0 && flux > 0.0)
	{
		residualNSC *= flux;
		if (rootExudateSugarCarbon[j][k]+residualNSC > 4.66*rootExudateAminoAcidCarbon[j][k])
		{
			residualNSC = 4.66*rootExudateAminoAcidCarbon[j][k] - rootExudateSugarCarbon[j][k];
		}
		rootExudateSugarCarbon[j][k] += residualNSC;
		rootNSC[j][k] -= residualNSC;
	}
}


//
//computeLeafAllocation()
//Update growth and leaf area and root-to-leaf area ratio if not in dormancy (root T < 5 C)
//What to do about LAI - need to modify lai and Al; lai_at_sat_kl is used only once per simulation
//Currently assumes allocation to new leaf growth is up to 45% of growth respiration, 
//    and leaf longevity is an input
//Assuming leaf biomass can be approximated by 86% (carbon in cellulose) of NSC use
//This computes potential lai at saturated kl, with actual determined by phenology
//NEED - LIVE LAI for canopy conductance since brown canopy LAI absorbs energy, but does not transpire
//
void BiogeochemicalCycles::computeLeafAllocation(trees_params& treesParams,
						 double newc[][MD],
						 double N_avail_rate_plant,
						 double kratio,
						 double nsc,
						 double psn,
						 double nscRatio,
						 double& rgrowth,
						 double& leafCfraction,
						 double lai,
						 int yday,
						 double& stressedLeafLifeSpan,
						 double SLA)
{
	double deltaLAI, unstressedLeafLifeSpan, SLA_instant;
	double reproduction, fromStem, CN;

        unstressedLeafLifeSpan = treesParams.leafLifeSpan*365.25*48.0;
	if (treesParams.usePhenology == 1 && treesParams.leafLifeSpan == 1.0 && lai < 0.2)
	{
		unstressedLeafLifeSpan *= 10.0;
	}
	if (treesParams.usePhenology == 0)
	{
		unstressedLeafLifeSpan *= 1000.0;
	}

//start with a base fraction of carbon allocation to leaf
        leafCfraction = 0.4*N_avail_rate_plant;
        if (leafCfraction < 0.0001)
        {
                leafCfraction = 0.0001;
        }
        else if (leafCfraction > 0.40)
        {
                leafCfraction = 0.40;
        }
//if available N is low then reduce allocation to leaf
	//CN = 50.0 - 26.0*N_avail_rate_plant;
//DSM - July 2020
	CN = treesParams.leaf_maxCN - (treesParams.leaf_maxCN-treesParams.leaf_minCN)*N_avail_rate_plant;

//---------------------------------------
//For Maize reproduction - DSM September 2019
	//if (yday > 150 && yday < 200)
	if (yday > 163 && yday < 215)
	{
		for (int k = 0; k < nRootOrders; k++)
		{
			rhizosphereNitrateNitrogen[0][k] += 0.05/48.0;
		}
	}
	//if (treesParams.usePhenology == 0 && treesParams.useLeafModule == 0 && 
	//	(treesParams.live_lai > treesParams.lai_at_full_canopy_height || yday > 200))
	//if (treesParams.usePhenology == 0 && treesParams.useLeafModule == 0 && 
	//	(treesParams.live_lai > treesParams.lai_at_full_canopy_height || yday > 210))
	if (treesParams.usePhenology == 0 && treesParams.useLeafModule == 0 && 
		(treesParams.live_lai > 0.99*treesParams.lai_at_full_canopy_height))
	{
		nsc = getLeafNSC()+getStemNSC();
		//if (yday > 200 || (yday > 193 && treesParams.live_lai > treesParams.lai_at_full_canopy_height))
		if (yday > 210 || (yday > 203 && treesParams.live_lai > 0.99*treesParams.lai_at_full_canopy_height))
		{
			//reproduction = rgrowth*N_avail_rate_plant*kratio;
			if (rgrowth > 0.8*0.95*nsc)
			{
				reproduction = 0.8*0.95*nsc;
			}
			else
			{
				reproduction = 0.8*rgrowth;
			}
		}
		else
		{
			reproduction = 0.0;
		}
		if (reproduction / CN > 0.95 * leafStoredNitrogen[0])
		{
			reproduction = CN * 0.95 * leafStoredNitrogen[0];
		}
		if (reproduction < 0.0)
		{
			reproduction = 0.0;
		}
		putChloroplastStarch(getChloroplastStarch() - 0.14*reproduction);
		fruitCarbon[0] += reproduction;
		fruitNitrogen[0] += reproduction / CN;
		fromStem = getStemNSC()/nsc;
		stemNSC[0] -= fromStem*reproduction;
		leafNSC[0] -= (1.0-fromStem)*reproduction;
		leafStoredNitrogen[0] -= reproduction / CN;
		leafCfraction *= 0.01;
		rgrowth *= 0.2;
	}
//--------------------------------------

//total NSC available
	nsc = getLeafNSC()+getStemNSC()+getRootNSC();

//prevent growth allocation of bringing NSC to zero
//increase or decrease growth by a nonlinear function of excess or deficit NSC
	if (rgrowth > 0.95*nsc)
	{
		rgrowth = 0.95*nsc;
	}

//lifespan in units of 30 minutes
//compute unstressed change in LAI
//Tardieu et al 1999 New Phytologist: SLA_instant = (dA/dt) / (An * p)
	double SLA_numerator, SLA_denominator;
        SLA_numerator = 0.4*rgrowth*treesParams.SLA_max/10000.0 - treesParams.lai_at_sat_kl/unstressedLeafLifeSpan;
	SLA_denominator = 0.0001 * leafCfraction * 0.5*(psn + getChloroplastStarch());
	SLA_instant = SLA_numerator / SLA_denominator;
        if (SLA_instant < treesParams.SLA_min)
	{
		 SLA_instant = treesParams.SLA_min;
	}
        if (SLA_instant > treesParams.SLA_max) 
	{
		SLA_instant = treesParams.SLA_max;
	}
	treesParams.SLA_instant = SLA_instant;

//compute unstressed change in LAI
        deltaLAI = leafCfraction*rgrowth*SLA_instant/10000.0 - treesParams.lai_at_sat_kl/unstressedLeafLifeSpan;
        treesParams.lai_at_sat_kl += deltaLAI;

	if (kratio < 0.001)
	{
		kratio = 0.001;
	}
//set the leaf stress level if it has K < 50% Ksat and stress is higher in leaf than in shallow root
        if (newc[1][1] > newc[3][1])
        {
                stressedLeafLifeSpan = kratio*unstressedLeafLifeSpan;
        }
        else
        {
                stressedLeafLifeSpan = unstressedLeafLifeSpan;
        }

//compute stressed change in LAI
        deltaLAI = leafCfraction*rgrowth*SLA_instant/10000.0 - treesParams.live_lai/stressedLeafLifeSpan;
        treesParams.live_lai += deltaLAI;

        if (treesParams.lai < 0.01)
        {
                treesParams.lai = 0.01;
        }
        if (treesParams.lai_at_sat_kl < 0.01)
        {
                treesParams.lai_at_sat_kl = 0.01;
        }
        if (treesParams.live_lai < 0.01)
        {
                treesParams.live_lai = 0.01;
        }
}

//
//updateLeafCarbonNitrogenPools()
//
void BiogeochemicalCycles::updateLeafCarbonNitrogenPools(int k,
							 trees_params& treesParams,
						 	 double delta_lai,
						 	 double RL,
							 double N_neg_fract,
							 double& N_neg_demand,
							 double& N_pos_demand)
{
	double CN, nsc, delta_nsc, deltaBiomassC, deltaBiomassN, GRcost, SLA_instant, deltaAreaL;
	double total_N_demand, residual_N_demand, NfromStorage;
	double N_avail_rate_plant;

	N_avail_rate_plant = plantNstatus[0];
	nsc = getLeafNSC();

//when leaves are expanding
//growth respiration is costs are taken from chloroplast sugar and/or starch
//when chloroplast sugar and starch are limiting expansion respiration is
//supported from stored reserves, but leaf growth is slower
//For leaf we convert from the kg N m-2 leaf * lai of m2 m-2 * 10^4 m2 ha-1
	if (delta_lai > 0.0)
	{
		delta_nsc = RL * delta_lai / treesParams.SLA * 10000.0 * (1.0/0.86);

//this function adjusts the leaf C/N ratio as a function of available N
//DSM - July 2020
		CN = treesParams.leaf_maxCN - (treesParams.leaf_maxCN-treesParams.leaf_minCN)*N_avail_rate_plant;

		if (delta_nsc > 0.99*nsc)
		{
			delta_nsc = 0.99*nsc;
		}
		GRcost = delta_nsc * 0.14;
		total_N_demand = (delta_nsc-GRcost) / CN;
		if (total_N_demand > 0.99*leafStoredNitrogen[0])
		{
			NfromStorage = 0.99*leafStoredNitrogen[0];
			delta_nsc = NfromStorage * CN * (1.0/0.86);
                }
		GRcost = delta_nsc * 0.14;
		residual_N_demand = total_N_demand;
		N_neg_demand += N_neg_fract*residual_N_demand;
		N_pos_demand += (1.0-N_neg_fract)*residual_N_demand;

//Update leaf and whole canopy carbon
		deltaBiomassC = delta_nsc-GRcost;
		leafBiomassCarbon[0] += deltaBiomassC;
		SingleLeafBiomassCarbon[k] += deltaBiomassC; // structural C, in kgC
		leafNSC[0] -= delta_nsc;

//Update leaf and whole canopy nitrogen
		deltaBiomassN = deltaBiomassC / CN;
		leafBiomassNitrogen[0] += deltaBiomassN;
		SingleLeafBiomassNitrogen[k] += deltaBiomassN; // structural N, in kgC
		leafStoredNitrogen[0] -= deltaBiomassN;

		SLA_instant = calcSLA(plantNstatus[0], treesParams.SLA_max, treesParams.SLA_min);
		delta_lai = deltaBiomassC / RL / 10000.0 * treesParams.SLA;
		deltaAreaL = delta_lai * treesParams.pot_size * SLA_instant / treesParams.SLA;
		SingleLeafArea[k] = SingleLeafArea[k] + deltaAreaL;
		treesParams.SLA = (treesParams.SLA * leafBiomassCarbon[0] + SLA_instant * deltaBiomassC) / (leafBiomassCarbon[0] + deltaBiomassC);
	}
}

void BiogeochemicalCycles::updateLeafCarbonNitrogenPools(trees_params& treesParams,
						 	 double delta_lai,
						 	 double RL,
							 double N_neg_fract,
							 double& N_neg_demand,
							 double& N_pos_demand)
{
	double CN, SLA_instant, minLeafC, baseLAI, nsc, delta_nsc, GRcost;
	double total_N_demand, residual_N_demand, NfromStorage;
	double N_avail_rate_plant = plantNstatus[0];

//when delta_lai is positive, then increase leaf biomass
//when delta_lai is negative, then retranslate 50% of leaf protein N and reduce biomass

//base LAI is what you set in the parameter file for lai
	baseLAI = treesParams.Al;
	nsc = getLeafNSC();

//set leaf carbon minimum for the leaf bud to equivalent of 0.1 m2 m-2 * input parameter LAI
	minLeafC = 0.1*baseLAI/treesParams.SLA*10000.0;
		
//
//when leaves are expanding
//growth respiration is costs are taken from chloroplast sugar and/or starch
//when chloroplast sugar and starch are limiting expansion respiration is
//supported from stored reserves, but leaf growth is slower
//For leaf we convert from the kg N m-2 leaf * lai of m2 m-2 * 10^4 m2 ha-1
//
	if (delta_lai > 0.0)
	{
//
//Instantaneous SLA for the current leaf expansion
//This should depend on factors that affect leaf expansion rate, photosynthesis, and carbon allocation to leaf
//
		SLA_instant = treesParams.SLA_instant;

		//delta_nsc = RL * delta_lai / treesParams.SLA * 10000.0 * (1.0/0.86);
		delta_nsc = delta_lai / SLA_instant * 10000.0 * (1.0/0.86);

//this function adjusts the leaf C/N ratio as a function of available N
//DSM - July 2020
		CN = treesParams.leaf_maxCN - (treesParams.leaf_maxCN-treesParams.leaf_minCN)*N_avail_rate_plant;

		if (delta_nsc > 0.99*nsc)
		{
			delta_nsc = 0.99*nsc;
		}
		GRcost = delta_nsc * 0.14;
		total_N_demand = (delta_nsc-GRcost) / CN;
		if (total_N_demand > 0.99*leafStoredNitrogen[0])
		{
			NfromStorage = 0.99*leafStoredNitrogen[0];
			delta_nsc = NfromStorage * CN * (1.0/0.86);
			GRcost = delta_nsc * 0.14;
			total_N_demand = (delta_nsc-GRcost) / CN;
                }
		residual_N_demand = total_N_demand;
		N_neg_demand += N_neg_fract*residual_N_demand;
		N_pos_demand += (1.0-N_neg_fract)*residual_N_demand;

		treesParams.SLA = (treesParams.SLA * leafBiomassCarbon[0] + SLA_instant * (delta_nsc-GRcost)) / (leafBiomassCarbon[0] + delta_nsc-GRcost);
		leafBiomassCarbon[0] += (delta_nsc-GRcost);
		leafNSC[0] -= delta_nsc;
		leafBiomassNitrogen[0] += (delta_nsc-GRcost)/CN;
		leafStoredNitrogen[0] -= (delta_nsc-GRcost)/CN;

	}
//when leaves are senescing and leaf carbon exceeds the leaf bud carbon target
	else if (delta_lai < 0.0 && leafBiomassCarbon[0] > minLeafC)
	{
		delta_nsc = delta_lai / treesParams.SLA * 10000.0;
		CN = leafBiomassCarbon[0]/leafBiomassNitrogen[0];
		leafBiomassCarbon[0] += delta_nsc;
//retain biomass C and N for next budburst
		if (leafBiomassCarbon[0] < minLeafC)
		{
			delta_nsc += (minLeafC-leafBiomassCarbon[0]);
			leafBiomassCarbon[0] = minLeafC;
		}
//0.5 denotes 50% retranslocation 
//apportioned to leaf bud and reserves
		leafBiomassNitrogen[0] += 0.5*delta_nsc/CN;
		leafBiomassNitrogen[0] += 0.25/baseLAI*delta_nsc/CN;
		leafStoredNitrogen[0] -= (0.5-0.25/baseLAI)*delta_nsc/CN;
		leafResidueCarbon[0] -= delta_nsc;
		leafResidueNitrogen[0] -= 0.5*delta_nsc/CN;
	}
}
	
//
//computeLeafNdemand()
//determine how much nitrogen is needed to support leaf growth of delta_lai
//assume that leaf expansion rates require nitrogen to be taken from storage
//any residual nitrogen needed will have to be supplied from belowground
//residual used to compute nitrogen demand
//
void BiogeochemicalCycles::computeLeafNdemand(trees_params treesParams,
						double& delta_lai,
						double N_neg_fract,
						double& N_neg_demand,
						double& N_pos_demand)
{
/*
	double leafCN, total_N_demand, residual_N_demand, NfromStorage;
	
	double N_avail_rate_plant = plantNstatus[0];

//this function adjusts the leaf C/N ratio as a function of available N
//allows C/N to vary from 35-146 (high-low N)
//	leafCN = 22.0/max(8.0,treesParams.SLA)*35.0 + (1.0-N_avail_rate_plant)*50.0;
//DSM - July 2020
	leafCN = treesParams.leaf_maxCN - (treesParams.leaf_maxCN - treesParams.leaf_minCN)*N_avail_rate_plant;
	total_N_demand = delta_lai / treesParams.SLA * 10000.0 / leafCN;
	NfromStorage = total_N_demand * plantNstatus[0];
	if (NfromStorage > 0.99*leafStoredNitrogen[0])
	{
		NfromStorage = 0.99*leafStoredNitrogen[0];
		delta_lai = NfromStorage * leafCN * treesParams.SLA / 10000.0;
	}
	residual_N_demand = total_N_demand;
	N_neg_demand += N_neg_fract*residual_N_demand;
	N_pos_demand += (1.0-N_neg_fract)*residual_N_demand;
*/
}

//
//updateStemCarbonPools()
//
void BiogeochemicalCycles::updateStemCarbonNitrogenPools(trees_params treesParams,
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
						 	 double SLA)
{
	double liveStemIncrement, deadStemIncrement, residueIncrement;
	double tgrowth, CN, stemCincrement, stem_Cfraction_upperLimit;

//increase the rate of live stem death linearly with temperature above 5 C until root warms to 25 C
//reduce root growth at temperatures higher than 25 C
	tgrowth = 0.0;
	if (t_canopy > 5.0)
	{
		tgrowth = (t_canopy-5.0)/20.0;
	}
//Live stem wood mortality
	CN = liveStemCarbon[0]/liveStemNitrogen[0];
	CN = max(CN, 20.0);
	liveStemIncrement = tgrowth*liveStemCarbon[0]/(5.0*lifeSpan);
	liveStemCarbon[0] -= liveStemIncrement;
	liveStemNitrogen[0] -= liveStemIncrement/CN;
	
//Dead stem wood
	deadStemIncrement = liveStemIncrement;
	deadStemCarbon[0] += deadStemIncrement;
	deadStemNitrogen[0] += deadStemIncrement/CN;

//Stem wood residue
	CN = deadStemCarbon[0]/deadStemNitrogen[0];
	CN = max(CN, 20.0);
	residueIncrement = deadStemCarbon[0]*0.01/365.25/48.0;
	stemResidueCarbon[0] += residueIncrement;
	stemResidueNitrogen[0] += residueIncrement/CN;
	deadStemCarbon[0] -= residueIncrement;
	deadStemNitrogen[0] -= residueIncrement/CN;

//increase the rate of growth linearly with temperature above 5 C until root warms to 25 C
//reduce root growth at temperatures higher than 25 C
	tgrowth = 0.0;
	if (t_canopy > 5.0)
	{
		tgrowth = (t_canopy-5.0)/20.0;
		if (tgrowth > 1.0)
		{
			tgrowth = 1.0/tgrowth;
		}
	}
//New stem growth
//when N availability is less than 50% then halt stem growth
	if (N_avail_rate_plant > 0.5 && nscRatio > 1.0)
	{
		if(treesParams.useLeafModule == 1)
		{
			stem_Cfraction_upperLimit = (1.0/(1.0+treesParams.root_to_shoot ))*(1.0/(1.0+treesParams.leaf_to_stem));
		}
		else
		{
//if leaf allocation is shutting down, then shut down stem allocation as well
			stem_Cfraction_upperLimit = 0.2 * leafCfraction / 0.4;
		}
		stemAllocation = stem_Cfraction_upperLimit * (N_avail_rate_plant-0.5)/0.5 * kratio * tgrowth;
	}
	else
	{
		stemAllocation = 0.0;
	}
	stemCincrement = stemAllocation*rgrowth;

//check to make sure there is enough NSC for stem growth
	if (stemCincrement > 0.99*stemNSC[0])
	{
		stemCincrement = 0.99*stemNSC[0];
	}
	stemNSC[0] -= stemCincrement;

	CN = 50.0 + (1.0-N_avail_rate_plant)*50.0;
	N_neg_demand += N_neg_fract*stemCincrement*0.86/CN;
        N_pos_demand += (1.0-N_neg_fract)*stemCincrement*0.86/CN;

	liveStemCarbon[0] += 0.86*stemCincrement;
	liveStemNitrogen[0] += 0.86*stemCincrement/CN;
	leafStoredNitrogen[0] -= 0.86*stemCincrement/CN;
	if (liveStemCarbon[0] < 0.01)
	{
		liveStemCarbon[0] = 0.01;
		liveStemNitrogen[0] = liveStemCarbon[0]/CN;
	}
	if (rgrowth > 0.0)
        {       
                stemAllocation = stemCincrement/rgrowth;
        }
	else
	{
		stemAllocation = 0.0;
	}
}

//
//updateRootCarbonNitrogenPools()
//set root allocation when N is not limiting to decline with root size;
//    start with 0.19, decrease by 0.02 per increment
//set root allocation when N is limiting to shift more C to finer roots
//    start with an average of 0.19 and 0.5, then 0.17 and 0.25, 0.15 and 0.125, 0.13 and 0.125
//    for 10 root classes this sums to 1.0
//DSM - July 2015
//
void BiogeochemicalCycles::updateRootCarbonNitrogenPools(trees_params treesParams,
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
						 	 double stemAllocation)
{
	double fineRootLow, fineRootHigh, rootLifeSpan, refLifeSpan;
	double referenceSLA, rootScalar, SLAscalar, rootAllocation, rootCincrement, rootCdecrement, residual;
	double CN, tgrowth, kratio_sum, root_relative_growth;

//assumption here is that plants with lower SLA have higher C:N ratios in roots
	referenceSLA = 22.0;
        SLAscalar = referenceSLA / treesParams.SLA;
	if (SLAscalar > 3.0)
	{
		SLAscalar = 3.0;
	}
	else if (SLAscalar < 1.0)
	{
		SLAscalar = 1.0;
	}
	refLifeSpan = treesParams.minRootLifespan*365.25*48.0; //units of 30 min

//sum of dimensionless root K values, each given by K / Ksat
	kratio_sum = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
		if (treesParams.drlat[j+3] < 0.00001)
		{
			kratio_vector[j] = 0.0;
		}
		kratio_sum += kratio_vector[j];
	}
//
//Cycle through all soil-root layers
//
        for (int j = 0; j < nRoots; j++)
        {
//allow for increased C allocation to finer roots at low N
//note that 0.19 + 0.17 +...0.01 = 1.0
//   and 0.5 + 0.25 +...0.00097656 ~= 1.0
        	fineRootLow = 0.19; //allocation to fine roots at optimal available N
        	fineRootHigh = 0.5;
//adjust root lifespan upward as the finest root C:N ratio increases
//root lifespan will also increase with diameter
//differences will be captured between root depths based on C:N
//e.g., McCormack et al 2012. Predicting fine root lifespan from plant functional traits in
//		temperate trees. New Phytologist, 195, 823-831.
		CN = fineRootBiomassCarbon[j][0] / (fineRootBiomassNitrogen[j][0] + 0.000001);
		if (CN > treesParams.fr_maxCN)
		{
			CN = treesParams.fr_maxCN;
		}
                rootLifeSpan = refLifeSpan*CN/20.0;
		assert(rootLifeSpan > 0.000001);
//increase the rate of growth linearly with root temperature above 5 C until root warms to 25 C
//reduce root growth at temperatures higher than 25 C
		tgrowth = 0.0;
		if (tempSoil[j] > 5.0)
		{
			tgrowth = (tempSoil[j]-5.0)/20.0;
			if (tgrowth > 1.0)
			{
				tgrowth = 1.0;
			}
		}
//compute relative allocation of carbon to the respective root-soil layers
		root_relative_growth = kratio_vector[j] / (kratio_sum+0.0000001);
//
//fine root dynamics
//
		rootScalar = 1.0;
                for (int k = 0; k < nFineRootOrders; k++)
                {
//compute root mortality
			CN = fineRootBiomassCarbon[j][k] / fineRootBiomassNitrogen[j][k];
			rootCdecrement = tgrowth*fineRootBiomassCarbon[j][k]/(rootLifeSpan+0.000001);
                        fineRootBiomassCarbon[j][k] -= rootCdecrement;
			fineRootBiomassNitrogen[j][k] -= rootCdecrement/CN;
//change in residue for updating rhizosphere additions (ADD)
			dCrootResidue[j][k] = rootCdecrement;
			dNrootResidue[j][k] = rootCdecrement/CN;
//add dead roots to residue, which will be released over time to the litter C and N
			rootResidueCarbon[j][k] += rootCdecrement;
			rootResidueNitrogen[j][k] += rootCdecrement/CN;

//compute new root growth
//rootAllocation is a weighting of optimal N and stress N allocation fractions
                	rootAllocation = N_avail_rate_plant*fineRootLow;
                        rootAllocation += (1.0-N_avail_rate_plant)*fineRootHigh;
//rootAllocation is further weighted based on hydraulic stress
			rootAllocation *= root_relative_growth;
//convert root NSC and N to root biomass
//CN ratio is increased when the plant has low N reserves
//DSM - July 2020
			CN = treesParams.fr_maxCN - (treesParams.fr_maxCN - treesParams.fr_minCN)*N_avail_rate_plant;
			CN *= rootScalar;

			rootCincrement = (1.0-leafCfraction-stemAllocation)*rgrowth*tgrowth*rootAllocation/0.86;
if (isnan(rootCincrement))
{
cout << "BiogeochemicalCycles::updateRootCarbonNitrogenPools():" << endl;
cout << leafCfraction << '\t' << stemAllocation << '\t' << rgrowth << '\t' << tgrowth << '\t' << N_avail_rate_plant << '\t' << rootAllocation << '\t' << CN << endl;
exit(1);
}
//fineRootBiomassCarbon[j][k] > maximumRootBiomassCarbon[j][k]

			if (fineRootBiomassCarbon[j][k]+0.86*rootCincrement > maximumRootBiomassCarbon[j][k])
			{
				rootCincrement = (maximumRootBiomassCarbon[j][k]-fineRootBiomassCarbon[j][k]);
			}
			if (rootCincrement/CN > 0.99 * rootMineralNitrogen[j][k])
			{
				residual = rootCincrement - 0.99 * CN * rootMineralNitrogen[j][k];
				residual *= min(1.0,leafStoredNitrogen[0]/leafBiomassNitrogen[0]*(1.0-leafCfraction));
				rootCincrement = 0.99 * CN * rootMineralNitrogen[j][k] + residual;
				rootCincrement /= 0.86;
				leafStoredNitrogen[0] -= residual/CN;
				rootMineralNitrogen[j][k] *= 0.01;
			}
			else
			{
				rootMineralNitrogen[j][k] -= rootCincrement/CN;
			}
			if (rootCincrement > 0.99*rootNSC[j][k])
			{
				rootCincrement = 0.99*rootNSC[j][k];
			}

			rootNSC[j][k] -= rootCincrement;
                        fineRootBiomassCarbon[j][k] += 0.86*rootCincrement;
                        fineRootBiomassNitrogen[j][k] += 0.86*rootCincrement/CN;
			N_neg_demand += 0.86*N_neg_fract*rootCincrement/CN;
			N_pos_demand += 0.86*(1.0-N_neg_fract)*rootCincrement/CN;
//increase root lifespan at sqrt(2) with each diameter doubling
                        rootLifeSpan *= treesParams.rootDiamMultiplier; 
                        fineRootLow -= 0.02;
                        fineRootHigh *= 0.5;
			rootScalar *= 1.25;
                }
//
//coarse root dynamics
//
                for (int k = nFineRootOrders; k < nRootOrders; k++)
                {
			CN = coarseRootBiomassCarbon[j][k]/coarseRootBiomassNitrogen[j][k];
//kill off roots
			rootCdecrement = tgrowth*coarseRootBiomassCarbon[j][k]/(rootLifeSpan+0.000001);
                        coarseRootBiomassCarbon[j][k] -= rootCdecrement;
                        coarseRootBiomassNitrogen[j][k] -= rootCdecrement/CN;
//change in residue for updating rhizosphere additions (ADD)
			dCrootResidue[j][k] = rootCdecrement;
			dNrootResidue[j][k] = rootCdecrement/CN;
//add dead roots to residue
			rootResidueCarbon[j][k] += rootCdecrement;
			rootResidueNitrogen[j][k] += rootCdecrement/CN;

//compute new root growth
//rootAllocation is a weighting of optimal N and stress N allocation fractions
                        rootAllocation = N_avail_rate_plant*fineRootLow;
                        rootAllocation += (1.0-N_avail_rate_plant)*fineRootHigh;
//rootAllocation is further weighted based on hydraulic stress
			rootAllocation *= root_relative_growth;
//convert root NSC and N to root biomass
//CN ratio is increased when the plant has low N reserves
//DSM - July 2020
			CN = treesParams.fr_maxCN - (treesParams.fr_maxCN - treesParams.fr_minCN)*N_avail_rate_plant;
			CN *= rootScalar;

			rootCincrement = (1.0-leafCfraction-stemAllocation)*rgrowth*tgrowth*rootAllocation;
			if (coarseRootBiomassCarbon[j][k]+0.86*rootCincrement > maximumRootBiomassCarbon[j][k])
			{
				rootCincrement = (maximumRootBiomassCarbon[j][k]-coarseRootBiomassCarbon[j][k])/0.86;
			}
			if (rootCincrement > 0.99*rootNSC[j][k])
			{
				rootCincrement = 0.99*rootNSC[j][k];
			}
			if (rootCincrement/CN > 0.99 * rootMineralNitrogen[j][k])
			{
				//rootCincrement = 0.99 * CN * rootMineralNitrogen[j][k];
				//rootCincrement /= 0.86;
				residual = rootCincrement - 0.99 * CN * rootMineralNitrogen[j][k];
				residual *= min(1.0,leafStoredNitrogen[0]/leafBiomassNitrogen[0]*(1.0-leafCfraction));
				rootCincrement = 0.99 * CN * rootMineralNitrogen[j][k] + residual;
				rootCincrement /= 0.86;
				leafStoredNitrogen[0] -= residual/CN;
				rootMineralNitrogen[j][k] *= 0.01;
			}
			else
			{
				rootMineralNitrogen[j][k] -= rootCincrement/CN;
			}
			rootNSC[j][k] -= rootCincrement;
                        coarseRootBiomassCarbon[j][k] += 0.86*rootCincrement;
                        coarseRootBiomassNitrogen[j][k] += 0.86*rootCincrement/CN;
			//rootMineralNitrogen[j][k] -= 0.86*rootCincrement/CN;
			N_neg_demand += 0.86*N_neg_fract*rootCincrement/CN;
			N_pos_demand += 0.86*(1.0-N_neg_fract)*rootCincrement/CN;
//roots double in lifespan with diameter doubling
                        rootLifeSpan *= treesParams.rootDiamMultiplier; 
                        fineRootLow -= 0.02;
                        fineRootHigh *= 0.5;
                }
	}
}

//
//The following is for root maintenance respiration
//   output units of kgC 30min-1
//
double BiogeochemicalCycles::root_respiration_rate(double resp_coef_root,
                                		   double resp_coefficient,
                                		   double t_soil,
                                		   double Croot,
                                		   double transport_rate)
{
        double rate;

transport_rate = 1.0;

        rate = transport_rate * resp_coef_root * exp(resp_coefficient*t_soil) * Croot;

        if (rate < 0.0)
        {
                rate = 0.0;
        }
        return (rate);
}

//
//The following is for stem maintenance respiration
//   output units of kgC 30min-1
//
double BiogeochemicalCycles::stem_respiration_rate(double resp_coef_stem,
                                		   double resp_coefficient,
                                		   double t_canopy,
                                		   double Cstem)
{
        double rate;

        rate = resp_coef_stem * exp(resp_coefficient*t_canopy) * exp(0.67*log(Cstem*10.0))/10.0;

        if (rate < 0.0)
        {
                rate = 0.0;
        }
        return (rate);
}

//
//The following is for leaf maintenance respiration
//   output units of kgC 30min-1
//
double BiogeochemicalCycles::leaf_respiration_rate(double resp_coef_leaf,
                                		   double resp_coefficient,
                                		   double tLeaf,
                                		   double lai,
                                		   double SLA)
{
        double rate;

        rate = resp_coef_leaf * exp(resp_coefficient*tLeaf) * lai / SLA;

        if (rate < 0.0)
        {
                rate = 0.0;
        }
        return (rate);
}

//
//getRhizosphereN()
//
void BiogeochemicalCycles::getRhizosphereN(double& N_neg,
				           double& N_pos)
{
	N_neg = N_pos = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        N_neg += getRhizosphereN_neg(j, k);
                        N_pos += getRhizosphereN_pos(j, k);
                }
        }
}
void BiogeochemicalCycles::getRhizosphereN(double& N_neg,
                                           double& N_pos,
					   int j)
{
        N_neg = N_pos = 0.0;
        for (int k = 0; k < nRootOrders; k++)
        {
                N_neg += getRhizosphereN_neg(j, k);
                N_pos += getRhizosphereN_pos(j, k);
        }
}

//
//getRhizosphereN_neg()
//
double BiogeochemicalCycles::getRhizosphereN_neg(int j,
                                        	 int k)
{
        double Nr;
        Nr = rhizosphereNitrateNitrogen[j][k];
        return(Nr);
}

//
//getRhizosphereN_pos()
//
double BiogeochemicalCycles::getRhizosphereN_pos(int j,
                                        	 int k)
{
        double Nr;
        Nr = rhizosphereAmmoniumNitrogen[j][k];
        return(Nr);
}

//
//getPlantN()
//
double BiogeochemicalCycles::getPlantN()
{
        double plantN;

	plantN = leafStoredNitrogen[0];
        for (int j = 0; j < nRoots; j++)
        {
                for (int k = 0; k < nRootOrders; k++)
                {
                        plantN += rootMineralNitrogen[j][k];
                }
        }
	return(plantN);
}

//
//getRootN()
//
double BiogeochemicalCycles::getRootN()
{
        double rootN = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                rootN += getRootN(j);
        }
	return(rootN);
}
double BiogeochemicalCycles::getRootN(int j)
{
        double rootN = 0.0;
        for (int k = 0; k < nRootOrders; k++)
        {
        	rootN += getRootN(j, k);
        }
        return(rootN);
}
double BiogeochemicalCycles::getRootN(int j,
				      int k)
{
        double rootN;
        rootN = rootMineralNitrogen[j][k];
        return(rootN);
}

//
//getFineRootN()
//
double BiogeochemicalCycles::getFineRootN()
{
        double rootN = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
		rootN += getFineRootN(j);
	}
	return(rootN);
}
double BiogeochemicalCycles::getFineRootN(int j)
{
        double rootN = 0.0;
        for (int k = 0; k < nFineRootOrders; k++)
        {
        	rootN += getFineRootN(j, k);
        }
        return(rootN);
}
double BiogeochemicalCycles::getFineRootN(int j,
					  int k)
{
        double rootN;
        rootN = rootMineralNitrogen[j][k];
        return(rootN);
}

//
//getRootBiomassN()
//
double BiogeochemicalCycles::getRootBiomassN()
{
        double rootN = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                rootN += getRootBiomassN(j);
        }
	return(rootN);
}
double BiogeochemicalCycles::getRootBiomassN(int j)
{
        double rootN = 0.0;
        for (int k = 0; k < nRootOrders; k++)
        {
                rootN += getRootBiomassN(j, k);
        }
        return(rootN);
}
double BiogeochemicalCycles::getRootBiomassN(int j,
					     int k)
{
        double rootN;
        rootN = fineRootBiomassNitrogen[j][k];
        rootN += coarseRootBiomassNitrogen[j][k];
        return(rootN);
}

//
//getFineRootBiomassN()
//
double BiogeochemicalCycles::getFineRootBiomassN()
{
        double rootN = 0.0;
        for (int j = 0; j < nRoots; j++)
        {
                rootN += getFineRootBiomassN(j);
        }
	return(rootN);
}
double BiogeochemicalCycles::getFineRootBiomassN(int j)
{
        double rootN = 0.0;
        for (int k = 0; k < nFineRootOrders; k++)
        {
                rootN += getFineRootBiomassN(j, k);
        }
        return(rootN);
}
double BiogeochemicalCycles::getFineRootBiomassN(int j,
						 int k)
{
        double rootN;
        rootN = fineRootBiomassNitrogen[j][k];
        return(rootN);
}

//
//computeRootNitrogenUptake()
//based on Porporato et al 2003 Advances in Water Resources
//nitrogen uptake from the rhizosphere volume by the root
//passive N uptake is proportional to rhizosphere flux into root
//active N uptake is proportion to ion gradient, for now assuming zero concentration inside root
//rhizopshere flux is already a proportion of canopy transpiration
//proportional to root area / total root area
// Pools are in kg ha-1
// Convert to g m-2 -> Pools * 1000 / 10000
// Water flux is m3 30 min-1
// Root area is per unit ground area
//
void BiogeochemicalCycles::computeRootNitrogenUptake(double UP_neg[][10],
						     double UP_pos[][10],
						     trees_params treesParams,
						     double* thetaSoil,
						     double* Rflux,
						     double ratot,
						     double field_capacity,
						     double DEM_neg,
						     double DEM_pos)
{
	double UP_neg_p, UP_neg_a, UP_pos_p, UP_pos_a, UPp, UPa, DEM_neg_A, DEM_pos_A;
	double F, DEMfract, rootLength, rhizVolume, rhizWaterVolume, kuN, concentration;
	double a, rhizWidth, rootRadius, fluxM, rootFract, absorbFract, totRoot;
	double rootC, totFineRootC, kC, kN, CostUP, CN;
	double DEM_tot, passiveUptake, activeUptake, uptakeRatio;
	double porosity, leafToRootCratio;
	double passiveNitrogenDemand, activeNitrogenDemand;
	double referenceSLA = 22.0;
        double refLifeSpan = treesParams.minRootLifespan; //years

	porosity = treesParams.porosity;

//
//Account for genes that regulate N uptake
//For active uptake, a sigmoidal function is used to represent the distribution of 
//  high and low affinity transporters. For passive uptake, assume a regulation of
//  aquiporin expression is negatively related to plant nitrogen status
//
	activeNitrogenDemand = 1.0+ 9.0/(1.0+exp(10.0*plantNstatus[0]-5.0));
	passiveNitrogenDemand = (1.0 - plantNstatus[0]);;
	
//S.A. Barber, 1995. Soil Nutrient Bioavailability: A mechanistic approach. i
//   John Wiley & Sons, New York
//   Table 4.5 on page 93
//   F = 5 x 10-7 cm2 s-1 * 1800 s 30min-1 * 0.0001 m2 cm-2 = 9 x 10-8 m2 30min-1
//   	F = 9.0e-8; // m2 30min-1;
//
//Jackson 1996, F = 2 x 10-6 cm2 s-1 * 1800 s 30 min-1 * 0.0001 m2 cm-2 = 3.6 x 10-7 m2 30min-1
   	F = 3.6e-7 * activeNitrogenDemand; // m2 30min-1;

//rhizosphere width
	rhizWidth = treesParams.rhizosphere_width/1000.0;

//For all soil-root zones
	UP_neg_p = UP_neg_a = UP_pos_p = UP_pos_a = 0.0;
	for (int j = 0; j < nRoots; j++)
	{
		totRoot = 0.0;
		for (int k = 0; k < nRootOrders; k++)
		{
			totRoot += rootArea[j][k];
		}

//convert rhizopshere flux from mmol m-2 s-1 units to m3 m-2 30min-1
		fluxM = 0.001*1800.0*EcConversion_StoT(Rflux[j], treesParams.lai);

//either assume nutrients are not allowed to be moved from root to rhizosphere
//or allow rhizodeposition of N
//pick your poison, commenting out this if{} will allow rhizodeposition
		if (fluxM < 0.0)
		{
			fluxM = 0.0;
		}
//For all rhizosphere zones
		rootRadius = 0.5 * treesParams.minRootDiam;
//assume 100% of finest root can absorb N, 50% of second finest, etc
		absorbFract = 1.0; 
		for (int k = 0; k < nRootOrders; k++)
		{
			rootFract = absorbFract*rootArea[j][k]/totRoot;
//rhizosphere water volume (m3 m-2)
			rhizVolume = computeRhizosphereVolume(j, k, porosity, rootRadius, rhizWidth);
			rootLength = rootArea[j][k]/(2.0*M_PI*rootRadius);

//how much liquid water is in the rhizosphere -- this is correct only when there is zero water uptake
			rhizWaterVolume = rhizVolume*thetaSoil[j]/porosity;
//
//NITRATE NITROGEN
//
//passive nitrate uptake
//conc = g m-2 / (m3 m-2) = g m-3 water
			a = 1.0; //for mobile nitrate ion
//N solubility declines at lower soil moisture
			if (thetaSoil[j] < field_capacity)
			{
				a *= sqrt(thetaSoil[j]/field_capacity+0.00000001);
			}
			a *= passiveNitrogenDemand;
			concentration = a * rhizosphereNitrateNitrogen[j][k] * 0.1 / (rhizWaterVolume+0.0000001);

//UPp = g m-3 x m3 30min-1 m-2 = g m-2 30min-1
			UPp = concentration * fluxM * rootFract;

//determine carbon cost of passive nitrate uptake
//downregulate uptake if NSC is not available
//From Cannell and Thornley, 2000, Annals of Botany
//Assume a respiratory cost of 1 mol NO3- N requires 0.4 mol C in glucose respired
//  In mass units, cost of NO3- uptake is 0.4 * (12/14) = 0.34g glucose per g NO3- N
//  For NH4+ uptake use 0.17 g glucose C per g NH4+ N taken up
//  These numbers assume use of NSC in form of C6H12O6
//Here is the cost of reducing nitrate to ammonium for assimilation
//  cost of nitrate reduction is 8 mol H (mol N)-1
//  glucose C cost of (8 x 6 x 12/24)/14 = 1.72 kg C kg-1 N reduced
			CostUP = 1.72*UPp*10.0;
			if (CostUP > 0.95*leafNSC[0] && CostUP > 0.0)
			{
				UPp *= 0.95*leafNSC[0]/CostUP;
				CostUP = 0.95*leafNSC[0];
				leafNSC[0] -= CostUP;
			}
			else if (CostUP > 0.0)
			{
				leafNSC[0] -= CostUP;
			}
			UP_neg_p += UPp*10.0;
			UP_neg[j][k] = UPp/rhizVolume;
			rootMineralNitrogen[j][k] += UPp*10.0;
//active nitrate uptake
//DEM_neg has units kg N ha-1
//convert to g m-2 30min-1
			DEM_neg_A = DEM_neg * 0.1;
			DEM_neg_A /= pow(plantNstatus[0]+0.0001, 2.0);
			DEMfract = DEM_neg_A*rootFract;
//diffusion coefficient F is m2 30min-1
//kuN = g m-3 x m-1 x m2 30min-1 = g m-2 30 min-1
//NEED CONCENTRATION GRADIENT (CONCENTRATION - 0)/HALF RHIZOSPHERE RADIUS
			kuN = concentration/(0.5*rhizWidth) * F * rootArea[j][k] * absorbFract;
			kuN *= pow(thetaSoil[j]/porosity, 3.0);
			UPa = 0.0;
			if (DEMfract < UPp)
			{
				UPa = 0.0;
			}
			else if (kuN < (DEMfract-UPp))
			{
				UPa = kuN;
			}
			else if ((kuN > (DEMfract-UPp)) && (DEMfract-UPp) > 0.0)
			{
				UPa = DEM_neg_A - UPp;
			}

//determine carbon cost of active nitrate uptake
//downregulate uptake if NSC is not available
//From Cannell and Thornley, 2000, Annals of Botany
//Assume a respiratory cost of 1 mol NO3- N requires 0.4 mol C in glucose respired
//  In mass units, cost of NO3- uptake is 0.4 * (12/14) = 0.34g glucose per g NO3- N
//  For NH4+ uptake use 0.17 g glucose C per g NH4+ N taken up
//  These numbers assume use of NSC in form of C6H12O6

			CostUP = 0.34*UPa*10.0;

//active N uptake is limited by available NSC in the root
//cout << "UPa = " << UPa << endl;
			if (CostUP > 0.95*rootNSC[j][k] && CostUP > 0.0) 
			{
				UPa *= rootNSC[j][k]/CostUP;
				CostUP = 0.95*rootNSC[j][k];
				rootNSC[j][k] -= CostUP;
			}
			else if (CostUP > 0.0)
			{
				rootNSC[j][k] -= CostUP;
			}
			UP_neg_a += UPa*10.0;
			UP_neg[j][k] += UPa/rhizVolume;

//Add in the cost of reducing nitrate to ammonium for assimilation
//  cost of nitrate reduction is 8 mol H (mol N)-1
//  glucose C cost of (8 x 6 x 12/24)/14 = 1.72 kg C kg-1 N reduced

			CostUP = 1.72*UPa*10.0;
			if (CostUP > 0.95*leafNSC[0] && CostUP > 0.0)
			{
				if ((CostUP - 0.95*leafNSC[0]) < 0.99*stemNSC[0])
				{
					CostUP -= 0.95*leafNSC[0];
					leafNSC[0] *= 0.05;
					stemNSC[0] -= CostUP;
				}
				else
				{
					leafNSC[0] *= 0.05;
					stemNSC[0] *= 0.01;
				}
			}
			else if (CostUP > 0.0)
			{
				leafNSC[0] -= CostUP;
			}

//move active N uptake to root and leaf in equal proportions
//convert from g m-3 to kg ha-1
			rootMineralNitrogen[j][k] += UPa*10.0; 

//
//AMMONIUM NITROGEN
//
//passive ammonium uptake
			a = 0.1; //assumes that most ammonium is absorbed in soil matrix
//N solubility declines at lower soil moisture
			if (thetaSoil[j] < field_capacity)
			{
				a *= sqrt(thetaSoil[j]/field_capacity+0.00000001);
			}
			a *= passiveNitrogenDemand;
			concentration = a * rhizosphereAmmoniumNitrogen[j][k] * 0.1 / (rhizWaterVolume+0.0000001);
			UPp = concentration * fluxM * rootFract;
			UP_pos_p += UPp*10.0;
			UP_pos[j][k] = UPp/rhizVolume;
			rootMineralNitrogen[j][k] += UPp*10.0;
//active ammonium uptake
			DEM_pos_A = DEM_pos * 0.1;
			DEM_pos_A /= pow(plantNstatus[0]+0.0001, 2.0);
			DEMfract = DEM_pos_A*rootFract;
//diffusion coefficient F is m2 30min-1
//kuN = g m-3 x m-1 x m2 30min-1 = g m-2 30 min-1
//NEED CONCENTRATION GRADIENT (CONCENTRATION - 0)/HALF RHIZOSPHERE RADIUS
			kuN = concentration/(0.5*rhizWidth) * F * rootArea[j][k] * absorbFract;
			kuN *= pow(thetaSoil[j]/porosity, 3.0);
			if (DEMfract < UPp)
			{
				UPa = 0.0;
			}
			else if (kuN < (DEMfract-UPp))
			{
				UPa = kuN;
			}
			else if ((kuN > (DEMfract-UPp)) && (DEMfract-UPp) > 0.0)
			{
				UPa = DEM_pos_A - UPp;
			}

//determine carbon cost of active ammonium uptake
//downregulate uptake if NSC is not available
//  For NH4+ uptake use 0.17 g glucose C per g NH4+ N taken up
//  These numbers assume use of NSC in form of C6H12O6

			CostUP = 0.17*UPa*10.0;

//active N uptake is limited by available NSC in the root
			if (CostUP > 0.95*rootNSC[j][k] && CostUP > 0.0) 
			{
				UPa *= 0.95*rootNSC[j][k]/CostUP;
				CostUP = 0.95*rootNSC[j][k];
				rootNSC[j][k] -= CostUP;
			}
			else if (CostUP > 0.0)
			{
				rootNSC[j][k] -= CostUP;
			}
			UP_pos_a += UPa*10.0;
			UP_pos[j][k] += UPa/rhizVolume;
//move active N uptake to root
//convert from g m-3 to kg ha-1
			rootMineralNitrogen[j][k] += UPa*10.0; 
			if (UP_neg[j][k] < 0.0)
			{
				UP_neg[j][k] = 0.0;
			}
			if (UP_pos[j][k] < 0)
			{
				UP_pos[j][k] = 0.0;
			}
//radius doubles with each root size class and fraction of root that takes up N halves
			rootRadius *= treesParams.rootDiamMultiplier;
			absorbFract *= 0.9;
		}
	}
//move passive N uptake to leaf
	passiveUptake = UP_neg_p + UP_pos_p;

//compute plant nitrogen status
	DEM_tot = DEM_neg + DEM_pos;

//proportion of leaf+root biomass in leaf
	totFineRootC = getFineRootCarbon();
	leafToRootCratio = leafBiomassCarbon[0] / (leafBiomassCarbon[0]+2.0*totFineRootC);

//
//Compute plantNstatus, with a range of 0 to 1.
//This assumes that available NSC determines an optimal amount of available N
	CN = getLeafBiomassCarbon() / getLeafBiomassN();
	//CN = treesParams.leaf_minCN;
	plantNstatus[0] = leafToRootCratio * CN * (leafStoredNitrogen[0] / (getLeafNSC()+0.000001));
	//plantNstatus[0] = min(leafToRootCratio, (leafStoredNitrogen[0] / getLeafNSC()) * treesParams.leaf_minCN);
	//plantNstatus[0] = min(leafToRootCratio, (leafStoredNitrogen[0] / getLeafNSC()) * CN);
	//plantNstatus[0] = leafToRootCratio * min(1.0, CN * (leafStoredNitrogen[0] / getLeafNSC()));
	CN = getFineRootCarbon() / getFineRootBiomassN();
	//CN = treesParams.fr_minCN;
	plantNstatus[0] += (1.0-leafToRootCratio) * CN * (getFineRootN() / (getFineRootNSC()+0.000001));
	//plantNstatus[0] += min((1.0-leafToRootCratio), (getFineRootN() / getFineRootNSC()) * treesParams.fr_minCN);
	//plantNstatus[0] += min((1.0-leafToRootCratio), (getFineRootN() / getFineRootNSC()) * CN);
	//plantNstatus[0] += (1.0-leafToRootCratio) * min(1.0, CN * (getFineRootN() / getFineRootNSC()));

	if (plantNstatus[0] > 1.0)
	{
		plantNstatus[0] = 1.0;
	}
	else if (plantNstatus[0] < 0.0001)
	{
		plantNstatus[0] = 0.0001;
	}
}

//
//computeLeaching()
//combine drainage proportion by rhizosphere and N concentration
//drain is in units of m 30 min-1
//
void BiogeochemicalCycles::computeLeaching(double LE_neg[][10],
		     			   double LE_pos[][10],
		     			   trees_params treesParams,
		     			   double* thetaSoil,
					   double* depthSoil,
		     			   double* bypassFlow)
{
	double rootRadius, rhizWidth, rhizVolume, totalRhizVolume[10], rhizWaterVolume, a, concentration;
	double porosity, leachate, rhizVolumeProportion, negLE[10], posLE[10];
	rhizWidth = 0.001*treesParams.rhizosphere_width;
	porosity = treesParams.porosity;
//For all soil-root layers
	for (int j = 0; j < nRoots; j++)
	{
		negLE[j] = 0.0;
		posLE[j] = 0.0;
	}
	for (int j = 0; j < nRoots; j++)
	{
		rootRadius = 0.5 * treesParams.minRootDiam;
		totalRhizVolume[j] = computeRhizosphereVolume(j, porosity, rootRadius, rhizWidth);
//For all rhizosphere zones
		for (int k = 0; k < nRootOrders; k++)
		{
//total radius is root radius plus rhizosphere radius
//rhizosphere water volume
			rhizVolume = computeRhizosphereVolume(j, k, porosity, rootRadius, rhizWidth);
			rhizVolumeProportion = rhizVolume / (totalRhizVolume[j]+0.0000001);
			rhizWaterVolume = rhizVolume*thetaSoil[j]/porosity;
//nitrate leaching
//units are g m-3
			a = 1.0; //for mobile nitrate ion
			concentration = 0.1 * a * rhizosphereNitrateNitrogen[j][k]/(rhizWaterVolume+0.0000001);
			LE_neg[j][k] = rhizVolumeProportion * bypassFlow[j] * concentration;
			negLE[j] += LE_neg[j][k];
//ammonium leaching
			a = 0.1; //assumes that most ammonium is absorbed in soil matrix
			concentration = 0.1 * a * rhizosphereAmmoniumNitrogen[j][k]/(rhizWaterVolume+0.0000001);
			LE_pos[j][k] = rhizVolumeProportion * bypassFlow[j] * concentration;
			posLE[j] += LE_pos[j][k];

			rootRadius *= treesParams.rootDiamMultiplier;
		}
//compute net advection of leachate
//negative results indicate net leachate gain; positive indicates net loss
//this assumes that advected leachate is well mixed
		leachate = 0.0;
		if (j > 0)
		{
			for (int k = 0; k < nRootOrders; k++)
			{
				LE_neg[j][k] -= 0.1*negLE[j-1]*totalRhizVolume[j]/depthSoil[j];
				LE_pos[j][k] -= 0.1*posLE[j-1]*totalRhizVolume[j]/depthSoil[j];
				leachate += LE_neg[j][k] + LE_pos[j][k];
			}
		}
//only net loss or no loss from the top root-soil zone
		else
		{
			for (int k = 0; k < nRootOrders; k++)
			{
				leachate += LE_neg[j][k] + LE_pos[j][k];
			}
		}
//we will output the net leaching in kgN ha-1 to be consistent with reporting
		nitrogenLeaching[j] = 10.0*leachate;
	}
}

//
//compute rhizosphere pore volume
//
double BiogeochemicalCycles::computeRhizosphereVolume(int rootNum,
						      double porosity,
						      double rootRadius,
						      double rhizosphereWidth)
{
	double rhizVolume;

	rhizVolume = 0.0;
	for (int k = 0; k < nRootOrders; k++)
	{
		rhizVolume += computeRhizosphereVolume(rootNum, k, porosity, rootRadius, rhizosphereWidth);
//treesParams.rootDiamMultiplier
		rootRadius *= 2.0;
	}
	return(rhizVolume);
}
double BiogeochemicalCycles::computeRhizosphereVolume(trees_params treesParams,
						      int rootNum)
{
	double rhizVolume, porosity, rootRadius, rhizosphereWidth, rootDiamMultiplier;

	rhizVolume = 0.0;
	porosity = treesParams.porosity;
	rhizosphereWidth = 0.001*treesParams.rhizosphere_width;
	rootDiamMultiplier = treesParams.rootDiamMultiplier;
	rootRadius = treesParams.minRootDiam;
	for (int k = 0; k < nRootOrders; k++)
	{
		rhizVolume += computeRhizosphereVolume(rootNum, k, porosity, rootRadius, rhizosphereWidth);
		rootRadius *= rootDiamMultiplier;
	}
	return (rhizVolume);
}
double BiogeochemicalCycles::computeRhizosphereVolume(trees_params treesParams,
						      int rootNum,
						      int rootOrder)
{
	double rhizVolume, porosity, rootRadius, rhizosphereWidth, rootDiamMultiplier;

	porosity = treesParams.porosity;
	rhizosphereWidth = 0.001*treesParams.rhizosphere_width;
	rootDiamMultiplier = treesParams.rootDiamMultiplier;
	rootRadius = treesParams.minRootDiam;
	for (int k = 0; k < nRootOrders; k++)
	{
		if (k == rootOrder)
		{
			rhizVolume = computeRhizosphereVolume(rootNum, k, porosity, rootRadius, rhizosphereWidth);
			break;
		}
		rootRadius *= rootDiamMultiplier;
	}
	return (rhizVolume);
}
double BiogeochemicalCycles::computeRhizosphereVolume(int rootNum,
						      int rootOrder,
						      double porosity,
						      double rootRadius,
						      double rhizosphereWidth)
{
	double totRadius, rootLength, rhizVolume;

	assert(rootNum >= 0);
	assert(rootNum < nRoots);
	assert(rootOrder >= 0);
	assert(rootOrder < nRootOrders);

	totRadius = rhizosphereWidth + rootRadius;
	rootLength = rootArea[rootNum][rootOrder]/(2.0*M_PI*rootRadius);
	rhizVolume = porosity*rootLength*M_PI*(totRadius*totRadius - rootRadius*rootRadius);
	return(rhizVolume);
}

//
//updateRhizospherePools()
//litter biomass pools
//  dCldt = ADD + BG - DECl
//  BD = kdCb
//  DEC = phi*fds*kl*Cl
//  dNldt = ADD/CNadd + BD/CNb - DECl/CNl
//humus pools
//  dChdt = rh*DECl - DECh
//  DECh = phi * fds * kh * Cb * Ch
//  dNhdt = rh * DECl/CNh - DECh/CNh
//biomass pool
//  dCbdt = (1-rh-rr)*DECl + (1-rr)*DECh - BD
//  dNbdt = (1-rh*CNl/CNh)*DECl/CNl + DECh/CNh - BD/CNb - PHI
//State variables have units of kg ha-1
//Internally these should be converted to g m-3 rhizosphere volume
//Rhizosphere volume given as m3 m-2
//AtoV: kg ha-1 * (1/10000) ha m-2 * 1000 g kg-1 / (m3 m-2) = g m-3
//
void BiogeochemicalCycles::updateRhizospherePools(trees_params treesParams,
						  double* thetaSoil,
						  double* tempSoil,
						  double* depthSoil,
						  double UP_neg[][10],
						  double UP_pos[][10],
						  double LE_neg[][10],
						  double LE_pos[][10])
{
	double PHI, phi, Cb, Ch, CNh, rr, CNb, Cl, CNl, rh, Cea, CNea, Ces;
	double ADD, CNadd, BD, DECl, DECh, DECea, DECes, dC, dN, dN2, dN3;
	double MIN, IMM_pos, IMM_neg;
	double NIT, fns, kd, kn;
	double theta, t_soil, depth, cumDepth, rootRadius, CN, bulkDensity;
	double dCldt, dNldt, dChdt, dNhdt, dCbdt, dNbdt, dNposdt, dNnegdt;
	double porosity, rhizVolume, rhizWidth, AtoV, AtoVbulk;
	double rootLifeSpan, refLifeSpan;
	double sumRespiration;
	double cumDist;
	double beta = 0.970;


	bulkDensity = treesParams.BD;
	rhizWidth = 0.001*treesParams.rhizosphere_width; //metres
	porosity = treesParams.porosity;
	cumDepth = 0.0;

//For all soil-root zones
	for (int j = 0; j < nRoots; j++)
	{
//
//adjust root lifespan upward as the finest root C:N ratio increases
//root lifespan will also increase with diameter
//differences will be captured between root depths based on C:N
//e.g., McCormack et al 2012. Predicting fine root lifespan from plant functional traits in
//              temperate trees. New Phytologist, 195, 823-831.
//
		refLifeSpan = treesParams.minRootLifespan*48.0*365.25;
                CN = fineRootBiomassCarbon[j][0] / fineRootBiomassNitrogen[j][0];
		CN = max(CN, 20.0);
                rootLifeSpan = refLifeSpan*CN/20.0;
                assert(rootLifeSpan > 0.000001);
		rootRadius = 0.5 * treesParams.minRootDiam; //meters
		theta = thetaSoil[j];
		t_soil = tempSoil[j];
		depth = depthSoil[j];
		cumDepth += depth;

//compute the dependency of decomposition rate on water
       		fns = nitrification_water_function(treesParams.theta_opt, theta);

//convert kg ha-1 basis to g m-3 bulk soil pore volume for layer i
		cumDist = 1.0-pow(beta, 100.0*(cumDepth-0.5*depth));
		if (cumDist > 0.95)
		{
			cumDist = 0.95;
		}
//cout << cumDist << endl;
		if (treesParams.usePhenology == true) //perennial plants
		{
			AtoVbulk = 0.1 / depth * bulkDensity / max(treesParams.ar[j+3], 1.0 - cumDist); 
		}
		else //annual plants
		{
			AtoVbulk = 0.1 / depth * bulkDensity;
		}
//humus C and C/N
//assume well-mixed humus throughout soil layer
		if (humusNitrogen[j] > 0.0)
		{
			Ch = humusCarbon[j];
			CNh = Ch / humusNitrogen[j];
			Ch *= AtoVbulk; //convert to g m-3 soil volume
		}
		else
		{
			Ch = 0.0;
			CNh = 1.0;
		}

		sumRespiration = 0.0;

//For all rhizosphere zones
		for (int k = 0; k < nRootOrders; k++)
		{
			rhizVolume = computeRhizosphereVolume(j, k, porosity, rootRadius, rhizWidth);
			if (rhizVolume < 0.000001)
			{
				rhizVolume = 0.000001;
			}
//Convert mass per area units to concentration by volume
//  	kg ha-1 * (1/10000) ha m-2 * 1000 g kg-1 / (m3 m-2) = g m-3
// 	AtoV = (0.1/depth/area) * (depth/area/rhizVolume);
			AtoV = 0.1/rhizVolume;


//Compute residue input C and C/N to the soil and update residue pools
			if (j == 0) //incorporate aboveground residue and root residue
			{
				if (treesParams.usePhenology == true) //perennial plants
				{
//add leaf residue 
					dC = 0.1/48.0*leafResidueCarbon[0]/(730.0/(fns+0.0001))*rhizVolume/depth;
					dN = dC * leafResidueNitrogen[0]/leafResidueCarbon[0];
					leafResidueCarbon[0] -= dC;
					leafResidueNitrogen[0] -= dN;
					ADD = dC;
//add stem residue
					dC = 0.1/48.0*stemResidueCarbon[0]/(36500.0/(fns+0.0001))*rhizVolume/depth;
					dN2 = dC * stemResidueNitrogen[0]/stemResidueCarbon[0];
					stemResidueCarbon[0] -= dC;
					stemResidueNitrogen[0] -= dN2;
					ADD += dC;
				}
				else //annual plants
				{
					ADD = 0.0;
					dN = dN2 = 0.0;
				}
//add root residue
				//dC = rootResidueCarbon[j][k]/(rootLifeSpan/(fns+0.0001));
				dC = rootResidueCarbon[j][k]/rootLifeSpan;
				dN3 = dC * rootResidueNitrogen[j][k]/rootResidueCarbon[j][k];
				rootResidueCarbon[j][k] -= dC;
				rootResidueNitrogen[j][k] -= dN3;
				ADD += dC;
				CNadd = ADD/(dN + dN2 + dN3);
//convert from kg ha-1 to g m-3 rhizosphere volume
				ADD *= AtoV; 

			}
			else //incorporate only root residue
			{
				ADD = rootResidueCarbon[j][k]/rootLifeSpan;
				dN = ADD * rootResidueNitrogen[j][k]/rootResidueCarbon[j][k];
				CNadd = ADD / dN;
				rootResidueCarbon[j][k] -= ADD;
				rootResidueNitrogen[j][k] -= dN;
				ADD *= AtoV; //convert to g m-3 rhizosphere volume
			}
//microbial C and C/N
			if (rhizosphereMicrobialNitrogen[j][k] > 0.0)
			{
				Cb = rhizosphereLiveMicrobialCarbon[j][k];
				CNb = Cb / rhizosphereMicrobialNitrogen[j][k];
				Cb *= AtoV; //convert to g m-3 rhizosphere volume
			}
			else
			{
				Cb = 0.0;
				CNb = 1.0;
			}

//litter C and C/N
			if (rhizosphereNl[j][k] > 0.0)
			{
				Cl = rhizosphereCl[j][k];
				CNl = Cl / rhizosphereNl[j][k];
				Cl *= AtoV; //convert to g m-3 rhizosphere volume
			}
			else
			{
				Cl = 0.0;
				CNl = 1.0;
			}
//exudate amino acid C and N
			if (rootExudateAminoAcidNitrogen[j][k] > 0.0)
			{
				Cea = rootExudateAminoAcidCarbon[j][k];
				CNea = Cea/rootExudateAminoAcidNitrogen[j][k];
				Cea *= AtoV;
			}
			else
			{
				Cea = 0.0;
				CNea = 1.0;
			}
//exudate sugar C
			Ces = rootExudateSugarCarbon[j][k];
			Ces *= AtoV;

//compute decomposition, mineralization-immobilization, and nitrification
//phi is the coefficient of nitrogen sufficiency for bacteria use (Eq. 23 in Porporato 2003 AWR
//
			phi = computeMineralizationImmobilization(MIN, PHI, IMM_pos, IMM_neg, 
							DECl, DECh, DECea, DECes, rh, rr, kd, kn, fns, j, k, 
							t_soil, theta, AtoV, AtoVbulk, depth, treesParams);
//amount of dead microbial biomass
//a linear function assumes no change in rate of predation/competition as the population gets crowded
			BD = kd * Cb * max(1.0,log10(Cb)-2.0);
			//BD = kd * Cb * log10(max(10.0,0.5*Cb));
			//BD = kd * Cb;

//
//update the carbon balance in the fast/litter pool (Eq. 4)
//assuming dead biomass going to labile pool declines as CNl declines towards 20
			dCldt = ADD + max(0.0,1.0-2.0*CNh/CNl)*BD - DECl;
			rhizosphereCl[j][k] += dCldt/AtoV; //convert to kg ha-1

//update the nitrogen balance in the fast/litter pool (Eq. 7)
			dNldt = ADD/CNadd + max(0.0,1.0-2.0*CNh/CNl)*BD/CNb - DECl/CNl;
			rhizosphereNl[j][k] += dNldt/AtoV; //convert to kg ha-1

//update the exudate amino acid C and N pools
			rootExudateAminoAcidCarbon[j][k] -= DECea/AtoV;
			rootExudateAminoAcidNitrogen[j][k] -= DECea/CNea/AtoV;
//update the exudate sugar C pool
			rootExudateSugarCarbon[j][k] -= DECes/AtoV;

//Update microbial biomass pools
//carbon balance in the biomass (Eq. 11)
//rr represents the fraction of carbon lost as CO2 in respiration
			dCbdt = (1.0 - rh -rr)*DECl+(1.0-rr)*(DECea+DECes) + (1.0 - rr)*DECh - BD;
			rhizosphereLiveMicrobialCarbon[j][k] += dCbdt/AtoV; //convert to kg ha-1
//adds microbial rain in to top layer only MCH 05082020
			if (j == 1)
			{
				CNb = rhizosphereLiveMicrobialCarbon[1][k]/rhizosphereMicrobialNitrogen[1][k];
				rhizosphereLiveMicrobialCarbon[1][k] += treesParams.microbialrainrate;
				rhizosphereMicrobialNitrogen[1][k] += treesParams.microbialrainrate / CNb;

                                //extra rain ins added MCH 23092020
                                rhizosphereAmmoniumNitrogen[1][k] += treesParams.raininAmmonium;
                                rhizosphereNitrateNitrogen[1][k] += treesParams.raininNitrate;
                                rhizosphereMineralNitrogen[1][k] += treesParams.raininMineralN;
                                rhizosphereLabileCarbon[1][k] += treesParams.raininLabileC;
			}			
//balance of nitrogen in the biomass (Eq. 12)
//PHI = MIN - (IMM_pos + IMM_neg);
//does not include rr because that is a carbon pathway only
			dNbdt = (1.0 - rh*CNl/CNh)*(DECl/CNl)+
					(DECea/CNea) + 
					DECh/CNh - BD/CNb - PHI;
			rhizosphereMicrobialNitrogen[j][k] += dNbdt/AtoV; //convert to kg ha-1

//Update ammonium in the rhizosphere
			dNposdt = MIN - IMM_pos - LE_pos[j][k] - UP_pos[j][k];
			rhizosphereAmmoniumNitrogen[j][k] += dNposdt/AtoV; //convert to kg ha-1
			if (rhizosphereAmmoniumNitrogen[j][k] < 0.00000001)
			{
				rhizosphereAmmoniumNitrogen[j][k] = 0.00000001;
			}
//then compute nitrification, NIT
			NIT = fns * kn * Cb * rhizosphereAmmoniumNitrogen[j][k]*AtoV;
			if (NIT > rhizosphereAmmoniumNitrogen[j][k]*AtoV)
			{
				NIT = 0.99*rhizosphereAmmoniumNitrogen[j][k]*AtoV;
			}
			dNposdt -= NIT;
			rhizosphereAmmoniumNitrogen[j][k] -= NIT/AtoV; //convert to kg ha-1
			if (rhizosphereAmmoniumNitrogen[j][k] < 0.00000001)
			{
				rhizosphereAmmoniumNitrogen[j][k] = 0.00000001;
			}
//Update nitrate in the rhizosphere
			dNnegdt = NIT - IMM_neg - LE_neg[j][k] - UP_neg[j][k];
			rhizosphereNitrateNitrogen[j][k] += dNnegdt/AtoV; //convert to kg ha-1
			if (rhizosphereNitrateNitrogen[j][k] < 0.00000001)
			{
				rhizosphereNitrateNitrogen[j][k] = 0.00000001;
			}
//Update total mineral nitrogen in the rhizosphere
			rhizosphereMineralNitrogen[j][k] += (dNposdt + dNnegdt)/AtoV; //convert to kg ha-1
			if (rhizosphereMineralNitrogen[j][k] < 0.00000001)
			{
				rhizosphereMineralNitrogen[j][k] = 0.00000001;
			}

//balance of carbon in the humus (Porporato Eq. 8)
			dChdt = rh * DECl + min(1.0,2.0*CNh/CNl)*BD - DECh;

//nitrogen balance equation for humus (Porporato Eq. 10)
			dNhdt = (rh * DECl - DECh)/CNh + min(1.0,2.0*CNh/CNl)*BD/CNb;
			humusCarbon[j] += dChdt/AtoV; //convert to kg ha-1
			humusNitrogen[j] += dNhdt/AtoV; //convert to kg ha-1

			rootRadius *= treesParams.rootDiamMultiplier;
			rootLifeSpan *= treesParams.rootDiamMultiplier;

//Note: CO2 release will be added here as rr * (DECl + DECh + DECea + DECes)
			sumRespiration += rr * (DECl + DECh + DECea + DECes) / AtoV;
			
		}
		heterotrophicRespiration[j] += sumRespiration;
	}
}

//
//computeMineralizationImmobilization()
//This function implements Eqs. 18 and 23 in Porporato et al 2003. Ad. Water Res, 26, 45-58.
//These equations regulate the dynamics of composition, mineralization, and immobilization.
//PHI = phi * fd(s)Cb{khCh[1/(C/N)h - (1-rr)/(C/N)b] + klCl[1/(C/N)l - rh/(C/N)h - (1-rh-rr)/(C/N)b}
//phi is a non-dimensional factor to take into account poor nitrogen and immobilization is not sufficient
//fd(s) nondimensional factor describing soil moisture effects on decomosition
// CNb = 8, requires CN of 24, 8 C metabolized for each N and 16 respired as CO2
//CNh varies from 10-12, CNb from 8-12, CNl from 20 to 50
//
double BiogeochemicalCycles::computeMineralizationImmobilization(double& MIN,
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
						   trees_params treesParams)
{
	double phi, phi_num, phi_den, Cb, kh, Ch, CNh, CNb, kl, Cl, CNl, kea, kes, Cea, CNea, Ces;
	double ki_pos, ki_neg, N_pos, N_neg, fds, fds_opt;
	double porosity;

//reaction coefficients, daily basis
	kd = treesParams.kd; //D-1
	kn = treesParams.kn; //m3 D-1 gC-1
	kea = treesParams.kea; //m3 D-1 gC-1
	kes = treesParams.kes; //m3 D-1 gC-1
	kl = treesParams.kl; //m3 D-1 gC-1
	kh = treesParams.kh; //m3 D-1 gC-1
//convert reaction coefficients to 30min-1 basis, time steps used in TREES
//D-1 -> 30min-1 = 1/48 = 0.02083333
	kd *= 0.02083333; //30-min-1
	kn *= 0.02083333; //m3 30-min-1 gC-1
	kea *= 0.02083333; //m3 30-min-1 gC-1
	kes *= 0.02083333; //m3 30-min-1 gC-1
	kl *= 0.02083333; //m3 30-min-1 gC-1
	kh *= 0.02083333; //m3 30-min-1 gN-1
	ki_pos = ki_neg = 1.0 * 0.02083333; //m3 30min-1 gC-1
	
	assert(root_num >= 0);
	assert(root_num < nRoots);
	assert(root_order >= 0);
	assert(root_order < 10);
//live microbial C and N
	if (rhizosphereMicrobialNitrogen[root_num][root_order] > 0.0)
	{
		Cb = rhizosphereLiveMicrobialCarbon[root_num][root_order];
		CNb = Cb / rhizosphereMicrobialNitrogen[root_num][root_order];
		Cb *= AtoV;
	}
	else
	{
		Cb = 0.0;
		CNb = 1.0;
	}
//humus C and N
	if (humusNitrogen[root_num] > 0.0)
	{
		Ch = humusCarbon[root_num];
		CNh = Ch / humusNitrogen[root_num];
		Ch *= AtoVbulk; //assumes well-mixed humus throughout soil layer
	}
	else
	{
		Ch = 0.0;
		CNh = 1.0;
	}
//litter C and N
	if (rhizosphereNl[root_num][root_order] > 0.0)
	{
		Cl = rhizosphereCl[root_num][root_order];
		CNl = Cl/rhizosphereNl[root_num][root_order];
		Cl *= AtoV;
	}
	else
	{
		Cl = 0.0;
		CNl = 1.0;
	}
//exudate amino acid C and N
	if (rootExudateAminoAcidNitrogen[root_num][root_order] > 0.0)
	{
		Cea = rootExudateAminoAcidCarbon[root_num][root_order];
		CNea = Cea/rootExudateAminoAcidNitrogen[root_num][root_order];
		Cea *= AtoV;
	}
	else
	{
		Cea = 0.0;
		CNea = 1.0;
	}
//exudate sugar C
	Ces = rootExudateSugarCarbon[root_num][root_order];
	Ces *= AtoV;

//mineral N
	N_pos = rhizosphereAmmoniumNitrogen[root_num][root_order]*AtoV;
	N_neg = rhizosphereNitrateNitrogen[root_num][root_order]*AtoV;

//cout << "Cb = " << Cb << "; Ch = " << Ch << "; Cl = " << Cl << endl;

//kl and kh are rates - see soil respiration rate equations
//rr (0 <= rr <= 1-rh is fraction of organic carbon that goes into CO2 production
//rr varies typically from 0.6 to 0.8, increases with labile carbon 
//rh is the isohumic coefficient in range 0.15 to 0.35 depending on plant residues
//assume rh increases with the litter C/N ratio

	rh = min(0.35, CNh/CNl);
	rr = 0.6;
/*
	if (rh < 0.15)
	{
		rh = 0.15;
	}
	else if (rh > 0.35)
	{
		rh = 0.35;
	}
	rr = 0.7 - 0.10 * rhizosphereLabileCarbon[root_num][root_order]*AtoV / 
				(rhizosphereLabileCarbon[root_num][root_order]*AtoV + Cl);
*/
	//rr = 0.6 + 0.2 * (kea*Cea+kes*Ces)/(kes*Cea+kes*Ces+kl*Cl);
	rr = 0.6 + 0.2 * (Cea+Ces)/(Cea+Ces+Cl);
	if (rr > (1.0-rh))
	{
		rr = 1.0-rh;
	}
//This is Eq. 23 in Porporato
	phi_num = ki_pos*N_pos + ki_neg*N_neg;
//	phi_den = kh*Ch*(1.0/CNh - (1.0-rr)/CNb) + kl*Cl*(1.0/CNl - rh/CNh - (1.0-rh-rr)/CNb);
//Eq. 23 modified to include amino acid and sugar exudates
	phi_den = kh*Ch*(1.0/CNh - (1.0-rr)/CNb) + 
			kl*Cl*(1.0/CNl - rh/CNh - (1.0-rh-rr)/CNb) +
			kea*Cea*(1.0/CNea - (1.0-rr)/CNb) +
			kes*Ces*(0.0 - (1.0-rr)/CNb);
	phi = -phi_num/phi_den;


	if (phi_den > 0.0) //net mineralization is occurring
	{
		phi = 1.0;
	}
	if (phi > 1.0) //may be at maximum immobilization
	{
		phi = 1.0;
	}

//In Porporato et al 2003, fds (decomp rate) depends on soil water content
//This is modified here to incorporate the role of temperature as well
//  and this computation is accomplished by using the DAMM (Davidson et al GCB) model
	fds_opt = DAMM_Cflux(treesParams, treesParams.optimal_soil_T, treesParams.theta_opt, 
								treesParams.porosity, depth, Cl);
	if (fds_opt > 0.0)
	{
		fds = DAMM_Cflux(treesParams, t_soil, theta, treesParams.porosity, depth, Cl);
		fds = fds/fds_opt;
	}
	else
	{
		fds = 0.0;
	}
	if (fds > 1.0)
	{
		fds = 1.0;
	}

	DECl = (phi*fds*kl*Cb)*Cl;
	if (DECl > 0.99*Cl)
	{
		DECl = 0.99*Cl;
	}
	DECh = (phi*fds*kh*Cb)*Ch;
	if (DECh > 0.99*Ch)
	{
		DECh = 0.99*Ch;
	}
	DECea = (phi*fds*kea*Cb)*Cea;
	if (DECea > 0.99*Cea)
	{
		DECea = 0.99*Cea;
	}
	DECes = (phi*fds*kes*Cb)*Ces;
	if (DECes > 0.99*Ces)
	{
		DECes = 0.99*Ces;
	}

//This is Eq. 18 in Porporato
// Modified to include amino acids - DSM 3/2/2018
// and sugars - DSM 3/3/2018
//PHI = DECh*(1.0/CNh - (1.0-rr)/CNb) + DECl*(1.0/CNl - rh/CNh - (1.0-rh-rr)/CNb);
//	PHI = phi*fds*Cb*(kh*Ch*(1.0/CNh - (1.0-rr)/CNb) + kl*Cl*(1.0/CNl - rh/CNh - (1.0-rh-rr)/CNb));
//Eq. 18 modified to include amino acid and sugar exudates
	PHI = phi*fds*Cb*(kh*Ch*(1.0/CNh - (1.0-rr)/CNb) + 
			kl*Cl*(1.0/CNl - rh/CNh - (1.0-rh-rr)/CNb) +
			kea*Cea*(1.0/CNea - (1.0-rr)/CNb) +
			kes*Ces*(0.0 - (1.0-rr)/CNb));
	if (PHI >= 0.0)
	{
		MIN = PHI;
		IMM_pos = IMM_neg = 0.0;
	}
	else
	{
		MIN = 0.0;
		IMM_pos = ki_pos*N_pos/(ki_pos*N_pos+ki_neg*N_neg) * (-1.0*PHI);
		IMM_neg = ki_neg*N_neg/(ki_pos*N_pos+ki_neg*N_neg) * (-1.0*PHI);
	}
	return(phi);
}


//
//nitrification_water_function()
//define the role of degree of staturation on nitrification rate
//Two parts:
//      (1) linear rise to an optimum (parameter theta_opt)
//      (2) less than linear decrease at soil_theta > theta_opt
//DSM November 12, 2008
//
double BiogeochemicalCycles::nitrification_water_function(double theta_opt,
                             			          double theta_soil)
{
        double ts1, ts2, increasing, decreasing, rate;

        if (theta_opt < 0.01)
	{
                 theta_opt = 0.01;
	}
        if (theta_opt > 1.0)
	{
                theta_opt = 1.0;
	}
        ts1 = theta_soil/theta_opt;
        ts2 = theta_opt/theta_soil;
        increasing = pow(ts1, 2.0);
        decreasing = pow(ts2, 2.0);
        rate = 1.0/(0.5*(increasing+decreasing));
        return rate;
}

//
//slow_temp_decomp_rate()
//
double BiogeochemicalCycles::slow_temp_decomp_rate(double t_soil,
                        			   double optimal_soil_T,
                        			   double slow_mineral_rate_parm)
{
        double rate;

        rate = pow(slow_mineral_rate_parm, (t_soil - optimal_soil_T)/10.0);

        return rate;

}

//
//fast_temp_decomp_rate()
//
double BiogeochemicalCycles::fast_temp_decomp_rate(double t_soil,
                        			   double optimal_soil_T,
                        			   double fast_mineral_rate_parm)
{
        double rate;

        rate = pow(fast_mineral_rate_parm, (t_soil - optimal_soil_T)/10.0);

        return rate;

}

//
//DAMM_Cflux()
//DAMM model for Rhet //
//DAMM model native output is mgC m-2 hr-1
//
double BiogeochemicalCycles::DAMM_Cflux(trees_params treesParams,
                  			double t_soil,
                  			double theta,
                  			double porosity,
                  			double depth,
                  			double Clitter)
{
	double Resp, areaCflux;

	Resp = DAMM_reactionVelocity(treesParams, t_soil, theta, porosity, depth, Clitter);

//areaCflux is now in mgC m-2 hr-1, convert to kgC ha-1 30min-1
    	areaCflux = 10000.0 * depth * Resp;
    	areaCflux *= (1.0/2000.0);  //   gC / m2 / 30min
    	areaCflux *= 10.0;        // to kg C / ha / 30min

    	return areaCflux;
}

//
//DAMM_reactionVelocity()
//DAMM model for Rhet //
//We assume this will yield a reaction velocity for each rhizosphere
//
double BiogeochemicalCycles::DAMM_reactionVelocity(trees_params treesParams,
                                        	   double t_soil,
                                        	   double theta,
                                        	   double porosity,
                                        	   double depth,
                                        	   double Clitter)
{
	double RkJ, EaSx, kMsx, xASx;
        double O2airfrac, Dliq, Dgas, kMO2, Clitter_gm2, Clitter_gm3, Clitter_gcm3, exponent;
        double Sx, O2, MMSx, MMO2, VmaxSx, R_Sx;

	RkJ = 8.314472 * pow(10.0,-3.0);

	EaSx = treesParams.EaSx;
	kMsx = treesParams.kMsx;
	xASx = treesParams.xASx;

        O2airfrac = 0.209; //L O2 per L air
        Dliq = 10.97;
        Dgas = 1.67;
        kMO2 = 0.0452;

// need C content in g/cm3
// Clitter is in units of kg ha-1, assume all fast C is in top 15cm
//UNIT CONVERSION HAS A PROBLEM
        Clitter_gm2 = Clitter * 0.1;
        Clitter_gm3 = Clitter_gm2 / depth;
        Clitter_gcm3 = Clitter_gm3 * pow(10.0,-6.0);
        exponent = 4.0/3.0;
        Sx = Clitter_gcm3 * Dliq * pow(theta, 3.0);
        O2 = Dgas * O2airfrac * pow((porosity - theta),exponent);
        MMSx = Sx / (kMsx+Sx);
	MMO2 = O2/(kMO2+O2);
        //VmaxSx = exp(xASx - EaSx/(RkJ * (t_soil + 273.15)));
        VmaxSx = xASx * exp(-EaSx/(RkJ * (t_soil + 273.15)));
        R_Sx = VmaxSx * MMSx * MMO2;
        return R_Sx;
}

//
// sequential leaf development
//

// "get" functions

// get biomass C of single leaf
double BiogeochemicalCycles::getSingleLeafBiomassCarbon(int k)
{
        return(SingleLeafBiomassCarbon[k]);
}

// get biomass N of single leaf
double BiogeochemicalCycles::getSingleLeafBiomassNitrogen(int k)
{
    	return(SingleLeafBiomassNitrogen[k]);
}

// get NSC of single leaf
double BiogeochemicalCycles::getSingleLeafNSC(int k)
{
    	return(SingleLeafNSC[k]);
}

// get chloroplast starch of single leaf
double BiogeochemicalCycles::getSingleLeafchloroplastStarch(int k)
{
    	return(SingleLeafchloroplastStarch[k]);
}

// get chloroplast sugar of single leaf
double BiogeochemicalCycles::getSingleLeafchloroplastSugar(int k)
{
    	return(SingleLeafchloroplastSugar[k]);
}

// get stored N of single leaf
double BiogeochemicalCycles::getSingleLeafStoredNitrogen(int k)
{
    	return(SingleLeafStoredNitrogen[k]);
}

// get Rubisco N of single leaf
double BiogeochemicalCycles::getSingleLeafRubiscoNitrogen(int k)
{
    	return(SingleLeafRubiscoNitrogen[k]);
}

// get leaf area of single leaf
double BiogeochemicalCycles::getSingleLeafArea(int k)
{
    	return(SingleLeafArea[k]);
}

// get leaf area potential of single leaf
double BiogeochemicalCycles::getSingleLeafAreaPotential(int k)
{
    	return(SingleLeafAreaPotential[k]);
}

// get thermal time of single leaf
double BiogeochemicalCycles::getSingleLeafThermTm(int k)
{
    	return(SingleLeafThermTm[k]);
}

//"put" functions

// put leaf area potential of single leaf for first time step

void BiogeochemicalCycles::putSingleLeafAreaPotential(double area_pot, 
						      int idx)
{
    	SingleLeafAreaPotential[idx] = area_pot;
}

//
// individual leaves are modeled on much smaller scales than
// leaf population of perenniel species.
// mass: grams
// area: cm2
// Below functions convert between units used for individual leaves
// and units used at whole plant level

// Area: cm2 to m2
double BiogeochemicalCycles::cm2_to_m2(double value)
{
    	return(value/10000.0);
}

// Area: cm2 to ha
double BiogeochemicalCycles::cm2_to_ha(double value)
{
    	return(value/100000000.0);
}

// Mass: g to kg
double BiogeochemicalCycles::g_to_kg(double value)
{
    	return(value/1000.0);
}

//
// relative growth rate functions
//
// leaf growth informed by genotype-specific parameters, K, N0, and r
// K = treesParams.leafAreaMax
// N0 = treesParams.initialLeafSize
// r = treesParams.leafArea_Rate
// To establish these parameters, leaf expansion time courses are recorded under
// non-limiting environmental conditions and fit to a logistic growth curve of
// three parameters (below). This set of parameters describe a genotype's
// "theoretical potential".
// Leaf growth as a function of thermal time follows
//                K
// A = ------------------------
//      1 + (K/(K-N0))e^(-rT)

// where K, N0 and r are described above and T is thermal time of plant
// Framework of for the thermal leaf development module is largely based on work
// described in:
// Granier and Tardieu 1998. Is thermal time adequate for expressing the effects of
// temperature on sunflower leaf develpment? Plant, Cell and Environment.
// relative rate of "division" i.e. area potential formation
// relative rate is the logarithmic derivative of the leaf growth logistic function
// this phase represents the formation of the potential area
// in practical terms, this is the phase of cell division, when cell NUMBER increases


vector<double> BiogeochemicalCycles::getQuantileVals(double alpha, double beta, double min_quant, double max_quant)
{
    vector<double> result_vals;
    double max_num;
    double min_num;
    double theta = 1/beta;
    
    //first estimate quantile values
    double sample_vec[10000] = {0};
    const int samples = 10000;  // number of samples
    
    default_random_engine generator;
    gamma_distribution<double> distribution(alpha, theta);
    
    for (int i=0; i< samples; ++i) {
        sample_vec[i] = distribution(generator);
    }
    
    // sort the array and get the quantile vals
    std::sort(sample_vec, sample_vec + samples);
 ///MCH changed to static_cast for best practices 24092020   
    int min_idx = static_cast<int>(min_quant*10000);
    int max_idx = static_cast<int>(max_quant*10000);
    min_num = sample_vec[min_idx];
    max_num = sample_vec[max_idx];

    result_vals.assign (2, 0); // two elements, initialized with value 0
    result_vals[0] = min_num;
    result_vals[1] = max_num;
    return result_vals;
}

double BiogeochemicalCycles::sampleTruncatedGamma(double alpha, double beta, double minlim, double maxlim)
{
    double sampledval;
    double theta = 1/beta;
    
    gamma_distribution<double> distribution(alpha, theta);
    random_device rd;
    default_random_engine generator(rd());
    
    double number = distribution(generator);
    while (number > maxlim || number < minlim)
    {
        number = distribution(generator);
    }
    sampledval = number;
    
    return(sampledval);
}



double BiogeochemicalCycles::calcRAR(double thermCD, 
				     trees_params& treesParams, double* Karray, double* Narray, double* rarray, int idx)
{
    	double A, f_prime, f;
    	double denom;
        double K = Karray[idx];
        double N = Narray[idx];
        double r = rarray[idx];
    
    	A = (K - N) / N;
    	denom = pow(1.0 + A * exp(-r * thermCD), 2.0);
    	f_prime = (K * A * r * exp(-r * thermCD))/denom;
    	f = K/(1.0 + ((K - N) / N) * exp(-r * thermCD));

    	return(f_prime/f);
}

// relative rate of expansion
// relative rate is the logarithmic derivative of the leaf growth logistic function
// this phase represents the realization of potential area, which is established during the cell division phase
// this phase is the expansion phase, when individual cell size increases

double BiogeochemicalCycles::calcRER(double thermEXP, 
				     trees_params& treesParams, double* Karray, double* Narray, double* rarray, int idx)
{
    	double A, f_prime, f;
    	double denom;
        double K = Karray[idx];
        double N = Narray[idx];
        double r = rarray[idx];
    
        A = (K - N) / N;
        denom = pow(1.0 + A * exp(-r * thermEXP), 2.0);
        f_prime = (K * A * r * exp(-r * thermEXP))/denom;
        f = K/(1.0 + ((K - N) / N) * exp(-r * thermEXP));
    
        return(f_prime/f);
}


//
// respiration cost associated with cell division
// based on construction cost of 1.33 g C substrate/g C product (proteins, nucleic acids, carbohydrates)
// linked to pseudo delta C based on change in potential area formation
// returns value in units of grams C substrate
// Cannell and Thornley 2000
double BiogeochemicalCycles::calcRespCD(double deltaAPot, 
					double pseudoC_coef)
{
    	return(1.33*pseudoC_coef*deltaAPot); // add explanation
}

//
// respiration cost associated with leaf expansion
// based on construction cost of 1.17 g C substrate/g C product (mostly carbohydrates, lignin associated with cell wall biosyn.)
// linked directly to delta C (kgC)
// returns value in units as input
// Cannell and Thornley 2000
double BiogeochemicalCycles::calcRespEXP(double deltaC)
{
    	return(1.17*deltaC);
}

// calculate SLA at this ts for a specific leaf increase
// SLA is in units m2/kg

double BiogeochemicalCycles::calcSLA(double plantNitrogenStatus, 
				     double SLA_max, 
				     double SLA_min)
{
// plantNstatus is from 0 to 1; 1 = happy
// happy = thicker leaves = lower SLA
// SLA_max and SLA_min can represent the limit for SLA range for species or genotype

    	double SLA_instant;

    	SLA_instant = SLA_max - (SLA_max - SLA_min)*plantNitrogenStatus;

    	return(SLA_instant);
}

// Compute delta C required for growth of leaf area and SLA:
// units: deltaAreaL in cm2
// SLA in m2/kg
// returns kg C

double BiogeochemicalCycles::calcdeltaC(double SLA, 
					double deltaAreaL)
{
    	double convertedArea;

    	convertedArea = cm2_to_m2(deltaAreaL);

    	return(convertedArea/SLA);
}

// Compute amount of delta N required for delta C (kgC
// adapted from computeLeafNdemand() from phenology module
// units: kgN
// 1 hectare = 1*10^8 cm2

double BiogeochemicalCycles::calcdeltaN(double deltaC, 
					trees_params treesParams)
{
    	double leafCN, SLA_instant, deltaN;
    	double N_avail_rate_plant = plantNstatus[0];

    	//SLA_instant = calcSLA(N_avail_rate_plant , treesParams.SLA_min*2, treesParams.SLA_min);
    	SLA_instant = calcSLA(N_avail_rate_plant , treesParams.SLA_max, treesParams.SLA_min);
    	//cout << "SLA at this ts of lf is " << SLA_instant << endl;

// this function adjusts the leaf C/N ratio as a function of available N
// allows C/N to vary from 35-146 (high-low N)
// adaptation of computeNdemand() function used in phenology module
    	//leafCN = 146.0 - 111.0*N_avail_rate_plant;
	leafCN = treesParams.leaf_maxCN - (treesParams.leaf_maxCN - treesParams.leaf_minCN)*N_avail_rate_plant;
    	deltaN = deltaC / leafCN;

    	return(deltaN);
}

// compute effect of hydraulic constraint on
// cell division phase of leaf growth
double BiogeochemicalCycles::calcHydRAR(double hydraulicStatus)
{
       return(1.0); // k_p_e.latK[1]/ksat[1][1];
}

// PLACEHOLDER: compute effect of nitrogen status on
// cell division phase of leaf growth
double BiogeochemicalCycles::calcNRAR(double NStatus)
{
    	//return min(NStatus*5.0, 1.0); // this means that if N status is < 0.2, relative growth will be limited\

// test new function, informed by Marc Brock's Timmer experiment data from 2016
    	//return(1.0- pow(26.5938, -1.0*NStatus));
    	//return min( (1.0- pow(26.5938, -1.0*NStatus))*2.0, 1.0 ); // try this -- less sensitivity
    //return(1.0);
//TESTING for fastplant from Marc and Diane's experiemental data. 
//FP ONLY! will update as parameters later. MCH 26102020 --working well! 
	return(0.4257+0.5811*(1-exp(-7.5352*NStatus)));
}

// PLACEHOLDER: compute effect of hydraulic constraint on
// expansion phase of leaf growth
double BiogeochemicalCycles::calcHydRER(double waterStat)
{
    
    double lim1, lim2;
    double coef;
    
    lim1 = 0.65;
    lim2 = 0.85;
    
    if(waterStat > lim2)
    {
        coef = 1.0;
    }
    
    if(waterStat <= lim2 && waterStat >= lim1)
    {
        coef = waterStat/(lim2-lim1) + lim1/(lim1-lim2);
    }
    
    if(waterStat < lim1)
    {
        coef = 0.0;
    }
    
    return(coef); // k_p_e.latK[1]/ksat[1][1];
}

// compute effect of NSC status on
// expansion phase of leaf growth
double BiogeochemicalCycles::calcNSCRER(double NSCStatus, 
					double deltaC)
{ // deltaC = C_needed (kgC/ha)
    	double adjustment, threshold;

    	threshold = 1.0*deltaC; // can play around w/ this threshold

    	if(NSCStatus >= threshold)
    	{
        	adjustment = 1.0;
    	} 
	else
    	{
        	adjustment = NSCStatus/threshold;
    	}

    	return(adjustment);
}

//
// main function to update growth for a single leaf
//

void BiogeochemicalCycles::updateLeaf(int k, 
				      double delta_thermTime, 
				      trees_params& treesParams, 
				      double ProjGrndArea_instant, 
				      double RL,
				      double N_neg_fract,
                                      double& N_neg_demand,
                                      double& N_pos_demand, double* Karray, double* Narray, double* rarray, int idx, double thetaRoot)
{
    	double thermCD, deltaApot, RAR, water_RAR_coef, N_RAR_coef;
    	double Rcd=0.0; // respiration associated with 'cell division' phase
    	double thermEXP, RER, deltaAreaL, SLA_instant, water_RER_coef, NSC_RER_coef, NSCStatus, C_needed;
    	double deltaC = 0.0; //kg
    	double deltaN = 0.0; // kg
    	double Rexp = 0.0; // respiration associated with 'expansion' phase
    	double duration_Leaf; // relationship with duration of expansion (input parameter)
    	double delta_nsc = 0.0; // respiration growth cost, kg/ha
    	double NfromStorage = 0.0; // N cost, kg/ha
        double waterStat = 1.0;
    
    	double K = Karray[idx];
    	double N = Narray[idx];
    	double r = rarray[idx];

//duration of leaf including CD and EXP
    	duration_Leaf = treesParams.dur_LeafExpansion/treesParams.proportion_CD;
    	SingleLeafThermTm[k] += delta_thermTime;

// if leaf is still developing
    	if ((SingleLeafThermTm[k] < duration_Leaf ) && 
			((SingleLeafArea[k] < SingleLeafAreaPotential[k] ) ) )
    	{

// --------------------
// in cell division(CD)
// --------------------
        	if ((SingleLeafThermTm[k] <= treesParams.proportion_CD * duration_Leaf) && 
				(SingleLeafArea[k] <= K ) )
        	{

// amt of thermal time in cell division
            		thermCD = SingleLeafThermTm[k]; 
// theoretical optimal RAR
            		RAR = calcRAR(thermCD, treesParams, Karray, Narray, rarray, idx);

// placeholder- hydraulic constraint on RAR
            		water_RAR_coef = calcHydRAR(1.0);
// placeholder- N constraint on RAR will go here
            		N_RAR_coef = calcNRAR(plantNstatus[0]);
            		deltaApot = RAR * SingleLeafAreaPotential[k] * delta_thermTime * water_RAR_coef * N_RAR_coef;

            		if ((SingleLeafAreaPotential[k] + deltaApot) > K)

// if the change will bring total greater than max
            		{ 
                		deltaApot =  K - SingleLeafAreaPotential[k];
            		}
            		SingleLeafAreaPotential[k]  += deltaApot;
            		Rcd = calcRespCD(deltaApot, 1.0*pow(10.0,-50.0)); // kgC
            		leafNSC[0] -= Rcd/cm2_to_ha(treesParams.pot_size);

            		//cout << "cost of cd, kgC, = " << Rcd << endl;
        	}

// ------------------------------
// and/or in leaf expansion (EXP)
// ------------------------------

        	if (SingleLeafThermTm[k] > ((1.0-treesParams.proportion_CD)* duration_Leaf))
        	{
            		thermEXP = SingleLeafThermTm[k] - ((1.0-treesParams.proportion_CD ) * duration_Leaf);
            		RER = calcRER(thermEXP, treesParams, Karray, Narray, rarray, idx); // theoretical optimal RER

// if current ts is first expansion step
            		if (SingleLeafArea[k] == 0.0) 
            		{
                		SingleLeafArea[k] = N; // cm2

// since intial size ought to be very close to 0, do not take carbon costs here
// bit of an abstraction

            		}

//optimal change in area
            		deltaAreaL = RER * SingleLeafArea[k] * delta_thermTime; 

//add limitations (H20, NSC)
                
                //if soil properties are well characterized, then denominator below set to porosity
                    	waterStat = thetaRoot/treesParams.porosity;
                    //waterStat = thetaRoot/0.32;
            		water_RER_coef = calcHydRER(waterStat);

            		deltaAreaL *= water_RER_coef;
            		SLA_instant = calcSLA(plantNstatus[0], treesParams.SLA_max, treesParams.SLA_min);
// returns in kgC/ha
            		C_needed = calcdeltaC(SLA_instant, deltaAreaL)/cm2_to_ha(treesParams.pot_size); 
            		NSCStatus = leafNSC[0];
            		NSC_RER_coef = calcNSCRER(NSCStatus, C_needed);
            		deltaAreaL *= NSC_RER_coef;
//Temporarily include the N constraint here 
// placeholder- N constraint on RAR will go here
//DSM DEC 2020 - Turn off N limitation during cell expansion stage
            	//	N_RAR_coef = calcNRAR(plantNstatus[0]);
            		deltaAreaL *= N_RAR_coef;

            		if ((SingleLeafArea[k]  + deltaAreaL) > SingleLeafAreaPotential[k] )
            		{
                		deltaAreaL = SingleLeafAreaPotential[k]  - SingleLeafArea[k];
            		}
//update area, cm2
//            		SingleLeafArea[k] = SingleLeafArea[k] + deltaAreaL;  

// calculate C and N requirements
            		deltaC = calcdeltaC(SLA_instant, deltaAreaL); // returns in kgC
 //           		SingleLeafBiomassCarbon[k] +=  deltaC; // structural C, in kgC
// structural N, in kgN, modulated by plant N status
            		deltaN= calcdeltaN(deltaC, treesParams);
                    	Rexp = calcRespEXP(deltaC);
                
                    	delta_nsc = Rexp/cm2_to_ha(treesParams.pot_size); //kgC/ha
                
//added - DSM
  //          		SingleLeafBiomassNitrogen[k] +=  deltaN; // structural N, in kgC
            		//cout << "deltaN at this ts is " << deltaN << endl;
            		Rexp = calcRespEXP(deltaC); // growth resp associated with expansion. kg C

// update costs - not used
            		delta_nsc = Rexp/cm2_to_ha(treesParams.pot_size); //kgC/ha

            		
//
// call original TREES function to update carbon pools
//The following treesParams.SLA/SLA_instant addresses a scaling between individual leaf SLA and whole plant SLA
//  so that the updateLeafCarbonNitrogenPools() uses SLA_instant of the current leaf
//DSM July 2020
//
			deltaAreaL = (deltaAreaL/treesParams.pot_size)*(treesParams.SLA/SLA_instant);

                	updateLeafCarbonNitrogenPools(k, treesParams, deltaAreaL, RL,
					N_neg_fract, N_neg_demand, N_pos_demand);
        	}
    	} 
	else
    	{
        // ----------------
        // else the leaf is mature
        // ----------------
	;

        // update resp and photo values? these can be calculated on the fly

    	}
}


// functions below to compute plant-level traits from
// leaf-level traits
// Unit conversions are done within functions in order to match
// units from original phenology module

//
// compute plant-level average SLA
// returns SLA in units m2/kg

double BiogeochemicalCycles::calcPlantSLA(int Lf_idx, 
					  trees_params treesParams)
{
    	double sum_area;
    	double sum_C;
	double SLA;

// biomass C and area at model initialization
    	double shoot_init, leaf_init;

// grams shoot mass total
    	//shoot_init = (1.0/treesParams.SLA_max*10.0)*treesParams.projectedArea_init; 
    	shoot_init = 1.0/(treesParams.SLA_max*8.0)*treesParams.projectedArea_init; 
// kg leaf mass total
    	leaf_init = 0.8*g_to_kg(shoot_init); 
//kg leaf biomass C
    	sum_C = leaf_init/(1.0+treesParams.leafNSCscalar); 
    	sum_area = treesParams.projectedArea_init; //cm2 leaf

// add biomass C and area of additional leaves
    	for (int i = 0; i < Lf_idx; i++)
    	{
        	sum_area += SingleLeafArea[i]; //cm2
        	sum_C += SingleLeafBiomassCarbon[i]; // kg
    	}

// unit conversion
// 1 m2 = 10,000 cm2
// SLA in m2/kg C
	SLA = sum_area/(sum_C*10000.0);
    	return(SLA); 
}

// Compute plant-level LAI from single leaves
// and projected ground area
// returns LAI in cm2/cm2 = m2/m2

double BiogeochemicalCycles::calcLAI(int Lf_idx, 
				     trees_params treesParams, double* Karray)
{
    	double sum_ProjArea;
    	double max_ProjArea, len, len_proj, lai;
        double avg_K;
        double sum_K = 0.0;
    
    	for (int k = 0; k < 500; k++)
    	{
        	sum_K += Karray[k];
    	}
    	avg_K = sum_K/500.0;
    //cout << avg_K << endl;
    
// maximum projected ground area (denominator of LAI)
    	max_ProjArea = 4.0*treesParams.leaf_len_to_width * avg_K *
				pow(sin(treesParams.leaf_insertAngle * PI/180.0), 2.0); //cm2

// projected area of leaves at model initialization (numerator of LAI)
    	sum_ProjArea = treesParams.projectedArea_init; //cm2

// add projected area of subsequent leaves
    	for (int i = 0; i < Lf_idx; i++)
    	{
        	len = 2.0*treesParams.leaf_len_to_width * 
				sqrt( SingleLeafArea[i]/(PI * treesParams.leaf_len_to_width));
        	len_proj = len * sin(treesParams.leaf_insertAngle * PI/180.0);
        	sum_ProjArea += PI*len_proj*sqrt(SingleLeafArea[i]/
				(PI * treesParams.leaf_len_to_width)); //cm2
    	}

    	//return max(1.0, sum_ProjArea/max_ProjArea);
    	return(sum_ProjArea/max_ProjArea); //cm2/cm2 = m2/m2
}


double BiogeochemicalCycles::calc_projGndArea(trees_params treesParams, double* Karray)
{
    	double max_ProjArea, len, len_proj, lai;
        double avg_K;
        double sum_K = 0.0;
    
    	for (int k = 0; k < 500; k++)
    	{
        	sum_K += Karray[k];
    	}
    	avg_K = sum_K/500;
    //cout << avg_K << endl;
    
// maximum projected ground area (denominator of LAI)
    	max_ProjArea = 4.0*treesParams.leaf_len_to_width * avg_K *
				pow(sin(treesParams.leaf_insertAngle * PI/180.0), 2.0); //cm2
	return(max_ProjArea);
}

double BiogeochemicalCycles::calc_projArea(int Lf_idx,
					   trees_params treesParams)
{
	double sum_ProjArea;
	double len, len_proj;

// projected area of leaves at model initialization (numerator of LAI)
        sum_ProjArea = treesParams.projectedArea_init; //cm2

// add projected area of subsequent leaves
        for (int i = 0; i < Lf_idx; i++)
        {
                len = 2.0*treesParams.leaf_len_to_width * 
                                sqrt( SingleLeafArea[i]/(PI * treesParams.leaf_len_to_width));
                len_proj = len * sin(treesParams.leaf_insertAngle * PI/180.0);
                sum_ProjArea += PI*len_proj*sqrt(SingleLeafArea[i]/
                                (PI * treesParams.leaf_len_to_width)); //cm2
        }     
	return(sum_ProjArea);
}


//
// compute plant-level leaf population biomass carbon
// from single leaves
// ProjGrndArea_instant in unit cm2
// returns value in kgC/ha
// does unit conversion within function
// 1 ha = 1x10^8 cm2

//NOT CURRENTLY USED

double BiogeochemicalCycles::getSumLeafBiomassCarbon(int Lf_idx, 
						     double ProjGrndArea_instant)
{
    	double leafBiomassCarbon;
    	double sum = 0.0;

    	for (int i = 0; i < Lf_idx; i++)
    	{
// check this... i think leaf state var are already in kg??
        	sum += SingleLeafBiomassCarbon[i]; // grams C 
    	}
    	leafBiomassCarbon = sum/ProjGrndArea_instant * 100000.0;

    	return(leafBiomassCarbon); //  kg/ha
}

//
// compute plant-level leaf population biomass N
// from single leaves
// ProjGrndArea_instant in unit cm2
// returns value in kgN/ha
// does unit conversion within function
// 1 ha = 1x10^8 cm2

double BiogeochemicalCycles::getSumLeafBiomassNitrogen(int Lf_idx, 
							double ProjGrndArea_instant)
{
    	double leafBiomassNitrogen;
    	double sum = 0.0;

    	for (int i = 0; i < Lf_idx; i++)
    	{
        	sum += SingleLeafBiomassNitrogen[i]; // grams N
    	}
    	leafBiomassNitrogen = sum/ProjGrndArea_instant * 100000.0;

    	return(leafBiomassNitrogen); //  kg/ha};
}

//
// compute total leaf growth respiration
// single leaf growth resp in units g substrate
// function returns total leaf growth respiration
// in units kg nsc/ha to match original function to update
// leaf carbon nitrogen pools
// 1 ha = 1x10^8 cm2

double BiogeochemicalCycles::getTotLeafGrowthRespiration(int Lf_idx, 
							 double ProjGrndArea_instant)
{
    	double TotLeafGrowthRespiration;
    	double sum = 0.0;

    	for(int i = 0; i < Lf_idx; i++)
    	{
        	sum += SingleLeafGrowthRespiration[i];  // grams substrate
    	}
    	TotLeafGrowthRespiration = sum/ProjGrndArea_instant * 100000.0;

    	return(sum);
}


//  functions below are placeholders
//  for leaf-level non-structural resources
//  ... for now, keep computing non-structural resources
//  at plant level


//
// compute plant-level leaf population NSC
// from single leaves
//
double BiogeochemicalCycles::getSumLeafNSC(int Lf_idx)
{
    	double sum= 0.0;

    	for (int i = 0; i < Lf_idx; i++)
    	{
        	sum += SingleLeafNSC[i];
    	}
    	return(sum);
}

//
// compute plant-level leaf population chl starch
// from single leaves
//
double BiogeochemicalCycles::getSumLeafchloroplastStarch(int Lf_idx)
{
    	double sum= 0.0;

    	for (int i = 0; i < Lf_idx; i++)
    	{
        	sum += SingleLeafchloroplastStarch[i];
    	}
    	return(sum);
}

//
// compute plant-level leaf population chl sugar
// from single leaves
//
double BiogeochemicalCycles::getSumLeafchloroplastSugar(int Lf_idx)
{
    	return(1.0);
}

//
// compute plant-level leaf population stored N
// from single leaves
//
double BiogeochemicalCycles::getSumLeafStoredNitrogen(int Lf_idx)
{
    	double sum= 0.0;

    	for (int i = 0; i < Lf_idx; i++)
    	{
        	sum += SingleLeafStoredNitrogen[i];
    	}
    	return(sum);
}

//
// compute plant-level leaf population Rubisco N
// from single leaves
//
double BiogeochemicalCycles::getSumLeafRubiscoNitrogen(int Lf_idx)
{
    	double sum= 0.0;

    	for (int i = 0; i < Lf_idx; i++)
    	{
        	sum += SingleLeafRubiscoNitrogen[i];
    	}
    	return(sum);
}

