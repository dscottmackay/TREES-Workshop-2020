2150.0	altitude, m@ChEAS
35.818	latitude
-106.305	longitude
6.0	z_ref, m 
2.35 20.0	lai, single sided
3.0 0.5 5.5	canopy_height, m
2.35 4.5 lai_at_full_canopy_height
1.0	l_angle, spherical, may sample
0.97	canopy emissivity
0.5	fPAR_beam
0.5	fPAR_diff
0.8	alpha_PAR
0.2	alpha_NIR
0.5 0.5 1.00	omega
2.0	p_crown
0.5	d_factor, C&Nfig5.5
0.1	zm_factor, C&Nfig5.5
0.2	zh_factor, 7.19
1 1 2   ps_model, photosynthesis model to use (1 = original C3, 2 = experimental C3, 3 = C4)
0.017 0.020 0.005      Rd_mult, Rd=Rd_mult*Vcmax
2.2 3.0    Jmax_mult, ratio of Jmax to Vcmax
0.8 0.7 0.9     thetaJ, curvature parameter
0.20 0.01    phiJ_sun, quantum yield, e-/photon
0.20 0.01    phiJ_shd, quantum yield, e-/photon
0.0017 0.0001 0.0040    Nleaf
1.0 0.785138 1.0 N_fixed_proportion
0.046 0.055 0.2  Nrubisco, proportion
40.4 35.0 55.0  Kc25, (Pa) MM const carboxylase, 25 deg C
2.1    q10Kc, (DIM) Q_10 for kc (default 2.1)
24800.0   Ko25, (Pa) MM const oxygenase, 25 deg C
1.2    q10Ko, (DIM) Q_10 for ko
3.6    act25, (umol/mgRubisco/min) Rubisco activity
2.4 2.0 3.0    q10act, (DIM) Q_10 for Rubisco activity (default was 2.4)
48.85 46.0 60.0 Vcmax25, maximum Rubisco activity at 25 C, umol m-2 s-1
460.0 58.0 100.0 Vpmax25, maximum PEP carbolylase activity at 25 C, umol m-2
175.0 400.0 400.0 Jmax25, maximum electron transport rate at 25 C, umol m-2 s-1
38.6 38.6 38.6 gammaStar25, compensation point at 25 C, umol
80.0 80.0 Kp25, Michaelis constant of PEP carboxylase for CO2 at 25 C, ubar
80.0 80.0 80.0 Vpr, PEP regeneration rate, umol m-2 s-1
0.0 0.15 0.15 f, correction for spectral quality of light
0.4 0.4 0.4 x, partitioning factor of electron transport rate
0.92 0.85 0.85 absorptance, fraction of irradiance absorbed
67.294 67.294 E_Vcmax, activation energy, maximum carboxylation rate, kJ mol-1
70.373 70.373 E_Vpmax, activation energy, maximum PEP rate, kJ mol-1
77.9 77.9 77.9 E_Jmax, activation energy, electron transport, kJ mol-1
36.3 36.3 E_Kp, activation energy, Michaelis reaction of PEP, kJ mol-1
59.36 59.36 E_kc, activation energy, Michaelis reaction of carboxylation, kJ mol-1
35.94 35.94 E_ko, activation energy, Michaelis reaction of oxygenation, kJ mol-1
66.3 66.3 66.3  E_Rd, activation energy, Michaelis reaction of mitochondrial respiration, kJ mol-1
23.4 23.4 E_gammaStar, activation energy, Michaelis reaction of compensation point, kJ mol-1
1.78 1.78 1.78 gm, mesophyll conductance to CO2, mol m-2 s-1
0.003 0.003 0.003 gbs, conductance of the bundle sheath, mol m-2 s-1
0.09 0.09 0.09 alphaGmax, fraction of glycolate carbon diverted to glycine during photorespiration
0.38 0.38 0.38 alphaSmax, fraction of glycolate carbon diverted to serine during photorespiration
1.5 1.5 1.5 Nmax, maximum rate of de novo nitrogen supply to the chloroplast, umol N m-2 s-1
0.024 0.024 0.0131 0 	Gsref
0.54 0.36 0.61  	m (proportion of Gsref)
0 1 0  isAmphistomatous, (1 or 0) has stomata on both sides of leaf
-3.584   Md
-2.13 0.1 midday_at_sat_kl(kolb and sperry)
1.5 0.54 e_at_saturated_kl
3.6 rhizosphere_width_(mm)
3 soilshells
0.04 0.0281 0.156 GMP_(mm)_geometric_mean_particle_diameter
22.0 8.8 GSD_geometric_standard_deviation_of_particle_size
1.4 BD_(MG/m3)_soil_bulk_density
0.43 0.4 0.5 porosity
0.20 silt_fraction
0.35 clay_fraction
0.109 residual, residual water content m3 m-3
1.0 frac_absorbing_length
0.01 Capacitance_(mol/Mpa*m2)_on_leaf_area_basis
1.0 axK:latKr_shoot_modules
1.0 axkr:latKr_root_modules
50.0 %total_R_in_root_system
-0.93 -1.06 0.05  pd_at_sat_kl
3.43 14.1 3.43 ax_Shoot-b_value_(weibull)
1.65 3.55 1.65 ax_Shoot-c_value_(weibull)
3.43 14.1 3.43 lat_Shoot-b_value_(weibull)
1.65 3.55 1.65 lat_Shoot-c_value_(weibull)
3.57 10.1 3.57 ax_Root-b_value_(weibull)
4.07 6.96 4.07 ax_Root-c_value_(weibull)
3.57 10.1 3.57 lat_Root-b_value_(weibull)
4.07 6.96 4.07 lat_Root-c_value_(weibull)
40.0 initial_conductivity(root)
0.01 decrement(root)- default 0.001
80.0 initial_conductivity(shoot)
0.02 decrement(shoot)
0.22 0.05 0.48 theta_opt
30.0 25.0 45.0 optimal_soil_T
1.00   growth_resp_proportion
0.0011 0.0010 0.0011 resp_coef_root, kg kg-1 day-1 deg 
0.0002 0.00018 0.0002 resp_coefficient_stem, kg kg-1 day-1 deg
0.0004 0.00036 0.0004 resp_coefficient_leaf, kg kg-1 day-1 deg
0.04 resp_coefficient (Q10), degC-1
72.26 71.22 73.30 EaSx, kjmol-1
0.000000995 0.000000877 0.00000111 kMsx, gCcm-3soil
0.000000000538 0.000000000347 0.000000000834 mgCcm-3soilh-1
0.0045 0.004256 0.0085 kd, d-1
0.15 0.6 0.6 kn, m3 d-1 gC-1
0.0 0.01 0.2 kea, m3 d-1 gC-1 (for exudates)
0.0 0.02 0.6 kes, m3 d-1 gC-1 (for exudates)
0.00002167 0.00001625 0.000065 kl, m3 d-1 gC-1
0.000000833 0.000000625 0.0000025 kh, m3 d-1 gC-1
58.0 40.0 fr_minCN, minimum fine root C:N ratio
116.0 22.0 fr_maxCN, maximum fine root C:N ratio
100.0 13.0 leaf_minCN, minimum leaf C:N ratio
150.0 22.0 leaf_maxCN, maximum leaf C:N ratio
114000.0 19782.4 44960 Cbelowground, kg ha-1
0.0015 0.00 0.986 Clitter_frac, dim
0.013 0.033 0.0416 Croot_frac, dim
10600.0 13710.0 29460.0 Cstem, kg ha-1
15300.0 6800.0 2.0 Csapwood, kg ha-1
0.20 0.001 0.1 Ccoarseroot_frac, dim
0.0000001 interception_per_leafArea, m m2 m-2
0.004 0.000001 0.05 litter_capacity, m
0.215 0.062 0.48 theta_deep0, initial
0.212 0.064 0.48 theta_mid0, initial
0.209 0.064 0.48 theta_shallow0, initial
0.004 0.05 litter_store, initial
7.6 3.5 3.5 SLA, m2 kgC-1 leaf
111.9 150.0 SRL1, m gC-1 specific root length at root diameter of 250 um
0.00039 0.0009 0.000125 minRootDiam, m diameter of thickest root
0.32 0.004 0.064 maxRootDiam, m diameter of thickest root
0.5 0.25 0.04 minRootLifespan, years, lifespan of finest root at lowest C:N ratio
1.06 0.1 1.0 LWP_spring_minimum, -MPa
4.5 1.5 2.5 LWP_stomatal_closure, -MPa
0 is_bryophyte (1 is yes, 0 is no)
0.0 0.000100 1.0 capRiseScalar (0 to 1)
1.0 0.0 1.0 precipReduction
1.0 drainScalar
0.12 0.12 0.1 leafNSCscalar (proportion of leaf structural carbon)
1 usePhenology
6.0 leafLifeSpan
10 max_iteration(the_max_number_of_iterations_to_achieve_convergence_Delta<THR
100.0 3.0 7.01 microbiomeScalar, unitless, multiplier for the initial nutrient status of microbiome
0.0 microbialrainrate
0.0 raininAmmonium
0.0 raininNitrate
0.0 raininMineralN
0.0 raininLabileC
0.0 snowpack_water_equivalent, m
0.0 snowpack_E_deficit_max, deg C
0.005 melt_Rcoef, m degC-1 30-min-1
1 0 1 useHydraulics
0 0 1 useInputStress
0 0 1 useInputWaterTable
213  dayToStopMaizeRefilling, used when both usePhenology and useLeafModule are set to zero (false)
0 useLeafModule
17.73 leafAreaMax // K, individual leaf area in cm2
0.00007596 initialLeafSize //A_pot_in in cm2
0.00090567 leafArea_Rate //r, shape of growth curve, dimensionless
16907.86 dur_LeafExpansion//d_exp deg C per hour, total thermal time available per leaf
7.6 7.6 400.0 SLA_max //SLA_max
9.6 7.6 60.0 SLA_min //SLA_min
60.0 leaf_insertAngle // leaf insertion angle
2.5 leaf_len_to_width // leaf length to width ratio
0.75 0.90 proportion_CD //a
1970.534 phyllochron //phyllochron C hour between leaf events
50000.0 floweringTime //TTF, thermal time at which bolting takes place
0.96 Tbase //Tb
4926.335 therm_plant_init
40.0 projectedArea_init // projected shoot area at initiation
289.0 40.0 pot_size //max projected area, cm2
0.5     root_to_shoot
1.59    leaf_to_stem
0 useLeafGamma
883.31 Kalpha
22.87 Kbeta
8.11 Nalpha
11.89 Nbeta
118.31 ralpha
131431.51 rbeta
3.41739e-05 0.0000001 0.0001 sd_err_Ec
5.47241 0.001 12.0 sd_err_NEE
0.5 0.0 1.0 sd_err_Ec_weight
