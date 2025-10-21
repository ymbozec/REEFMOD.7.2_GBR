# REEFMOD.7.2_GBR

This repository contains the scripts of ReefMod-GBR, a coral individual-based model that simulates coral populations across 2,300 km of Australia's Great Barrier Reef (GBR).

The model reconstructs coral trajectories across the GBR between 2008-2024 and forecasts possible coral futures (2024-2100) based on temporally- and spatially-explicit forcing of water quality, cyclones, heat stress (mass coral bleaching) and the simulated population dynamics of the coral-eating crown-of-thorns starfish (CoTS). Heat stress projections for the GBR are available from multiple climate models of the CMIP-5 and CMIP-6 suites, and multiple scenarios of greenhouse gas emission (RCPs/SSPs). Following version 7.0, the model simulates mechanisms of coral adaptation to heat stress (adaptation through selection). This version has been used for simulating the FY2025 Counterfactual scenarios (SSP1-2.6, SSP2-4.5, SSP3-7.0) and is being used for the exploration of new management interventions within the Reef Restoration and Adaptation Program (RRAP: https://gbrrestoration.org/).

Version 7.2 has refined the implementation of interventions with mechanisms of coral adaptation introduced in version 7.0:
   - the outplanting of corals of specified species group, size and heat tolerance at a specified density
   - the enrichment of coral larvae of specified species group with a specified density
    
Version 7.2 also integrates the last implementation (June 2024, Christina Skinner @UQ) of the CoTS control module used for ReefMod smulations in the COTS Control Innovation Program (https://barrierreef.org/cots-control-innovation-program), the extension of hindcast forcing up to 2024 (allowing to simulate the impacts of the 2024 GBR mass bleaching), new connectivity for coral larvae (1 km resolution) and the first implementation of Allee effects on coral fertilisation (Mumby et al. 2024).

See MAIN_REEFMOD_GBR.m for the full list of new implementations.

## Note on Allee effects (Oct 2025)
The current version integrates the refinement of coral fertilisation success (FS) when the option of Allee effects is turned OFF. 

Previously FS = 1 and now FS = 0.4. This change enables a fair comparison of coral dynamics with vs. without Allee effects (as the maximum fertilisation success *with* Allee effects is around 0.3-0.4). As a consequence, the shape parameter of the larval-stock recruitment relationship had to be re-calibrated (coral recovery dynamics consistent with LTMP observations + simulated number of coral juveniles consistent with GBR observations, as in Bozec et al. 2022).

These modifications were performed *after* running the RRAP Counterfactual FY2025 study (Feb 2025). The FY2025 Counterfactuals were obtained with the option Allee effects turned OFF (MAIN_REEFMOD_GBR), the previous parametrisation of the shape parameter of the larval-stock recruitment (CORAL.BH_beta = 5x1e6xones(6,1) in f_multiple_reef) and fertilisation success (FS) set to 1 (f_runmodel). Running the Counterfactuals with the new parameterisation (and the option Allee effects still OFF) should give very similar (if not the exact same) results.

## Extra requirements

Running the code requires:
1. inclusion into the folder /data/Climatology of the CMIP6 projections of heat stress. These projections are in the folder /Climatology/Future/CMIP6 of the repository [REEFMOD.7.0_GBR](https://github.com/ymbozec/REEFMOD.7.0_GBR)
2. inclusion into the folder /data of the file GBR1_CONNECT.mat (new coral connectivity at 1 km resolution).
This file exceeds GitHub’s size limit and is therefore not included here. Please contact me (y.bozec@uq.edu.au) to obtain a copy.

## Instructions

The code is written in MATLAB (2023b or earlier versions).
To execute the model:
1. Download all the necessary scripts and folders.
2. Add them to your current MATLAB path.
3. In the Command Window, type:
    > run('MAIN_REEFMOD_GBR.m')

This will start the simulation. The current settings run a projection of one climate change scenario (ie, one CMIP-5 or CMIP-6 climate model under a specific scenario of carbon emission RCP/SSP - as specified by the user) for the period 2008-2100, with CoTS control as the only management intervention (counterfactual simulation). See REEFMOD.6.8_GBR/MAIN_REEFMOD_GBR.m for using CMIP-5 models as input.

The number of repeat simulations can be set with 'NB_SIMULATIONS' (currently set to 20). Simulations are then executed sequentially, each identified by the iterator "simul" (eg, from 1 to 20), which sets set a specific seed for the MATLAB random number generator, ensuring reproducibility of the results. Each simulation is stochastic, incorporating several randomised components, including the timing of future heat stress within each decade, the selection of a specific scenario of future cyclones, the initialisation of coral cover and Crown-of-Thorns Starfish (CoTS) density, the magnitude of coral mortality events, the forcing scheme of water quality. Because the runtime of one complete simulation (ie, from year 2008 to year 2100) is about 2 hours, the use of HPC resources is recommended. Shorter simulations can be obtained by setting a lower number of 6-month time steps ('NB_TIME_STEPS').


## Citation
Bozec, Y.-M., A. A. Adam, B. Arellano-Nava, A. K. Cresswell, V. Haller-Bull, T. Iwanaga, L. Lachs, S. A. Matthews, J. K. McWhorter, K. R. N. Anthony, S. A. Condie, P. R. Halloran, J. C. Ortiz, C. Riginos, and P. J. Mumby. 2025. A rapidly closing window for coral persistence under global warming. bioRxiv. https://www.biorxiv.org/content/10.1101/2025.01.23.634487v1.full

## Earlier model versions for the GBR
Bozec, Y.-M., K. Hock, R. A. Mason, M. E. Baird, C. Castro-Sanguino, S. A. Condie, M. Puotinen, A. Thompson, and P. J. Mumby. 2022. Cumulative impacts across Australia’s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs 92(1), e01494
https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecm.1494

Castro-Sanguino, C., Y.-M. Bozec, S. A. Condie, C. S. Fletcher, K. Hock, C. Roelfsema, D. A. Westcott, and P. J. Mumby. 2023. Control efforts of crown‐of‐thorns starfish outbreaks to limit future coral decline across the Great Barrier Reef. Ecosphere 14:e4580. https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.4580


## Contact
Yves-Marie Bozec, The University of Queensland (y.bozec@uq.edu.au)
