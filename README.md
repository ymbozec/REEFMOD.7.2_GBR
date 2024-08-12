REEFMOD.7.2_GBR

This repository contains the scripts of ReefMod-GBR, a coral individual-based model that simulates coral populations across 2,300 km of Australia's Great Barrier Reef (GBR).

The model reconstructs coral trajectories across the GBR between 2008-2024 and forecasts possible coral futures (2024-2100) based on temporally- and spatially-explicit forcing of water quality, cyclones, heat stress (mass coral bleaching) and the simulated population dynamics of the coral-eating crown-of-thorns starfish (CoTS). Heat stress projections for the GBR are available from multiple climate models of the CMIP-5 and CMIP-6 suites, and multiple scenarios of greenhouse gas emission (RCPs/SSPs). Following version 7.0, the model simulates mechanisms of coral adaptation to heat stress (adaptation through selection). This version is being used for the exploration of new management interventions within the Reef Restoration and Adaptation Program (RRAP: https://gbrrestoration.org/).

Version 7.2 has refined the implementation of interventions with mechanisms of coral adaptation introduced in version 7.0:
   - the outplanting of corals of specified species group, size and heat tolerance at a specified density
   - the enrichment of coral larvae of specified species group with a specified density
    
Version 7.2 also integrates the last implementation (Jun 24, from UQ Christina Skinner) of the CoTS control module used for ReefMod smulations within the COTS Control Innovation Program (https://barrierreef.org/cots-control-innovation-program)

WARNING: To be run the code requires inclusion of the folder REEFMOD.7.0_GBR/data into the Matlab working directory.

Either CMIP-5 or CMIP-6 models can be selected from the front script MAIN_REEFMOD_GBR.m (but see REEFMOD.6.8_GBR/MAIN_REEFMOD_GBR.m for using CMIP-5 models as input). The script is designed to run the model under one climate change scenario (ie, one CMIP-5 or CMIP-6 model with the available scenario of carbon emission RCP/SSP) chosen by the user. The number of repeat simulations, ie, the stochastic simulation of the same warming scenario under different projections of cyclones and other random events (including the initialisation of coral cover and CoTS density, the magnitude of coral mortality, the forcing scheme of water quality,...) can be set with 'NB_SIMULATIONS' (set to 20 in RRAP). Because the runtime of one model run is about 2 hours, use of HPC resources is recommended.


Yves-Marie Bozec, The University of Queensland (y.bozec@uq.edu.au)

Citation: Bozec, Y.-M., K. Hock, R. A. Mason, M. E. Baird, C. Castro-Sanguino, S. A. Condie, M. Puotinen, A. Thompson, and P. J. Mumby. 2022. Cumulative impacts across Australia’s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs 92(1), e01494 https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecm.1494# REEFMOD.7.2_GBR
