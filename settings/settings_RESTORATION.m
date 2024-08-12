% REEFMOD-GBR restoration settings
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 09/2019 (RRAP investment case)
%
% Re-worked in 08/2021 for greater flexibility + extended options
% Feb-Apr 2022: specific to deployment in the Moore Reef cluster
%__________________________________________________________________________

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Options for coral outplanting (aquaculture corals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timing of effort
META.doing_coral_outplanting = RESTORATION.doing_coral_outplanting ; % timing of coral deployment

% Restoration effort in terms of total number of outplants available at each time step
META.total_nb_outplants = RESTORATION.total_nb_outplants;
% Restoration effort in terms of number of reefs where coral outplanting is undertaken at each time step
META.nb_reefs_outplanted = [] ; % ? Set to Inf if unlimited OR if specific reefs are restored (listed in SETTINGS_RESTORATION
% If 0, outplanting cannot happen. For the counterfactual, set to Inf with outplanted_density = 0 for ghost deployment
META.outplant_density_variable = 1; % 1: deployed density on a reef is variable and deployment area (as a proportion of reef area) is fixed;
% 0: density is fixed (META.outplanted_density), and deployment area is variable (NOT IMPLEMENTED YET)

% Set the proportion of the different coral types to be outplanted
META.outplant_species_prop = [0.02 0.14 0.14 0 0.7 0]; % must sum to 1
META.outplanted_density = 6.8; % Density (m-2) of deployed outplants on a reef (6.8 corals per m-2 corresponds to the density of 1-yr old corals
% observed in the Philippines after experimenting larval re-seeding). RRAP 2022 intervention study was assuming density for outplants = density after
% larval re-seeding for comparing the efficiency of the two techniques (focusing on other parameters)

% META.outplant_density_aquaculture.slope = -0.72 ; % function to adjust outplant density (NOT IMPLEMENTED YET)
% META.outplant_density_aquaculture.intercept = 11.8 ; % function to adjust outplant density (NOT IMPLEMENTED YET)

% size of coral outplant (as 1yr old corals)
META.outplant_coral_diameter_mean = [2.56 2.56 2.56 1.41 1.41 1.41] ; % mean diameter (in cm) of outplants for each deployed group
META.outplant_coral_diameter_sd = [0.26 0.26 0.26 0.14 0.14 0.14] ; % sd diameter (in cm) of outplants for each deployed group

% Mean heat tolerance of outplants relative to the contemporary mean of the group
META.MEAN_HT_outplants = RESTORATION.MEAN_HT_outplants;  % as +°C-week (DHW)
% Variance of heat tolerance of outplants
META.VAR_HT_outplants = CORAL.VAR_HT;  % as +°C-week (DHW) - same as native corals?
% Maximum heat tolerance of outplants
META.MAX_HT_outplants = CORAL.MAX_HT;  % as +°C-week (DHW) - same as native corals?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Options for moving corals (larval enrichment)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timing of effort
META.doing_larval_enrichment = RESTORATION.doing_larval_enrichment ; % timing of rubble stab

% Set the total amount of larvae delivered
META.total_nb_larvae = RESTORATION.total_nb_larvae ; % defined in MAIN. Set to Inf if using density as input
% Restoration effort: number of reefs where larval enrichment is undertaken at each time step
META.nb_reefs_enriched = RESTORATION.nb_reefs_enriched ;
META.enrichment_density_variable = 0 ; % 1: (NOT IMPLEMENTED YET) deployed density on a reef is variable and deployment area (as a proportion of reef area) is fixed;
% 0: density is fixed (META.enriched_density), and deployment area is fixed (number of reefs is variable)

% Set the proportion of the different coral types in the pool of larvae delivered
% META.enriched_species_prop = [ 0.05  0.25  0.25  0.15  0.15  0.15]; % ~average props on midshelf/outer reefs (Jez) ;
% META.enriched_species_prop = CORAL.prop_settlers ; % use default proportion of settlers
META.enriched_species_prop = META.outplant_species_prop; % for comparison
% META.enriched_species_prop = [0.01 0.1 0.1 0.07 0.49 0.23]; % proportion of juveniles at Moore Reef
META.enriched_density = 6.8; % Density (m-2) after larval re-seeding on a reef (6.8 corals per m-2 corresponds to the density of 1-yr old corals
% observed in the Philippines after experimenting larval re-seeding).

% size of coral larvae (as 1yr old corals) - use the same as outplants for now
META.enriched_coral_diameter_mean = META.outplant_coral_diameter_mean;
META.enriched_coral_diameter_sd = META.outplant_coral_diameter_sd;

% Mean heat tolerance of seeded larvae relative to the contemporary mean of the group
META.MEAN_HT_larvae = RESTORATION.MEAN_HT_larvae;  % as +°C-week (DHW)
% Variance of heat tolerance of larvae
META.VAR_HT_larvae = CORAL.VAR_HT;  % as +°C-week (DHW) - same as native corals?
% Maximum heat tolerance of larvae
META.MAX_HT_larvae = CORAL.MAX_HT;  % as +°C-week (DHW) - same as native corals?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) CORAL DEPLOYMENT OPTIONS (works both for coral outplanting and larval enrichment)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.deploy_random_density = 1 ; % if 1, nb of outplants in a cell determined by Poisson distribution; if 0, nb=ceil(density)
META.threshold_for_deploying.min_cover = 0 ; % Minimum percent coral cover (all corals) for deploying on a reef
META.threshold_for_deploying.max_cover = 100 ; % Maximum percent coral cover (all corals) for deploying on a reef
META.deployment_area_prop = 0.1 ; % Maximum proportion of reef area under coral deployment 

MinDeploymentArea_km2 = 0.01; % Logistically, it's not worth deploying on a reef less than 0.01 km2 in size
% Setting deployment area to a minimum 0.01 km2 would only exclude 25 reefs in the Cairns region and 183 GBR-wide.
MinDeploymentCells = 3 ;  % minimum number of treated cells (=site within a reef) to avoid potential artefacts
% Also, this gives enough replicates for deployment (random density)

% See version 6.8 for more complex determination of number of cells to restore on a grid as a function of reef area
% Here we simply assume that what is restored is a given proportion of the grid (META.deployment_area_prop)
% All reef grids have the same number of cells (typically 400)
% And we exclude reefs that have an area less than the MinDeploymentArea_km2
GridSize = META.grid_x_count*META.grid_y_count;
Reef_ID = META.reef_ID;
Reef_area_km2 = META.area_habitat;
Reef_name = GBR_REEFS.GBR_NAME(META.reef_ID);
NumberCellsTreated = uint16(GridSize*META.deployment_area_prop*ones(META.nb_reefs,1));
DeploymentArea_km2 = META.deployment_area_prop*Reef_area_km2;

% Create table for deployment characteristics
TMP = table(Reef_ID, Reef_name, Reef_area_km2, DeploymentArea_km2, NumberCellsTreated);
I = find(TMP.Reef_area_km2 < MinDeploymentArea_km2); % find reefs than are smaller than the min deployment area -> will be excluded
TMP.DeploymentArea_km2(I) = 0;
TMP.NumberCellsTreated(I) = 0;

% %%%% From version 6.8
% GridSize = META.grid_x_count*META.grid_y_count;
% NumberCellsTreated = uint16(GridSize*MinDeploymentArea_km2./META.area_habitat);
% DeploymentArea_km2 = MinDeploymentArea_km2*ones(size(NumberCellsTreated));
% Reef_ID = META.reef_ID;
% Reef_area_km2 = META.area_habitat;
% Reef_name = GBR_REEFS.GBR_NAME(META.reef_ID);
% 
% % Create table for deployment characteristics
% TMP = table(Reef_ID, Reef_name, Reef_area_km2, DeploymentArea_km2, NumberCellsTreated);
% I = find(TMP.Reef_area_km2 < MinDeploymentArea_km2); % find reefs than are smaller than the min deployment area -> will be excluded
% J = find(TMP.NumberCellsTreated < MinDeploymentCells); % find reefs than are too big so that nb of deployed cells is too small
% TMP.NumberCellsTreated(J) = MinDeploymentCells;
% TMP.DeploymentArea_km2(J) = (MinDeploymentCells/GridSize)*TMP.Reef_area_km2(J);
% TMP.DeploymentArea_km2(I) = 0;
% TMP.NumberCellsTreated(I) = 0;


% %%%%% Intervention study 2022: PORTFOLIO Moore cluster (3D deployment areas for 1M corals at 6.8 corals per m2):
% select_ID = [9 11 12 177 178]; % ONLY WORKS WITH META.reef_ID defined for the 190 reefs around CAIRNS
% MyReefCluster = table(META.reef_ID(select_ID), Reef_name(select_ID), (2*[22059  29412 0 7352 14796]/1e6)', ...
%     'VariableNames',{'Reef_ID', 'Reef_name', 'DeploymentAreaKm2'});
% % Deployment areas for Moore, Elford, Briggs, Milln and Thetford, based on IPMF optimal site deployment
% 
% select_cluster = ismember( TMP.Reef_ID , MyReefCluster.Reef_ID);
% TMP.DeploymentArea_km2(select_cluster==0)=0;
% TMP.NumberCellsTreated(select_cluster==0)=0;
% TMP.DeploymentArea_km2(select_ID)= MyReefCluster.DeploymentAreaKm2;
% TMP.NumberCellsTreated(select_ID)= GridSize.*MyReefCluster.DeploymentAreaKm2./TMP.Reef_area_km2(select_ID);

META.coral_deployment = TMP; % Store the final table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) Options for rubble stabilisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timing of effort
META.doing_rubble_stabilisation = RESTORATION.doing_rubble_stabilisation ; % timing of rubble stab

% Restoration effort: number of reefs where rubble is stabilised at each time step
META.nb_reefs_stabilised = RESTORATION.nb_reefs_stabilised ;
META.proportion_rubble_stabilised = 1 ; % proportion of rubble that will be stabilised on a reef
META.threshold_for_stabilisation = 0 ; % Minimum percent rubble cover above which reef is selected for rubble stabilisation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) Options for Solar Radiation Management (fogging & cooling)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fogging
META.doing_fogging = RESTORATION.doing_fogging; % timing of fogging (no fogging if sum(META.doing_fogging)=0)
META.nb_reefs_fogged = RESTORATION.nb_reefs_fogged; %depreciated?
META.fogged_reef_ID = RESTORATION.fogged_reef_ID;

META.bleaching_mortality_under_fogging = 0.8 ; % coefficient applied to the calculated bleaching mortality

% If no targeted reef IDs, set priority list
META.priority_option_Fogging.region = 0 ;%[1 2 3];
META.priority_option_Fogging.shelf = 0 ;%[2 3 1];
META.priority_option_Fogging.link_strength = 0 ; %1 (increasing, min external supply first) or %2 (decreasing, max first)
META.priority_option_Fogging.link_number = 0 ; %1 (increasing, min external supply first) or %2 (decreasing, max first)
META.priority_option_Fogging.reef_area = 0 ;
META.priority_option_Fogging.focal_reef = []; % 695-Moore; 697-Elford; 698-Briggs ; 969-Milln; 970-Thetford

% Cooling (not impleented yet)
META.doing_cooling = RESTORATION.doing_cooling ; % timing of cooling (no cooling if sum(META.doing_cooling)=0)
META.cooling_factor = 0; % Set 0 if no cooling, otherwise consider the different levels (Bozec & Mumby 2019) [-0.3 ; -0.7 ; -1.3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE PRIORITY LISTS FOR RESTORATION (1 for each technique)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Currently five different criteria to prioritise reefs for restoration.

% 1) Regional
META.priority_option_Outplant.region = 0;%[1 2 3];
% Takes either 0(no regional priority) or a vector with 1(North) 2(Central) 3(South) in any priority order

% 2) Shelf position
META.priority_option_Outplant.shelf = 0;%[2 3 1];
% Takes either 0(no shelf priority) or a vector with 1(inshore) 2(midshelf) 3(outer) in any priority order

% 3) Likelihood of larval sink as the sum of all inbound link strengths
META.priority_option_Outplant.link_strength = 0 ; %1 (increasing, min external supply first) or %2 (decreasing, max first)
% Takes either 0(no priority), 1(decreasing, highest potential first) or 2(increasing, lowest potential first)

% 4) Likelihood of larval sink as the total number of inbound links
META.priority_option_Outplant.link_number = 0 ; %1 (increasing, min external supply first) or %2 (decreasing, max first)
% Takes either 0(no priority), 1(decreasing, highest potential first) or 2(increasing, lowest potential first)

% 5) Reef size (area)
META.priority_option_Outplant.reef_area = 0 ;
% Takes either 0(no priority), 
% Using GBRMPA reef outline as habitat area: 1(increasing, smallest area first) or 2(decreasing, largest area first)
% Using Gemorphic map as habitat area: 3(increasing, smallest area first) or 4(decreasing, largest area first)

% 5) Focused on a group of reef (from a focal reef, visit all reefs at increasing distance)
META.priority_option_Outplant.focal_reef = 695 ; %Moore Reef

%% Populate same options for Rubble Stabilisation...
META.priority_option_RubbleStab = META.priority_option_Outplant;
%... And refine if necessary:
% META.priority_option_RubbleStab.region = 0 ; % [1 2 3]
% META.priority_option_RubbleStab.shelf = 0 ; % [1 2 3]
% META.priority_option_RubbleStab.link_strength = 0 ;
% META.priority_option_RubbleStab.link_number = 0 ;
% META.priority_option_RubbleStab.reef_area = 0 ;

%% Populate same options for Larval enrichment...
META.priority_option_LarvalEnrich = META.priority_option_Outplant;
%... And refine if necessary:
% META.priority_option_LarvalEnrich.region = 0 ; % [1 2 3]
% META.priority_option_LarvalEnrich.shelf = 0 ; % [1 2 3]
% META.priority_option_LarvalEnrich.link_strength = 0 ;
% META.priority_option_LarvalEnrich.link_number = 0 ; 
% META.priority_option_LarvalEnrich.reef_area = 0 ; 

