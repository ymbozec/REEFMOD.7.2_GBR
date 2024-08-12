% REEFMOD-GBR model run script
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 06/2023
%__________________________________________________________________________

function [META, REEF, RESULT, RECORD] = f_multiple_reef(OPTIONS, RESTORATION, nb_time_steps, simul)

PARAMETERS
META.nb_time_steps = nb_time_steps;

%% Reef areas
load('GBR_REEF_POLYGONS_2024.mat') % with updated estimates of ungrazable substratum from benthic maps 10/2023
% New reef areas based on the 3D surface areas of geomorphic classes down to 20m depth (Roelfsema et al. 2021 Rem. Sens.).
% Available for all geomorphic classes (Geom_total) and coral geomorphic classes (Geom_CH). 
% Also contains reef specific estimates of UNGRAZABLE (from benthic maps, Roelfsema et al. 2021),
% and corrected assignment of shelf position (Caro) and zoning status (Tina)

%% Reef selection
META.reef_ID = [1:3806]'; % Entire GBR
% META.reef_ID = GBR_REEFS.Reef_ID(GBR_REEFS.LAT<-15.4 & GBR_REEFS.LAT >-15.7); % Region around Moore reef
% META.reef_ID = GBR_REEFS.Reef_ID(GBR_REEFS.LAT<-15.76 & GBR_REEFS.LAT >-17.34); % Cairns region reduced (190 reefs) for restoration

META.nb_reefs = length(META.reef_ID);
META.outside_reef_ID = GBR_REEFS.Reef_ID(ismember(GBR_REEFS.Reef_ID,META.reef_ID)==0); % empty when all 3,806 reefs are included
META.reef_lat = GBR_REEFS.LAT(META.reef_ID); % required for CoTS control
META.reef_lon = GBR_REEFS.LON(META.reef_ID); % required for CoTS control

% Define which habitat area to use 
META.area_habitat = GBR_REEFS.Geom_CH_km2(META.reef_ID); % 3D area of coral habitat (CH) geomorphic classes 
% META.area_habitat = GBR_REEFS.Geom_total_km2(META.reef_ID); % 3D area of all geomorphic classes
% META.area_habitat = GBR_REEFS.Reference_Area_km2(META.reef_ID); % 2D area based on GBRMPA reference reef outline 

%% Set up bleaching
META.doing_bleaching = OPTIONS.doing_bleaching ;

% Adjust bleaching to align with the shallow (-2m) mortality recorded by Hughes et al. (2018) or deep (-7m) as in Baird et al. (2018)
% CORAL.bleaching_depth = 1 % shallow bleaching as simulated in Bozec et al. (2022)
% Now calculated at depth from empirical data (Baird et al. 2018 MEPS)
DEPTH = 7; % representative depth (m) for bleaching-induced mortality calculations (use DEPTH = 2 to simulate the equivalent of Hughes' mortalities)
CORAL.bleaching_depth = (0.420 + 0.272*DEPTH)^-1;

% Parameters of heat tolerance (assume for now same value for all taxonomic groups)
CORAL.MAX_HT = 8 * [1 ; 1 ; 1 ; 1 ; 1 ; 1 ]; % limits of heat tolerance that can be achieved in simulations (+/- 8 degC-week)
CORAL.heritability_HT = OPTIONS.heritability_HT * [1 ; 1 ; 1 ; 1 ; 1 ; 1 ]; % as proportion of phenotypic variance due to additive genetic variance

%% Set up cyclones
META.doing_hurricanes = OPTIONS.doing_cyclones ;

% 03/2024 - randomise the incidence of a predicted cylone (based on orbital velocity thresholds):
META.randomize_hurricane_strike = 0; 
% if no (=0), REEF.hurricane_strike_incidence = 1 (= strike probability for a reef)
% if yes (=1), REEF.hurricane_strike_incidence gets reef specific proba (STRIKE_INCIDENCE in settings_GBR) 

%% Set up CoTS
% NOTE: don't refine COTS parameters here because they will be erased by settings_COTS
META.doing_COTS = OPTIONS.doing_COTS ;
META.doing_COTS_control = OPTIONS.doing_COTS_control ;
META.report_COTS_control = 0 ; % Put 1 to record the detailed results of CoTS control (will create RESULT.COTS_records)

REEF_COTS =[]; % (required, even if not doing CoTS)

%% Set up WQ
META.doing_water_quality = OPTIONS.doing_WQ;
META.doing_Chl_forcing = 1; % with (1) or without (0) CoTS larval survival driven by Chlorophyll concentration

%% Track colony size distributions
META.doing_size_frequency = OPTIONS.doing_size_frequency; % Track coral colonies sizes (SLOW!)

%% Connectivity
META.doing_coral_connectivity = 1 ;
% META.coral_immigration = 1e6*ones(1,META.nb_time_steps) ; % Forced larval input for a 400m2 area - works only if connectivity is OFF

META.recruitment_type = 1; % turn into 0 for fixed recruitment (but then connect and genetics won't work)
% See ReefMod.6.8 for an example of forced larval input for Moore Reef cluster

% Parameter a of the B-H function (same for all reefs), calibrated with
CORAL.BH_alpha = 15*CORAL.prop_settlers; % per m2
CORAL.BH_beta = 5*1e6*ones(6,1); % for a 400m2 reef

% Force self-seeding of coral larvae
META.coral_min_selfseed = 0.28 ; % relative proportion of produced larvae that stay on the reef (Bozec et al. 2022)

%% Rubble (standard)
META.tracking_rubble = 1;
META.rubble_decay_rate = 0.128 ; % 2/3 stabilised after 4 years
META.convert_rubble = 1; % Conversion factor from coral loss to rubble cover
META.convert_rubble_lag = 6; % 3 years delay for converting coral loss due to bleaching and CoTS (delayed structural loss).

%% Grazing
REEF.herbivory = 1; % full grazing
ALGAL.nb_step_algal_dynamics = 1 ; %%%%% ONLY TO SPEED-UP THE CODE WHEN FULL GRAZING (otherwise set to 6)

%% Restoration
META.doing_restoration = OPTIONS.doing_restoration;

%% Genetic adaptation as in Bozec & Mumby 2019 (needs revision)
% not to be confused with the inheritance model of heat tolerance
META.doing_genetics = OPTIONS.doing_genetics ;

if META.doing_genetics==1
    
    settings_GENETICS;
    
    META.genetics.SIGMA_COLD = OPTIONS.genetic_parms(1) + [ 0 0 0 0 0 0 ]; % for the cold side (when temp<Topt)
    META.genetics.SIGMA_HOT = OPTIONS.genetic_parms(2) + [ 0 0 0 0 0 0 ]; % for the hot side (when temp>Topt)
    META.genetics.esd = OPTIONS.genetic_parms(3) + [ 0 0 0 0 0 0 ]; % SD of environmental effect on fitness (mean=0), on top of genetics.
    META.genetics.enhanced_tolerance = OPTIONS.thermal_tolerance_outplants ;
    
    % Load the pre-adapted pool of QTL
    load(['QTL_pool_sigma_c' num2str(OPTIONS.genetic_parms(1)) '_h' num2str(OPTIONS.genetic_parms(2))...
        '_esd' num2str(OPTIONS.genetic_parms(3)) '.mat']);
    
    CORAL.growth_rate = CORAL.growth_rate/mean_fitness; % average fitness across the GBR after burn-in
    % this allows adjusting fitness so that mean individual growth rate is close to ReefMod's default.
    % Note this average fitness is relative to the local environment corals are adapted to, ie we are not 
    % modelling latitudinal differences in growth rates (100% fitness in the far South gives the same growth rate than 
    % in the far North, while in reality coral growth declines at with latitude
    
end

%% INITIALISATION
rng(simul); % to get a repeatable scheme of random number generation in RAND, RANDI, RANDN

INITIALISATION

settings_GBR

if META.doing_restoration == 1
    
    settings_RESTORATION;
    
    % Then generate priority lists for each technique = list of reef ID sorted from highest to lowest priority
    % Note the list is set only once and remains the same throughout the simulation
    MY_REEFS = GBR_REEFS(META.reef_ID,:);
%     META.priority_list_Outplant = f_generate_priority_list_NEW(META.priority_option_Outplant, MY_REEFS, CONNECT);
%     META.priority_list_RubbleStab = f_generate_priority_list_NEW(META.priority_option_RubbleStab, META.reef_ID, MY_REEFS, CONNECT);
%     META.priority_list_LarvalEnrich  = f_generate_priority_list_NEW(META.priority_option_LarvalEnrich, MY_REEFS, CONNECT);
%     META.priority_list_Fogging = f_generate_priority_list_NEW(META.priority_option_Fogging, MY_REEFS, CONNECT);

    % ONLY FOR THE BCA: focus on the two reef clusters and deploy fixed density of outplants
    % Deployment areas are set in settings_RESTORATION
    META.priority_list_Outplant = find(META.coral_deployment.DeploymentArea_km2>0);
    META.priority_list_LarvalEnrich = find(META.coral_deployment.DeploymentArea_km2>0);
    META.priority_list_RubbleStab = 1:length(META.reef_ID); % do it everywhere for the moment
    
    reef_fogging_list = ismember(META.coral_deployment.Reef_ID,META.fogged_reef_ID);
    META.priority_list_Fogging = find(reef_fogging_list == 1);
    
else % This is needed in MULTIPLE_REEF_SETUP
    
    META.doing_cooling=0;
    META.cooling_factor = 0; 
    
end

% Adjustments in case of coral outplanting
if META.nb_coral_types > 6
    
    % Add initial cover for outplants (0%)
    init_coral_cover = [init_coral_cover zeros(META.nb_reefs, META.nb_coral_types-6)];
    
    if META.doing_COTS == 1
        % Extend the vector of CoTS preferences
        X = META.COTS_feeding_prefs;
        X = [X ; X];
        META.COTS_feeding_prefs = X/sum(X); % feeding prefs sum to 1
%         META.COTS_pref_corals = [META.COTS_pref_corals META.COTS_pref_corals(META.outplanted_species)+5]; % don't need to change this
    end
end

%% REFINE IF SPECIFIC FORCING CONDITIONS (eg, for building look-up table for the ReefMod Engine)
% (only works with single reef simulations)
if isempty(OPTIONS.init_coral_cover)==0
    init_coral_cover = OPTIONS.init_coral_cover;  
end

if isempty(OPTIONS.init_sand_cover)==0
    init_sand_cover = OPTIONS.init_sand_cover;  
end

if isempty(OPTIONS.init_rubble_cover)==0
    init_rubble_cover = OPTIONS.init_rubble_cover;  
end    
    
if isempty(OPTIONS.ssc)==0 
    CORAL_recruit_survival = (1 - 1.88*0.001*OPTIONS.ssc)^(180/40);
    CORAL_juvenile_growth = 1 - 0.176*log(OPTIONS.ssc+1);
    
    FERT = exp(4.579 - 0.010*OPTIONS.ssc)/exp(4.579); % fertilization success
    SETT = (99.571 - 10.637*log(OPTIONS.ssc+1))/99.571; % settlement success
    CORAL_larvae_production =  FERT*SETT;
    
    % Store in REEF_POP (same value for all years)
    for y=1:size(REEF_POP,2)
        
        REEF_POP(y).CORAL_recruit_survival = CORAL_recruit_survival;
        REEF_POP(y).CORAL_juvenile_growth = CORAL_juvenile_growth;
        REEF_POP(y).CORAL_larvae_production = CORAL_larvae_production;    
    end
end

MULTIPLE_REEF_SETUP

%% REFINE BACKGROUND COTS DENSITY ON SPECIFIC REEFS
if META.doing_COTS == 1
    for n=1:length(META.reef_COTS_INIT_BOX)
        id_reef = find(META.reef_ID==META.reef_COTS_INIT_BOX(n));
        
        if isempty(id_reef) == 0
            REEF(id_reef).COTS_background_density = META.COTS_background_density_INIT_BOX;
        end
    end
end

%% APPLY CHLOROPHYLL FORCING OF COTS LARVAL SURVIVAL ONLY ON INSHORE REEFS (ONLY WORKS IF doing_Chl_forcing=1)
if META.doing_Chl_forcing == 1 && META.doing_COTS == 1
    for n=1:META.nb_reefs
        if GBR_REEFS.Shelf_position(META.reef_ID(n)) > 1
            for yr=1:size(REEF_POP,2)
                REEF_POP(yr).COTS_larvae_survival(n) = META.COTS_min_larval_survival;
            end
        end
    end
end

%% If doing COTS control 
if META.doing_COTS_control == 1
    settings_COTS_CONTROL % %Tina: load in new parameter file with all control info
end
             
clearvars -except META REEF CORAL ALGAL CONNECT REEF_POP REEF_COTS

[RESULT, RECORD] = f_runmodel(META, REEF, CORAL, ALGAL, CONNECT, REEF_POP, REEF_COTS) ;
