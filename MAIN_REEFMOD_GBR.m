%__________________________________________________________________________
%
% REEFMOD-GBR MAIN SCRIPT, version 7.2
% Yves-Marie Bozec, y.bozec@uq.edu.au, Mar 2024

% This is an incremental version of 7.0 (which was used for simulating FY2024 counterfactuals - Feb 2024)
% and version 7.1 (first attempt at implementing wave-energy based cyclone refugia - March 2024).
% This new version takes the last developments of CoTS control modelling from Tina, including the reduced mortality
% in non protected zones, new list for targeted/priority reefs, new manta tow CoTS observations (excluding RHIS)
% Code for restoration intervention used in 2022 has been tidied up.

%__________________________________________________________________________
clear

SaveDir ='';

NB_SIMULATIONS = 1; % Number of repeated runs

% Always run the hindcast before future projections. Initialisation = winter 2007
% NB_TIME_STEPS has to be an even number, eg. 32 (hindcast 2008-2023) + 154 (forecast 2024-2100)

% NB_TIME_STEPS = 32; % HINDCAST: summer 2008 to winter 2023
NB_TIME_STEPS = 32+154; % HINDCAST+FORECAST summer 2008 - winter 2100

OutputName = 'R0_GBR.7.2'; options = [1 1 1 1 0 1 0.3]; % see options below

%% select the Global Circulation Model for climate change projection (CMIP-6)
GCM = 1; 
% 1=CNRM-ESM2-1, 2=EC-Earth3-Veg, 3=IPSL-CM6A-LR, 4=MRI-ESM2-0, 5=UKESM1-0-LL, ...
% 6=GFDL-ESM4, 7=MIROC-ES2L, 8=MPI-ESM1-2-HR, 9=MIROC6, 10=NorESM2-LM

SSP = 3; % 1=SSP1-1.9, 2=SSP1-2.6, 3=SSP2-4.5, 4=SSP3-7.0, 5=SSP5-8.5
% Note SSP1-1.9 is not available for GFDL-ESM4, MPI-ESM1-2-HR and NorESM2-LM

%% --------------------------------------------------------------------------------
% Climate forecasts from CMIP6 global circulation models
All_GCMs = ["CNRM-ESM2-1" ; "EC-Earth3-Veg" ; "IPSL-CM6A-LR" ; "MRI-ESM2-0" ; "UKESM1-0-LL" ; ...
        "GFDL-ESM4" ; "MIROC-ES2L" ; "MPI-ESM1-2-HR" ; "MIROC6" ; "NorESM2-LM" ];
All_SSPs = ["119" ; "126" ; "245" ; "370" ; "585" ];

OPTIONS.GCM = All_GCMs(GCM);
OPTIONS.SSP = All_SSPs(SSP);

%% --------------------------------------------------------------------------------
% Stressor options: yes(1)/no(0)
OPTIONS.doing_cyclones = options(1);
OPTIONS.doing_bleaching = options(2) ; 
OPTIONS.doing_COTS = options(3);
OPTIONS.doing_WQ = options(4);
OPTIONS.doing_restoration = options(5) ;
OPTIONS.doing_COTS_control= options(6);
OPTIONS.heritability_HT = options(7); % heritability of heat tolerance for all groups (inheritance model of HT phenotypes)

OPTIONS.doing_size_frequency = 1; % for tracking population size structure (incl. juveniles)
OPTIONS.doing_genetics = 0 ; % for running the genetic model (Bozec & Mumby 2019) - needs revision
OPTIONS.genetic_parms = [ 1 1 1.5 ]; % sigma cold/sigma hot/esd (Bozec & Mumby 2019) - useless if doing_genetics = 0

%% --------------------------------------------------------------------------------
% Below options are used to force simulations with specific starting conditions - keep empty if of no use
OPTIONS.init_coral_cover = []; %0.01*ones(1,6); % as proportional cover (vector of 6 values)
OPTIONS.init_sand_cover = []; % as proportional cover 
OPTIONS.init_rubble_cover = []; % as proportional cover 
OPTIONS.ssc = []; %0.1; %in mg/L

%% --------------------------------------------------------------------------------
% RESTORATION INTERVENTIONS (high level controls; set parameters in settings_RESTORATION)

% 1) Outplanting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timing of outplanting (set to 0 for not doing this intervention)
RESTORATION.doing_coral_outplanting = uint8(zeros(1,NB_TIME_STEPS)); % if zero, no coral deployment at time step for all reefs
% RESTORATION.doing_coral_outplanting(1,38:2:46) = 1 ; % set to 1 to indicate outplanting: first in summer 2026 and last in summer 2030 (BCA)
% RESTORATION.doing_coral_outplanting(1,[38:2:46 78:2:86 118:2:126 158:2:166] ) = 1 ; % outplanting starts 2026-2030, 2046-2050, 2066-2070, 2086-2090

RESTORATION.total_nb_outplants = 100e6; % Nb of outplants available for the GBR at each deployment (= time step). Set to Inf if density/nb_reefs are the drivers
RESTORATION.outplanted_density = 6.8;  % only in the case of fixed density of outplants (ignores RESTORATION.total_nb_outplants?)
RESTORATION.MEAN_HT_outplants = 2*[1 1 1 1 1 1];% as +°C-week (DHW) compared to the contemporary mean of the group

% 2) Larval enrichment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timing of larval enrichment (set to 0 for not doing this intervention)
RESTORATION.doing_larval_enrichment = uint8(zeros(1,NB_TIME_STEPS));
% RESTORATION.doing_larval_enrichment(1,38:2:46) = 1 ; % set to 1 to do outplanting
% RESTORATION.doing_larval_enrichment(1,2:2:30) = 1 ;

RESTORATION.nb_reefs_enriched = 10 ; % max number of reefs where larval enrichment is undertaken at each time step
% Will be picked from a priority list generated in f_generate_priority_list OR specified in SETTINGS_RESTORATION)
% If 0, enrichment cannot happen. For the counterfactual, set to Inf with total_nb_larvae = 0 for ghost deployment
RESTORATION.total_nb_larvae = Inf; % Max number of 'larvae' (ie, 1 yr old corals) available at each time step. Set to Inf if unlimited.
RESTORATION.MEAN_HT_larvae = 2*[1 1 1 1 1 1];% as +°C-week (DHW) compared to the contemporary mean of the group

% 3) Rubble stabilisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the restoration effort: number of reefs where rubble is stabilised at each time step
RESTORATION.nb_reefs_stabilised = 0 ; % (if 0 rubble stabilisation cannot happen) 
% Set the timing of intervention (if 1 intervention is deployed at step t, if 0 no intervention at t)
RESTORATION.doing_rubble_stabilisation = uint8(zeros(1,NB_TIME_STEPS));
% RESTORATION.doing_rubble_stabilisation(1,38:1:end) = 1 ; % set to 1 to do outplanting: first in summer 2026 and last in summer 2030

% 4) Solar Radiation Management  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fogging for RRAP 2022 intervention simulations
RESTORATION.doing_fogging = uint8(zeros(1,NB_TIME_STEPS)); %if zero don't do fogging at time step
% RESTORATION.doing_fogging(1,32:2:end) = 1 ; % set to 1 to do fogging (only in summer)

RESTORATION.nb_reefs_fogged = []; % number of reefs if random deployment, otherwise =[] but need to listthe  reef IDs
RESTORATION.fogged_reef_ID = [695 697 969 970]; % using 5 fogging units (22.5 km2) covering 23 km2 (equivalent 2D reef areas)
% 695-Moore; 697-Elford; 698-Briggs ; 969-Milln; 970-Thetford
% Leave empty otherwise []

% Cloud brightening (=cooling) from Bozec & Mumby (2019) (needs revision)
RESTORATION.doing_cooling = uint8(zeros(1,NB_TIME_STEPS)); %if zero don't do cooling at time step
% RESTORATION.doing_fogging(1,32:2:end) = 1 ; % set to 1 to do cooling (only in summer)


%% --------------------------------------------------------------------------------
%% Saving options
if NB_TIME_STEPS > 32
    
    OPTIONS.OutputFileName = [SaveDir OutputName '_herit' num2str(OPTIONS.heritability_HT) '_SSP' char(OPTIONS.SSP) '_' char(All_GCMs(GCM)) '.mat'];
    % Display the running scenario
    ['Running scenario ' OutputName '_herit' num2str(OPTIONS.heritability_HT) '_SSP' char(OPTIONS.SSP) '_' char(All_GCMs(GCM)) ' .....']
else
    OPTIONS.OutputFileName = [SaveDir OutputName '_herit' num2str(OPTIONS.heritability_HT) '_HINDCAST.mat' ];
    % Display the running scenario
    ['Running scenario ' OutputName '_herit' num2str(OPTIONS.heritability_HT) '_HINDCAST'  ' .....' ]

end


%% --------------------------------------------------------------------------------
%% RUN REEFMOD
OUTPUTS = struct('REEF', [],'RESULT', [],'RECORD', []);
TEMP_META = struct('META', []);

% parfor run_id = 1:NB_SIMULATIONS
for run_id = 1:NB_SIMULATIONS
run_id
    [meta, REEF, RESULT, RECORD] = f_multiple_reef(OPTIONS, RESTORATION, NB_TIME_STEPS, run_id);
      
    OUTPUTS(run_id).RESULT = RESULT ;
    OUTPUTS(run_id).RECORD = RECORD ;
    OUTPUTS(run_id).REEF = REEF ;
    TEMP_META(run_id).META = meta ;
end

META = TEMP_META(1).META; % Keep only one META because common to all simulations
clear TEMP_META ADAPT run_id meta REEF RESULT RECORD GCM GCM_list options SaveDir OutputName RCP RCP_list


%% --------------------------------------------------------------------------------
%% FORMAT, OPTIMISE & SAVE
%% --------------------------------------------------------------------------------
%% Memory allocation
% 1) coral outputs
coral_cover_per_taxa = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');
coral_larval_supply = coral_cover_per_taxa;
nb_coral_offspring = coral_cover_per_taxa;
nb_coral_recruit = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'uint16');
total_shelter_volume_per_taxa = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');
coral_HT_mean = nan(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');
coral_HT_var = nan(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');
selection_diff = nan(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');

if OPTIONS.doing_size_frequency == 1   
    % Bin edges for juvenile bins (cm): 1  3  5 (2 size classes)
    % Adolescents (cm):  5  9  13  17 (3 size classes)
    % Adults (cm): 17 31 45 59 73 87 101 (6 size classes)
    nb_coral_juv = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 2, 'uint16') ;
    nb_coral_adol = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 3, 'uint16') ;
    nb_coral_adult = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 6, 'uint16') ; 
end

% 2) Coral cover loss following stressors
coral_cover_lost_bleaching = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps, META.nb_coral_types, 'single');
coral_cover_lost_cyclones = coral_cover_lost_bleaching;
coral_cover_lost_COTS = coral_cover_lost_bleaching;

% 3) Stress records
record_applied_DHWs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'single');
record_applied_cyclones = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'uint16');

if exist('RECORD.applied_bleaching_mortality')==1
    record_applied_bleaching_mortality  = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'single');
end

% 4) Restoration records
if OPTIONS.doing_restoration==1
    
    % Total nb of outplants per reef per time step
    record_total_outplants_deployed = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps, META.nb_coral_types,'uint16');
    record_outplanted_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'uint16');
    
    record_total_larvae_deployed = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps, META.nb_coral_types,'uint16');
    record_enriched_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'uint16');  
    
    record_rubble_pct2D_stabilised = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'single');
    record_stabilised_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'uint16');
    
    record_fogged_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'uint16');
 
    coral_cover_per_taxa_restored_sites = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');
    
else
    clear RESTORATION
end

% 5) Other variables
rubble = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');
nongrazable = zeros(NB_SIMULATIONS, META.nb_reefs,'single');
macroTurf = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');
macroEncrustFleshy = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');
macroUprightFleshy = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');

% 6) CoTS outputs
if OPTIONS.doing_COTS == 1
    COTS_mantatow = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 'single');
    COTS_densities = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 16, 'single'); % 16 age classes
    COTS_predicted_densities = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 16, 'single'); % 16 age classes
    COTS_predicted_mantatow = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 'single');
    COTS_settler_density = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 'uint16');
    COTS_larval_supply = COTS_mantatow;
    COTS_larval_output = COTS_mantatow;
end

if OPTIONS.doing_COTS_control == 1 && META.nb_time_steps > 23
    COTS_CONTROL_culled_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps); % gives which reefs were culled or not (1/0)
    COTS_CONTROL_remaining_dives = zeros(NB_SIMULATIONS, META.nb_time_steps); %remaining number of control dives available after intervention
    COTS_CONTROL_culled_density = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps); % density extracted by the intervention (adults)
end

%% Populate outputs
for simul = 1:NB_SIMULATIONS
    
    coral_cover_per_taxa(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_pct2D));
    coral_larval_supply(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_larval_supply)); % nb of incoming larvae per unit of reef area (400m2)
    nb_coral_recruit(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_settler_count));
    nb_coral_offspring(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_total_fecundity)); % nb of larvae produced per unit of reef area (400m2)
    total_shelter_volume_per_taxa(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.total_shelter_volume_per_taxa));
    
    coral_HT_mean(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_HT_mean));
    coral_HT_var(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_HT_var));
    selection_diff(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.selection_diff));

    if OPTIONS.doing_size_frequency == 1
        nb_coral_juv(simul,:,:,:,:)= squeeze(cat(4,OUTPUTS(simul).RESULT.coral_juv_count(:,:,:,:)));
        nb_coral_adol(simul,:,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_adol_count(:,:,:,:)));
        nb_coral_adult(simul,:,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_adult_count(:,:,:,:)));
    end
    
    coral_cover_lost_bleaching(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_bleaching);
    coral_cover_lost_cyclones(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_cyclones);
    coral_cover_lost_COTS(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_COTS);
    
    record_applied_cyclones(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.hurricane_events);
    record_applied_DHWs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.applied_DHWs);
    
    if exist('RECORD.applied_bleaching_mortality')==1
        record_applied_bleaching_mortality(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.applied_bleaching_mortality);
    end
    
    if  OPTIONS.doing_restoration==1
        record_total_outplants_deployed(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.total_outplanted);
        record_outplanted_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.outplanted_reefs);
        record_total_larvae_deployed(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.total_enriched);
        record_enriched_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.enriched_reefs);
        record_rubble_pct2D_stabilised(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.rubble_cover_pct2D_stabilised(:,1:end));
        record_stabilised_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.stabilised_reefs);
        record_fogged_reefs(simul,:,:) = RECORD.fogged_reefs;
        
        coral_cover_per_taxa_restored_sites(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_pct2D_restored_sites));
    end
    
    rubble(simul,:,:)= squeeze(OUTPUTS(simul).RESULT.rubble_cover_pct2D);
    nongrazable(simul,:) = squeeze(cat(4,OUTPUTS(simul).REEF.nongrazable_substratum));
    
    macroTurf(simul,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,4)));
    macroEncrustFleshy(simul,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,2)));
    macroUprightFleshy(simul,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,3)));
      
    if OPTIONS.doing_COTS == 1
        % Assuming 0.6 CoTS per grid ~ 0.22 CoTS per tow
        % (0.22 per tow is equivalent to 1500 COTS per km2 (Moran & De'ath 92), so that 1 COTS per grid (400m2) is equivalent to 0.22*2500/1500
        COTS_mantatow(simul,:,:) = (0.22/0.6)*squeeze(cat(4,OUTPUTS(simul).RESULT.COTS_total_perceived_density));
        COTS_densities(simul,:,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_all_densities); % Density for 400m2
        COTS_predicted_densities(simul,:,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_all_densities_predicted); % Density for 400m2
        COTS_predicted_mantatow(simul,:,:) = (0.22/0.6)*squeeze(OUTPUTS(simul).RESULT.COTS_total_perceived_density_predicted);
        
        COTS_settler_density(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_settler_densities); % Density for 400m2
        COTS_larval_supply(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_larval_supply); % Density for 400m2
        COTS_larval_output(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_larval_output); % Density for 400m2
    end
    
    if OPTIONS.doing_COTS_control == 1
        COTS_CONTROL_culled_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_culled_reefs);
        COTS_CONTROL_remaining_dives(simul,:) = squeeze(OUTPUTS(simul).RESULT.COTS_control_remaining_dives);
        COTS_CONTROL_culled_density(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_culled_density);
        
        if META.report_COTS_control==1
            CONTROL_RECORD(simul).COTS_control_records = OUTPUTS(simul).RESULT.COTS_records;
        end
    end
end

if META.doing_COTS_control == 1 && META.report_COTS_control==0
    META=rmfield(META,{'boatProperties','COTS_control_sites','cntrl_reefID','cntrl_sites'});
end

% Express shelter volume as total across all coral group in dm3 per m2
reef_shelter_volume_absolute = sum(total_shelter_volume_per_taxa,4)/(META.total_area_cm2/1e4);

% Express shelter volume relative to maximum possible value (ie, 1 tabular coral of the size of the cell in every cell)
total_shelter_volume_max = META.grid_x_count*META.grid_y_count*exp(-8.32 + 1.50*log(META.cell_area_cm2));
reef_shelter_volume_relative = sum(total_shelter_volume_per_taxa,4)/total_shelter_volume_max;
reef_shelter_volume_relative(reef_shelter_volume_relative>1)=1;
reef_shelter_volume_relative(reef_shelter_volume_relative<0)=0;

% New (08/2021): only record COTS densities by yearly classes to reduce output size
% COTS_densities = COTS_densities0(:,:,:,1:2:end)+COTS_densities0(:,:,:,2:2:end);
% 09/2021: now just sum across all juveniles and across all adults
if OPTIONS.doing_COTS == 1
    COTS_juv_densities = sum(COTS_densities(:,:,:,1:(META.COTS_adult_min_age-1)),4);
    COTS_adult_densities = sum(COTS_densities(:,:,:,META.COTS_adult_min_age:end),4);
end

clear OUTPUTS simul s NB_TIME_STEPS

save (OPTIONS.OutputFileName)
