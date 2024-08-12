%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tina Skinner, MSEL, June 2023.
%
% Parametrisation of COTS control settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

META.do_CSIRO_COTSctrl= 1; %whether we do CSIRO-like algorithms for COTS control (though not aware of any other algorithms)

META.COTS_global_trigger = 0; %check global trigger for CSIRO COTS control - not sure what this refers to

META.COTS_reef_trigger = 0; %check reef-level trigger for CSIRO COTS control - not sure what this refers to

META.COTS_control_start = 22; %timestep in 6-month intervals when control should start at 22, e.g. to allow for burning in of the model. This way control starts at start of 2019 which is 23.

%%% CHOOSING REEFS TO CONTROL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.COTS_reefs2cull_strat = 1;  %choose method of selecting reefs to manage from f_makeReefList_TS

META.min_control_cover = 20; %minimum coral cover needed for control to still happen when high COTS in f_make_ReefList_TS

% Cull sites: Tina April 2023 new based on Geom_CH_2D_km2. For each reef, 2D coral habitat area in first column, no. cull sites in second. 
CS=load('COTS_sites_new.mat'); %Note, all cull sites integers or model crashes. All 0's rounded up to 1 as number of sites per reef should not be 0, given it's calculated on the area of 2D coral habitat.
META.COTS_control_sites=CS.COTS_sites;

Newlist=[META.reef_ID META.COTS_control_sites(META.reef_ID,2)]; %%Add index to track other metrics %% YM correction for sizeable code
META.cntrl_sites=rmmissing(Newlist); %%Remove missing reefs
META.cntrl_reefID=[META.cntrl_sites(:,1)];%% get the reefID for reefs within MPA

META.COTS_fixed_list = 1; %Specifies whether boat order for reef visitation is fixed or not

META.max_COTS = 3; %Specifies max COTS per tow above which control wouldn't happen as too many.

%META.targetreef_n = 300; July 23 - now have fixed target list. Can update.

%%% CONTROL PROGRAM EFFORT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

META.COTS_cull_boats = 5; % Specify the available effort in number of boats and days per boat; develop surveys later
META.COTS_cull_days = 90;  %Number of boat days at sea - 100 days/6months per AMPTO; some also lost to travel
META.COTS_cull_voyages = 9;  %number of discrete voyages per boat per 6 months - 13 per AMPTO - reduced to 9 to align with effort from actual Control Program
META.divers = 8; % number of divers per vessel
META.diven = 4; % number of dives per day

%META.COTS_pref_coral_groups=1:4;%which coral groups are taken into account for ET
META.COTS_postcontrol_proportions=[ 1 1 (1-META.COTS_detectability(3:end))./(sum(1-META.COTS_detectability(3:end)))];%this is the population structure at ET after control; based on detectability, i.e. survivng population=1-detectability

%ecological threshold above which a reef must be to be treated - not currently used
%META.COTS_ecological_threshold=0.22;

%%% GIVING BOATS PROPERTIES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assign properties to each boat

boats = META.COTS_cull_boats; %Specify boats using number of cull boats

for i = 1:boats %here you could remove the loop, if numb is a scalar integer
    %identifyTheBoat
    META.boatProperties.boatID = 1:boats ;
    %alotted boat days per 6 months
    META.boatProperties.boatDays = repmat(META.COTS_cull_days,1,boats);%
    %alotted voyages per 6 months
    META.boatProperties.voyages = repmat(META.COTS_cull_voyages,1,boats);%
    %number of divers onboard 
    META.boatProperties.divers = repmat(META.divers,1,boats); %
    %boat visit reefs in specific fixed order (1) of their ranking, or choose randomly (0) from the top X/reefs2cull reefs where to go first
    META.boatProperties.fixedOrder=repmat(META.COTS_fixed_list,1,boats);
end

%%% CULLING EFFORT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating culling effort for each vessel

META.control_effort_allocation = 0.9; %Proportion of control effort to go to culling after some allocated to mantatowing/RHIS. Here 10% to manta towing/RHIS

for i=1:boats %for each boat
        
META.boatProperties.totalTeamDives(i)=META.boatProperties.boatDays(i)*META.diven*META.control_effort_allocation;%total dives a team can make; each should be on a new site
META.boatProperties.totalInidvDives(i)=META.boatProperties.boatDays(i)*META.diven*META.boatProperties.divers(i)*META.control_effort_allocation;%tdives of individual divers; note that this assumes divers fromt eh same boat can all be on different sites on the same reef durign the same dive which is probably unrealistic

end

META.max_dives_per_site = 300; %For high effort reefs, specify stopping rule threshold number of dives at cull site level - for hours convert * (40/60) = 200 hours

META.max_dives_per_reef = 3000; % For high effort reefs, specify stopping rule threshold number of dives at reef level - for hours convert * (40/60) = 2000 hours

%%% NOT IMPLEMENTED YET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify whether full state of the system is known; implement partial knowledge later - not currently developed
%META.COTS_state_known=1; 

% cull surrounding reefs - not currently developed
%META.COTS_cull_surrounding = 0; 

%Number of survey boats - not developed yet
%META.COTS_survey_boats = 0;  

%Number of survey voyages - not developed yet
%META.COTS_survey_voyages = 0;  

%Is this vessel a cull boat (1) or a survey vessel (0) currently not developed
%META.boatProperties.mission = repmat(1,1,boats);%randi([0 1],1,boats);

%the location on the coast that the vessel is stationed: 1 = Port Douglas, 2 = Cairns, 3 = Townsville, 4 = Mackay
%homeports=[2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4]; % extend this vector if your wish to have more than 8 vessels
%homeports=[2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4]; % extended for modelling up to 40 vessels
%META.boatProperties.homePort = homeports(1:boats);%randi([1 4],1,boats); Doesn't seem to be used?

%META.boatProperties.maxDays_atSea = repmat(13,1,boats);%   %maximum days the vessel can be at sea at one time

%%% OBSOLETE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No record of these parameters throughout the model, have commented them and left them here for now in case 

%META.listFileName = 1;   

%choose control strategy - obsolete? Not used anymore it seems.
%META.COTS_control_strat = 3;   

%picking the highest priority reef
%META.top_reef_picks=1;  %What is this?
%META.COTS_top_reef_picks = 1;  

% META.calculate_effort = 2;   %Not sure what this is ??

%effort quota, or area cleaned by boat per 6 months; 11.52km2 per 6
%months * proportion effort available
%META.boatProperties.effortQuota = zeros(1, length(META.boatProperties.boatID));
%META.boatProperties.effortQuota(i) = area_day*boat_days*META.control_effort_allocation;

%average distance of 20 min swim in m; calculated from timed swims, assume AMPTO moves at same speed, although culls probably slower
%swimd=480;
%average width covered during swim in m; manual reach; no changes due to habitat complexity etc, no slowdown due to high densities
%swimw=3;
%average area covered by a diver per hour in m2; 4320m2, or 66x66m
%swima=swimd*swimw*3;
%number of hours dived per dive; AMPTO
%divet=2/3;%was 2/3 originally, maybe still is?
%area covered by boat per day; 115200m2, or 340x340m, or 0.1152km2
%area_day=swima*META.divers*divet*diven;
%META.boatProperties.areaPerDay(i) = area_day;
