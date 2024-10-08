%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Populate the reef grid randomly with coral colonies of random sizes
% Number of colonies within a cell for each species is random and bounded by META.max_colonies
% Now integrates clade C/D, which replaces the bleaching history of each colony

% 09/2022: now calling REEFn. (== REEF(n).) the structure array of reef ID = n for this function and
% f_derivecoralcover to avoid confusion

function [coral, algal] = f_initialise_population(META, REEFn, CORAL)    

% Dimensions of the grid
m = META.grid_x_count ;
n = META.grid_y_count ;

% Initialise metrics of coral colonies
coral(META.nb_coral_types).cover_cm2 = sparse(zeros(m*n, 1)) ; % planar area in cm2 (~size)

% Initialise the cover of each algal type
algal(META.nb_algal_types).cover_cm2 = sparse(zeros(m*n, 1)) ;

% Initialise temporary variables
coral_cm2 = zeros(m*n, META.nb_coral_types,META.max_colonies) ; % area (in cm2) of all colonies
coral_ID = coral_cm2 ; % ID number of all colonies
coral_HT = coral_cm2 ; % heat tolerance effect size
algal_cm2 = zeros(m*n, META.nb_algal_types) ;

% List of ID of grazable cells (exclude sand in loops)
list_cell = uint16(1:(m*n)) ;
list_cell = list_cell(REEFn.grazable_cell==1) ;

%%%%% ADD CORALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If initial coral cover is non-null, lets fill the grid with coral colonies
if sum([REEFn.initial_coral_cover])~=0
    
    all_colony_sizes=[];
    all_colony_types=[];
    all_colony_ID=[];
    all_colony_HT=[];
    
    for s = 1:META.nb_coral_types
        
        if REEFn.initial_coral_cover(s)>0
            
            % PULL IN FREQUENCY FOR EACH SIZE CLASS
            colony_size = f_derivecoralcover(s, REEFn, CORAL); % generate colonies of different sizes (cm2)
            
            % STORE THE CORRESPONDING CORAL TYPE
            colony_type = s*ones(size(colony_size));
            
            % CREATE an ID for each colony
            colony_ID = 1:1:length(colony_type);
            
            % ASSIGN heat tolerance at random within specified range
            colony_HT = normrnd(0,sqrt(CORAL.VAR_HT(s)), 1, length(colony_type));
            colony_HT(abs(colony_HT)>=CORAL.MAX_HT(s))=0 ; % replace extreme values (unlikely) by the mean

            % CONCATENATE
            all_colony_sizes = [all_colony_sizes colony_size]; % list of sizes of all colonies, all species combined
            all_colony_types = [all_colony_types colony_type]; % corresponding coral type
            all_colony_ID = [all_colony_ID colony_ID]; %corresponding ID (to be combined later with colony type for full identification)
            all_colony_HT = [all_colony_HT colony_HT];
        end
    end

    [all_colony_sizes,rank] = sort(all_colony_sizes,'descend'); % sort colony size by decreasing value (largest colonies get settled first)
    all_colony_types = all_colony_types(rank); % sort colony type accordingly
    all_colony_ID = all_colony_ID(rank); % sort colony type accordingly
    all_colony_HT = all_colony_HT(rank); % sort heat tolerance effect sizes accordingly
    
    % THEN FILL THE GRAZABLE CELLS WITH CORAL COLONIES (largest get settled first)  
    i = 0; % iterator for visiting the cells
    r = randperm(length(list_cell)) ;
    list_cell2 = list_cell(r) ;
    
    list_cell2 = repmat(list_cell2, 1,50);

    current_total_cover = zeros(m*n,1); % initialise the counter of total coral cover per cell
    current_total_colonies = zeros(m*n,META.nb_coral_types); % counter of number of colonies per cell
    
    
    for c = 1:length(all_colony_sizes) % pick up every colony, one by one
    
        accept = 1;  % decide whether the colony has been processed or not

        while accept == 1           

            i = i+1; % iterator of cell visiting -> goes to the next cell until accept = 2   
%             if i>length(list_cell2)
%                 i
%                 REEF.initial_coral_cover
%                 REEF.nongrazable_substratum
%             end
            cell = list_cell2(i); 
            colony_count = current_total_colonies(cell,all_colony_types(c))+1;
                           
            if (colony_count <= META.max_colonies) ...
                    && ((all_colony_sizes(c)+current_total_cover(cell,1))<=REEFn.substrate_SA_cm2(cell))

                coral_cm2(cell, all_colony_types(c), colony_count) = all_colony_sizes(c) ; % allocate it to the cell
                coral_ID(cell, all_colony_types(c), colony_count) = all_colony_ID(c) ;
                coral_HT(cell, all_colony_types(c), colony_count) = all_colony_HT(c) ;
                
                current_total_cover(cell,1) = current_total_cover(cell,1) + all_colony_sizes(c) ; %update total coral cover in that cell
                current_total_colonies(cell,all_colony_types(c)) = colony_count;
                accept = 2; % process the next colony
   
            end
        end
    end
end

%%%%% ADD MACROALGAE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total amount of corals to restrict allocation of macroalgae
total_coral_cm2 = round(sum(sum(coral_cm2,3),2));

if sum(REEFn.initial_algal_cover(2:3))>0
    
    % 25/02/2017: from now on, macroalgae are generated by patches of 200cm2
    % this contrasts with previous versions whereby all cells could only contain one type of algae
    % and where macroaglae was present, together with the coral cover the cell was 100% covered
    total_macroalgal_cm2 = REEFn.initial_algal_cover*sum(REEFn.substrate_SA_cm2);
    
    Algal_patch_cm2 = 200*ones(1,10000);
    sum_Algal_patch = cumsum(Algal_patch_cm2) ;
    
    I2 = find(sum_Algal_patch < (total_macroalgal_cm2(2)-mean(Algal_patch_cm2)));
    patch(2).selection = Algal_patch_cm2(1:(max(I2)+1));
    
    I3 = find(sum_Algal_patch < (total_macroalgal_cm2(3)-mean(Algal_patch_cm2)));
    patch(3).selection = Algal_patch_cm2(1:(max(I3)+1));
    
    % The last patch matches the desired total cover
    patch(2).selection(end)=patch(2).selection(end)+total_macroalgal_cm2(2)-sum(patch(2).selection);
    patch(3).selection(end)=patch(3).selection(end)+total_macroalgal_cm2(3)-sum(patch(3).selection);
    
    i = 0 ;
    r = randperm(length(list_cell)) ;
    list_cell2 = list_cell(r) ;
    list_cell2 = [list_cell2 list_cell2 list_cell2 list_cell2 list_cell2 list_cell2];
    list_cell2 = [list_cell2 list_cell2 list_cell2 list_cell2 list_cell2 list_cell2];
    
    for a = 2:3 % process Dictyota first
        
        while total_macroalgal_cm2(a) > 0
            
            for p = 1:length(patch(a).selection) % pick up every patch, one by one
                
                accept = 1;  % decide whether the patch has been processed or not
                
                while accept == 1
                    
                    i = i+1; % iterator of cell visiting -> goes to the next cell until accept = 2
                    cell = list_cell2(i);
                    
                    if (patch(a).selection(p)+total_coral_cm2(cell)+sum(algal_cm2(cell,2:3))<=REEFn.substrate_SA_cm2(cell))
                        
                        algal_cm2(cell, a) = algal_cm2(cell, a) + patch(a).selection(p) ; % allocate it to the cell (and add to previous lob if revisiting)
                        total_macroalgal_cm2(a) = total_macroalgal_cm2(a) - patch(a).selection(p) ; % update amount of Lob still to dispatch
                        accept = 2; % process the next colony
                        
                    end
                end
            end
        end
    end
end

% Where lob and dict colonize, the rest is EAM
algal_cm2(:,1)= REEFn.substrate_SA_cm2 - total_coral_cm2 - sum(algal_cm2(:,2:3),2);
% algal_cm2(sum(algal_cm2(:,2:3),2)==0,4)=0; 
algal_cm2(REEFn.grazable_cell==0,1) = 0 ; % no EAM on sand

%%%%% FINAL ALLOCATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:META.nb_coral_types
    
    temp_coral_cm2 = squeeze(coral_cm2(:,s,:)) ;
    id_col = spones(temp_coral_cm2);
    id_sum=sum(id_col,1);
    
    coral(s).cover_cm2 = temp_coral_cm2(:,id_sum~=0) ;
    coral(s).colony_ID = squeeze(coral_ID(:,s,id_sum~=0)) ;
    coral(s).heat_tolerance = squeeze(coral_HT(:,s,id_sum~=0)) ;
    
    if META.doing_clades == 1
        % Assign a clade (affects mortality to thermal stress and growth rate)
        coral(s).clade = spones(coral(s).cover_cm2); %default is 1 (thermally sensitive, clade C)
        rand_clade = sprand(coral(s).cover_cm2) ;
        coral(s).clade(rand_clade > CORAL.clade_prop) = 2 ; % alternative is 2 (thermally tolerant, clade D)
    else
        coral(s).clade = 0;
    end
    
    if META.doing_3D==1    %Initialise 3D variables        
        coral(s).surface_cm2 = zeros(size(coral(s).cover_cm2)) ;
        coral(s).volume_cm3 = coral(s).surface_cm2 ;        
    else
        coral(s).surface_cm2 = 0 ;
        coral(s).volume_cm3 = 0 ;        
    end
end

for a = 2:META.nb_algal_types    
    tmp_a = zeros(m*n,1) ; 
    tmp_a(:) = round(algal_cm2(:,a)) ;
    algal(a).cover_cm2 = sparse(tmp_a) ;   
end

% Fill the rest of each cell with turf
algal(1).cover_cm2 = REEFn.substrate_SA_cm2 ...
	- sum([algal(2:META. nb_algal_types).cover_cm2],2) - sum([coral.cover_cm2],2) ;
algal(1).cover_cm2(REEFn.grazable_cell==0) = 0 ;
algal(1).cover_cm2 = sparse(algal(1).cover_cm2);

% Final check
total = sum(cat(2,algal.cover_cm2),2) + sum(cat(2,coral.cover_cm2),2);
total(REEFn.grazable_cell==0) = REEFn.substrate_SA_cm2(REEFn.grazable_cell==0);

if sum(total) ~= sum(REEFn.substrate_SA_cm2)
    error('inconsistent filling across the grid')
end