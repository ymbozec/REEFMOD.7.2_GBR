%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [fecundity_adol, fecundity_adult] = f_estimate_fecundity (coral_cm2, CORAL)
% 
% % Calculate fecundity output for a reef in terms of total number of larvae produced over 
% % the grid for a coral species (will be transformed to larval input for other reefs using
% % the connectivity matrix).
% id1 = zeros(size(coral_cm2)) ;
% id2 = id1 ;
% 
% % Calculate fecundity for adults
% id1(coral_cm2 >= CORAL.adult_size) = 1 ;
% fecundity_adult = sum(sum(2*pi*((sqrt(coral_cm2.*id1/pi)).^2)*216)) ;
% 
% % calculate fecundity for adol
% id2(coral_cm2 >= CORAL.adol_size) = 1 ;
% id2 = id2 - id1;
% fecundity_adol = sum(sum(0.25*2*pi*((sqrt(coral_cm2.*id2/pi)).^2)*216)) ;

function [total_fecundity, select_HT] = f_estimate_fecundity (coral_cm2, F_list, fecund_min_size, a, b, HT_pop_size)

% Select colony sizes with 100% gravid
I=find(coral_cm2>= fecund_min_size);

if isempty(I)==0
    adult_coral_sizes = coral_cm2(I);
    adult_relfitness = F_list(I);
    
    % Allometric relationship based on Hall and Hughes 1996:
    % Egg volume = exp(a+b*log(size))
    
    % all_egg_volumes = exp(a + b*log(adult_relfitness.*adult_coral_sizes)) ; % mm3 of eggs produced by each colony
    % Correction 29/09/2022: logs in Hall&Hughes were log10, not natural logs!!
    % Also, makes more sense to apply relative fitness on number of eggs rather than colony size (log...)
    all_egg_volumes = adult_relfitness.*(10.^(a + b*log10(adult_coral_sizes))) ; % mm3 of eggs produced by each colony
    
    total_fecundity = floor(sum(sum(all_egg_volumes))/0.1) ; %0.1 mm3 is the average volume of an egg
    
    
%% 10/2023: generate N representative HT phenotypes of the parents (using individual fecundity as weight)
    if length(I)==1 % if only one colony spawning then only select this one for representative HT
        % (necessary, otherwise randsample is messed up below)
        select_HT = I*ones(HT_pop_size,1);
        
    else %otherwise sample the ID among the represented ones using fecundity as weight
        
        select_HT = randsample(I, HT_pop_size, 'true', all_egg_volumes/sum(all_egg_volumes));
    end

else % if no fecund corals, then no larvae, and no HT
    total_fecundity = 0;
    select_HT = [];
end