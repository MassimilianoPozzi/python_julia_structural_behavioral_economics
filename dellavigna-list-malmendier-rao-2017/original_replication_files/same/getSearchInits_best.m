% =================================================================
% Updated by Alex Steiny 9/8/2016
% This program returns start points for the analysis
% Use top 20 estimates from estimation results below
% =================================================================

function [searchInits,lower_start,upper_start] = getSearchInits_best(noSearchInits,...
    search_init_file,sim_voters,rand_set, votesimFile)

noOfParams=15;

if search_init_file == 'fixsige490'
lower_start = [0.2, 0.2, eps, -50, eps, eps, eps, -20, -30, eps, eps, -1, eps, -30, 50]';
upper_start = [0.4, 0.4, 0.5, eps, 50,  10,  100,  20,  10,  30,  30,  1, 20,  100, 200]';
load ('fixsige490.mat'); 
end

if search_init_file == 'bench_best'
lower_start = [0.2, 0.2, eps, -50, eps, eps, eps, -20, -30, eps, eps, -1, eps, -30, 50]';
upper_start = [0.4, 0.4, 0.5, eps, 50,  10,  100,  20,  10,  30,  30,  1, 20,  100, 200]';
load ('bench_best.mat'); 
end

%take best estimates as new start points
searchInits_new = convergedParamHats(1:noSearchInits,:);
searchInits = searchInits_new';
end