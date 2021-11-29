% =================================================================
% This program returns start points for the analysis
% don't restrict turnout since we have 100% or 0%
% =================================================================


function [searchInits,lower_start,upper_start] = getSearchInits_vnv_sep(noSearchInits,...
    L_setting,sim_voters,rand_set,votesimFile)
noOfParams=15;

if L_setting == 'cons'
    lower_start = [0.2, 0.2, eps, -50, eps, eps, eps, -20, -30, eps, eps, -1, eps, -30, 50]';
    upper_start = [0.4, 0.4, 0.5, eps, 50,  10,  100,  20,  10,  30,  30,  1, 20,  100, 200]';
end

if L_setting == 'dist' 
    lower_start = [0.2, 0.2, eps, -50, eps, eps, eps, -20, -30, eps, eps, -1, eps, -30, 50]';
    upper_start = [0.4, 0.4, 0.5, eps, 50,  10,  100,  20,  10,  30,  30,  1, 1,   100, 200]';
end

lb=repmat(lower_start,1,noSearchInits);
ub=repmat(upper_start,1,noSearchInits);
searchInits=lb+(ub-lb).*rand(noOfParams,noSearchInits);
end