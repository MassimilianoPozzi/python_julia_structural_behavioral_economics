% =================================================================
% This program returns start points for the analysis, using only start
%   points that give turnout between 50-70% ("reasonable start points")
% Allow non-selection param to vary across voter and non-voter
% Sigma Epsilon is fixed
% =================================================================

function [searchInits,lower_start,upper_start] = getSearchInits_fix_sigep(fixed_sigma_eps,...
    noSearchInits,L_setting,sim_voters,rand_set,votesimFile)
noOfParams=22;

if L_setting == 'cons' 
    lower_start = [0.2, 0.2, 0.2, 0.2, eps, eps, -50, -50, eps, eps, eps, eps, eps, eps, -20, -30, eps, eps, -1, eps, -30, 50]';
    upper_start = [0.4, 0.4, 0.4, 0.4, 0.5, 0.5, eps, eps,  50,  50,  10,  10, 100, 100,  20,  10,  30,  30,  1,  20,  100, 200]';
end

if L_setting == 'dist' 
    lower_start = [0.2, 0.2, 0.2, 0.2, eps, eps, -50, -50, eps, eps, eps, eps, eps, eps, -20, -30, eps, eps, -1, eps, -30, 50]';
    upper_start = [0.4, 0.4, 0.4, 0.4, 0.5, 0.5, eps, eps,  50,  50,  10,  10, 100, 100,  20,  10,  30,  30,  1,  1,  100, 200]';
end

% draw X times as many and then only keep first noSearchInits with turnout
% in the right range
X=10;
lb=repmat(lower_start,1,noSearchInits*X);
ub=repmat(upper_start,1,noSearchInits*X);
searchInits_all=lb+(ub-lb).*rand(noOfParams,noSearchInits*X);

% capture turnout at each of the start points (moment 101)
turnout=repmat(0,1,noSearchInits*X);
for i=1:noSearchInits*X
    %replace fixed sigma epsilon as the last parameter
    searchInits_all(noOfParams,i) = fixed_sigma_eps;
    %record turnout
    temp = votesimFile(searchInits_all(:,i),rand_set);
    turnout(i) = temp(101);
end 

%keep only observations with turnout in the range
searchInits_valid = searchInits_all(:,turnout>0.5 & turnout<0.7);
searchInits = searchInits_valid(:,1:noSearchInits);

end

