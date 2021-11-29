% =================================================================
% This program returns start points for the analysis, using only start
%   points that give turnout between 40-80% ("reasonable start points")
% Inputs: number of start points to return, L-setting (cons or dist),
%   number of simulated voters, random set, and name of the votesim file
% Ouputs: matrix of start points, range from which start points were
%   selected
% Same parameters for voters and non-voters
% =================================================================

function [searchInits,lower_start,upper_start] = getSearchInits(noSearchInits,...
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

% draw X times as many and then only keep first noSearchInits with turnout
% in the right range
X=4;
lb=repmat(lower_start,1,noSearchInits*X);
ub=repmat(upper_start,1,noSearchInits*X);
searchInits_all=lb+(ub-lb).*rand(noOfParams,noSearchInits*X);

% capture turnout at each of the start points (moment 101)
turnout=repmat(0,1,noSearchInits*X);
for i=1:noSearchInits*X
    temp = votesimFile(searchInits_all(:,i),rand_set);
    turnout(i) = temp(101);
end 

%keep only observations with turnout in the range
searchInits_valid = searchInits_all(:,turnout>0.4 & turnout<0.8);
searchInits = searchInits_valid(:,1:noSearchInits);
end