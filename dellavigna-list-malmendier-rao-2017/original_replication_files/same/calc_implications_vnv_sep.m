% =================================================================
% This file replicates the beginning of the relevant voteSim file and then 
% outputs list of model implications for the tables
% Running voters and non-voters separately - exogenous voter status
% =================================================================

function [implications]= calc_implications_vnv_sep(parameters,sim_voters,L_setting,rand_set)

% =================================================================
% Set Parameters
% ================================================================

% =========================================
% =========================================
simSize=sim_voters; 
N=5.4; % Number of times asked in Congressional
N_P=10.1; % times asked in presidential

h0	=	parameters(1);
r	=	parameters(2);
eta	=	parameters(3);
mu_s	=	parameters(4);
sigma_s	=	parameters(5);
S_svy	=	parameters(6);
timeval	=	parameters(7);

mu_sv	=	parameters(8);
mu_sn	=	parameters(9);
sigma_sv	=	parameters(10);
sigma_sn	=	parameters(11);
rho	=	parameters(12);

if L_setting == 'cons' 
    L=parameters(13);
end

if L_setting == 'dist' 
    L_lambda=parameters(13);
end

mu_eps = parameters(14);
sigma_eps = parameters(15);

% Values of survey incentives (relative to 0d10m)
% XdYm = X dollars and Y min
D_0d10m = 0;
D_0d5m = timeval*5/60;
D_10d5m = 10+timeval*5/60;
D_10d10m = 10;

% Extra incentive if say "not vote"
%  5m survey: +1m + $5 
%  10m survey: -8m
I_5d1m = 5-timeval*1/60;
I_8m = timeval*8/60;

% ==========================================================================
% Draw from random variables (distributions of s, sv, sn, eps)
% transform already set random components
% =======================================================================

% order of rand_set: s eps sv sn

%%% Simulate utility of doing a 0d10m survey 
s = mu_s + sigma_s.*rand_set(:,1);

%%% Simulate epsilons ("other" reasons to vote)
eps = mu_eps + sigma_eps.*rand_set(:,2);

%%% Simulate sv and sn
sv = mu_sv + sigma_sv.*rand_set(:,3);
sn = mu_sn + sigma_sn.*rand_set(:,4);

if L_setting == 'dist' 
    %%% Simulate lying cost
    % This is an expoential distribution 
    L = -log(rand_set(:,5)) / L_lambda;
end

% ==========================================================================
% Vote or Not
% =======================================================================

% Net expected utility of voting (relative to not voting)
sigVal = (max(sv,sn-L) - max(sn,sv-L)); % net sig utility voter - non-voter (1 ask)
sigVal_x_N = N.*sigVal; % asked N times
utilityVoting = sigVal_x_N + eps; % net utility
voted = utilityVoting>0; % Vote if net utility > 0

% =================================================================
% Lying if Asked Section - utility from being asked once
% ================================================================

% utilVotingQuestion= utility you get from being asked one 
%     time (max of lie or not lie)
wouldLieIfAsked = voted.*(sn-L>sv) + (1-voted).*(sv-L>sn);
utilVotingQuestion = voted.*max(sn-L,sv) + (1-voted).*max(sv-L,sn);

% ========================================================================
% Output vector: 
% signal value * N with L= 0, 2, 5, 10
% Utility from being asked about voting

% Mean Value of saying voted (sV for voters, sV-L for nonvoters)
% Mean Value of saying didn't vote (sN-L for voters, sN for nonvoters)
value_say_vote = voted.*(sv) + (1-voted).*(sv-L);
value_say_notvote = voted.*(sn-L) + (1-voted).*(sn);


implications = zeros(11,1); 
    % Mean Value of saying voted
    implications(1) = mean(value_say_vote);
    % Mean Value of saying didn't vote
    implications(2) = mean(value_say_notvote);
    
    % signal value * N by L
    implications(3) = mean(N.*(max(sv,sn-0) - max(sn,sv-0)));
    implications(4) = mean(N.*(max(sv,sn-2) - max(sn,sv-2)));
    implications(5) = mean(N.*(max(sv,sn-5) - max(sn,sv-5)));
    implications(6) = mean(N.*(max(sv,sn-10) - max(sn,sv-10)));
    % utility asked once
    implications(7) = mean([utilVotingQuestion]);
    implications(8) = std([utilVotingQuestion]);
    implications(9) = prctile([utilVotingQuestion],25);
    implications(10) = prctile([utilVotingQuestion],50);
    implications(11) = prctile([utilVotingQuestion],75);

end

