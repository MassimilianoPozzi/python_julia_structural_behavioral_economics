% =================================================================
% Takes as inputs the parameter values, simulation N, type of lying cost, 
%   and set of random components. Outputs the simulation moments
% Same parameters for voters and non-voters
% Measurement Error - 10% of apparent voters are true non-voters
% no measurement error in turnout
% =================================================================

function [simMoments] = voteSimEndog_meas_v(parameters,rand_set)

% =================================================================
% Set Parameters
% ================================================================

% Set simulation paramaters

simSize=length(rand_set); 
N=5.4; % Number of times asked in Congressional
N_P=10.1; % times asked in presidential


% Parameters:
% h0 = baseline probability of being at home
% r = prob of seeing flyer
% eta = elasticity of response to sorting in and out of the house
% s= how much you like doing a survey
% S= social pressure cost of doing a survey
% sv= value of saying you voted (appearing as a voter)
% sn= value of saying you didn't vote (appearing as a non-voter)
% rho= correlation between sv and sn
% we will typically set sigma_sv = sigma_sn and rho = 0
% eps= all other reasons to vote

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

L=parameters(13);

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


% ==========================================================================
% Vote or Not
% =======================================================================

% Net expected utility of voting (relative to not voting)
sigVal = (max(sv,sn-L) - max(sn,sv-L)); % net sig utility voter - non-voter (1 ask)
sigVal_x_N = N.*sigVal; % asked N times
utilityVoting = sigVal_x_N + eps; % net utility
voted = utilityVoting>0; % Vote if net utility > 0

% return a list of indices (e.g., 1, 3, 4) identifying voters and non-voters
voterIndex = find(voted==1);
nonvoterIndex = find(voted==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measurement Error
% 10% of apparent voters are non-voters
% 0% of apparent non-voters are voters

appear_voted = voted;

% use random component in column 6 to determine who is mismeasured
prop_mismeas = rand_set(:,6);
prop_mismeas_v = prop_mismeas(voterIndex);
prop_mismeas_nv = prop_mismeas(nonvoterIndex);

% determine final number of apparent voters and non-voters
N_apparent_v = round(length(voterIndex) / 0.9);
N_apparent_nv = simSize - N_apparent_v;

% determine how many of each we will need to change
N_mismeas_as_v = round(0.1 * N_apparent_v);
%N_mismeas_as_nv = round(0.1 * N_apparent_nv);

% need final number of apparent voters and non-voters and number of each
% that are reassigned to be > 0
if N_apparent_v>0 && N_apparent_nv > 0 && N_mismeas_as_v > 0 % && N_mismeas_as_nv > 0
    % Identify true non-voters that will appear as voters 
    % we will change N_mismeas_v of them
    mismeas_temp = flipud(unique(sort(prop_mismeas_nv))); %prop for nv to be mismeasured
    mismeas_topval=mismeas_temp(N_mismeas_as_v); % cutoff value
    nv_appear_v = nonvoterIndex(find(prop_mismeas_nv>=mismeas_topval));

    % Identify true voters that will appear as non-voters 
    % we will change N_mismeas_nv of them
    %mismeas_temp = flipud(unique(sort(prop_mismeas_v))); %prop for v to be mismeasured
    %mismeas_topval=mismeas_temp(N_mismeas_as_nv); % cutoff value
    %v_appear_nv = voterIndex(find(prop_mismeas_v>=mismeas_topval));
   
    %reassign true non-voters who appear as voters
    for i=1:length(nv_appear_v)
      appear_voted(nv_appear_v(i)) = 1;
    end 
    
    %reassign true voters who appear as non-voters
    %for i=1:length(v_appear_nv)
    %  appear_voted(v_appear_nv(i)) = 0;
    %end
end

appear_voterIndex = find(appear_voted==1);
appear_nonvoterIndex = find(appear_voted==0);
mismeasIndex = find(appear_voted~=voted);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
% turnout calculated based on true turnout
%        (no measurement error in turnout)
% share who turn out in the control group (no GOTV intervention)
Turnout_control = mean(voted);  

% GOTV
% Assumes intervention = expect to be asked h0 times more

% net utility of voting if seen the GOTV flyer
sigVal_x_N_GOTV = (N+h0).*sigVal;
utilityVoting_GOTV = sigVal_x_N_GOTV + eps; 
voted_GOTV = utilityVoting_GOTV>0; 

% share who turn out with the GOTV intervention
% assume everyone sees the flyer, counts as N+h0 times asked
Turnout_GOTV = mean(voted_GOTV);

%PRESIDENTIAL
sigVal_x_N_P = N_P.*sigVal; % asked N_P times
sigVal_x_N_P_GOTV = (N_P+h0).*sigVal; % asked N_P times

utilityVoting_P = sigVal_x_N_P + eps;
utilityVoting_P_GOTV = sigVal_x_N_P_GOTV + eps; 

voted_P = utilityVoting_P>0; 
voted_P_GOTV = utilityVoting_P_GOTV>0; 

% turnout 
Turnout_P_control = mean(voted_P); 
Turnout_P_GOTV = mean(voted_P_GOTV);


% =================================================================
% Lying if Asked - with measurement error, opposite if you appear as opposite
% ================================================================

% utilVotingQuestion= utility you get from being asked one 
%     time (max of lie or not lie)
true_wouldLieIfAsked = voted.*(sn-L>sv) + (1-voted).*(sv-L>sn);
utilVotingQuestion = voted.*max(sn-L,sv) + (1-voted).*max(sv-L,sn);

% Response to incentives to say "not vote"
true_wouldLieIfAsked_5d1m = voted.*(sn-L+I_5d1m>sv) + (1-voted).*(sv-L>sn+I_5d1m);
true_wouldLieIfAsked_8m = voted.*(sn-L+I_8m>sv) + (1-voted).*(sv-L>sn+I_8m);

%Reverse if you are mismeasured
wouldLieIfAsked = true_wouldLieIfAsked;
wouldLieIfAsked_5d1m = true_wouldLieIfAsked_5d1m;
wouldLieIfAsked_8m = true_wouldLieIfAsked_8m;
for i=1:length(mismeasIndex)
  wouldLieIfAsked(mismeasIndex(i)) = 1-true_wouldLieIfAsked(mismeasIndex(i)); 
  wouldLieIfAsked_5d1m(mismeasIndex(i)) = 1-true_wouldLieIfAsked_5d1m(mismeasIndex(i)); 
  wouldLieIfAsked_8m(mismeasIndex(i)) = 1-true_wouldLieIfAsked_8m(mismeasIndex(i));       
end 


% ==========================================================================
% ==========================================================================
% Calculate Moments
% =======================================================================
% ==========================================================================

% =================================================================
% Utility from Doing Survey - use true voter index
% ================================================================

% NF= no flyer, F = flyer, FV= flyer + voting, 
% OO = opt-out, OOV = opt-out + voting

% anticipated utility from doing survey and voting survey (VF or I)
util_svyOnly_0d5m = s+D_0d5m;
util_svyPlusVotingQues_0d5m = util_svyOnly_0d5m + utilVotingQuestion;

util_svyOnly_10d10m = s+D_10d10m;
util_svyPlusVotingQues_10d10m = util_svyOnly_10d10m + utilVotingQuestion;

util_svyOnly_10d5m = s+D_10d5m;
util_svyPlusVotingQues_10d5m = util_svyOnly_10d5m + utilVotingQuestion;


% If asked, do survey if greater than the social pressure cost
% NI = not informed that survey is about voting, I=informed (VF or I)
doesSvyIfAsked_NI_0d5m = util_svyOnly_0d5m > -S_svy;
doesSvyIfAsked_I_0d5m = util_svyPlusVotingQues_0d5m > -S_svy;

doesSvyIfAsked_NI_10d10m = util_svyOnly_10d10m > -S_svy;
doesSvyIfAsked_I_10d10m = util_svyPlusVotingQues_10d10m > -S_svy;

doesSvyIfAsked_NI_10d5m = util_svyOnly_10d5m > -S_svy;
doesSvyIfAsked_I_10d5m = util_svyPlusVotingQues_10d5m > -S_svy;


% anticipated utility given you are asked to do the survey
anticipatedUtil_Svy_NI_0d5m = max(util_svyOnly_0d5m,-S_svy); 
anticipatedUtil_Svy_I_0d5m = max(util_svyPlusVotingQues_0d5m,-S_svy);

anticipatedUtil_Svy_NI_10d10m = max(util_svyOnly_10d10m,-S_svy); 
anticipatedUtil_Svy_I_10d10m = max(util_svyPlusVotingQues_10d10m,-S_svy);

anticipatedUtil_Svy_NI_10d5m = max(util_svyOnly_10d5m,-S_svy); 
anticipatedUtil_Svy_I_10d5m = max(util_svyPlusVotingQues_10d5m,-S_svy);


% opt-out if anticipated utility is negative
optsOutIfSees_OO_0d5m = anticipatedUtil_Svy_NI_0d5m < 0;
optsOutIfSees_OOV_0d5m= anticipatedUtil_Svy_I_0d5m < 0;

optsOutIfSees_OO_10d10m = anticipatedUtil_Svy_NI_10d10m < 0;
optsOutIfSees_OOV_10d10m= anticipatedUtil_Svy_I_10d10m < 0;

optsOutIfSees_OO_10d5m = anticipatedUtil_Svy_NI_10d5m < 0;
optsOutIfSees_OOV_10d5m= anticipatedUtil_Svy_I_10d5m < 0;


% choosing probability of being at home is bounded between 0 and 1
hStar_F_0d5m = max(0,min(1,h0+eta.*anticipatedUtil_Svy_NI_0d5m));
hStar_FV_0d5m = max(0,min(1,h0+eta.*anticipatedUtil_Svy_I_0d5m));

hStar_F_10d10m = max(0,min(1,h0+eta.*anticipatedUtil_Svy_NI_10d10m));
hStar_FV_10d10m = max(0,min(1,h0+eta.*anticipatedUtil_Svy_I_10d10m));

hStar_F_10d5m = max(0,min(1,h0+eta.*anticipatedUtil_Svy_NI_10d5m));
hStar_FV_10d5m = max(0,min(1,h0+eta.*anticipatedUtil_Svy_I_10d5m));


% =============================================================
% Separate Voters and Nonvoters - use apparent voter index
% =============================================================

% split into separate vectors of voters and non-voter vectors 
% note: they will be of different length

%voters
hStar_F_0d5m_v = hStar_F_0d5m(appear_voterIndex);
hStar_FV_0d5m_v = hStar_FV_0d5m(appear_voterIndex);
doesSvyIfAsked_NI_0d5m_v = doesSvyIfAsked_NI_0d5m(appear_voterIndex);
doesSvyIfAsked_I_0d5m_v = doesSvyIfAsked_I_0d5m(appear_voterIndex);
optsOutIfSees_OO_0d5m_v=optsOutIfSees_OO_0d5m(appear_voterIndex);
optsOutIfSees_OOV_0d5m_v=optsOutIfSees_OOV_0d5m(appear_voterIndex);

hStar_F_10d10m_v = hStar_F_10d10m(appear_voterIndex);
hStar_FV_10d10m_v = hStar_FV_10d10m(appear_voterIndex);
doesSvyIfAsked_NI_10d10m_v = doesSvyIfAsked_NI_10d10m(appear_voterIndex);
doesSvyIfAsked_I_10d10m_v = doesSvyIfAsked_I_10d10m(appear_voterIndex);
optsOutIfSees_OO_10d10m_v=optsOutIfSees_OO_10d10m(appear_voterIndex);
optsOutIfSees_OOV_10d10m_v=optsOutIfSees_OOV_10d10m(appear_voterIndex);

hStar_F_10d5m_v = hStar_F_10d5m(appear_voterIndex);
hStar_FV_10d5m_v = hStar_FV_10d5m(appear_voterIndex);
doesSvyIfAsked_NI_10d5m_v = doesSvyIfAsked_NI_10d5m(appear_voterIndex);
doesSvyIfAsked_I_10d5m_v = doesSvyIfAsked_I_10d5m(appear_voterIndex);
optsOutIfSees_OO_10d5m_v=optsOutIfSees_OO_10d5m(appear_voterIndex);
optsOutIfSees_OOV_10d5m_v=optsOutIfSees_OOV_10d5m(appear_voterIndex);

wouldLieIfAsked_v = wouldLieIfAsked(appear_voterIndex);
wouldLieIfAsked_5d1m_v = wouldLieIfAsked_5d1m(appear_voterIndex);
wouldLieIfAsked_8m_v=wouldLieIfAsked_8m(appear_voterIndex);

% non-voters
hStar_F_0d5m_nv = hStar_F_0d5m(appear_nonvoterIndex);
hStar_FV_0d5m_nv = hStar_FV_0d5m(appear_nonvoterIndex);
doesSvyIfAsked_NI_0d5m_nv = doesSvyIfAsked_NI_0d5m(appear_nonvoterIndex);
doesSvyIfAsked_I_0d5m_nv = doesSvyIfAsked_I_0d5m(appear_nonvoterIndex);
optsOutIfSees_OO_0d5m_nv=optsOutIfSees_OO_0d5m(appear_nonvoterIndex);
optsOutIfSees_OOV_0d5m_nv=optsOutIfSees_OOV_0d5m(appear_nonvoterIndex);

hStar_F_10d10m_nv = hStar_F_10d10m(appear_nonvoterIndex);
hStar_FV_10d10m_nv = hStar_FV_10d10m(appear_nonvoterIndex);
doesSvyIfAsked_NI_10d10m_nv = doesSvyIfAsked_NI_10d10m(appear_nonvoterIndex);
doesSvyIfAsked_I_10d10m_nv = doesSvyIfAsked_I_10d10m(appear_nonvoterIndex);
optsOutIfSees_OO_10d10m_nv=optsOutIfSees_OO_10d10m(appear_nonvoterIndex);
optsOutIfSees_OOV_10d10m_nv=optsOutIfSees_OOV_10d10m(appear_nonvoterIndex);

hStar_F_10d5m_nv = hStar_F_10d5m(appear_nonvoterIndex);
hStar_FV_10d5m_nv = hStar_FV_10d5m(appear_nonvoterIndex);
doesSvyIfAsked_NI_10d5m_nv = doesSvyIfAsked_NI_10d5m(appear_nonvoterIndex);
doesSvyIfAsked_I_10d5m_nv = doesSvyIfAsked_I_10d5m(appear_nonvoterIndex);
optsOutIfSees_OO_10d5m_nv=optsOutIfSees_OO_10d5m(appear_nonvoterIndex);
optsOutIfSees_OOV_10d5m_nv=optsOutIfSees_OOV_10d5m(appear_nonvoterIndex);

wouldLieIfAsked_nv = wouldLieIfAsked(appear_nonvoterIndex);
wouldLieIfAsked_5d1m_nv = wouldLieIfAsked_5d1m(appear_nonvoterIndex);
wouldLieIfAsked_8m_nv=wouldLieIfAsked_8m(appear_nonvoterIndex);

% =============================================================
% Disaggregated Moments - this section doesn't change with same parameters
% =============================================================

% voters
% ===========================

%%% PH = probability of being at home
PH_NF_0d5m_v = h0;
PH_F_0d5m_v = (1-r)*h0 + r.*mean(hStar_F_0d5m_v);
PH_FV_0d5m_v = (1-r)*h0 + r.*mean(hStar_FV_0d5m_v);
PH_OO_0d5m_v = (1-r)*h0 + r.*mean((1-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v); 
PH_OOV_0d5m_v = (1-r)*h0 + r.*mean((1-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v);

PH_NF_10d10m_v = h0;
PH_F_10d10m_v = (1-r)*h0 + r.*mean(hStar_F_10d10m_v);
PH_FV_10d10m_v = (1-r)*h0 + r.*mean(hStar_FV_10d10m_v);
PH_OO_10d10m_v = (1-r)*h0 + r.*mean((1-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v); 
PH_OOV_10d10m_v = (1-r)*h0 + r.*mean((1-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v);

PH_NF_10d5m_v = h0;
PH_F_10d5m_v = (1-r)*h0 + r.*mean(hStar_F_10d5m_v);
PH_FV_10d5m_v = (1-r)*h0 + r.*mean(hStar_FV_10d5m_v);
PH_OO_10d5m_v = (1-r)*h0 + r.*mean((1-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v); 
PH_OOV_10d5m_v = (1-r)*h0 + r.*mean((1-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v);


% PSV=unconditional prob of doing the survey (not cond on opening door)  
% PSV < PH mechanically

% 0d5m
PSV_NF_NI_0d5m_v = h0.*mean(doesSvyIfAsked_NI_0d5m_v);
PSV_NF_I_0d5m_v = h0.*mean(doesSvyIfAsked_I_0d5m_v);

PSV_F_NI_0d5m_v = (1-r)*PSV_NF_NI_0d5m_v + r*mean(hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v);
PSV_F_I_0d5m_v = (1-r)*PSV_NF_I_0d5m_v + r*mean(hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v);

PSV_FV_NI_0d5m_v = (1-r)*PSV_NF_NI_0d5m_v + r*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v);
PSV_FV_I_0d5m_v = (1-r)*PSV_NF_I_0d5m_v + r*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v);

PSV_OO_NI_0d5m_v = (1-r)*PSV_NF_NI_0d5m_v + r*mean((1-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v);
PSV_OO_I_0d5m_v = (1-r)*PSV_NF_I_0d5m_v + r*mean((1-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v);

PSV_OOV_NI_0d5m_v = (1-r)*PSV_NF_NI_0d5m_v + r*mean((1-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v);
PSV_OOV_I_0d5m_v = (1-r)*PSV_NF_I_0d5m_v + r*mean((1-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v);

% 10d10m
PSV_NF_NI_10d10m_v = h0.*mean(doesSvyIfAsked_NI_10d10m_v);
PSV_NF_I_10d10m_v = h0.*mean(doesSvyIfAsked_I_10d10m_v);
 
PSV_F_NI_10d10m_v = (1-r)*PSV_NF_NI_10d10m_v + r*mean(hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v);
PSV_F_I_10d10m_v = (1-r)*PSV_NF_I_10d10m_v + r*mean(hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v);
 
PSV_FV_NI_10d10m_v = (1-r)*PSV_NF_NI_10d10m_v + r*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v);
PSV_FV_I_10d10m_v = (1-r)*PSV_NF_I_10d10m_v + r*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v);
 
PSV_OO_NI_10d10m_v = (1-r)*PSV_NF_NI_10d10m_v + r*mean((1-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v);
PSV_OO_I_10d10m_v = (1-r)*PSV_NF_I_10d10m_v + r*mean((1-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v);
 
PSV_OOV_NI_10d10m_v = (1-r)*PSV_NF_NI_10d10m_v + r*mean((1-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v);
PSV_OOV_I_10d10m_v = (1-r)*PSV_NF_I_10d10m_v + r*mean((1-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v);

% 10d5m
PSV_NF_NI_10d5m_v = h0.*mean(doesSvyIfAsked_NI_10d5m_v);
PSV_NF_I_10d5m_v = h0.*mean(doesSvyIfAsked_I_10d5m_v);
 
PSV_F_NI_10d5m_v = (1-r)*PSV_NF_NI_10d5m_v + r*mean(hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v);
PSV_F_I_10d5m_v = (1-r)*PSV_NF_I_10d5m_v + r*mean(hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v);
 
PSV_FV_NI_10d5m_v = (1-r)*PSV_NF_NI_10d5m_v + r*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v);
PSV_FV_I_10d5m_v = (1-r)*PSV_NF_I_10d5m_v + r*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v);
 
PSV_OO_NI_10d5m_v = (1-r)*PSV_NF_NI_10d5m_v + r*mean((1-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v);
PSV_OO_I_10d5m_v = (1-r)*PSV_NF_I_10d5m_v + r*mean((1-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v);
 
PSV_OOV_NI_10d5m_v = (1-r)*PSV_NF_NI_10d5m_v + r*mean((1-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v);
PSV_OOV_I_10d5m_v = (1-r)*PSV_NF_I_10d5m_v + r*mean((1-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v);


% POO=prob of opting out (not conditional on seeing flyer)
% Scaled by baseline likelihood of being at home
POO_OO_0d5m_v =  h0*r*mean(optsOutIfSees_OO_0d5m_v);
POO_OOV_0d5m_v = h0*r*mean(optsOutIfSees_OOV_0d5m_v);

POO_OO_10d10m_v =  h0*r*mean(optsOutIfSees_OO_10d10m_v);
POO_OOV_10d10m_v = h0*r*mean(optsOutIfSees_OOV_10d10m_v);

POO_OO_10d5m_v =  h0*r*mean(optsOutIfSees_OO_10d5m_v);
POO_OOV_10d5m_v = h0*r*mean(optsOutIfSees_OOV_10d5m_v);


% empirical moments are total lying in treatments / total doing
% survey in treatments

% PSVL = unconditional percent who do survey and lie 
% No flyer treatment only, simplifies later code

% PL=cond on agreeing to do the survey, did you lie?
% incentive to lie is a surprise later (doesn't affect PH or PSV)

% 0d5m, 5d1m incentive
PSVL_NF_NI_0d5m_v = mean(h0.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_v);
PSVL_NF_I_0d5m_v = mean(h0.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_NF_NI_0d5m_5d1m_v = mean(h0.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_NF_I_0d5m_5d1m_v = mean(h0.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);

PSVL_F_NI_0d5m_v = (1-r)*PSVL_NF_NI_0d5m_v + r*mean(hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_v);
PSVL_F_I_0d5m_v = (1-r)*PSVL_NF_I_0d5m_v + r*mean(hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_F_NI_0d5m_5d1m_v = (1-r)*PSVL_NF_NI_0d5m_5d1m_v + r*mean(hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_F_I_0d5m_5d1m_v = (1-r)*PSVL_NF_I_0d5m_5d1m_v + r*mean(hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);

PSVL_FV_NI_0d5m_v = (1-r)*PSVL_NF_NI_0d5m_v + r*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_FV_I_0d5m_v = (1-r)*PSVL_NF_I_0d5m_v + r*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_FV_NI_0d5m_5d1m_v = (1-r)*PSVL_NF_NI_0d5m_5d1m_v + r*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_FV_I_0d5m_5d1m_v = (1-r)*PSVL_NF_I_0d5m_5d1m_v + r*mean(hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);

PSVL_OO_NI_0d5m_v = (1-r)*PSVL_NF_NI_0d5m_v + r*mean((1-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_v);
PSVL_OO_I_0d5m_v = (1-r)*PSVL_NF_I_0d5m_v + r*mean((1-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_OO_NI_0d5m_5d1m_v = (1-r)*PSVL_NF_NI_0d5m_5d1m_v + r*mean((1-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_NI_0d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_OO_I_0d5m_5d1m_v = (1-r)*PSVL_NF_I_0d5m_5d1m_v + r*mean((1-optsOutIfSees_OO_0d5m_v).*hStar_F_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);

PSVL_OOV_NI_0d5m_v = (1-r)*PSVL_NF_NI_0d5m_v + r*mean((1-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_OOV_I_0d5m_v = (1-r)*PSVL_NF_I_0d5m_v + r*mean((1-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_v);
PSVL_OOV_NI_0d5m_5d1m_v = (1-r)*PSVL_NF_NI_0d5m_5d1m_v + r*mean((1-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_OOV_I_0d5m_5d1m_v = (1-r)*PSVL_NF_I_0d5m_5d1m_v + r*mean((1-optsOutIfSees_OOV_0d5m_v).*hStar_FV_0d5m_v.*doesSvyIfAsked_I_0d5m_v.*wouldLieIfAsked_5d1m_v);


% 10d10m, 8m incentive
PSVL_NF_NI_10d10m_v = mean(h0.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_v);
PSVL_NF_I_10d10m_v = mean(h0.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_NF_NI_10d10m_8m_v = mean(h0.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_8m_v);
PSVL_NF_I_10d10m_8m_v = mean(h0.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
 
PSVL_F_NI_10d10m_v = (1-r)*PSVL_NF_NI_10d10m_v + r*mean(hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_v);
PSVL_F_I_10d10m_v = (1-r)*PSVL_NF_I_10d10m_v + r*mean(hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_F_NI_10d10m_8m_v = (1-r)*PSVL_NF_NI_10d10m_8m_v + r*mean(hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_8m_v);
PSVL_F_I_10d10m_8m_v = (1-r)*PSVL_NF_I_10d10m_8m_v + r*mean(hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
 
PSVL_FV_NI_10d10m_v = (1-r)*PSVL_NF_NI_10d10m_v + r*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_FV_I_10d10m_v = (1-r)*PSVL_NF_I_10d10m_v + r*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_FV_NI_10d10m_8m_v = (1-r)*PSVL_NF_NI_10d10m_8m_v + r*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
PSVL_FV_I_10d10m_8m_v = (1-r)*PSVL_NF_I_10d10m_8m_v + r*mean(hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
 
PSVL_OO_NI_10d10m_v = (1-r)*PSVL_NF_NI_10d10m_v + r*mean((1-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_v);
PSVL_OO_I_10d10m_v = (1-r)*PSVL_NF_I_10d10m_v + r*mean((1-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_OO_NI_10d10m_8m_v = (1-r)*PSVL_NF_NI_10d10m_8m_v + r*mean((1-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_NI_10d10m_v.*wouldLieIfAsked_8m_v);
PSVL_OO_I_10d10m_8m_v = (1-r)*PSVL_NF_I_10d10m_8m_v + r*mean((1-optsOutIfSees_OO_10d10m_v).*hStar_F_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
 
PSVL_OOV_NI_10d10m_v = (1-r)*PSVL_NF_NI_10d10m_v + r*mean((1-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_OOV_I_10d10m_v = (1-r)*PSVL_NF_I_10d10m_v + r*mean((1-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_v);
PSVL_OOV_NI_10d10m_8m_v = (1-r)*PSVL_NF_NI_10d10m_8m_v + r*mean((1-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);
PSVL_OOV_I_10d10m_8m_v = (1-r)*PSVL_NF_I_10d10m_8m_v + r*mean((1-optsOutIfSees_OOV_10d10m_v).*hStar_FV_10d10m_v.*doesSvyIfAsked_I_10d10m_v.*wouldLieIfAsked_8m_v);



% 10d5m, 5d1m incentive
PSVL_NF_NI_10d5m_v = mean(h0.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_v);
PSVL_NF_I_10d5m_v = mean(h0.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_NF_NI_10d5m_5d1m_v = mean(h0.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_NF_I_10d5m_5d1m_v = mean(h0.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
 
PSVL_F_NI_10d5m_v = (1-r)*PSVL_NF_NI_10d5m_v + r*mean(hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_v);
PSVL_F_I_10d5m_v = (1-r)*PSVL_NF_I_10d5m_v + r*mean(hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_F_NI_10d5m_5d1m_v = (1-r)*PSVL_NF_NI_10d5m_5d1m_v + r*mean(hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_F_I_10d5m_5d1m_v = (1-r)*PSVL_NF_I_10d5m_5d1m_v + r*mean(hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
 
PSVL_FV_NI_10d5m_v = (1-r)*PSVL_NF_NI_10d5m_v + r*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_FV_I_10d5m_v = (1-r)*PSVL_NF_I_10d5m_v + r*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_FV_NI_10d5m_5d1m_v = (1-r)*PSVL_NF_NI_10d5m_5d1m_v + r*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_FV_I_10d5m_5d1m_v = (1-r)*PSVL_NF_I_10d5m_5d1m_v + r*mean(hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
 
PSVL_OO_NI_10d5m_v = (1-r)*PSVL_NF_NI_10d5m_v + r*mean((1-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_v);
PSVL_OO_I_10d5m_v = (1-r)*PSVL_NF_I_10d5m_v + r*mean((1-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_OO_NI_10d5m_5d1m_v = (1-r)*PSVL_NF_NI_10d5m_5d1m_v + r*mean((1-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_NI_10d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_OO_I_10d5m_5d1m_v = (1-r)*PSVL_NF_I_10d5m_5d1m_v + r*mean((1-optsOutIfSees_OO_10d5m_v).*hStar_F_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
 
PSVL_OOV_NI_10d5m_v = (1-r)*PSVL_NF_NI_10d5m_v + r*mean((1-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_OOV_I_10d5m_v = (1-r)*PSVL_NF_I_10d5m_v + r*mean((1-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_v);
PSVL_OOV_NI_10d5m_5d1m_v = (1-r)*PSVL_NF_NI_10d5m_5d1m_v + r*mean((1-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);
PSVL_OOV_I_10d5m_5d1m_v = (1-r)*PSVL_NF_I_10d5m_5d1m_v + r*mean((1-optsOutIfSees_OOV_10d5m_v).*hStar_FV_10d5m_v.*doesSvyIfAsked_I_10d5m_v.*wouldLieIfAsked_5d1m_v);




% non-voters
% ===========================

%%% PH = probability of being at home
PH_NF_0d5m_nv = h0;
PH_F_0d5m_nv = (1-r)*h0 + r.*mean(hStar_F_0d5m_nv);
PH_FV_0d5m_nv = (1-r)*h0 + r.*mean(hStar_FV_0d5m_nv);
PH_OO_0d5m_nv = (1-r)*h0 + r.*mean((1-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv); 
PH_OOV_0d5m_nv = (1-r)*h0 + r.*mean((1-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv);
 
PH_NF_10d10m_nv = h0;
PH_F_10d10m_nv = (1-r)*h0 + r.*mean(hStar_F_10d10m_nv);
PH_FV_10d10m_nv = (1-r)*h0 + r.*mean(hStar_FV_10d10m_nv);
PH_OO_10d10m_nv = (1-r)*h0 + r.*mean((1-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv); 
PH_OOV_10d10m_nv = (1-r)*h0 + r.*mean((1-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv);
 
PH_NF_10d5m_nv = h0;
PH_F_10d5m_nv = (1-r)*h0 + r.*mean(hStar_F_10d5m_nv);
PH_FV_10d5m_nv = (1-r)*h0 + r.*mean(hStar_FV_10d5m_nv);
PH_OO_10d5m_nv = (1-r)*h0 + r.*mean((1-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv); 
PH_OOV_10d5m_nv = (1-r)*h0 + r.*mean((1-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv);
 
 
% PSV=unconditional prob of doing the survey (not cond on opening door)  
% PSV < PH mechanically
 
% 0d5m
PSV_NF_NI_0d5m_nv = h0.*mean(doesSvyIfAsked_NI_0d5m_nv);
PSV_NF_I_0d5m_nv = h0.*mean(doesSvyIfAsked_I_0d5m_nv);
 
PSV_F_NI_0d5m_nv = (1-r)*PSV_NF_NI_0d5m_nv + r*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv);
PSV_F_I_0d5m_nv = (1-r)*PSV_NF_I_0d5m_nv + r*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
 
PSV_FV_NI_0d5m_nv = (1-r)*PSV_NF_NI_0d5m_nv + r*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
PSV_FV_I_0d5m_nv = (1-r)*PSV_NF_I_0d5m_nv + r*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
 
PSV_OO_NI_0d5m_nv = (1-r)*PSV_NF_NI_0d5m_nv + r*mean((1-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv);
PSV_OO_I_0d5m_nv = (1-r)*PSV_NF_I_0d5m_nv + r*mean((1-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
 
PSV_OOV_NI_0d5m_nv = (1-r)*PSV_NF_NI_0d5m_nv + r*mean((1-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
PSV_OOV_I_0d5m_nv = (1-r)*PSV_NF_I_0d5m_nv + r*mean((1-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv);
 
% 10d10m
PSV_NF_NI_10d10m_nv = h0.*mean(doesSvyIfAsked_NI_10d10m_nv);
PSV_NF_I_10d10m_nv = h0.*mean(doesSvyIfAsked_I_10d10m_nv);
 
PSV_F_NI_10d10m_nv = (1-r)*PSV_NF_NI_10d10m_nv + r*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv);
PSV_F_I_10d10m_nv = (1-r)*PSV_NF_I_10d10m_nv + r*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
 
PSV_FV_NI_10d10m_nv = (1-r)*PSV_NF_NI_10d10m_nv + r*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
PSV_FV_I_10d10m_nv = (1-r)*PSV_NF_I_10d10m_nv + r*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
 
PSV_OO_NI_10d10m_nv = (1-r)*PSV_NF_NI_10d10m_nv + r*mean((1-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv);
PSV_OO_I_10d10m_nv = (1-r)*PSV_NF_I_10d10m_nv + r*mean((1-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
 
PSV_OOV_NI_10d10m_nv = (1-r)*PSV_NF_NI_10d10m_nv + r*mean((1-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
PSV_OOV_I_10d10m_nv = (1-r)*PSV_NF_I_10d10m_nv + r*mean((1-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv);
 
% 10d5m
PSV_NF_NI_10d5m_nv = h0.*mean(doesSvyIfAsked_NI_10d5m_nv);
PSV_NF_I_10d5m_nv = h0.*mean(doesSvyIfAsked_I_10d5m_nv);
 
PSV_F_NI_10d5m_nv = (1-r)*PSV_NF_NI_10d5m_nv + r*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv);
PSV_F_I_10d5m_nv = (1-r)*PSV_NF_I_10d5m_nv + r*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
 
PSV_FV_NI_10d5m_nv = (1-r)*PSV_NF_NI_10d5m_nv + r*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
PSV_FV_I_10d5m_nv = (1-r)*PSV_NF_I_10d5m_nv + r*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
 
PSV_OO_NI_10d5m_nv = (1-r)*PSV_NF_NI_10d5m_nv + r*mean((1-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv);
PSV_OO_I_10d5m_nv = (1-r)*PSV_NF_I_10d5m_nv + r*mean((1-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
 
PSV_OOV_NI_10d5m_nv = (1-r)*PSV_NF_NI_10d5m_nv + r*mean((1-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
PSV_OOV_I_10d5m_nv = (1-r)*PSV_NF_I_10d5m_nv + r*mean((1-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv);
 
 
% POO=prob of opting out (not conditional on seeing flyer)
% Scaled by baseline likelihood of being at home
POO_OO_0d5m_nv =  h0*r*mean(optsOutIfSees_OO_0d5m_nv);
POO_OOV_0d5m_nv = h0*r*mean(optsOutIfSees_OOV_0d5m_nv);
 
POO_OO_10d10m_nv =  h0*r*mean(optsOutIfSees_OO_10d10m_nv);
POO_OOV_10d10m_nv = h0*r*mean(optsOutIfSees_OOV_10d10m_nv);
 
POO_OO_10d5m_nv =  h0*r*mean(optsOutIfSees_OO_10d5m_nv);
POO_OOV_10d5m_nv = h0*r*mean(optsOutIfSees_OOV_10d5m_nv);
 
% PSVL = unconditional percent who do survey and lie 
% No flyer treatment only, simplifies later code

% PL=cond on agreeing to do the survey, did you lie?
% incentive to lie is a surprise later (doesn't affect PH or PSV)

% 0d5m, 5d1m incentive
PSVL_NF_NI_0d5m_nv = mean(h0.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_NF_I_0d5m_nv = mean(h0.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_NF_NI_0d5m_5d1m_nv = mean(h0.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_NF_I_0d5m_5d1m_nv = mean(h0.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_F_NI_0d5m_nv = (1-r)*PSVL_NF_NI_0d5m_nv + r*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_F_I_0d5m_nv = (1-r)*PSVL_NF_I_0d5m_nv + r*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_F_NI_0d5m_5d1m_nv = (1-r)*PSVL_NF_NI_0d5m_5d1m_nv + r*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_F_I_0d5m_5d1m_nv = (1-r)*PSVL_NF_I_0d5m_5d1m_nv + r*mean(hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_FV_NI_0d5m_nv = (1-r)*PSVL_NF_NI_0d5m_nv + r*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_FV_I_0d5m_nv = (1-r)*PSVL_NF_I_0d5m_nv + r*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_FV_NI_0d5m_5d1m_nv = (1-r)*PSVL_NF_NI_0d5m_5d1m_nv + r*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_FV_I_0d5m_5d1m_nv = (1-r)*PSVL_NF_I_0d5m_5d1m_nv + r*mean(hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_OO_NI_0d5m_nv = (1-r)*PSVL_NF_NI_0d5m_nv + r*mean((1-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_OO_I_0d5m_nv = (1-r)*PSVL_NF_I_0d5m_nv + r*mean((1-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_OO_NI_0d5m_5d1m_nv = (1-r)*PSVL_NF_NI_0d5m_5d1m_nv + r*mean((1-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_NI_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_OO_I_0d5m_5d1m_nv = (1-r)*PSVL_NF_I_0d5m_5d1m_nv + r*mean((1-optsOutIfSees_OO_0d5m_nv).*hStar_F_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_OOV_NI_0d5m_nv = (1-r)*PSVL_NF_NI_0d5m_nv + r*mean((1-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_I_0d5m_nv = (1-r)*PSVL_NF_I_0d5m_nv + r*mean((1-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_NI_0d5m_5d1m_nv = (1-r)*PSVL_NF_NI_0d5m_5d1m_nv + r*mean((1-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_OOV_I_0d5m_5d1m_nv = (1-r)*PSVL_NF_I_0d5m_5d1m_nv + r*mean((1-optsOutIfSees_OOV_0d5m_nv).*hStar_FV_0d5m_nv.*doesSvyIfAsked_I_0d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
 
% 10d10m, 8m incentive
PSVL_NF_NI_10d10m_nv = mean(h0.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_NF_I_10d10m_nv = mean(h0.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_NF_NI_10d10m_8m_nv = mean(h0.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_8m_nv);
PSVL_NF_I_10d10m_8m_nv = mean(h0.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
 
PSVL_F_NI_10d10m_nv = (1-r)*PSVL_NF_NI_10d10m_nv + r*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_F_I_10d10m_nv = (1-r)*PSVL_NF_I_10d10m_nv + r*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_F_NI_10d10m_8m_nv = (1-r)*PSVL_NF_NI_10d10m_8m_nv + r*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_8m_nv);
PSVL_F_I_10d10m_8m_nv = (1-r)*PSVL_NF_I_10d10m_8m_nv + r*mean(hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
 
PSVL_FV_NI_10d10m_nv = (1-r)*PSVL_NF_NI_10d10m_nv + r*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_FV_I_10d10m_nv = (1-r)*PSVL_NF_I_10d10m_nv + r*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_FV_NI_10d10m_8m_nv = (1-r)*PSVL_NF_NI_10d10m_8m_nv + r*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
PSVL_FV_I_10d10m_8m_nv = (1-r)*PSVL_NF_I_10d10m_8m_nv + r*mean(hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
 
PSVL_OO_NI_10d10m_nv = (1-r)*PSVL_NF_NI_10d10m_nv + r*mean((1-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_OO_I_10d10m_nv = (1-r)*PSVL_NF_I_10d10m_nv + r*mean((1-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_OO_NI_10d10m_8m_nv = (1-r)*PSVL_NF_NI_10d10m_8m_nv + r*mean((1-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_NI_10d10m_nv.*wouldLieIfAsked_8m_nv);
PSVL_OO_I_10d10m_8m_nv = (1-r)*PSVL_NF_I_10d10m_8m_nv + r*mean((1-optsOutIfSees_OO_10d10m_nv).*hStar_F_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
 
PSVL_OOV_NI_10d10m_nv = (1-r)*PSVL_NF_NI_10d10m_nv + r*mean((1-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_I_10d10m_nv = (1-r)*PSVL_NF_I_10d10m_nv + r*mean((1-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_NI_10d10m_8m_nv = (1-r)*PSVL_NF_NI_10d10m_8m_nv + r*mean((1-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
PSVL_OOV_I_10d10m_8m_nv = (1-r)*PSVL_NF_I_10d10m_8m_nv + r*mean((1-optsOutIfSees_OOV_10d10m_nv).*hStar_FV_10d10m_nv.*doesSvyIfAsked_I_10d10m_nv.*wouldLieIfAsked_8m_nv);
 
 
 
% 10d5m, 5d1m incentive
PSVL_NF_NI_10d5m_nv = mean(h0.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_NF_I_10d5m_nv = mean(h0.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_NF_NI_10d5m_5d1m_nv = mean(h0.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_NF_I_10d5m_5d1m_nv = mean(h0.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_F_NI_10d5m_nv = (1-r)*PSVL_NF_NI_10d5m_nv + r*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_F_I_10d5m_nv = (1-r)*PSVL_NF_I_10d5m_nv + r*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_F_NI_10d5m_5d1m_nv = (1-r)*PSVL_NF_NI_10d5m_5d1m_nv + r*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_F_I_10d5m_5d1m_nv = (1-r)*PSVL_NF_I_10d5m_5d1m_nv + r*mean(hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_FV_NI_10d5m_nv = (1-r)*PSVL_NF_NI_10d5m_nv + r*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_FV_I_10d5m_nv = (1-r)*PSVL_NF_I_10d5m_nv + r*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_FV_NI_10d5m_5d1m_nv = (1-r)*PSVL_NF_NI_10d5m_5d1m_nv + r*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_FV_I_10d5m_5d1m_nv = (1-r)*PSVL_NF_I_10d5m_5d1m_nv + r*mean(hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_OO_NI_10d5m_nv = (1-r)*PSVL_NF_NI_10d5m_nv + r*mean((1-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_OO_I_10d5m_nv = (1-r)*PSVL_NF_I_10d5m_nv + r*mean((1-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_OO_NI_10d5m_5d1m_nv = (1-r)*PSVL_NF_NI_10d5m_5d1m_nv + r*mean((1-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_NI_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_OO_I_10d5m_5d1m_nv = (1-r)*PSVL_NF_I_10d5m_5d1m_nv + r*mean((1-optsOutIfSees_OO_10d5m_nv).*hStar_F_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
 
PSVL_OOV_NI_10d5m_nv = (1-r)*PSVL_NF_NI_10d5m_nv + r*mean((1-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_I_10d5m_nv = (1-r)*PSVL_NF_I_10d5m_nv + r*mean((1-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_nv);
PSVL_OOV_NI_10d5m_5d1m_nv = (1-r)*PSVL_NF_NI_10d5m_5d1m_nv + r*mean((1-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);
PSVL_OOV_I_10d5m_5d1m_nv = (1-r)*PSVL_NF_I_10d5m_5d1m_nv + r*mean((1-optsOutIfSees_OOV_10d5m_nv).*hStar_FV_10d5m_nv.*doesSvyIfAsked_I_10d5m_nv.*wouldLieIfAsked_5d1m_nv);


% =================================================================
% Form Moments - this section doesn't change
% ================================================================

ph=NaN(30,1); 
psv=NaN(30,1);
poo=NaN(12,1);
psvini=NaN(20,1);
pl=NaN(8,1);
pv=NaN(4,1);

% 30 PH moments (5 tx * 3 length * 2 v/nv)
ph(1)	=		PH_NF_0d5m_v	;		
ph(2)	=		PH_NF_10d10m_v	;		
ph(3)	=		PH_NF_10d5m_v	;		
ph(4)	=		PH_F_0d5m_v	;		
ph(5)	=		PH_F_10d10m_v	;		
ph(6)	=		PH_F_10d5m_v	;		
ph(7)	=		PH_FV_0d5m_v	;		
ph(8)	=		PH_FV_10d10m_v	;		
ph(9)	=		PH_FV_10d5m_v	;		
ph(10)	=		PH_OO_0d5m_v	;		
ph(11)	=		PH_OO_10d10m_v	;		
ph(12)	=		PH_OO_10d5m_v	;		
ph(13)	=		PH_OOV_0d5m_v	;		
ph(14)	=		PH_OOV_10d10m_v	;		
ph(15)	=		PH_OOV_10d5m_v	;		
ph(16)	=		PH_NF_0d5m_nv	;		
ph(17)	=		PH_NF_10d10m_nv	;		
ph(18)	=		PH_NF_10d5m_nv	;		
ph(19)	=		PH_F_0d5m_nv	;		
ph(20)	=		PH_F_10d10m_nv	;		
ph(21)	=		PH_F_10d5m_nv	;		
ph(22)	=		PH_FV_0d5m_nv	;		
ph(23)	=		PH_FV_10d10m_nv	;		
ph(24)	=		PH_FV_10d5m_nv	;		
ph(25)	=		PH_OO_0d5m_nv	;		
ph(26)	=		PH_OO_10d10m_nv	;		
ph(27)	=		PH_OO_10d5m_nv	;		
ph(28)	=		PH_OOV_0d5m_nv	;		
ph(29)	=		PH_OOV_10d10m_nv	;		
ph(30)	=		PH_OOV_10d5m_nv	;		

% 30 PSV moments (5 tx * 3 length * 2 v/nv)
% 50/50 average across NI and I treatments

psv(1)  =   mean([  PSV_NF_I_0d5m_v     PSV_NF_NI_0d5m_v    ]); 
psv(2)  =   mean([  PSV_NF_I_10d10m_v   PSV_NF_NI_10d10m_v  ]); 
psv(3)  =   mean([  PSV_NF_I_10d5m_v    PSV_NF_NI_10d5m_v   ]); 
psv(4)  =   mean([  PSV_F_I_0d5m_v      PSV_F_NI_0d5m_v     ]); 
psv(5)  =   mean([  PSV_F_I_10d10m_v    PSV_F_NI_10d10m_v   ]); 
psv(6)  =   mean([  PSV_F_I_10d5m_v     PSV_F_NI_10d5m_v    ]); 
psv(7)  =   mean([  PSV_FV_I_0d5m_v     PSV_FV_NI_0d5m_v    ]); 
psv(8)  =   mean([  PSV_FV_I_10d10m_v   PSV_FV_NI_10d10m_v  ]); 
psv(9)  =   mean([  PSV_FV_I_10d5m_v    PSV_FV_NI_10d5m_v   ]); 
psv(10) =   mean([  PSV_OO_I_0d5m_v     PSV_OO_NI_0d5m_v    ]); 
psv(11) =   mean([  PSV_OO_I_10d10m_v   PSV_OO_NI_10d10m_v  ]); 
psv(12) =   mean([  PSV_OO_I_10d5m_v    PSV_OO_NI_10d5m_v   ]); 
psv(13) =   mean([  PSV_OOV_I_0d5m_v    PSV_OOV_NI_0d5m_v   ]); 
psv(14) =   mean([  PSV_OOV_I_10d10m_v  PSV_OOV_NI_10d10m_v ]); 
psv(15) =   mean([  PSV_OOV_I_10d5m_v   PSV_OOV_NI_10d5m_v  ]); 
psv(16) =   mean([  PSV_NF_I_0d5m_nv    PSV_NF_NI_0d5m_nv   ]); 
psv(17) =   mean([  PSV_NF_I_10d10m_nv  PSV_NF_NI_10d10m_nv ]); 
psv(18) =   mean([  PSV_NF_I_10d5m_nv   PSV_NF_NI_10d5m_nv  ]); 
psv(19) =   mean([  PSV_F_I_0d5m_nv     PSV_F_NI_0d5m_nv    ]); 
psv(20) =   mean([  PSV_F_I_10d10m_nv   PSV_F_NI_10d10m_nv  ]); 
psv(21) =   mean([  PSV_F_I_10d5m_nv    PSV_F_NI_10d5m_nv   ]); 
psv(22) =   mean([  PSV_FV_I_0d5m_nv    PSV_FV_NI_0d5m_nv   ]); 
psv(23) =   mean([  PSV_FV_I_10d10m_nv  PSV_FV_NI_10d10m_nv ]); 
psv(24) =   mean([  PSV_FV_I_10d5m_nv   PSV_FV_NI_10d5m_nv  ]); 
psv(25) =   mean([  PSV_OO_I_0d5m_nv    PSV_OO_NI_0d5m_nv   ]); 
psv(26) =   mean([  PSV_OO_I_10d10m_nv  PSV_OO_NI_10d10m_nv ]); 
psv(27) =   mean([  PSV_OO_I_10d5m_nv   PSV_OO_NI_10d5m_nv  ]); 
psv(28) =   mean([  PSV_OOV_I_0d5m_nv   PSV_OOV_NI_0d5m_nv  ]); 
psv(29) =   mean([  PSV_OOV_I_10d10m_nv PSV_OOV_NI_10d10m_nv    ]); 
psv(30) =   mean([  PSV_OOV_I_10d5m_nv  PSV_OOV_NI_10d5m_nv ]); 

% 12 POO moments (2 tx * 3 length * 2 v/nv)
poo(1)	=		POO_OO_0d5m_v	;		
poo(2)	=		POO_OO_10d10m_v	;		
poo(3)	=		POO_OO_10d5m_v	;		
poo(4)	=		POO_OOV_0d5m_v	;		
poo(5)	=		POO_OOV_10d10m_v	;		
poo(6)	=		POO_OOV_10d5m_v	;		
poo(7)	=		POO_OO_0d5m_nv	;		
poo(8)	=		POO_OO_10d10m_nv	;		
poo(9)	=		POO_OO_10d5m_nv	;		
poo(10)	=		POO_OOV_0d5m_nv	;		
poo(11)	=		POO_OOV_10d10m_nv	;		
poo(12)	=		POO_OOV_10d5m_nv	;	

% 20 PSV by info moments (5 tx * 2 I/NI * 2 v/nv)
psvini(1)   =   mean([  PSV_NF_NI_0d5m_v    PSV_NF_NI_10d10m_v      PSV_NF_NI_10d5m_v   ]);
psvini(2)   =   mean([  PSV_NF_I_0d5m_v     PSV_NF_I_10d10m_v       PSV_NF_I_10d5m_v  ]);
psvini(3)   =   mean([  PSV_F_NI_0d5m_v     PSV_F_NI_10d10m_v       PSV_F_NI_10d5m_v    ]);
psvini(4)   =   mean([  PSV_F_I_0d5m_v      PSV_F_I_10d10m_v        PSV_F_I_10d5m_v   ]);
psvini(5)   =   mean([  PSV_FV_NI_0d5m_v    PSV_FV_NI_10d10m_v      PSV_FV_NI_10d5m_v   ]);
psvini(6)   =   mean([  PSV_FV_I_0d5m_v     PSV_FV_I_10d10m_v       PSV_FV_I_10d5m_v  ]);
psvini(7)   =   mean([  PSV_OO_NI_0d5m_v    PSV_OO_NI_10d10m_v      PSV_OO_NI_10d5m_v   ]);
psvini(8)   =   mean([  PSV_OO_I_0d5m_v     PSV_OO_I_10d10m_v       PSV_OO_I_10d5m_v  ]);
psvini(9)   =   mean([  PSV_OOV_NI_0d5m_v   PSV_OOV_NI_10d10m_v     PSV_OOV_NI_10d5m_v  ]);
psvini(10)  =   mean([  PSV_OOV_I_0d5m_v    PSV_OOV_I_10d10m_v      PSV_OOV_I_10d5m_v ]);
psvini(11)  =   mean([  PSV_NF_NI_0d5m_nv   PSV_NF_NI_10d10m_nv     PSV_NF_NI_10d5m_nv  ]);
psvini(12)  =   mean([  PSV_NF_I_0d5m_nv    PSV_NF_I_10d10m_nv      PSV_NF_I_10d5m_nv ]);
psvini(13)  =   mean([  PSV_F_NI_0d5m_nv    PSV_F_NI_10d10m_nv      PSV_F_NI_10d5m_nv   ]);
psvini(14)  =   mean([  PSV_F_I_0d5m_nv     PSV_F_I_10d10m_nv       PSV_F_I_10d5m_nv  ]);
psvini(15)  =   mean([  PSV_FV_NI_0d5m_nv   PSV_FV_NI_10d10m_nv     PSV_FV_NI_10d5m_nv  ]);
psvini(16)  =   mean([  PSV_FV_I_0d5m_nv    PSV_FV_I_10d10m_nv      PSV_FV_I_10d5m_nv ]);
psvini(17)  =   mean([  PSV_OO_NI_0d5m_nv   PSV_OO_NI_10d10m_nv     PSV_OO_NI_10d5m_nv  ]);
psvini(18)  =   mean([  PSV_OO_I_0d5m_nv    PSV_OO_I_10d10m_nv      PSV_OO_I_10d5m_nv ]);
psvini(19)  =   mean([  PSV_OOV_NI_0d5m_nv  PSV_OOV_NI_10d10m_nv    PSV_OOV_NI_10d5m_nv ]);
psvini(20)  =   mean([  PSV_OOV_I_0d5m_nv   PSV_OOV_I_10d10m_nv     PSV_OOV_I_10d5m_nv ]);

% 8 PL moments (1 tx * 2 10m/5m * 2 incentives)

% Empirical moments are sum of people lying in relevant tx
% divided by the sum of people answering the survey in relevant tx.

pl(1)   =   mean([  PSVL_NF_NI_0d5m_v PSVL_NF_I_0d5m_v PSVL_NF_NI_10d5m_v PSVL_NF_I_10d5m_v ...
                    PSVL_F_NI_0d5m_v PSVL_F_I_0d5m_v PSVL_F_NI_10d5m_v PSVL_F_I_10d5m_v ...
                    PSVL_FV_NI_0d5m_v PSVL_FV_I_0d5m_v PSVL_FV_NI_10d5m_v PSVL_FV_I_10d5m_v ...
                    PSVL_OO_NI_0d5m_v PSVL_OO_I_0d5m_v PSVL_OO_NI_10d5m_v PSVL_OO_I_10d5m_v ...
                    PSVL_OOV_NI_0d5m_v PSVL_OOV_I_0d5m_v PSVL_OOV_NI_10d5m_v PSVL_OOV_I_10d5m_v ]) / ...
            mean([  PSV_NF_NI_0d5m_v PSV_NF_I_0d5m_v PSV_NF_NI_10d5m_v PSV_NF_I_10d5m_v ...
                    PSV_F_NI_0d5m_v PSV_F_I_0d5m_v PSV_F_NI_10d5m_v PSV_F_I_10d5m_v ...
                    PSV_FV_NI_0d5m_v PSV_FV_I_0d5m_v PSV_FV_NI_10d5m_v PSV_FV_I_10d5m_v ...
                    PSV_OO_NI_0d5m_v PSV_OO_I_0d5m_v  PSV_OO_NI_10d5m_v PSV_OO_I_10d5m_v ...
                    PSV_OOV_NI_0d5m_v PSV_OOV_I_0d5m_v PSV_OOV_NI_10d5m_v PSV_OOV_I_10d5m_v ]) ;    % v, 5m, no incentive

pl(2)   =   mean([  PSVL_NF_NI_0d5m_5d1m_v PSVL_NF_I_0d5m_5d1m_v PSVL_NF_NI_10d5m_5d1m_v PSVL_NF_I_10d5m_5d1m_v ...
                    PSVL_F_NI_0d5m_5d1m_v PSVL_F_I_0d5m_5d1m_v PSVL_F_NI_10d5m_5d1m_v PSVL_F_I_10d5m_5d1m_v ...
                    PSVL_FV_NI_0d5m_5d1m_v PSVL_FV_I_0d5m_5d1m_v PSVL_FV_NI_10d5m_5d1m_v PSVL_FV_I_10d5m_5d1m_v ...
                    PSVL_OO_NI_0d5m_5d1m_v PSVL_OO_I_0d5m_5d1m_v PSVL_OO_NI_10d5m_5d1m_v PSVL_OO_I_10d5m_5d1m_v ...
                    PSVL_OOV_NI_0d5m_5d1m_v PSVL_OOV_I_0d5m_5d1m_v PSVL_OOV_NI_10d5m_5d1m_v PSVL_OOV_I_10d5m_5d1m_v ]) / ...
            mean([  PSV_NF_NI_0d5m_v PSV_NF_I_0d5m_v PSV_NF_NI_10d5m_v PSV_NF_I_10d5m_v ...
                    PSV_F_NI_0d5m_v PSV_F_I_0d5m_v PSV_F_NI_10d5m_v PSV_F_I_10d5m_v ...
                    PSV_FV_NI_0d5m_v PSV_FV_I_0d5m_v PSV_FV_NI_10d5m_v PSV_FV_I_10d5m_v ...
                    PSV_OO_NI_0d5m_v PSV_OO_I_0d5m_v  PSV_OO_NI_10d5m_v PSV_OO_I_10d5m_v ...
                    PSV_OOV_NI_0d5m_v PSV_OOV_I_0d5m_v PSV_OOV_NI_10d5m_v PSV_OOV_I_10d5m_v ]) ;    % v, 5m, incentive

pl(3)   =   mean([  PSVL_NF_NI_10d10m_v PSVL_NF_I_10d10m_v ...
                    PSVL_F_NI_10d10m_v PSVL_F_I_10d10m_v ...
                    PSVL_FV_NI_10d10m_v PSVL_FV_I_10d10m_v ...
                    PSVL_OO_NI_10d10m_v PSVL_OO_I_10d10m_v ...
                    PSVL_OOV_NI_10d10m_v PSVL_OOV_I_10d10m_v  ]) / ...
            mean([  PSV_NF_NI_10d10m_v PSV_NF_I_10d10m_v ...
                    PSV_F_NI_10d10m_v PSV_F_I_10d10m_v ...
                    PSV_FV_NI_10d10m_v PSV_FV_I_10d10m_v ...
                    PSV_OO_NI_10d10m_v PSV_OO_I_10d10m_v ...
                    PSV_OOV_NI_10d10m_v PSV_OOV_I_10d10m_v    ]) ;    % v, 10m, no incentive
                
pl(4)   =   mean([  PSVL_NF_NI_10d10m_8m_v PSVL_NF_I_10d10m_8m_v ...
                    PSVL_F_NI_10d10m_8m_v PSVL_F_I_10d10m_8m_v ...
                    PSVL_FV_NI_10d10m_8m_v PSVL_FV_I_10d10m_8m_v ...
                    PSVL_OO_NI_10d10m_8m_v PSVL_OO_I_10d10m_8m_v ...
                    PSVL_OOV_NI_10d10m_8m_v PSVL_OOV_I_10d10m_8m_v    ]) / ...
            mean([  PSV_NF_NI_10d10m_v PSV_NF_I_10d10m_v ...
                    PSV_F_NI_10d10m_v PSV_F_I_10d10m_v ...
                    PSV_FV_NI_10d10m_v PSV_FV_I_10d10m_v ...
                    PSV_OO_NI_10d10m_v PSV_OO_I_10d10m_v ...
                    PSV_OOV_NI_10d10m_v PSV_OOV_I_10d10m_v    ]) ; % v, 10m, incentive
                
pl(5)   =   mean([  PSVL_NF_NI_0d5m_nv PSVL_NF_I_0d5m_nv PSVL_NF_NI_10d5m_nv PSVL_NF_I_10d5m_nv ...
                    PSVL_F_NI_0d5m_nv PSVL_F_I_0d5m_nv PSVL_F_NI_10d5m_nv PSVL_F_I_10d5m_nv ...
                    PSVL_FV_NI_0d5m_nv PSVL_FV_I_0d5m_nv PSVL_FV_NI_10d5m_nv PSVL_FV_I_10d5m_nv ...
                    PSVL_OO_NI_0d5m_nv PSVL_OO_I_0d5m_nv PSVL_OO_NI_10d5m_nv PSVL_OO_I_10d5m_nv ...
                    PSVL_OOV_NI_0d5m_nv PSVL_OOV_I_0d5m_nv PSVL_OOV_NI_10d5m_nv PSVL_OOV_I_10d5m_nv ]) / ...
            mean([  PSV_NF_NI_0d5m_nv PSV_NF_I_0d5m_nv PSV_NF_NI_10d5m_nv PSV_NF_I_10d5m_nv ...
                    PSV_F_NI_0d5m_nv PSV_F_I_0d5m_nv PSV_F_NI_10d5m_nv PSV_F_I_10d5m_nv ...
                    PSV_FV_NI_0d5m_nv PSV_FV_I_0d5m_nv PSV_FV_NI_10d5m_nv PSV_FV_I_10d5m_nv ...
                    PSV_OO_NI_0d5m_nv PSV_OO_I_0d5m_nv  PSV_OO_NI_10d5m_nv PSV_OO_I_10d5m_nv ...
                    PSV_OOV_NI_0d5m_nv PSV_OOV_I_0d5m_nv PSV_OOV_NI_10d5m_nv PSV_OOV_I_10d5m_nv ]) ;    % nv, 5m, no incentive
 
pl(6)   =   mean([  PSVL_NF_NI_0d5m_5d1m_nv PSVL_NF_I_0d5m_5d1m_nv PSVL_NF_NI_10d5m_5d1m_nv PSVL_NF_I_10d5m_5d1m_nv ...
                    PSVL_F_NI_0d5m_5d1m_nv PSVL_F_I_0d5m_5d1m_nv PSVL_F_NI_10d5m_5d1m_nv PSVL_F_I_10d5m_5d1m_nv ...
                    PSVL_FV_NI_0d5m_5d1m_nv PSVL_FV_I_0d5m_5d1m_nv PSVL_FV_NI_10d5m_5d1m_nv PSVL_FV_I_10d5m_5d1m_nv ...
                    PSVL_OO_NI_0d5m_5d1m_nv PSVL_OO_I_0d5m_5d1m_nv PSVL_OO_NI_10d5m_5d1m_nv PSVL_OO_I_10d5m_5d1m_nv ...
                    PSVL_OOV_NI_0d5m_5d1m_nv PSVL_OOV_I_0d5m_5d1m_nv PSVL_OOV_NI_10d5m_5d1m_nv PSVL_OOV_I_10d5m_5d1m_nv ]) / ...
            mean([  PSV_NF_NI_0d5m_nv PSV_NF_I_0d5m_nv PSV_NF_NI_10d5m_nv PSV_NF_I_10d5m_nv ...
                    PSV_F_NI_0d5m_nv PSV_F_I_0d5m_nv PSV_F_NI_10d5m_nv PSV_F_I_10d5m_nv ...
                    PSV_FV_NI_0d5m_nv PSV_FV_I_0d5m_nv PSV_FV_NI_10d5m_nv PSV_FV_I_10d5m_nv ...
                    PSV_OO_NI_0d5m_nv PSV_OO_I_0d5m_nv  PSV_OO_NI_10d5m_nv PSV_OO_I_10d5m_nv ...
                    PSV_OOV_NI_0d5m_nv PSV_OOV_I_0d5m_nv PSV_OOV_NI_10d5m_nv PSV_OOV_I_10d5m_nv ]) ;    % nv, 5m, incentive

pl(7)	=	mean([	PSVL_NF_NI_10d10m_nv PSVL_NF_I_10d10m_nv ...
                    PSVL_F_NI_10d10m_nv PSVL_F_I_10d10m_nv ...
                    PSVL_FV_NI_10d10m_nv PSVL_FV_I_10d10m_nv ...
                    PSVL_OO_NI_10d10m_nv PSVL_OO_I_10d10m_nv ...
                    PSVL_OOV_NI_10d10m_nv PSVL_OOV_I_10d10m_nv	]) / ...
            mean([	PSV_NF_NI_10d10m_nv PSV_NF_I_10d10m_nv ...
                    PSV_F_NI_10d10m_nv PSV_F_I_10d10m_nv ...
                    PSV_FV_NI_10d10m_nv PSV_FV_I_10d10m_nv ...
                    PSV_OO_NI_10d10m_nv PSV_OO_I_10d10m_nv ...
                    PSV_OOV_NI_10d10m_nv PSV_OOV_I_10d10m_nv	]) ;	% nv, 10m, no incentive
                
pl(8)	=	mean([	PSVL_NF_NI_10d10m_8m_nv PSVL_NF_I_10d10m_8m_nv ...
                    PSVL_F_NI_10d10m_8m_nv PSVL_F_I_10d10m_8m_nv ...
                    PSVL_FV_NI_10d10m_8m_nv PSVL_FV_I_10d10m_8m_nv ...
                    PSVL_OO_NI_10d10m_8m_nv PSVL_OO_I_10d10m_8m_nv ...
                    PSVL_OOV_NI_10d10m_8m_nv PSVL_OOV_I_10d10m_8m_nv	]) / ...
            mean([	PSV_NF_NI_10d10m_nv PSV_NF_I_10d10m_nv ...
                    PSV_F_NI_10d10m_nv PSV_F_I_10d10m_nv ...
                    PSV_FV_NI_10d10m_nv PSV_FV_I_10d10m_nv ...
                    PSV_OO_NI_10d10m_nv PSV_OO_I_10d10m_nv ...
                    PSV_OOV_NI_10d10m_nv PSV_OOV_I_10d10m_nv	]) ; % nv, 10m, incentive


                
% Percent of people who vote

pv(1) = Turnout_control ;	
pv(2) = Turnout_GOTV ;
pv(3) = Turnout_P_control ;	
pv(4) = Turnout_P_GOTV ;

% return the simulated moment or 9999 if NaN
simMoments = min(9999,[ph ; psv ; poo ; psvini ; pl ; pv]);

end
