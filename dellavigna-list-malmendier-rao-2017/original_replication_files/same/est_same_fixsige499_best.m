% =================================================================
% Estimate the model
% Same parameters for V/NV, fix sigma epsilon to 499.9
% Use the best estimate from the benchmark model (est_same_bench_best)
% as the start point
% =================================================================

clc; clear; 
warning off all

curDir = pwd;
mkdir('../output','same');
outputDir=('../output/same/');

dataFile = strcat('same-fixsige499-best-',datestr(now,'dd-mmm-yyyy-HH_MM'));
outputFile = strcat(dataFile,'.txt');
diary(strcat(outputDir,outputFile));

% set votesim file that generates moments
votesimFile = @voteSimEndogenousVoting
calc_implicationsFile = @calc_implications

% choose L: cons (fixed L across people) or dist (L_lambda expon dist)
L_setting = 'cons'

% Choose L value if constant
cons_fixedL = 10;

% Choose mu and sigma epsilon (for simulation or fixed value)
fixed_mu_eps = 50;
fixed_sigma_eps = 499.9;

% Estimation options:
% Number of Simulated Voters in each simulation
sim_voters=500000
%Number of search start points
% NOTE: this will pull the top estimate from the results of the
% benchmark model estimated using as start points the top 20 start points 
% only use noSearchInits=1 or else fixed sigma-epsilon is not fixed
noSearchInits=1
search_init_file = 'bench_best'

% Empirical moments and var-cov matrix: 1 x 100 moments and matrix 100 x 100
% Baseline turnout and GOTV added to Empirical Moments
load '../Moments.mat';    
format short g
Moments_org = Moments;
VCcontrol_org = VCcontrol;
Moments = [Moments_org; 0.6000; 0.6102]; 
VCcontrol = blkdiag(VCcontrol_org,0.0109^2,0.01376^2);

%Different estimations will use different sets of moments, defined here
allMoments=1:101;
noTurnoutMoments = 1:100;
noLyingMoments = [1:92 101];
onlyLyingMoments = 93:100;
baseLyingMoments = [1:92 93:2:99 101]; % no lying incentives
voterMoments=[1:15 31:45 61:66 73:82 93:96]; % only voter moments
nonvoterMoments=[16:30 46:60 67:72 83:92 97:100]; % only non-voter moments
C_GOTVMoments=1:102; % congressional GOTV
max_empirical=1:102;

noIncentLyingMoments = [1:93 95 97 99 101];
noINIMoments = [1:72 93:101];

% includes out of sample moments in the predictions
noOfMoments = 104;

%Different estimations estimate different sets of parameters, defined here.
paramLabels={'h0', 'r', 'eta', 'mu_s', 'sigma_s', 'S_svy', 'timeval', ...
'mu_sv','mu_sn', 'sigma_sv', 'sigma_sn', 'rho', 'L', 'mu_eps', 'sigma_eps'};
      
noOfParams=15;
allParams = 1:noOfParams;
fixedL_Rho=[1:11 14:15]; % NOTE: fix L and rho (exclude)
fixedRho_SigEps=[1:11 13:14]; % NOTE: fix rho and sig_eps
fixedL_RhoSigEps=[1:11 14]; % NOTE: fix L, rho, and sigma_eps
fixed_LRhoEps=[1:11]; % NOTE: fix L, rho, mu_eps and sigma_eps
fixedRho=[1:11 13:15]; 

% changed language to "master" and "copy".
% this code sets sigma_sn = sigma_sv. 

%We also want to fix certain parameters to be equal to each other.
% Copy is forced equal to Matched (i.e. Matched is the "master" param)
%MasterParams must be in ascending order! (unless it's a particularly simple case...)
%Master Params are the lower index parameters that later parameters are set equal to.
Master_sigmas = [10]; % sigma_sv
Copy_sigmas = [11]; %sigma_sn

%set tolerance parameters for numerical minimum distance estimation
% Y and X tolerance
TolFun=1e-10;
TolX=1e-10;
maxSSE=1e-9;


% =========================================================================
% *************************************************************************
% =========================================================================

% Sets value to use when parameters are designated as fixed
if L_setting == 'cons' 
svyEstimates  = [0.37, 0.34, 0.1, -28, 30, 1, 42, -4, -13, 15, 15, 0, cons_fixedL, fixed_mu_eps, fixed_sigma_eps]';
end
if L_setting == 'dist' 
svyEstimates  = [0.37, 0.34, 0.1, -28, 30, 1, 42, -4, -13, 15, 15, 0, 0.05, fixed_mu_eps, fixed_sigma_eps]';
end

paramsToSim=[  svyEstimates ];

%Setting estimation parameters  
% Note: this code does not allow varying rho (must be fixed)
paramsToUse=fixedRho_SigEps;  % OR: fixedL_Rho fixedRho
display 'fixedRho_SigEps';
momentsToUse=allMoments;  %allMoments C_GOTVMoments
display 'allMoments';

fixedParams=paramsToSim; % fixing L and rho to the paramstoSim file (starting point)
MasterParams=Master_sigmas ;
CopyParams=Copy_sigmas ;

% weight moments by the inv variance (VCcontrol matrix)
W=inv(diag(diag(VCcontrol(momentsToUse,momentsToUse))));
% W=eye(length(momentsToUse));
display 'W: inv VC';

% shuffle to new random setting if desired - currently set to same seed 
% before rand set and get search inits 
%rng shuffle; rng

%Draw random set for simulations
% order of rand_set: s eps sv sn L
% NOTE: This code forces rho=0
rng(686)
rho	=	0;
SIGMA_svsn = [1 rho ; rho 1]; 
SIGMA_rand = blkdiag(eye(2),SIGMA_svsn);
rand_set = mvnrnd(zeros(1,4),SIGMA_rand,sim_voters);
rand_set = [rand_set,rand(sim_voters,1)];


% =========================================================================
% *************************************************************************
% =========================================================================

%Correct the paramsToUse vector so that any matched parameters to fixed
% params are also treated as unused
%(This should still allow for fixing one parameter to another fixed one
%even if they aren't the same values in the fixedParams vector. Ie, 
%paramsToUse is established first, and then MasterParams is applied.)
%Again, be careful that paramsToUse and MasterParams and CopyParams
%lead to the expected behavior...
changedParams=1;
while (changedParams~=0)
   changedParams=0;
   for i=1:length(MasterParams)
      if ismember(MasterParams(i),paramsToUse)==0 && ismember(CopyParams(i),paramsToUse)==1
         index=find(paramsToUse==CopyParams(i));
         paramsToUse(index)=[];
	 changedParams=1;
      end
   end
end

%really minimal error checking...
if (length(MasterParams) ~= length(CopyParams))
    error('Master parameter vectors must be same length');
end


%Names of fixed parameters
fixedParamNames='';
for kk=1:noOfParams 
    if(ismember(kk,paramsToUse)==0)
        fixedParamNames=horzcat(fixedParamNames, sprintf('\t') ,paramLabels{kk}); 
    end
end

%Names of parameters that are set equal
MasterParamNames='';
for kk=1:noOfParams 
    if(ismember(kk,MasterParams)==1)
        MasterParamNames=horzcat(MasterParamNames, sprintf('\t') ,paramLabels{kk}); 
    end
end
CopyParamNames='';
for kk=1:noOfParams 
    if(ismember(kk,CopyParams)==1)
        CopyParamNames=horzcat(CopyParamNames, sprintf('\t') ,paramLabels{kk}); 
    end
end

%Build matrix of initial search values in the parameter space
% Fix random seed, get back start space bounds for display
rng(747)
[searchInits,lower_start,upper_start] =  getSearchInits_best(noSearchInits,...
    search_init_file,sim_voters,rand_set,votesimFile);
for j=1:noOfParams
    if (ismember(j,paramsToUse)==0)
        searchInits(j,:)=fixedParams(j);
    end
end
for j=1:length(MasterParams)
    if (CopyParams(j)<=MasterParams(j))
	error('Master Parameters must come before the parameters that are matched to them.');
    end
    searchInits(CopyParams(j),:)=searchInits(MasterParams(j),:);
end


%These will store the results of the minimum distance estimations
paramHats=zeros(noOfParams,noSearchInits);
fval= repmat(-9999,1,noSearchInits); % Value of the min dist (weighted SSE)
sse= repmat(-9999,1,noSearchInits);
converged=repmat(-9999,1,noSearchInits); % indicates whether search converged

%Define the search space - depends on L or L_lambda
if L_setting == 'cons' 
lowerBound=[eps, eps, eps, -100, eps, eps, eps, -100, -100, eps, eps, -1, eps, -500, eps];
upperBound=[1,    1,  0.5,  100, 100, 100, 200,  100,  100, 100, 100,  1,  50,  500, 2000];
end
if L_setting == 'dist' 
lowerBound=[eps, eps, eps, -100, eps, eps, eps, -100, -100, eps, eps, -1, eps, -500, eps];
upperBound=[ 1,    1, 0.5,  100, 100, 100, 200,  100,  100, 100, 100,  1,  1,   500, 2000];
end

disp('=================================================================');
disp(' [Lower_Start  Upper_Start Lower_Bound Upper_Bound]');
disp([ lower_start  upper_start lowerBound'  upperBound']);
disp('=================================================================');

% Make sure all search inits are in bounds (should only be relevant when
% using getsearchinits_best)
for j=1:noSearchInits
searchInits(:,j) = min(max(searchInits(:,j),lowerBound'),upperBound');
end

% MINSEARCH FOR NEW PARAMETER BOUNDS AND TO CALL NEW ENDOG VOTING
% SIMULATED MOMENTS
% Note: This loop performs the actual estimations. Passes minsearch the
% empirical moments (moments.mat), weighting matrix (only for moments
% used), search initialization points... Returns estimates, min value 0/1
% indication for convergence, weighted sse
% Uses the fminsearch routine allowing for bounds on the parameter space
iter = 0;
funcCount= 0;

parfor j=1:noSearchInits
    if mod(j,10)==0
        disp(j);
    end
   [paramHats(:,j), fval(j), converged(j), sse(j), temp_iter, temp_funcCount] ...
       = minSearch(Moments,W,...
       searchInits(:,j), paramsToUse, MasterParams, CopyParams, ...
       momentsToUse, sim_voters, lowerBound, upperBound, ...
       L_setting, rand_set, votesimFile, TolFun, TolX); 
   iter = iter + temp_iter;
   funcCount = funcCount + temp_funcCount;
end

searchInits=searchInits';
paramHats=paramHats';

%display some information about the estimates
%for converged moments, want to keep even out of sample (use noOfMoments)
noOfAcceptableSolutions=sum(converged); % converged estimates
disp('Converged solutions:  '); disp(noOfAcceptableSolutions);
convergedParamHats=repmat(-9999,noOfAcceptableSolutions,noOfParams);
convergedInits=repmat(-9999,noOfAcceptableSolutions,noOfParams);
convergedSse=repmat(-9999,noOfAcceptableSolutions,1);
convergedMoments=repmat(-9999,noOfAcceptableSolutions,noOfMoments);

%We only care about the estimates that converged
j=0;
for i=1:noSearchInits
    if converged(i)==1 
         j=j+1;
         convergedParamHats(j,:)=paramHats(i,:);
         convergedInits(j,:)=searchInits(i,:);
         convergedSse(j)=sse(i);
         convergedMoments(j,:)=votesimFile(convergedParamHats(j,:)',rand_set)';
    end
end
    
sseParams=[convergedSse convergedParamHats convergedInits convergedMoments];
sseParamsSorted=sortrows(sseParams); % note: this is weighted sse
convergedSse=sseParamsSorted(:,1);
convergedParamHats=sseParamsSorted(:,2:noOfParams+1);    
convergedInits=sseParamsSorted(:,noOfParams+2:noOfParams+1+noOfParams);
convergedMoments=sseParamsSorted(:,noOfParams+noOfParams+2:end);


%Number of estimates within the desired SSE range
noSmallSSE=sum(convergedSse<maxSSE);

%Estimate with the smallest SSE
bestEstimate=convergedParamHats(1,:);
predictedMoments=convergedMoments(1,:)';
m=Moments(momentsToUse)-predictedMoments(momentsToUse);
unwt_sse=m'*m;

%Implications for best estimate
implications=calc_implicationsFile(bestEstimate,sim_voters,L_setting,rand_set);

%save data
save(strcat(outputDir,dataFile,'-prelim.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate standard errors using jacobianest 
% (adjustment for empirical N / simulation N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% empirical N (13197) / simulation N
Jm_Js = 13197 / sim_voters;

%JACBOIANEST 
[DFDY_CSD_jacest, error_jacest] = ...
    jacobianest(@(parameters)votesimFile(parameters,rand_set),...
    bestEstimate);

% Note: here we need to take account for the fact that some params move
% together
%Master parameters are effectively the same, so the total derivative wrt
%that parameter should be the sum of the partials wrt each one
%(need to use while loop in case there is a chain of matched parameters)
changed=1;
while changed~=0
   changed=0;
   for i=1:length(MasterParams)
      if sum(DFDY_CSD_jacest(:,CopyParams(i))) ~= 0
         changed=1;
         totalderiv = DFDY_CSD_jacest(:,MasterParams(i))+DFDY_CSD_jacest(:,CopyParams(i));
         DFDY_CSD_jacest(:,MasterParams(i))=totalderiv;
         DFDY_CSD_jacest(:,CopyParams(i))=0;
      end
   end
end


%weed out the parameters that are duplicates so the matrix A will be invertible...
paramsToEstimate=paramsToUse;
for i=1:length(CopyParams)
   index=find(paramsToEstimate==CopyParams(i));
   paramsToEstimate(index)=[];
end

%JACBOIANEST
%this should be correct with the corrections to paramsToUse/Estimate above
% simulation adjustment
sim_adjust = 1 + Jm_Js;
DFDY_CSD_jacest =-DFDY_CSD_jacest(momentsToUse,paramsToEstimate);
A=DFDY_CSD_jacest'*W*DFDY_CSD_jacest;
B=DFDY_CSD_jacest'*W*(sim_adjust.*VCcontrol(momentsToUse,momentsToUse))*W*DFDY_CSD_jacest;
VC = A\B/A;
sesTemp=diag(VC).^0.5;
ses_jacest = zeros(noOfParams,1);
for i=1:length(paramsToEstimate)
    ses_jacest(paramsToEstimate(i))=sesTemp(i);
end

% Display best estimates and jacobianest SE results
disp('=================================================================');
disp(' [ Best_Estimate  ses_jacest ]');
disp([ bestEstimate'  ses_jacest ]);
disp('=================================================================');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Errors for Implications
fun = @(x) calc_implicationsFile(x,sim_voters,L_setting,rand_set);
output=fun(bestEstimate');
DFDY_sigVal=jacobianest(fun,bestEstimate');
DFDY_sigVal =-DFDY_sigVal(:,paramsToEstimate);
VC_imp=DFDY_sigVal*VC*DFDY_sigVal';
SES_imp=diag(VC_imp).^(0.5);
disp('================================================================');
disp('Implications Output')
disp([output SES_imp]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Display
m=Moments(momentsToUse)-predictedMoments(momentsToUse);
wsse_used=m'*W*m;

m_all=Moments(allMoments)-predictedMoments(allMoments);
W_all=inv(diag(diag(VCcontrol(allMoments,allMoments))));
wsse_all=m_all'*W_all*m_all;


% SHOW ALL MOMENTS, BEST PARAMETERS AND ALL ESTIMATED PARAMETERS
disp('     Emp_Moments  Simulated_Moments ');
disp([ Moments(max_empirical) predictedMoments(max_empirical)]);
disp('=================================================================');
fprintf(' WSSE_used = %g   WSSE_all_noGOTV = %g  \n', wsse_used, wsse_all );          
disp('=================================================================');

fprintf('\n\n\n');
fprintf('DETAILED RESULTS:\n\n');
disp(sprintf('%d out of %d estimates converged \n',noOfAcceptableSolutions,noSearchInits));
disp('[Index sseMoments] and [Estimated Parameters] ='); disp([ [1:noOfAcceptableSolutions]' convergedSse convergedParamHats]);

%transpose for transfer to excel
excel_bounds = [lower_start upper_start lowerBound' upperBound'];
excel_ParamHats = convergedParamHats';
excel_Moments = convergedMoments';

save(strcat(outputDir,dataFile,'.mat'));

p = gcp;
delete(p)

diary off
