% set votesim file that generates moments
votesimFile = @voteSimEndog_2N
calc_implicationsFile = @calc_implications_2N

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


implications = output;
save(strcat(dataFile,'.mat'));
