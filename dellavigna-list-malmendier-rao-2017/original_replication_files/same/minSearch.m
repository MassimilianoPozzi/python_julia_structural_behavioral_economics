% =================================================================
% Find the minimum distance estimator
% This code uses the non-derivative based fminsearchbnd routine

% Inputs: moments, weighting matrix, start point, paramters to minimize 
%   over, master and copy params to set equal, moments to use for objective 
%   function, number of simulated potential voters, parameter search space, 
%   L-setting (cons or dist), fixed random components, votesim file name, 
%   and tolerance levels
% Outputs: solution (parameter vals), weighted SSE for minimzation, 
%   indicator for whether search converged, total weighted sse, number of
%   iterations, and number of function evaluations
% =================================================================

function [theta,fval, converged, sse, iter, funcCount] = minSearch(moments,W,parsinit,...
    paramsToUse, MasterParams, CopyParams, momentsToUse, sim_voters, ...
    lowerBound, upperBound, L_setting, rand_set, votesimFile, varargin)

noOfParams=size(parsinit,1);

if size(varargin,2)==2 
    TolFun=varargin{1};
    TolX=varargin{2};
elseif size(varargin,2)==1 
    TolFun=varargin{1};
    TolX=1e-10;
else 
    TolFun=1e-14;
    TolX=1e-10; 
end

noOfOptimizationParams=length(paramsToUse);  % no of params that the algorithm will optimize over

%also subtract off the parameters that are fixed
for i=1:length(MasterParams)
  %make sure the pegged parameter is not separately fixed
  if ismember(CopyParams(i),paramsToUse)
    noOfOptimizationParams = noOfOptimizationParams-1;
  end
end

paramsToOptimize=-9999*ones(noOfOptimizationParams,1);
lb=-9999*ones(noOfOptimizationParams,1);
ub=9999*ones(noOfOptimizationParams,1);

% NOTE: make sure bounds are consistent for matched parameters?
j=0;
for i=1:noOfParams
    if ismember(i,paramsToUse)==1 & ismember(i,CopyParams)==0
        j=j+1;
        paramsToOptimize(j)=parsinit(i);
        lb(j)=lowerBound(i);
        ub(j)=upperBound(i);
    end
end


%Set options for search   

%find minimum distance estimate
converged=-9999;

% FMINSEARCHBND
options=optimset('Display','none', 'MaxFunEvals',20000,'MaxIter',10000, ...
    'LargeScale','off', 'TolFun', TolFun, 'TolX', TolX);
%note TolX applies to transformed variable

[thetaTemp,fval, exitFlag,output] = fminsearchbnd(@(p)minimand(p,L_setting,rand_set), ...
        paramsToOptimize,lb,ub,options);
iter = output.iterations;
funcCount = output.funcCount;

% NOTE: identify whether converged
if (exitFlag <= 0) 
    converged=0;
end

if (exitFlag > 0) 
    converged=1;
end


% NOTE: generate theta - update paramsToOptimize to add back the 
% fixed values and making sure that the copied values are correct
theta=-99999*ones(noOfParams,1);

%set fixed params to the fixed values
m=0;
for n=1:noOfParams
    if ismember(n,paramsToUse)==0 | ismember(n,CopyParams)==1 
        theta(n)=parsinit(n);
    else
        m=m+1;
        theta(n)=thetaTemp(m);
    end
end
%and set Copy params to the Master values
for n=1:length(MasterParams)
  theta(CopyParams(n))=theta(MasterParams(n));
end


%report the estimate and weighted SSE
% caculate moments using the votesim file
mSimOpt=votesimFile(theta,rand_set);
mTemp=moments(momentsToUse)-mSimOpt(momentsToUse);
sse=mTemp'*W*mTemp;



%This is the objective function to minimize in fmincon above
function  [y] = minimand(paramsToOptimizeArg,L_setting,rand_set)    

    j=0;
    parameters=-9999*ones(noOfParams,1);
    for ii=1:noOfParams
        if ismember(ii,paramsToUse)==0 || ismember(ii,CopyParams)==1
            parameters(ii)=parsinit(ii);
        else 
            j=j+1;
            parameters(ii)=paramsToOptimizeArg(j);
        end
    end
    for ii=1:length(MasterParams)
      parameters(CopyParams(ii))=parameters(MasterParams(ii));
    end


    % calculate simulated moments
    mSim=votesimFile(parameters,rand_set);
    m=moments(momentsToUse)-mSim(momentsToUse);
    
    y=m'*W*m;

    end

end

