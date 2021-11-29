%Find the minimum distance estimate in the benchmark model

function [theta,fval, converged, sse] = minSearch(moments,W,parsinit,paramsToUse, MatchedParams, ToMatchParams,  momentsToUse, varargin)

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

%Search space
lowerBound= [0.001 0.001 0.001 -999 0.001 0.001 0.001 0.001 0.001 0.001 -999 0.001 0.001 0.001 0.001 0.001 0.001 -999 0.001 0.001 0.001 0.001 0.001 0.001 -999 0.001 0.001 0.001 -999 -999 0.001 0.001 -1 0.001   -999 -999 0.001 0.001 -1 0.001 -999 -999 0.001 0.001 -1 0.001 -999 -999 0.001 0.001 -1 0.001];
upperBound=[1 1 999 999 999 999 999 1 1 10 999 999 999 999 1 1 10 999 999 999 999 1 1 10 999 999 999 999 999 999 999 999 1 999 999 999 999 999 1 999 999 999 999 999 1 999 999 999 999 999 1 999];

noOfOptimizationParams=length(paramsToUse);  % no of params that the algorithm will optimize over

%also subtract off the parameters that are fixed between the genders
for i=1:length(MatchedParams)
  %make sure the pegged parameter is not separately fixed
  if ismember(ToMatchParams(i),paramsToUse)
    noOfOptimizationParams = noOfOptimizationParams-1;
  end
end
%for i=17:31
%    if ismember(i,paramsToUse)==1 & ismember(i-16,genderlessParams)==1
%        noOfOptimizationParams = noOfOptimizationParams-1;
%    end
%end

paramsToOptimize=-9999*ones(noOfOptimizationParams,1);
lb=-9999*ones(noOfOptimizationParams,1);
ub=9999*ones(noOfOptimizationParams,1);

j=0;
for i=1:noOfParams
    if ismember(i,paramsToUse)==1 & ismember(i,ToMatchParams)==0
        j=j+1;
        paramsToOptimize(j)=parsinit(i);
        lb(j)=lowerBound(i);
        ub(j)=upperBound(i);
    end
end


%Set options for search   
 options=optimset('MaxFunEvals',20000,'MaxIter',10000, 'Display', 'none', 'LargeScale','off', 'TolFun', TolFun, 'TolX', TolX);
 options.Algorithm='interior-point';
 options.AlwaysHonorConstraints='bounds';

%find minimum distance estimate
 converged=-9999;
[thetaTemp,fval, exitFlag] = fmincon(@minimand, paramsToOptimize,[],[],[],[],lb,ub,[],options);
if (exitFlag <= 0) 
    converged=0;
end

if (exitFlag > 0) 
    converged=1;
end



theta=-99999*ones(noOfParams,1);

%set fixed params to the fixed values
m=0;
for n=1:noOfParams
    if ismember(n,paramsToUse)==0 | ismember(n,ToMatchParams)==1 
        theta(n)=parsinit(n);
    else
        m=m+1;
        theta(n)=thetaTemp(m);
    end
end
%and set matched params to the matched values
for n=1:length(MatchedParams)
  theta(ToMatchParams(n))=theta(MatchedParams(n));
end

%report the estimate and SSE
mSimOpt=voteSimWrapper(theta);
mTemp=moments-mSimOpt;
mTemp = mTemp(momentsToUse);
sse=mTemp'*W*mTemp;



%This is the objective function to minimize in fmincon above
function  [y] = minimand(paramsToOptimizeArg)    

j=0;
parameters=-9999*ones(noOfParams,1);
for ii=1:noOfParams
    if ismember(ii,paramsToUse)==0 || ismember(ii,ToMatchParams)==1
        parameters(ii)=parsinit(ii);
    else 
        j=j+1;
        parameters(ii)=paramsToOptimizeArg(j);
    end
end
for ii=1:length(MatchedParams)
  parameters(ToMatchParams(ii))=parameters(MatchedParams(ii));
end

mSim=voteSimWrapper(parameters);
mTemp=moments-mSim;
m=mTemp(momentsToUse);

y=m'*W*m;

end

end

