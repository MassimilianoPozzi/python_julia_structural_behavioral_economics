clc
clear all

% W=eye(4);
VCcontrol=diag([0.0022 0.0025 0.0022 0.0022 0.0023]);
moments = [0.33 0.49 0.65 0.37 0.5229]';
% VCcontrol=diag([0.0022 0.0025 0.0022 0.0022 ]);
% moments = [0.33 0.49 0.65 0.37 ]';
W = inv(VCcontrol);
noOfParams=3;
lowerBound= [-10 -999 0.001];
upperBound= [10 999 10000];
inits=[0.3 7 20]';


%Set options for search   
TolFun=1e-14;
TolX=1e-10; 
options=optimset('MaxFunEvals',20000,'MaxIter',10000, 'Display', 'none', 'LargeScale','off', 'TolFun', TolFun, 'TolX', TolX);
options.Algorithm='interior-point';
options.AlwaysHonorConstraints='bounds';

%find minimum distance estimate
converged=-9999;
[theta,fval, exitFlag] = fmincon(@minimandGneezyLieAlt, inits,[],[],[],[],lowerBound,upperBound,[],options);

if (exitFlag <= 0) 
    converged=0;
end

if (exitFlag > 0) 
    converged=1;
end

% % 
invVC=inv(VCcontrol);
[DFDY, FAC] =jacobianest(@simGneezyLieAlt,theta);
A=DFDY'*W*DFDY;
B=DFDY'*W*VCcontrol*W*DFDY;
VC = A\B/A;
ses=diag(VC).^0.5;

%report the estimate and SSE
mSimOpt=simGneezyLieAlt(theta);
mTemp=moments-mSimOpt;
sse=mTemp'*W*mTemp;


disp(' ### ESTIMATES AND STANDARD ERRORS ###');
disp('=============================================================================================================================');
disp('    Estimates  SES');    
disp(' [          diag(VC)         ');
disp([ theta ses]);
%disp([ bestEstimate']);
disp('============================================================================================================================');
[moments mSimOpt]
sse



