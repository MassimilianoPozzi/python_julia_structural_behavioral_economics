clear all

% Use exponential (exponential==1) or power function estimates (exponential==0)?
exponential = 1;
% Use parametric (parametric==1) or non-parametric version (parametric==0)?
parametric = 0;

if exponential == 1
    % exponential NLLS estimates
    s = 3.7225938E-6;
    k = 1.70926702E-16;
    gamma = 0.015641071;
    alpha = 1;
    RMSE = 653.578104;
else
    % power NLLS estimates estimates
    s = 3.16824120E-6;
    k = 5.11863E-70;
    gamma = 20.5461432;
    alpha = 1;
    RMSE = 0.63605161;
end



noRandomDraws = 40*540;
maxEffort = 3500;
effort1000 = 1000;
effort2000 = 2000;
effort3000 = 3000;


effortNoPay = xlsread('dtafiles/buttonpressesNoPay.xls','Sheet1');
effortLong = repmat(effortNoPay,noRandomDraws/540,1);
effortShocks = rand(size(effortLong));
effortShocks(1:size(effortNoPay,1)) = zeros(size(effortNoPay));
effortLong = sort(effortLong + 0.3*(effortShocks-0.5).*effortLong);

effortLowPay = xlsread('dtafiles/buttonpressesLowPay.xls','Sheet1');
effortLowPay = [effortLowPay;mean(effortLowPay);mean(effortLowPay)];
effortLowLong = repmat(effortLowPay,noRandomDraws/540,1);
effortLowShocks = rand(size(effortLowLong));
effortLowShocks(1:size(effortLowPay,1)) = zeros(size(effortLowPay));
effortLowLong = sort(effortLowLong + 0.05*(effortLowShocks-0.5).*effortLowLong);

ki = sort(s./exp(effortLong*gamma), 'descend');
if parametric == 1
    epsilon = RMSE*randn(size(effortLong));
    ki = sort(k*exp(gamma*epsilon), 'descend');
end

if exponential == 0
    effortLong = log(effortLong);
    effortLowLong = log(effortLowLong);
    effort1000 = log(effort1000);
    effort2000 = log(effort2000);
    effort3000 = log(effort3000);
end



effort = 1:maxEffort;
piecerate = zeros(maxEffort,1);
for i = 1:3
    piecerate(effort>=i*1000) = 0.001*i;
end    

%% Calculate marginal benefit (mb) and marginal costs (mc)
mbNoPay = (s+0.00)*effort;
mbPayOnce = (s+0.001)*effort;
mbPayTwice = (s+2*0.001)*effort;
mbPay = (s+piecerate').*effort;
mc = ki*exp(gamma*effort);

% Calculate utility for different scenarios
utilityNoPay = repmat(mbNoPay,noRandomDraws,1) - mc;
utilityPayOnce = repmat(mbPayOnce,noRandomDraws,1) - mc;
utilityPayTwice = repmat(mbPayTwice,noRandomDraws,1) - mc;
utilityPay = repmat(mbPay,noRandomDraws,1) - mc;

% Define indicators whether staying is preferred to jumping to next 1000
indicatorNoPay = ((utilityNoPay - repmat(utilityPay(:,1000),1,maxEffort)>0 & repmat(utilityPay(:,1000),1,maxEffort)>=0) | repmat(utilityPay(:,1000),1,maxEffort)<0);
indicatorPayOnce = ((utilityPayOnce - repmat(utilityPay(:,2000),1,maxEffort)>0 & repmat(utilityPay(:,2000),1,maxEffort)>=0) | repmat(utilityPay(:,2000),1,maxEffort)<0);
indicatorPayTwice = ((utilityPayTwice - repmat(utilityPay(:,3000),1,maxEffort)>0 & repmat(utilityPay(:,3000),1,maxEffort)>=0) | repmat(utilityPay(:,3000),1,maxEffort)<0);

%% Identify cutoff types
% first cutoff
MNoPay = max(indicatorNoPay,[],2); % indicates whether it is preferable for some type not to switch
indexSwitchingPointNoPay = find(1-MNoPay,1)-1; % find index of the switching type
marginalTypeNoPay = ki(indexSwitchingPointNoPay); % return cost scaling factor ki for marginal type
marginalTypeEffortNoPay = (1/gamma)*log(s/ki(indexSwitchingPointNoPay)); % return provided effort of the marginal type

% second cutoff
MPayOnce = max(indicatorPayOnce,[],2); % indicates whether it is preferable for some type not to switch
indexSwitchingPointPayOnce = find(1-MPayOnce,1)-1;
marginalTypePayOnce = ki(indexSwitchingPointPayOnce);
marginalTypeEffortPayOnce = (1/gamma)*log(s/ki(indexSwitchingPointPayOnce));

% third cutoff
MPayTwice = max(indicatorPayTwice,[],2); % indicates whether it is preferable for some type not to switch
indexSwitchingPointPayTwice = find(1-MPayTwice,1)-1;
marginalTypePayTwice = ki(indexSwitchingPointPayTwice);
marginalTypeEffortPayTwice = (1/gamma)*log(s/ki(indexSwitchingPointPayTwice));

% find position of effort levels of 1000, 2000, and 3000
index1000 = find((effortLong>=effort1000),1);
index2000 = find((effortLong>=effort2000),1);
index3000 = find((effortLong>=effort3000),1);

% Simulate effort with 1c/1000 piecerate
simulatedEffort2 = effortLong;
simulatedEffort2(marginalTypeEffortPayTwice<effortLong & effortLong<effort3000) = effort3000;
simulatedEffort2(marginalTypeEffortPayOnce<effortLong & effortLong<effort2000) = effort2000;
simulatedEffort2(marginalTypeEffortNoPay<effortLong & effortLong<effort1000) = effort1000;

% Plot effort observations
hist([effortLong simulatedEffort2  effortLowLong ],50)
g=findobj(gca,'Type','patch');
set(g(1),'FaceColor',[0 0 0],'EdgeColor',[0 0 0])
set(g(2),'FaceColor',[.1 .3 1],'EdgeColor',[.1 .3 1])
set(g(3),'FaceColor',[.6 .6 .6],'EdgeColor',[.6 .6 .6])
title('Simulated Effort in Low Pay Treatment')
line(mean(effortLong), ylim)

fprintf('Cutoff for jumping to 1000 is %6.2f\n', marginalTypeEffortNoPay)
fprintf('Share of jumpers is %6.2f\n',(index1000-indexSwitchingPointNoPay)/noRandomDraws)
fprintf('Cutoff for jumping to 2000 is %6.2f\n', marginalTypeEffortPayOnce)
fprintf('Share of jumpers is %6.2f\n',(index2000-indexSwitchingPointPayOnce)/noRandomDraws)
fprintf('Cutoff for jumping to 3000 is %6.2f\n', marginalTypeEffortPayTwice)
fprintf('Share of jumpers is %6.2f\n',(index3000-indexSwitchingPointPayTwice)/noRandomDraws)

fprintf('The mean effort of NoPay is %6.2f\n', mean(effortLong))
fprintf('The mean effort of SimLowPay is %6.2f\n', mean(simulatedEffort2))
fprintf('The mean effort of LowPay is %6.2f\n', mean(effortLowLong))