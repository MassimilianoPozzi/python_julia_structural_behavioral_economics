function [simMoments] = simGneezyLie(parameters)

% warning off MATLAB:quad:MinStepSize

mu=parameters(1);
sigma=parameters(2);
alpha=parameters(3);


if sum(parameters([2])<0)>0
    simMoments=NaN(5,1);
    return;
end


% =============================================================
% Survey Stuff
% =============================================================

shareLieMin1930 = normcdf(-1+alpha*10,mu,sigma);
shareLie2130 = normcdf(1+alpha*10,mu,sigma);
shareLie3030 = normcdf(10+alpha*10,mu,sigma);
shareLie2115 = normcdf(1-alpha*5,mu,sigma);
shareLie3020 = normcdf(10+alpha*0,mu,sigma);

simMoments = [shareLieMin1930 shareLie2130 shareLie3030 shareLie2115 shareLie3020]';


end



