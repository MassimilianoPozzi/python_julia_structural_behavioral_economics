function [simMoments] = simGneezyLieAlt(parameters)

% warning off MATLAB:quad:MinStepSize

a=parameters(1);
L=parameters(2);
sigma_eps = parameters(3);


if sum(parameters([3])<0)>0
    simMoments=NaN(5,1);
    return;
end


% =============================================================
% Survey Stuff
% =============================================================

shareLie1930 = normcdf(-1+10*a-L,0,sigma_eps);
shareLie2130 = normcdf(+1+10*a-L,0,sigma_eps);
shareLie3030 = normcdf(10+10*a-L,0,sigma_eps);
shareLie2115 = normcdf(+1-5*a-L,0,sigma_eps);
shareLie3020 = normcdf(10+0*a-L,0,sigma_eps);

simMoments = [shareLie1930 shareLie2130 shareLie3030 shareLie2115 shareLie3020]';

% simMoments = [shareLieMin1930 shareLie2130 shareLie3030 shareLie2115 ]';

end



