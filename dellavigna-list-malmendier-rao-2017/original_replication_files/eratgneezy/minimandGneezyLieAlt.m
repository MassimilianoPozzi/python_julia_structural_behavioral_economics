function  [y] = minimandGneezyLieAlt(parameters)    
VCcontrol=diag([0.0022 0.0025 0.0022 0.0022 0.0023]);
moments = [0.33 0.49 0.65 0.37  0.5229]';
% VCcontrol=diag([0.0022 0.0025 0.0022 0.0022]);
% moments = [0.33 0.49 0.65 0.37  ]';
W = inv(VCcontrol);
mSimOpt=simGneezyLieAlt(parameters);
mTemp=moments-mSimOpt;
y=mTemp'*W*mTemp;

end