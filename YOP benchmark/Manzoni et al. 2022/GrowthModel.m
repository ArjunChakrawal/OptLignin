function [dx, y] = GrowthModel(time, state, control,param)
% global param

beta=param(1);
emax =param(2);
r =param(3);
mu=param(4);
gamma = param(5);
k=param(6);

Cs = state;
uptake = control;
dx =  -uptake - Cs * gamma + mu*(emax.*beta.*(uptake - r)/(uptake + beta)-k*uptake) ;

y.Cs = Cs;
y.uptake = uptake;

end