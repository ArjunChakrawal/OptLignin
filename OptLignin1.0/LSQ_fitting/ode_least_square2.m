function dydt = ode_least_square2(y, est_par,param)
vhmax= est_par(1);
vomax = est_par(2);
mo= est_par(3);
ro= est_par(4);
CT = y(1);
Co = y(2);
a=param.a;b=param.b;g=param.g;
L = Co /CT;
Dh = vhmax * g(L, a, b) * (CT-Co);
Do = vomax * Co;
e = param.emax - ro * vomax;
G = e * (Dh + Do);
dydt = [-Dh-Do +  G; ...
   -Do + mo * G];
end