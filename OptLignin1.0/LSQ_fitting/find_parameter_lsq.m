function [par, sol] = find_parameter_lsq(obs_data, param, ...
   init_guess, ysim_state_space_lsq)

y1norm = obs_data.Ct_obs ./ obs_data.Ct_obs(1);
y2norm = obs_data.Co_obs ./ obs_data.Co_obs(1);
ydata = [y1norm; y2norm];

fun = @(par, t) ysim_state_space_lsq(par, param, obs_data);
options1 = optimoptions('lsqcurvefit', 'Display', 'final-detailed', ...
   ... %         'PlotFcn', {@optimplotx, @optimplotfirstorderopt, @optimplotresnorm}, ...
   'MaxFunEvals', 2000, 'MaxIter', 1000, 'TolFun', 1e-6, 'TolX', 1e-6, ...
   'FiniteDifferenceStepSize', 1e-6, 'FiniteDifferenceType', 'forward');
try
   lb = param.lb;
   ub = param.ub;
catch
   lb = [0.0001,0.0001, param.mo, 0.01]; % [vh_max, mo,ro]
   ub = [0.04,0.04, param.mo, 400];
end

problem = createOptimProblem('lsqcurvefit', ...
   'objective', fun, ...
   'xdata', obs_data.tobs, 'ydata', ydata, ...
   'x0', init_guess, ...
   'lb', lb, 'ub', ub, ...
   'options', options1);

par = lsqcurvefit(problem);
[ysim,sol]= ysim_state_space_lsq(par, param, obs_data);
% err=ysim-ydata;
% rmse=sqrt(mean(err.^2));
% % rmseN=sqrt(mean(err.^2))/mean(ydata);
% dev= mean(ydata) -ydata;
% SST = sum(dev.^2);
% SSR = sum(err.^2);
% r2 = 1 - SSR / SST;

end