function [par, sol, rmse, r2] = find_parameter(obs_data, param, ...
    init_guess, g, fig, makeplot_state_space, ysim_state_space, n)

[ocp, ~] = opt_con(param, g, init_guess, obs_data.tobs(end)*n);
y1norm = obs_data.Ct_obs ./ obs_data.Ct_obs(1);
y2norm = obs_data.Co_obs ./ obs_data.Co_obs(1);
ydata = [y1norm; y2norm];

IC = [param.CT_0, param.CO_0];
fun = @(par, t) ysim_state_space(par, ocp, obs_data, IC);
options1 = optimoptions('lsqcurvefit', 'Display', 'final-detailed', ...
    ... %         'PlotFcn', {@optimplotx, @optimplotfirstorderopt, @optimplotresnorm}, ...
    'MaxFunEvals', 2000, 'MaxIter', 1000, 'TolFun', 1e-6, 'TolX', 1e-6, ...
    'FiniteDifferenceStepSize', 1e-9, 'FiniteDifferenceType', 'central');
try
    lb = param.lb;
    ub = param.ub;
catch
    lb = [0.0001, 0.001, 1]; % [vh_max, mo,ro]
    ub = [0.005, 0.2, 400];
end

problem = createOptimProblem('lsqcurvefit', ...
    'objective', fun, ...
    'xdata', obs_data.tobs, 'ydata', ydata, ...
    'x0', init_guess, ...
    'lb', lb, 'ub', ub, ...
    'options', options1);

% 
[par, ~, ~, ~, ~, ~, ~] = lsqcurvefit(problem); 
% tic;
% ms=MultiStart("UseParallel",false);
% [par,fval,exitflag,outpt,solutions] = run(ms,problem,50);
% toc;
[~, sol] = opt_con(param, g, par, obs_data.tobs(end)*n);

[r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT]= est_r2_rmse(obs_data,sol);


end