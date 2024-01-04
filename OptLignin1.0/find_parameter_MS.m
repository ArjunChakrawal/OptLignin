function [par,sol, rmse,rsquare] = find_parameter_MS(ocp, obs_data,param, ...
    init_guess,g,fig,makeplot_state_space,output_var_)
ydata = [obs_data.Ct_obs./obs_data.Ct_obs(1);...
obs_data.Co_obs./obs_data.Co_obs(1)];
fun = @(par,t) output_var_(par, ocp, obs_data, [param.CT_0, param.CO_0]);
options1 = optimoptions('lsqcurvefit', 'Display', 'final-detailed', ...
    ...%     'PlotFcn', {@optimplotx, @optimplotfirstorderopt, @optimplotresnorm}, ...
    'MaxFunEvals', 2000, 'MaxIter', 1000, 'TolFun', 1e-6, 'TolX', 1e-6, ...
    'FiniteDifferenceStepSize', 1e-9, 'FiniteDifferenceType', 'forward');
lb = [0.0001, 0.001,1]; % [vh_max, mo,ro]
ub = [0.01, 0.15,200];
problem = createOptimProblem('lsqcurvefit',...
    'objective',fun,...
    'xdata',obs_data.tobs,'ydata',ydata,...
    'x0',init_guess,...
    'lb',lb,'ub',ub,...
    'options',options1);
resn=[];parest=[];
fac=[0.01,0.1,2,20];
for i =1:4
    problem.x0(1)=init_guess(1)*fac(i);
    [par, resnorm,~, ~, ~, ~, ~] = lsqcurvefit(problem);
    [ysim, sol] = fun(par, obs_data.tobs);
    rmse = sqrt(mean(ysim-ydata).^2);
    resn=[resn;rmse];
    parest=[parest;par];
end
par=parest(resn==min(resn),:);
% ms=MultiStart("UseParallel",true);
% [par,fval,exitflag,outpt,solutions] = run(ms,problem,5);

[ysim, sol] = fun(par, obs_data.tobs);
makeplot_state_space(sol, '',param,obs_data,g,fig)
rmse = sqrt(mean(ysim-ydata).^2);
SSR = 0;
SST = sum((mean([obs_data.Ct_obs; obs_data.Co_obs])...
    -[obs_data.Ct_obs; obs_data.Co_obs]).^2);
rsquare= 1-SSR/SST;
disp([par,rmse])
end