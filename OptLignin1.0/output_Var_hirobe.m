function [ydata, solf]= output_Var_hirobe(par, ocp, obs_data, IC)

ocp.Control.setUpper(par(1))
 ocp.Parameter(1:3).setBound(par(2:4))
ocp.State(1:2).setInitial([1-IC;IC]) % set initial value of total c ,lignin
ocp.parameterize();
solf = ocp.optimize(... %     'initialGuess', initialGuess, ...
    'ipopt', struct('max_iter', 2000,"print_level",1,"max_cpu_time",50) ...
    );
tobs=obs_data.tobs;
tobs_Co=obs_data.tobs_Co;
time = solf.NumericalResults.Independent;
Ch = solf.NumericalResults.State(1, :);
Co = solf.NumericalResults.State(2, :);
time_normalized=time;
Ctotal =  interp1(time_normalized, Ch, tobs')...
    +interp1(time_normalized, Co, tobs');

Cosim = interp1(time_normalized, Co, tobs_Co');
ydata = [Ctotal, Cosim ./ Cosim(1)]';
% disp("final time is " +time(end))
end