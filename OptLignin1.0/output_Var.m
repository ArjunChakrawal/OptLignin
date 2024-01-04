function [ysim, solf]= output_Var(par, ocp, obs_data)
fCO_0 = obs_data.Co_obs(1); % initial fraction of CO in litter[-]
fCH_0 = obs_data.Ct_obs(1) - fCO_0; % initial fraction of CH in litter[-]

ocp.Control.setUpper(par(1))
 ocp.Parameter(1:2).setBound(par(2:4))
ocp.State(1:2).setInitial([fCH_0;fCO_0]) % set initial value of total c ,lignin
ocp.parameterize();
solf = ocp.optimize(... %     'initialGuess', initialGuess, ...
    'ipopt', struct('max_iter', 2000,"print_level",1,"max_cpu_time",50) ...
    );
time = solf.NumericalResults.Independent;
Ch = solf.NumericalResults.State(1, :);
Co = solf.NumericalResults.State(2, :);
try
   Ctotal =  interp1(time, Ch, obs_data.tobs')...
    +interp1(time, Co, obs_data.tobs');
    Cosim = interp1(time, Co, obs_data.tobs_Co');
    ysim = [Ctotal, Cosim ./ Cosim(1)]';
catch
   Ctotal =  interp1(time, Ch, obs_data.tobs')...
    +interp1(time, Co, obs_data.tobs');
    Cosim = interp1(time, Co, obs_data.tobs');
    ysim = [Ctotal, Cosim ./ Cosim(1)]';
end

% disp("final time is " +time(end))
end