function [ysim, solf] = ysim_state_space(par, ocp, obs_data, IC)
% ocp.Control.setUpper(par(1))
ocp.Parameter(1:3).setBound(par) % fixes parameter value of vhmax by making lower=upper
ocp.State(1:2).setInitial(IC) % set initial value of total c ,lignin
ocp.parameterize();
solf = ocp.optimize( ... %     'initialGuess', initialGuess, ...
    'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
    );
time = solf.NumericalResults.Independent;
CT = solf.NumericalResults.State(1, :);
Co = solf.NumericalResults.State(2, :);
try
    Ctotal = interp1(time, CT, obs_data.tobs','pchip');
    Cosim = interp1(time, Co, obs_data.tobs_Co','pchip');
catch
    Ctotal = interp1(time, CT, obs_data.tobs','pchip');
    Cosim = interp1(time, Co, obs_data.tobs','pchip');
end

ysim = [Ctotal./Ctotal(1), Cosim./Cosim(1)]';
end