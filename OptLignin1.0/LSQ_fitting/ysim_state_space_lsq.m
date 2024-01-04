function [ysim,sol]= ysim_state_space_lsq(par, param, obs_data)


IC=[obs_data.Ct_obs(1), obs_data.Co_obs(1)];

sol = ode23s(@(t,y)ode_least_square2(y, par, param), ...
   [0,obs_data.tobs(end)] ,IC, odeset('RelTol',1e-8,'AbsTol',1e-10));

CTsim = sol.y(1,:);
Cosim = sol.y(2,:);

try
   Ctotal = interp1(sol.x, CTsim, obs_data.tobs','pchip');
   Cosim = interp1(sol.x, Cosim, obs_data.tobs_Co','pchip');
catch
   Ctotal = interp1(sol.x, CTsim, obs_data.tobs','pchip');
   Cosim = interp1(sol.x, Cosim, obs_data.tobs','pchip');
end

ysim = [Ctotal./Ctotal(1), Cosim./Cosim(1)]';

end