function [r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol)
y1norm = obs_data.Ct_obs ./ obs_data.Ct_obs(1);
y2norm = obs_data.Co_obs ./ obs_data.Co_obs(1);
ydata = [y1norm; y2norm];

try
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    Co = sol.NumericalResults.State(2, :);
    p=2;
catch 
   time = sol.x;
   CT = sol.y(1, :);
   Co = sol.y(2, :);
   p=3;
end

try
    Ctotal = interp1(time, CT, obs_data.tobs','pchip');
    Cosim = interp1(time, Co, obs_data.tobs_Co','pchip');
catch
    Ctotal = interp1(time, CT, obs_data.tobs','pchip');
    Cosim = interp1(time, Co, obs_data.tobs','pchip');
end
ysim = [Ctotal./Ctotal(1), Cosim./Cosim(1)]';

CtotalNorm = Ctotal'./Ctotal(1);
CosimNorm= Cosim'./Cosim(1);
[r2, rmse,AIC] = r2_rmse_aic(ydata, ysim,p);
[r2_co, rmse_co,AIC_Co]=r2_rmse_aic(y2norm, CosimNorm,p);
[r2_CT, rmse_CT,AIC_CT]=r2_rmse_aic(y1norm, CtotalNorm,p);

% p
end
