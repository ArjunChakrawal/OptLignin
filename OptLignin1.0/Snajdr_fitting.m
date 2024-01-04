clear all
clc
close all
%%

set_up
study="Snajdr2011";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";

fig_pathCombined = "fig\"+scenario_name+"\combined\";
fig_vo=figure;tiledlayout('flow');
Ligdecomposition_Starts = zeros(1, 1);
C_remain_lig_dec_start=Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo= avgvo;
max_vo = init_vo; 
vo_at_Ligdecomposition_Starts=init_vo;
%% Read_data
init_CN = 49;
t= table2array(readtable('data/Snajdr et al 2011.xlsx','Range','P4:P12')); %months
C_remaining=table2array(readtable('data/Snajdr et al 2011.xlsx','Range','M4:M12')); %gC litter
lignin_C_remaining=table2array(readtable('data/Snajdr et al 2011.xlsx','Range','F22:F26')); % g C lignin
tobs_Co = table2array(readtable('data/Snajdr et al 2011.xlsx','Range','L22:L26'));  %months
aromaticC = aromatic_fraction_inAIS(lignin_C_remaining);
param.emax = emax_fun(init_CN);

obs_data=[];
obs_data.tobs  =  t;
obs_data.Ct_obs = C_remaining;
obs_data.tobs_Co  = tobs_Co;
obs_data.Co_obs  =aromaticC;

tnorm = obs_data.tobs;
y=log(obs_data.Ct_obs./obs_data.Ct_obs(1));
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=(y(end)-y(1))./tnorm(end);
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);
tt=0:tnorm(end)*n;

figure;
subplot(121)
plot(obs_data.tobs,obs_data.Ct_obs,'o-','DisplayName','TC'); hold on
plot(tt,obs_data.Ct_obs(1).*exp(k.*tt))
scatter(terminalTime,obs_data.Ct_obs(1).*exp(k.*terminalTime),100,'Marker','*', ...
    'MarkerEdgeColor','red')
plot(obs_data.tobs_Co,obs_data.Co_obs,'o-','DisplayName','ligninC'); hold on
xlabel("time"); ylabel('gC/bag')

for j = 1:length(tobs_Co)
    idtemp(j) = find(t == tobs_Co(j));
end
subplot(122)
plot(1-obs_data.Ct_obs(idtemp)./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1),'-o'); hold on
xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')

%%
i=1;
param.emax = emax_fun(init_CN);
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);
param.lb=[0.0005, param.mo,1];
param.ub=[0.04, param.mo,400];
init_guess = [-k, 0.05,10]; % [vh_max, mo,ro];

fig=figure;
[par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
    init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);

param.rsquare = rsquare;param.rmse = rmse;
makeplot_state_space(sol, '',param,obs_data,g,fig)
exportgraphics(fig, fig_path_all+"Snajdr2011.png",'Resolution',100)

LC0=lignin_C_remaining(1)/C_remaining(1); % AISC/totalC

param.numiter=200;
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
par_ = [1, par,LC0,init_CN,LC0*init_CN,[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT,AIC,AIC_Co,AIC_CT] ];


save_fname = excel_path+ "fitted_par_Snajdr.txt";
save(save_fname,"par_",'-ascii','-double','-tabs')

vo = sol.NumericalResults.Control;
time = sol.NumericalResults.Independent;
CT = sol.NumericalResults.State(1, :);
voN = vo ./ max(vo);
id = find(voN > 0.001);
Ligdecomposition_Starts(i) = time(id(1)) ;
C_remain_lig_dec_start(i) = CT(id(1))/CT(1);
avgvo(i) = mean(vo(time<obs_data.tobs(end)));
init_vo(i) = vo(1);
max_vo(i) = max(vo);
vo_at_Ligdecomposition_Starts(i)=vo(id(1));
figure(fig_vo);
nexttile
plot(time, voN,'linewidth', 2); hold on
plot([1, 1]*time(id(1)), [0, 1],'--k','linewidth', 2)
xlabel('day'); ylabel('vo/max(vo)')
title("i=" +i)

exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_, Ligdecomposition_Starts, avgvo,init_vo,max_vo,...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT','AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%%
par_est = load("est_params\" +scenario_name+"\fitted_par_Snajdr.txt");

