clear all
clc
close all
%%

set_up
study="Fioretto2005";

disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";

fig_pathCombined = "fig\"+scenario_name+"\combined\";
fig_vo=figure;tiledlayout('flow');
Ligdecomposition_Starts = zeros(3, 1);
C_remain_lig_dec_start=Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo= avgvo;
max_vo=init_vo;
vo_at_Ligdecomposition_Starts=init_vo;

param.lb=[0.0005, param.mo,1];
param.ub=[0.04, param.mo,400];
%% Read_data
i=1;
CIncanus= readtable('data/Fioretto et al 2005.xlsx','Range','B10:F21');
CN0 = 48;
param.emax = emax_fun(CN0);

obs_data=[];
obs_data.tobs  =  CIncanus.year.*365;
obs_data.Ct_obs = CIncanus.CRemainingG;
obs_data.Co_obs  = aromatic_fraction_inAIS(CIncanus.ligninCRemainingG);
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);

% figure;
% scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
% xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')

tnorm = obs_data.tobs;
f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=f.b ;
tt=0:tnorm(end)*2;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);

par_Fioretto=[];

fig=figure;
init_guess = [-k, 0.05,50]; % [vh_max, mo,ro];
[par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
    init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
param.rmse=rmse;param.rsquare=rsquare;
makeplot_state_space(sol, '',param,obs_data,g,fig)
exportgraphics(fig, fig_path_all+"Fioretto_C_Incanus.png",'Resolution',100)

LC0=CIncanus.ligninCRemainingG(1)/obs_data.Ct_obs(1);
param.numiter=200;
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);

par_Fioretto = [par_Fioretto;i, par,LC0,CN0,LC0*CN0,[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT,AIC,AIC_Co,AIC_CT]];

vo = sol.NumericalResults.Control;
time = sol.NumericalResults.Independent;
CT = sol.NumericalResults.State(1, :);
voN = vo ./ max(vo);
id = find(voN > vo_thres);
Ligdecomposition_Starts(i) = time(id(1));
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
%%
i=2;
MCommunis= readtable('data/Fioretto et al 2005.xlsx','Range','B22:F33');
CN0 = 62;
param.emax = emax_fun(CN0);

obs_data=[];
obs_data.tobs  =  MCommunis.year.*365;
obs_data.Ct_obs = MCommunis.CRemainingG;
obs_data.Co_obs  = aromatic_fraction_inAIS(MCommunis.ligninCRemainingG);
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);
% figure;
% scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
% xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')

tnorm = obs_data.tobs;
f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=f.b ;
tt=0:tnorm(end)*2;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);
fig=figure;
param.numiter=50;

init_guess = [-k, 0.05,50]; % [vh_max, mo,ro];
[par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
    init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
param.rmse=rmse;param.rsquare=rsquare;
makeplot_state_space(sol, '',param,obs_data,g,fig)
exportgraphics(fig, fig_path_all+"Fioretto_MCommunis.png",'Resolution',100)

LC0=MCommunis.ligninCRemainingG(1)/obs_data.Ct_obs(1);
param.numiter=200;
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);

par_Fioretto = [par_Fioretto;i, par,LC0,CN0,LC0*CN0,[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT,AIC,AIC_Co,AIC_CT]];

% par_Fioretto = [par_Fioretto;i, par,LC0,CN0,LC0*CN0,rmse,rsquare];
vo = sol.NumericalResults.Control;
time = sol.NumericalResults.Independent;
CT = sol.NumericalResults.State(1, :);
voN = vo ./ max(vo);
id = find(voN > vo_thres);
Ligdecomposition_Starts(i) = time(id(1));
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
%%
i=3;
Qliex= readtable('data/Fioretto et al 2005.xlsx','Range','B34:F44');
CN0 = 46;
param.emax = emax_fun(CN0);

obs_data=[];
obs_data.tobs  =  Qliex.year.*365;
obs_data.Ct_obs = Qliex.CRemainingG;
obs_data.Co_obs  = aromatic_fraction_inAIS(Qliex.ligninCRemainingG);
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);
% figure;
% scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
% xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')

tnorm = obs_data.tobs;
f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=f.b ;
tt=0:tnorm(end)*2;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);
fig=figure;
param.numiter=50;

init_guess = [0.01, 0.02,0.300]; % [vh_max, mo,ro];
[par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
    init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
param.rmse=rmse;param.rsquare=rsquare;
makeplot_state_space(sol, '',param,obs_data,g,fig)

exportgraphics(fig, fig_path_all+"Fioretto_Qliex.png",'Resolution',100)

LC0=Qliex.ligninCRemainingG(1)/obs_data.Ct_obs(1);
param.numiter=200;
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);

par_Fioretto = [par_Fioretto;i, par,LC0,CN0,LC0*CN0,[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT,AIC,AIC_Co,AIC_CT]];

% par_Fioretto = [par_Fioretto;i, par,LC0,CN0,LC0*CN0,rmse,rsquare];
vo = sol.NumericalResults.Control;
time = sol.NumericalResults.Independent;
CT = sol.NumericalResults.State(1, :);
voN = vo ./ max(vo);
id = find(voN > vo_thres);
Ligdecomposition_Starts(i) = time(id(1));
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
%%
save_fname = excel_path+ "fitted_par_Fioretto.txt";
save(save_fname,"par_Fioretto",'-ascii','-double','-tabs')

exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_Fioretto, Ligdecomposition_Starts, avgvo,init_vo,max_vo,...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i','vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT','AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts','C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")
