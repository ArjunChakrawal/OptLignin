clear all
clc
close all
%%

set_up
study="Fioretto2005";

disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
par_=[];
fig = figure; fig.Position=[239 607 784 262];
tiledlayout('flow');

%% Read_data
i=1;
CIncanus= readtable('../data/Fioretto et al 2005.xlsx','Range','B10:F21');
CN0 = 48;
param.emax = emax_fun(CN0);

obs_data=[];
obs_data.tobs  =  CIncanus.year.*365;
obs_data.Ct_obs = CIncanus.CRemainingG;
obs_data.Co_obs  = aromatic_fraction_inAIS(CIncanus.ligninCRemainingG);

init_guess = [0.003	0.003 0.2	50];
[par, sol] = find_parameter_lsq(obs_data, param, ...
   init_guess, @ysim_state_space_lsq)  ;
[r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);

nexttile
fit_plotting(sol, obs_data);lh = legend('show');lh.Box = "off";
title("i=" +i)
CT=sol.y(1,:);

dcodt = diff(sol.y(2,:))./diff(sol.x);
neg_idx  = find(dcodt<0, 1);
if(isempty(neg_idx))
   neg_idx=length(sol.x);
end
tau=sol.x(neg_idx);
CTN=CT./CT(1); Ctau=CTN(neg_idx);

LC0=CIncanus.ligninCRemainingG(1)/obs_data.Ct_obs(1);
par_ = [par_;par,LC0,CN0,LC0*CN0,...
         r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];


%%
i=2;
MCommunis= readtable('..\data/Fioretto et al 2005.xlsx','Range','B22:F33');
CN0 = 62;
param.emax = emax_fun(CN0);

obs_data=[];
obs_data.tobs  =  MCommunis.year.*365;
obs_data.Ct_obs = MCommunis.CRemainingG;
obs_data.Co_obs  = aromatic_fraction_inAIS(MCommunis.ligninCRemainingG);

init_guess = [0.003	0.003 0.2	50];
[par, sol] = find_parameter_lsq(obs_data, param, ...
   init_guess, @ysim_state_space_lsq)  ;
[r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);

nexttile
fit_plotting(sol, obs_data);lh = legend('show');lh.Box = "off";
title("i=" +i)
CT=sol.y(1,:);

dcodt = diff(sol.y(2,:))./diff(sol.x);
neg_idx  = find(dcodt<0, 1);
if(isempty(neg_idx))
   neg_idx=length(sol.x);
end
tau=sol.x(neg_idx);
CTN=CT./CT(1); Ctau=CTN(neg_idx);

LC0=MCommunis.ligninCRemainingG(1)/obs_data.Ct_obs(1);
par_ = [par_;par,LC0,CN0,LC0*CN0,...
         r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];


%%
i=3;
Qliex= readtable('..\data/Fioretto et al 2005.xlsx','Range','B34:F44');
CN0 = 46;
param.emax = emax_fun(CN0);

obs_data=[];
obs_data.tobs  =  Qliex.year.*365;
obs_data.Ct_obs = Qliex.CRemainingG;
obs_data.Co_obs  = aromatic_fraction_inAIS(Qliex.ligninCRemainingG);
init_guess = [0.003	0.003 0.2	50];
[par, sol] = find_parameter_lsq(obs_data, param, ...
   init_guess, @ysim_state_space_lsq)  ;
[r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);

nexttile
fit_plotting(sol, obs_data);lh = legend('show');lh.Box = "off";
title("i=" +i)
CT=sol.y(1,:);

dcodt = diff(sol.y(2,:))./diff(sol.x);
neg_idx  = find(dcodt<0, 1);
if(isempty(neg_idx))
   neg_idx=length(sol.x);
end
tau=sol.x(neg_idx);
CTN=CT./CT(1); Ctau=CTN(neg_idx);

LC0=Qliex.ligninCRemainingG(1)/obs_data.Ct_obs(1);
par_ = [par_;par,LC0,CN0,LC0*CN0,...
         r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];

%%
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)



