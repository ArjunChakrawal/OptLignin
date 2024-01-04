clear all
clc
close all
%%

set_up
study="He2016";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
%% Read_data
data= readtable('..\data/He et al 2016.xlsx','Range','A21:G32'); 
CN0=[35.29,25.93]; % Fargesia nitida, Salix paraplesia
fig = figure;
tiledlayout('flow');
par_=[];
%%
i=1;
obs_data=[];
obs_data.tobs  =  data.day;
obs_data.Ct_obs = data.C_g_FN;
obs_data.Co_obs  = aromatic_fraction_inAIS(data.LC_g_FN);

param.emax = emax_fun(CN0(i));

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

LC0=data.LC_g_FN(1)/data.C_g_FN(1);
par_ = [par_; par,LC0,CN0(i),LC0*CN0(i),...
         r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];

%%
i=2;
obs_data=[];
obs_data.tobs  =  data.day;
obs_data.Ct_obs = data.C_g_SP;
obs_data.Co_obs  = aromatic_fraction_inAIS(data.LC_g_SP);
param.emax = emax_fun(CN0(i));

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

LC0=data.LC_g_FN(1)/data.C_g_FN(1);
par_ = [par_; par,LC0,CN0(i),LC0*CN0(i),...
         r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];


%%
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)


