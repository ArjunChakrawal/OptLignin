clear all
clc
close all
%%

set_up
close all
set_up
study ="He2019";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
param.lb = [0.0005,0.0001, param.mo, 0.01];
param.ub = [0.04,0.04, param.mo, 400];
%%

day=table2array(readtable('..\data/He et al 2019.xlsx','Range','H3:H9'));  %gC litter
C_remaining= table2array(readtable('..\data/He et al 2019.xlsx','Range','M3:P9'));
lignin_C_remaining= table2array(readtable('..\data/He et al 2019.xlsx','Range','M12:P18'));
init_CN=[37.10 37.1 32.7 32.7];
SP = readtable('..\data/He et al 2019.xlsx','Range','I2:L2').Properties.VariableNames;
%%

par_=[];
fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];
for i =1:length(SP)
   param.emax = emax_fun(init_CN(i));
   obs_data=[];
   obs_data.tobs  = day;
   obs_data.Ct_obs = C_remaining(:,i);
   obs_data.Co_obs  = lignin_C_remaining(:,i);
   
   init_guess = [0.003	0.003 0.2	50];
   [par, sol] = find_parameter_lsq(obs_data, param, ...
      init_guess, @ysim_state_space_lsq)  ;
   [r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);
   
   nexttile
   fit_plotting(sol, obs_data);lh = legend('show');lh.Box = "off";
   title("i=" +i)
   CT=sol.y(1,:);Co=sol.y(2,:);
   
   dcodt = diff(sol.y(2,:))./diff(sol.x);
   neg_idx  = find(dcodt<0, 1);
   tau=sol.x(neg_idx);
   CTN=CT./CT(1); Ctau=CTN(neg_idx); 
   
   LC0=lignin_C_remaining(1,i)/obs_data.Ct_obs(1);
   par_ = [par_; par,LC0,init_CN(i),LC0*init_CN(i),...
      r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];
   
end

final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)

