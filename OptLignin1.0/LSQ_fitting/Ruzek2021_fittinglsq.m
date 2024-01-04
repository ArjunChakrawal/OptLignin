clear all
clc
close all
%%
close all
set_up
study ="Ruzek2021";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
param.lb = [0.0005,0.0001, param.mo, 0.01];
param.ub = [0.04,0.04, param.mo, 400];
%% Rooibos tea

Cday=table2array(readtable('..\data/Ruzek et al 2021.xlsx','sheet', ...
   'processed_data','Range','B6:B10'));
C_remaining= table2array(readtable('..\data/Ruzek et al 2021.xlsx','sheet', ...
   'processed_data','Range','G6:J10'));
Lday = table2array(readtable('..\data/Ruzek et al 2021.xlsx','sheet', ...
   'processed_data','Range','L6:L9'));
lignin_C_remaining= table2array(readtable('..\data/Ruzek et al 2021.xlsx','sheet', ...
   'processed_data','Range','Q6:T9'));
SP=readtable('..\data/Ruzek et al 2021.xlsx','sheet', ...
   'processed_data','Range','M11:P11').Properties.VariableNames;
init_CN=53.2;
%%
par_=[];
fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];
for i =1:size(C_remaining,2)
   param.emax = emax_fun(init_CN);
   obs_data=[];
   obs_data.tobs  = Cday;
   obs_data.Ct_obs = C_remaining(:,i);
   obs_data.tobs_Co =Lday;
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
   par_ = [par_;par,LC0,init_CN,LC0*init_CN,...
      r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];
   
end
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)


%% green tea
% Cday=table2array(readtable('..\data/Ruzek et al 2021.xlsx','sheet', ...
%     'processed_data','Range','B6:B10'));  %gC litter
%
% C_remaining= table2array(readtable('..\data/Ruzek et al 2021.xlsx','sheet', ...
%     'processed_data','Range','C6:F10'));
% Lday = table2array(readtable('..\data/Ruzek et al 2021.xlsx','sheet', ...
%     'processed_data','Range','L6:L9'));
% lignin_C_remaining= table2array(readtable('..\data/Ruzek et al 2021.xlsx','sheet', ...
%     'processed_data','Range','M6:P9'));
% init_CN=13.4;
%
% par_Ruzek=[];
% for i =1:size(C_remaining,2)
% param.emax = emax_fun(init_CN);
% obs_data=[];
% obs_data.tobs  = Cday;
% obs_data.Ct_obs = C_remaining(:,i);
% obs_data.tobs_Co =Lday;
% obs_data.Co_obs  = lignin_C_remaining(:,i);
% for j = 1:length(Lday)
%     idtemp(j) = find(Cday== Lday(j));
% end
% tnorm = obs_data.tobs;
% f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
% final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
% k=f.b ;
% terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
% % n=terminalTime./obs_data.tobs(end);
% n=1.5;
% tt=0:tnorm(end)*n;
%
% param.CO_0=obs_data.Co_obs(1);
% param.CT_0=obs_data.Ct_obs(1);
%
% figure;
% subplot(131)
% plot(obs_data.tobs,obs_data.Ct_obs,'o-','DisplayName','TC'); hold on
% plot(tt,obs_data.Ct_obs(1)*exp(k.*tt))
% scatter(terminalTime,exp(k.*terminalTime),100,'Marker','*', ...
%     'MarkerEdgeColor','red')
% plot(obs_data.tobs_Co,obs_data.Co_obs,'o-','DisplayName','ligninC'); hold on
% xlabel("time"); ylabel('gC/bag')
% subplot(132)
% scatter(1-obs_data.Ct_obs(idtemp)./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
% xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
% subplot(133)
% scatter(obs_data.tobs_Co,obs_data.Co_obs./obs_data.Ct_obs(idtemp));
% ylabel('L/C');xlabel("time");
%
% end
% %
% fig=figure;
% init_guess = [0.01, 0.0001,0.005]; % [vh_max, mo,ro];
% [ocp,sol] =  opt_con(param,g,init_guess,obs_data.tobs(end)*n);
% makeplot_state_space(sol, '',param,obs_data,g,fig)
%
%
% init_guess = [-k, 0.05,50]; % [vh_max, mo,ro];
%
% [par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
%     init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
% makeplot_state_space(sol, '',param,obs_data,g,fig)
% exportgraphics(fig, fig_path_all+"Yue_"+SP{i}+".png",'Resolution',100)
%
% LC0=lignin_C_remaining(1,i)/obs_data.Ct_obs(1);
% par_Ruzek = [par_Ruzek;i, par,init_CN(i),LC0*init_CN(i),rmse,rsquare];
%
% end