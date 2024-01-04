clear all
clc
close all
%%

set_up
study="Snajdr2011";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";

%% Read_data
init_CN = 49;
t= table2array(readtable('..\data/Snajdr et al 2011.xlsx','Range','P4:P12')); %months
C_remaining=table2array(readtable('..\data/Snajdr et al 2011.xlsx','Range','M4:M12')); %gC litter
lignin_C_remaining=table2array(readtable('..\data/Snajdr et al 2011.xlsx','Range','F22:F26')); % g C lignin
tobs_Co = table2array(readtable('..\data/Snajdr et al 2011.xlsx','Range','L22:L26'));  %months
aromaticC = aromatic_fraction_inAIS(lignin_C_remaining);
param.emax = emax_fun(init_CN);

obs_data=[];
obs_data.tobs  =  t;
obs_data.Ct_obs = C_remaining;
obs_data.tobs_Co  = tobs_Co;
obs_data.Co_obs  =aromaticC;


%%
i=1;
param.emax = emax_fun(init_CN);
init_guess = [0.003	0.003 0.2	50];
[par, sol] = find_parameter_lsq(obs_data, param, ...
   init_guess, @ysim_state_space_lsq)  ;
[r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);
par_=[];
fig = figure;
nexttile
fit_plotting(sol, obs_data);lh = legend('show');lh.Box = "off";
title("i=" +i)
CT=sol.y(1,:);Co=sol.y(2,:);

dcodt = diff(sol.y(2,:))./diff(sol.x);
neg_idx  = find(dcodt<0, 1);
if(isempty(neg_idx))
   neg_idx=length(sol.x);
end
tau=sol.x(neg_idx);
CTN=CT./CT(1); Ctau=CTN(neg_idx);

LC0=lignin_C_remaining(1)/C_remaining(1); % AISC/totalC
par_ = [par,LC0,init_CN,LC0*init_CN,...
   r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau ];

final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)




