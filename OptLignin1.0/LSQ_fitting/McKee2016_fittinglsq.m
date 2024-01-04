clear all
clc
close all
%%

set_up

study = "McKee2016";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";

%%
init_CN=44.3/1.47;
param.emax = emax_fun(init_CN);
fieldmass = [18.4 12.6313 6.21871 5.0279 4.076 ]'; % in g
obs_data=[];
obs_data.tobs  =[0 180.00 364.00 546.00 730.00]' ;
obs_data.Ct_obs = fieldmass*0.443; % in g
AUR =  [0.723 1.2673 1.0148 0.818 0.6913 ]' ;
obs_data.Co_obs  =aromatic_fraction_inAIS( AUR*fraction_of_C_in_AUR);

init_guess = [0.003	0.003 0.2	50];
[par, sol] = find_parameter_lsq(obs_data, param, ...
   init_guess, @ysim_state_space_lsq)  ;
[r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);

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

LC0=obs_data.Co_obs(1)/obs_data.Ct_obs(1);
par_ = [ par,LC0,init_CN,LC0*init_CN,...
   r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau ];

final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)


