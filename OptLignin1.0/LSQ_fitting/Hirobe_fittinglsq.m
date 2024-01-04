clear all
clc
close all
%%

set_up
study ="Hirobe2004";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
% param.lb = [0.0005,0.0001, param.mo, 0.01];
% param.ub = [0.04,0.04, param.mo, 400];
%% Read_data
initial_fraction_of_C = 0.5;
Hirobe04_data = readtable('../data/Hirobe04.xls', 'Sheet', 'C and Mineral');
organic_hirobe = readtable('../data/Hirobe04.xls', 'Sheet', 'Organic');

%% Hirobe fitting
close all
datacode = unique(organic_hirobe.Species);

fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];

par_=[];
for i = 1:length(datacode)
   sample_hirobe
   CN0=Ct_obs(1)/Nt_obs(1);
   param.CO_0=obs_data.Co_obs(1);
   param.CT_0=obs_data.Ct_obs(1);
   param.emax = emax_fun(CN0);
   init_guess = [0.003	0.003 0.2	50];
   
   [par, sol] = find_parameter_lsq(obs_data, param, ...
      init_guess, @ysim_state_space_lsq)  ;
   [r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);
   
   nexttile
   fit_plotting(sol, obs_data);lh = legend('show');lh.Box = "off";
   title("i=" +i+"(" +Hirobe04_data.SpeciesName{i}+")")
   CT=sol.y(1,:);Co=sol.y(2,:);
   
   dcodt = diff(sol.y(2,:))./diff(sol.x);
   neg_idx  = find(dcodt<0, 1);
   tau=sol.x(neg_idx);
   CTN=CT./CT(1); Ctau=CTN(neg_idx); 
   
   AISC0=AIS_C(1)/Ct_obs(1);
   LN0=AIS_C(1)/Nt_obs(1);
   
   par_ = [par_; par, AISC0, CN0, LN0, r2,...
      rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];
   
   disp("current simulation is i ="+ i )
   
end
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)

