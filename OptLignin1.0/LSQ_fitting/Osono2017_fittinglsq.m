clear all
clc
close all
%%

set_up
close all
set_up
study ="Osono17";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
% param.lb = [0.0005,0.0001, param.mo, 0.01];
% param.ub = [0.04,0.04, param.mo, 400];
%% Read_data
t = readtable('../data\Osono17.xlsx', "Sheet",'Sampling_Dates',"Range", 'B2:B12');
Osono2017 = readtable('../data\Osono17.xlsx', 'Sheet', 'Origianal_data');
iddata = find(Osono2017.Collection == 0);
day = [0, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36]'*30;
treename_Osono2017 = unique(Osono2017.Tree);


%% Osono2017
close all

par_=[];
fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];

numSpecies = length(treename_Osono2017);
for i = 1:numSpecies
   sample_Osono17
   CN0=amount_C(tsrt)/amount_N(tsrt);
   param.emax = emax_fun(CN0);

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
   
   AISC0=amount_AUR_C(tsrt)/amount_C(tsrt);
   LN0=amount_AUR_C(tsrt)/amount_N(tsrt);
   
   par_ = [par_; par,AISC0,CN0,LN0,...
      r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];
   
   disp("current simulation is i ="+ i )
end
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);

writetable(final_table,excel_path+"data summary_LSQ.xlsx")

set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)



