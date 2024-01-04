clear all
clc
close all
%%
close all
set_up
study ="Osono2004Upper";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
% param.lb = [0.0005,0.0001, param.mo, 0.01];
% param.ub = [0.04,0.04, param.mo, 400];
% Upper fitting
treename = {'FC', 'QC', 'AM', 'AR', 'SC', 'BG', 'MJ', 'PH', 'CC', 'CL', 'MO', 'PR', 'AT', 'CJ'};
osono_t = table2array(readtable('..\data\Osono04.xls', "Sheet", 'Sampling_Dates',"Range", 'C2:C12'));
g_bag = table2array(readtable('..\data\Osono04.xls', "Sheet",'1)Content',"Range", 'B4:O14'));
mgAUR_gdLitter =  table2array(readtable('..\data\Osono04.xls',"Sheet", '1)Content',"Range",'B32:O42'));
mgN_gdLitter_U =table2array(readtable('..\data\Osono04.xls',"Sheet", '1)Content',"Range",'B144:O154'));

amount_Lig_C0 = (0.001 * mgAUR_gdLitter(1, :))'; % initial fraction of lignin C in litter C
amount_C = g_bag.*initial_fraction_of_C ; %gC/ bag
amount_N0 = mgN_gdLitter_U(1, :) .* g_bag(1, :)*0.001; % g N/ bag
CN0_osono  = (initial_fraction_of_C * g_bag(1, :))./ (amount_N0);
amount_AUR_C_osono04 = fraction_of_C_in_AUR.*mgAUR_gdLitter .* g_bag*0.001; % g lignin / bag
LN0_osono = amount_AUR_C_osono04(1,:)./amount_N0;% g lignin/g N


sz =size(amount_C);
par_=[];
fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];

for i = 1:sz(2)
   sample_Osono04
   param.CO_0=obs_data.Co_obs(1);
   param.CT_0=obs_data.Ct_obs(1);
   param.emax = emax_fun(CN0_osono(i));
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
   
   par_ = [par_; par,AISC0,CN0_osono(i),LN0_osono(i),r2,...
      rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];
   
   disp("current simulation is i ="+ i )
end

final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)

%% Lower fitting
osono_t = table2array(readtable('..\data\Osono04.xls', "Sheet", 'Sampling_Dates',"Range", 'C2:C12'));
g_bag = table2array(readtable('..\data\Osono04.xls', "Sheet",'1)Content',"Range", 'B17:O27'));
mgAUR_gdLitter =  table2array(readtable('..\data\Osono04.xls',"Sheet", '1)Content',"Range",'B45:O55'));
mgN_gdLitter_U =table2array(readtable('..\data\Osono04.xls',"Sheet", '1)Content',"Range",'B157:O167'));
amount_Lig_C0 = (0.001 * mgAUR_gdLitter(1, :))'; % initial fraction of lignin C in litter C

amount_C = g_bag.*initial_fraction_of_C ; %gC/ bag
amount_N0 = mgN_gdLitter_U(1, :) .* g_bag(1, :)*0.001; % g N/ bag
CN0_osono  = (initial_fraction_of_C * g_bag(1, :))./ (amount_N0);
amount_AUR_C_osono04 = fraction_of_C_in_AUR.*mgAUR_gdLitter .* g_bag*0.001; % g lignin / bag
LN0_osono = amount_AUR_C_osono04(1,:)./amount_N0;% g lignin/g N

% Ososno et al. 2004
close all
study ="Osono2004Lower";

sz =size(amount_C);
par_=[];
fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];
for i = 1:sz(2)
   sample_Osono04
   param.CO_0=obs_data.Co_obs(1);
   param.CT_0=obs_data.Ct_obs(1);
   param.emax = emax_fun(CN0_osono(i));
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
   
   par_ = [par_; par,AISC0,CN0_osono(i),LN0_osono(i),r2,...
      rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau];
   
   disp("current simulation is i ="+ i )
   
end
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);
writetable(final_table,excel_path+"data summary_LSQ.xlsx")
set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)



