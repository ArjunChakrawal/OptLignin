clear all
clc
close all
%%

set_up

study="He2016";

disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";
fig_pathCombined = "fig\"+scenario_name+"\combined\";

fig_vo=figure;tiledlayout('flow');
Ligdecomposition_Starts = zeros(2, 1);
C_remain_lig_dec_start=Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo= avgvo;
max_vo=init_vo;
vo_at_Ligdecomposition_Starts=init_vo;

%% Read_data
data= readtable('data/He et al 2016.xlsx','Range','A21:G32'); 
CN0=[35.29,25.93]; % Fargesia nitida, Salix paraplesia
%%
i=1;
obs_data=[];
obs_data.tobs  =  data.day;
obs_data.Ct_obs = data.C_g_FN;
obs_data.Co_obs  = aromatic_fraction_inAIS(data.LC_g_FN);
tnorm = obs_data.tobs;
f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=f.b ;
tt=0:tnorm(end)*2;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);

figure;
scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
par_=[];

fig=figure;
param.emax = emax_fun(CN0(i));
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);
param.lb=[0.0005, param.mo,1];
param.ub=[0.004, param.mo,400];

init_guess = [0.002, 0.025,150]; % [vh_max, mo,ro];
[par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
    init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
param.rmse=rmse;param.rsquare=rsquare;
makeplot_state_space(sol, '',param,obs_data,g,fig)
exportgraphics(fig, fig_path_all+"He2016_"+"Fargesia nitida"+".png",'Resolution',100)
LC0=data.LC_g_FN(1)/data.C_g_FN(1);
param.numiter=200;
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
par_ = [par_;i, par,LC0,CN0(i),LC0*CN0(i),[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT,AIC,AIC_Co,AIC_CT]];

vo = sol.NumericalResults.Control;
time = sol.NumericalResults.Independent;
CT = sol.NumericalResults.State(1, :);
voN = vo ./ max(vo);
id = find(voN > vo_thres);
Ligdecomposition_Starts(i) = time(id(1)) ;
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
obs_data=[];
obs_data.tobs  =  data.day;
obs_data.Ct_obs = data.C_g_SP;
obs_data.Co_obs  = aromatic_fraction_inAIS(data.LC_g_SP);
param.emax = emax_fun(CN0(i));
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);

tnorm = obs_data.tobs;
f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=f.b ;
tt=0:tnorm(end)*2;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);

figure;
scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
param.numiter=50;
fig=figure;
init_guess = [0.002, 0.025,150]; % [vh_max, mo,ro];
[par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
    init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
param.rmse=rmse;param.rsquare=rsquare;
makeplot_state_space(sol, '',param,obs_data,g,fig)
exportgraphics(fig, fig_path_all+"He2016_"+"Salix paraplesia"+".png",'Resolution',100)

LC0=data.LC_g_FN(1)/data.C_g_FN(1);
param.numiter=200;
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
par_ = [par_;i, par,LC0,CN0(i),LC0*CN0(i),[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT,AIC,AIC_Co,AIC_CT]];

save_fname = excel_path+ "fitted_par_He2016.txt";
save(save_fname,"par_",'-ascii','-double','-tabs')

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
exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_, Ligdecomposition_Starts, avgvo,init_vo,max_vo,...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i','vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT','AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts','C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")



