clear all
clc
close all
%%

set_up
% l=0:0.001:0.05;
% figure;
% plot(l, g(l, a*0.1, b*0.1), 'LineWidth', 2, 'DisplayName', ...
%     "sigmoid $g = 1-\frac{1}{1+exp \left(-\frac{L-0.125}{0.0125}\right)}$"); hold on
% ylabel("Factor reducing carbohydrate uptake")
% xlabel('Lignocellulose index = Co/(Ch+Co)');
study = "McKee2016";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";
fig_pathCombined = "fig\"+scenario_name+"\combined\";
fig_vo=figure;tiledlayout('flow');
Ligdecomposition_Starts = zeros(1, 1);
C_remain_lig_dec_start=Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo= avgvo;
max_vo=init_vo;
vo_at_Ligdecomposition_Starts=init_vo;


param.lb=[0.0005, param.mo,1];
param.ub=[0.04, param.mo,400];

%%
init_CN=44.3/1.47;
param.emax = emax_fun(init_CN);
fieldmass = [18.4 12.6313 6.21871 5.0279 4.076 ]'; % in g
obs_data=[];
obs_data.tobs  =[0 180.00 364.00 546.00 730.00]' ;
obs_data.Ct_obs = fieldmass*0.443; % in g
AUR =  [0.723 1.2673 1.0148 0.818 0.6913 ]' ;
obs_data.Co_obs  =aromatic_fraction_inAIS( AUR*fraction_of_C_in_AUR);

tnorm = obs_data.tobs;
f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=f.b ;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);
tt=0:tnorm(end)*n;

param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);

figure;
subplot(121)
plot(obs_data.tobs,obs_data.Ct_obs,'o-','DisplayName','TC'); hold on
plot(tt,obs_data.Ct_obs(1)*exp(k.*tt))
scatter(terminalTime,exp(k.*terminalTime),100,'Marker','*', ...
    'MarkerEdgeColor','red')
plot(obs_data.tobs,obs_data.Co_obs,'o-','DisplayName','ligninC'); hold on
xlabel("time"); ylabel('gC/bag')
subplot(122)
scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')

%%
fig=figure;
init_guess = [0.0075, 0.1,50]; % [vh_max, mo,ro];
[par,sol]  = find_parameter(obs_data,param,...
    init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);

makeplot_state_space(sol, '',param,obs_data,g,fig)
exportgraphics(fig, fig_path_all+"McKee_grassLitter.png",'Resolution',100)

LC0=obs_data.Co_obs(1)/obs_data.Ct_obs(1);

param.numiter=200;
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);

par_ = [1, par,LC0,init_CN,LC0*init_CN,[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT,AIC,AIC_Co,AIC_CT] ];

% par_ = [1,par,LC0,init_CN,LC0*init_CN,rmse,rsquare];

save_fname = excel_path+ "fitted_par_McKee2016.txt";
save(save_fname,"par_",'-ascii','-double','-tabs')

i=1;
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

exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_, Ligdecomposition_Starts, avgvo,init_vo,max_vo,...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i','vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT','AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts','C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

