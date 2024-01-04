clear all
clc
close all
%%

set_up
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";

param.lb=[0.0005, param.mo,1];
param.ub=[0.04, param.mo,400];
%% Read_data
close all
day=table2array(readtable('data/Yue et al 2016.xlsx','Range','A30:A35'));  %gC litter
C_remaining=table2array(readtable('data/Yue et al 2016.xlsx','Range','J30:M35'));  %gC litter
lignin_C_remaining=table2array(readtable('data/Yue et al 2016.xlsx','Range','N30:Q35'));  %gC lignin

init_C = table2array(readtable('data/Yue et al 2016.xlsx','Range','B23:E23'));
init_CN = table2array(readtable('data/Yue et al 2016.xlsx','Range','B25:E25'));
SP = readtable('data/Yue et al 2016.xlsx','Range','B22:E22').Properties.VariableNames;

%%
par_Yue=[];
for i =1:length(SP)
param.emax = emax_fun(init_CN(i));
obs_data=[];
obs_data.tobs  = day;
obs_data.Ct_obs = C_remaining(:,i);
obs_data.Co_obs  = aromatic_fraction_inAIS(lignin_C_remaining(:,i));

tnorm = obs_data.tobs;
f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=f.b ;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);
tt=0:tnorm(end)*n;

param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);

figure(1);
subplot(121)
plot(obs_data.tobs,obs_data.Ct_obs./obs_data.Ct_obs(1),'o-','DisplayName','TC'); hold on
plot(tt,1*exp(k.*tt))
scatter(terminalTime,exp(k.*terminalTime),100,'Marker','*', ...
    'MarkerEdgeColor','red')
plot(obs_data.tobs,obs_data.Co_obs./obs_data.Co_obs(1),'o-','DisplayName','ligninC'); hold on
xlabel("time"); ylabel('gC/bag')
subplot(122)
scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
if(i==1)
    n=1;
end

fig=figure(2);
init_guess = [-k, 0.05,50]; % [vh_max, mo,ro];
[par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
    init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
param.rmse=rmse;param.rsquare=rsquare;
makeplot_state_space(sol, '',param,obs_data,g,fig)
exportgraphics(fig, fig_path_all+"Yue_"+SP{i}+".png",'Resolution',100)

LC0=lignin_C_remaining(1,i)/obs_data.Ct_obs(1);
par_Yue = [par_Yue;i, par,LC0,init_CN(i),LC0*init_CN(i),rmse,rsquare];

end

save_fname = excel_path+ "fitted_par_Yue.txt";
save(save_fname,"par_Yue",'-ascii','-double','-tabs')

