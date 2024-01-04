clear all
clc
close all
%%

set_up
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";

%% Read_data
close all
t= table2array(readtable('data/Peng et al 2022/Peng et al 2022.xlsx','Range','A3:A16')); %months
gCt_obs=readtable('data/Peng et al 2022/Peng et al 2022.xlsx','Range','H2:M16');  %gC litter
gCo_obs=readtable('data/Peng et al 2022/Peng et al 2022.xlsx','Range','H19:M33'); % g C lignin
init_N = table2array(readtable('data/Peng et al 2022/Peng et al 2022.xlsx','Range','Y2:Y3'));
spcode = gCt_obs.Properties.VariableNames;

%%
close all
fig=figure;ax=gca;
tend = length(t)-1;

% tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
for i= 1:length(spcode)
    sample_Peng
    plot(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1),'-o'); hold on
    xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
end

%%
% n=1;
% close all
tend = length(t)-0;
i=1;
sample_Peng
fig=figure;
param.vomax=vomax;
param.emax = emax;
param.f = 1;
param.a= a;
param.b = b;
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);
% param.Tmax=obs_data.tobs(end)*5;
init_guess = [0.002, 0.025,150]; % [vh_max, mo,ro];
[ocp,sol] =  opt_con(param,g,init_guess,obs_data.tobs(end)*n);

makeplot_state_space(sol, '',param,obs_data,g,fig)

ysim=ysim_state_space(init_guess, ocp, obs_data, [param.CT_0, param.CO_0]);
y1norm=obs_data.Ct_obs./obs_data.Ct_obs(1);
y2norm= obs_data.Co_obs./obs_data.Co_obs(1);
ydata = [y1norm;y2norm];

rmse = sqrt(mean(ysim-ydata).^2);
SSR =sum((ysim-ydata).^2);
SST = sum((mean([y1norm; y2norm])-[y1norm; y2norm]).^2);
rsquare= 1-SSR/SST;
nexttile;scatter(ysim, ydata)

par_Peng = [[];i, init_guess,table2array(gCo_obs(i,1))/param.CT_0(1),param.CT_0(1)/init_N(1),...
    table2array(gCo_obs(i,1))/init_N(1),rmse,rsquare];

%% fitting
% close all
fig=figure;

par_Peng=[];
i=6;
% for i = 1:length(spcode)
sample_Peng
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);
param.lb=[0.0005, 0.001,1];
param.ub=[0.004, 0.1,400];

[par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
    init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
makeplot_state_space(sol, '',param,obs_data,g,fig)
par_Peng = [par_Peng;i, par,table2array(gCo_obs(i,1))/param.CT_0(1),param.CT_0(1)/init_N(1),...
    table2array(gCo_obs(i,1))/init_N(1),rmse,rsquare];
disp("current simulation is i ="+ i )

% end
save_fname = excel_path+ "fitted_par_Peng.txt";
save(save_fname,"par_Peng",'-ascii','-double','-tabs')
