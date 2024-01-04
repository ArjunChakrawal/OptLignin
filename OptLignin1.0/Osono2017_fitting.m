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
t = readtable('data\Osono17.xlsx', "Sheet",'Sampling_Dates',"Range", 'B2:B12');
Osono2017 = readtable('data\Osono17.xlsx', 'Sheet', 'Origianal_data');
iddata = find(Osono2017.Collection == 0);
day = [0, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36]'*30;
treename_Osono2017 = unique(Osono2017.Tree);

%% test data from berg with Fixed parameters
% close all
% i=1;
% sample_Osono17
% CN0=amount_C(tsrt)/amount_N(tsrt);
% 
% fig=figure;
% param.emax = emax_fun(CN0);
% param.CO_0=obs_data.Co_obs(1);
% param.CT_0=obs_data.Ct_obs(1);
% param.Tmax=obs_data.tobs(end)*5;
% init_guess = [0.03, 0.05,100]; % [vh_max, mo,ro]; 
% % init_guess = [2.9240125131500642e-03	   5.0000000000375307e-02	   1.0627020499413811e+02]; %exp
% [ocp,sol] =  opt_con(param,g,init_guess,obs_data.tobs(end)*n);
% % [ocp,sol] =  opt_con_free_terminal_Time(param,g,init_guess);
% makeplot_state_space(sol, '',param,obs_data,g,fig)

%% Osono2017
close all
fig=figure(2);
par_Osono17=[];
numSpecies = length(treename_Osono2017);
for i = 1:numSpecies
    sample_Osono17
    CN0=amount_C(tsrt)/amount_N(tsrt);
    param.emax = emax_fun(CN0);
    param.CO_0=obs_data.Co_obs(1);
    param.CT_0=obs_data.Ct_obs(1);
    init_guess = [0.03, 0.05,100]; % [vh_max, mo,ro]; 
    [par,sol,rmse,rsquare]  = find_parameter( obs_data, param,...
        init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
    param.rmse=rmse;param.rsquare=rsquare;makeplot_state_space(sol, '',param,obs_data,g,fig)
    exportgraphics(fig, fig_path_all+"Osono17_"+i+".png",'Resolution',100)

    AISC0=amount_AUR_C(tsrt)/amount_C(tsrt);
    LN0=amount_AUR_C(tsrt)/amount_N(tsrt);
    par_Osono17 = [par_Osono17;i, par,AISC0,CN0,LN0,rmse,rsquare];
    disp("current simulation is i ="+ i )

end
save_fname = excel_path+ "fitted_par_Osono17.txt";
save(save_fname,"par_Osono17",'-ascii','-double','-tabs')
