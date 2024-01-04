clear all
clc
close all
%%

set_up

disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";
param.lb=[0.0005, param.mo,0.5];
param.ub=[0.04, param.mo,400];
%% Read_data
initial_fraction_of_C = 0.5;
Hirobe04_data = readtable('data/Hirobe04.xls', 'Sheet', 'C and Mineral');
organic_hirobe = readtable('data/Hirobe04.xls', 'Sheet', 'Organic');

%% test data from berg with Fixed parameters
% close all
% datacode = unique(organic_hirobe.Species);
% i=1;
% sample_hirobe
% CN0=Ct_obs(1)/Nt_obs(1);
% param.emax = emax_fun(CN0);
% 
% fig=figure;
% param.CO_0=obs_data.Co_obs(1);
% param.CT_0=obs_data.Ct_obs(1);
% param.Tmax=obs_data.tobs(end)*5;
% init_guess = [0.03, 0.05,10]; % [vh_max, mo,ro];
% % init_guess = [2.9240125131500642e-03	   5.0000000000375307e-02	   1.0627020499413811e+02]; %exp
% [ocp,sol] =  opt_con(param,g,init_guess,obs_data.tobs(end)*n);
% % [ocp,sol] =  opt_con_free_terminal_Time(param,g,init_guess);
% % makeplot_state_space(sol, '',param,obs_data,g,fig)


%% Hirobe fitting
close all
datacode = unique(organic_hirobe.Species);

fig=figure;

par_Hirobe=[];
for i = 1:length(datacode)
    sample_hirobe
    CN0=Ct_obs(1)/Nt_obs(1);
    param.CO_0=obs_data.Co_obs(1);
    param.CT_0=obs_data.Ct_obs(1);
    param.emax = emax_fun(CN0);
    init_guess = [0.03, 0.05,10]; % [vh_max, mo,ro];
    [par,sol,rmse,rsquare]  = find_parameter(obs_data,param,...
        init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
    param.rmse=rmse;param.rsquare=rsquare;makeplot_state_space(sol, '',param,obs_data,g,fig)
    exportgraphics(fig, fig_path_all+"Hirobe_"+i+"_"+Hirobe04_data.SpeciesName{i}+".png",'Resolution',100)
    AISC0=AIS_C(1)/Ct_obs(1);
    LN0=AIS_C(1)/Nt_obs(1);
    par_Hirobe = [par_Hirobe;i, par,AISC0,CN0,LN0,rmse,rsquare];
    disp("current simulation is i ="+ i )

end
save_fname = excel_path+ "fitted_par_Hirobe.txt";
save(save_fname,"par_Hirobe",'-ascii','-double','-tabs')
