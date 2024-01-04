clear all
clc
close all
%%

set_up
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";

initial_fraction_of_C = 0.5;
param.lb=[0.0005, param.mo,1];
param.ub=[0.04, param.mo,400];
%% Upper fitting
treename = {'FC', 'QC', 'AM', 'AR', 'SC', 'BG', 'MJ', 'PH', 'CC', 'CL', 'MO', 'PR', 'AT', 'CJ'};
osono_t = table2array(readtable('data\Osono04.xls', "Sheet", 'Sampling_Dates',"Range", 'C2:C12'));
g_bag = table2array(readtable('data\Osono04.xls', "Sheet",'1)Content',"Range", 'B4:O14'));
mgAUR_gdLitter =  table2array(readtable('data\Osono04.xls',"Sheet", '1)Content',"Range",'B32:O42')); 
mgN_gdLitter_U =table2array(readtable('data\Osono04.xls',"Sheet", '1)Content',"Range",'B144:O154'));

amount_Lig_C0 = (0.001 * mgAUR_gdLitter(1, :))'; % initial fraction of lignin C in litter C
amount_C = g_bag.*initial_fraction_of_C ; %gC/ bag
amount_N0 = mgN_gdLitter_U(1, :) .* g_bag(1, :)*0.001; % g N/ bag
CN0_osono  = (initial_fraction_of_C * g_bag(1, :))./ (amount_N0);
amount_AUR_C_osono04 = fraction_of_C_in_AUR.*mgAUR_gdLitter .* g_bag*0.001; % g lignin / bag
LN0_osono = amount_AUR_C_osono04(1,:)./amount_N0;% g lignin/g N
% scatter(CN0_osono,LN0_osono)
% scatter(amount_Lig_C0,CN0_osono)
% test data from berg with Fixed parameters
% close all
% i=1;
% sample_Osono04
% init_guess = [0.0013, 0.001,100]; % [vh_max, mo,ro];  
% param.emax = emax_fun(CN0_osono(i));
% param.CO_0=obs_data.Co_obs(1);
% param.CT_0=obs_data.Ct_obs(1);
% param.Tmax=obs_data.tobs(end)*5;
% [ocp,sol] =  opt_con(param,g,init_guess,obs_data.tobs(end)*n);
% [ocp,sol] =  opt_con_free_terminal_Time(param,g,init_guess);
% makeplot_state_space(sol, '',param,obs_data,g,fig)

sz =size(amount_C);
par_Osono=[];
fig=figure;
for i = 1:sz(2)
    sample_Osono04
    init_guess = [0.0013, 0.001,100]; % [vh_max, mo,ro];  
    param.CO_0=obs_data.Co_obs(1);
    param.CT_0=obs_data.Ct_obs(1);
    param.emax = emax_fun(CN0_osono(i));
    [par,sol,rmse,rsquare]  = find_parameter( obs_data, param,...
        init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
    param.rmse=rmse;param.rsquare=rsquare;makeplot_state_space(sol, '',param,obs_data,g,fig)
    exportgraphics(fig, fig_path_all+"Osono04Upper_"+i+".png",'Resolution',100)
    par_Osono = [par_Osono;i, par,AISC0,CN0_osono(i),LN0_osono(i),rmse,rsquare];
    disp("current simulation is i ="+ i )
end
save_fname = excel_path+ "fitted_par_Osono04U.txt";
save(save_fname,"par_Osono",'-ascii','-double','-tabs')


%% Lower fitting
osono_t = table2array(readtable('data\Osono04.xls', "Sheet", 'Sampling_Dates',"Range", 'C2:C12'));
g_bag = table2array(readtable('data\Osono04.xls', "Sheet",'1)Content',"Range", 'B17:O27'));
mgAUR_gdLitter =  table2array(readtable('data\Osono04.xls',"Sheet", '1)Content',"Range",'B45:O55')); 
mgN_gdLitter_U =table2array(readtable('data\Osono04.xls',"Sheet", '1)Content',"Range",'B157:O167'));
amount_Lig_C0 = (0.001 * mgAUR_gdLitter(1, :))'; % initial fraction of lignin C in litter C

amount_C = g_bag.*initial_fraction_of_C ; %gC/ bag
amount_N0 = mgN_gdLitter_U(1, :) .* g_bag(1, :)*0.001; % g N/ bag
CN0_osono  = (initial_fraction_of_C * g_bag(1, :))./ (amount_N0);
amount_AUR_C_osono04 = fraction_of_C_in_AUR.*mgAUR_gdLitter .* g_bag*0.001; % g lignin / bag
LN0_osono = amount_AUR_C_osono04(1,:)./amount_N0;% g lignin/g N
% scatter(CN0_osono,LN0_osono)
% scatter(amount_Lig_C0,CN0_osono)
% test data from berg with Fixed parameters
% close all
% i=1;
% sample_Osono04
% fig=figure;
% param.emax = emax_fun(CN0_osono(i));
% param.CO_0=obs_data.Co_obs(1);
% param.CT_0=obs_data.Ct_obs(1);
% param.Tmax=obs_data.tobs(end)*5;
% init_guess = [0.0013, 0.001,100]; % [vh_max, mo,ro]; with parabolic g for berg
% init_guess = [2.9240125131500642e-03	   5.0000000000375307e-02	   1.0627020499413811e+02]; %exp
% [ocp,sol] =  opt_con(param,g,init_guess,obs_data.tobs(end)*n);
% [ocp,sol] =  opt_con_free_terminal_Time(param,g,init_guess);
% makeplot_state_space(sol, '',param,obs_data,g,fig)

% Ososno et al. 2004
close all
sz =size(amount_C);
par_Osono=[];
fig=figure;
for i = 1:sz(2)
    sample_Osono04
    param.CO_0=obs_data.Co_obs(1);
    param.CT_0=obs_data.Ct_obs(1);
    param.emax = emax_fun(CN0_osono(i));
    init_guess = [0.0013, 0.001,100]; 
    [par,sol,rmse,rsquare]  = find_parameter( obs_data, param,...
        init_guess,g,fig,@makeplot_state_space,@ysim_state_space,n);
    param.rmse=rmse;param.rsquare=rsquare;makeplot_state_space(sol, '',param,obs_data,g,fig)
    exportgraphics(fig, fig_path_all+"Osono04Lower_"+i+".png",'Resolution',100)
    AISC0=  obs_data.Co_obs(1)/obs_data.Ct_obs(1);
    par_Osono = [par_Osono;i, par,AISC0,CN0_osono(i),LN0_osono(i),rmse,rsquare];
    disp("current simulation is i ="+ i )

end
save_fname = excel_path+ "fitted_par_Osono04L.txt";
save(save_fname,"par_Osono",'-ascii','-double','-tabs')



