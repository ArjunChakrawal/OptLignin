clear all
clc
close all

%%
set_up
disp("this is the g-function with " +scenario_name+" scenario")
excel_path = "est_params\" +scenario_name + "\";
fig_path_all = "fig\" +scenario_name + "\";

%% Read_data
initial_fraction_of_C = 0.5;
berg_data = readtable('data\Berg and McClaugherty 1989 - Suppl Mat.xlsx');
datacode = unique(berg_data.DatasetCode);

%% test data from berg with Fixed parameters
close all
test_data_from_berg
param.emax = emax_fun(CN0);
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);
init_guess = [0.00109	0.1	16.4018];

[ocp,sol] =  opt_con(param,g,init_guess,obs_data.tobs(end)*n);
[r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);

makeplot_state_space(sol, '', param, obs_data, g, figure)
%% Berg fitting
close all
ix = 0;
fig = figure;
par_berg = [];
for i = 1:length(datacode)
    id = find(berg_data.DatasetCode == i);
    dataset = berg_data(id(1):id(end), :);
    spcode = unique(dataset.SpeciesCode);
    for j = 1:length(spcode)
        ix = ix + 1;
        sample_berg
        param.CO_0 = obs_data.Co_obs(1);
        param.CT_0 = obs_data.Ct_obs(1);
        param.emax = emax_fun(CN0);

        param.lb = [0.0005, param.mo, 1];
        param.ub = [0.04, param.mo, 400];
        init_guess = [0.0009, 0.1, 10];
        [par, sol, rmse, rsquare] = find_parameter(obs_data, ...
            param, init_guess, g, fig, @makeplot_state_space, @ysim_state_space, n);
        param.rmse = rmse;
        param.rsquare = rsquare;
        makeplot_state_space(sol, '', param, obs_data, g, fig)

        exportgraphics(fig, fig_path_all+"Berg" +decdata2.Dataset{1}+"_j" +j+".png", 'Resolution', 100)
        par_berg = [par_berg; i, j, par, AISC0, CN0, LN0, rmse, rsquare];
        disp("current simulation is i =" +i+" j=" +j)
    end
end
%%
save_fname = excel_path + "fitted_par_Berg.txt";
save(save_fname, "par_berg", '-ascii', '-double', '-tabs')



