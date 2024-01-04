clear all
clc
close all

%%

set_up
disp("this is the g-function with " +scenario_name+" scenario")
excel_path = "est_params\" +scenario_name + "\";
fig_path_all = "fig\" +scenario_name + "\";


%% Read_data MAR site data
PA_initial = readtable('data/preston_2009.xlsx', 'Sheet', 'PA_initial');
PAdata = readtable('data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'A1:J60');
% PAnmr = readtable('data/preston_2009.xlsx', 'Sheet', 'NMR75', 'Range', 'A2:K35');
sites = {'MAR'};
PA = PAdata(strcmpi(PAdata.site, sites{1}), :);

massLoss = readtable('data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'M2:X9');
massLoss = removevars(massLoss, {'Cdc'}); % because for cdc lignin at t=2 N/A
massLoss = removevars(massLoss, {'Gfh'}); % excluding grass

%%
% close all
spcode = massLoss.Properties.VariableNames(2:end);

% figure;
% tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
% for i= 1:length(spcode)
%     sample_Preston2009
% end

%%
% i=9;
% sample_Preston2009
% fig=figure;
% param.vomax=vomax;
% param.emax = emax_fun(init_CN);
% param.CO_0=obs_data.Co_obs(1);
% param.CT_0=obs_data.Ct_obs(1);
% init_guess = [0.001, 0.5,100]; % [vh_max, mo,ro];
% [ocp,sol] =  opt_con(param,g,init_guess,obs_data.tobs(end)*n);
% makeplot_state_space(sol, '',param,obs_data,g,fig)

%% Preston09 fitting
close all
fig = figure;

par_Preston09 = [];
for i = 1:length(spcode)
    sample_Preston2009
    param.emax = emax_fun(init_CN);
    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    CO_N= max(obs_data.Co_obs./obs_data.Co_obs(1));

    if(CO_N>1.8)
        mo_f=3;
        disp("mo_f=2")
    else
        mo_f = 1;
    end

    param.lb = [0.0001, param.mo*mo_f, 0.1];
    param.ub = [0.04, param.mo*mo_f, 400];

    init_guess = [0.001, 0.5, 100]; % [vh_max, mo,ro];
    [par, sol, ~, ~] = find_parameter(obs_data, param, ...
        init_guess, g, fig, @makeplot_state_space, @ysim_state_space, n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    param.rsquare = r2;
    param.rmse = rmse;
    makeplot_state_space(sol, '', param, obs_data, g, fig)
    exportgraphics(fig, fig_path_all+"Preston09_" +temp{i}+".png", 'Resolution', 100)
    LC0 = AIS_C(1) / Ct_obs(1);
    par_Preston09 = [par_Preston09; 0, i, par, LC0, init_CN,...
        LC0 * init_CN, rmse, r2]; % 0=MAR site
    disp("current simulation is i =" +i+"  "+temp(i))

end
save_fname = excel_path + "fitted_par_Preston09.txt";
save(save_fname, "par_Preston09", '-ascii', '-double', '-tabs')

%% 'CBR','NH1','TER'
close all
sites = {'CBR', 'NH1', 'TER'};
massLoss = readtable('data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'AC2:AF23');

Initial_litter_mass = 1; % grams
tobs = massLoss.t_year_(1:7);
par_ = [];
for i = 1:3
    PA = PAdata(strcmpi(PAdata.site, sites{i}), :);
    mltemp = massLoss(strcmpi(massLoss.site, sites{i}), :);
    temp = massLoss.Properties.VariableNames(3:end);
    for j = 1:2
        initial_fraction_of_C = 0.001 * PA_initial.TotalC_g_kg_1_(strcmpi(PA_initial.CIDETCode, temp{j}));
        init_CN = PA_initial.C_N_atomic_(strcmpi(PA_initial.CIDETCode, temp{j}));

        id = strcmpi(PA.spCode, temp{j});
        aur_C_conc = PA.AURC__gKg_1_(id); % g/kg litter or mg/g litter

        aur_C_conc0 = PA_initial.AUR_C_gKg_1_(strcmpi(PA_initial.CIDETCode, temp{j}));

        ml = table2array(mltemp(:, temp{1}));
        Ct_obs = Initial_litter_mass .* 0.01 .* ml .* initial_fraction_of_C; % g litter/ bag * fraction of litter*fraction of C = (gC/g bag)

        tobs_Co = PA.ElapsedYears(id);

        for jj = 1:length(tobs_Co)
            idtemp(jj) = find(tobs == tobs_Co(jj));
        end
        AIS_C0 = Initial_litter_mass .* (0.01 * ml(1)) .* aur_C_conc0 .* 0.001;
        AIS_C = [AIS_C0; Initial_litter_mass .* 0.01 * ml(idtemp) .* aur_C_conc .* 0.001];

        aromaticC = aromatic_fraction_inAIS(AIS_C); % true lignin


        obs_data = [];
        obs_data.tobs = tobs * 365;
        obs_data.Ct_obs = Ct_obs;
        obs_data.tobs_Co = [0; PA.ElapsedYears(id) * 365];
        obs_data.Co_obs = aromaticC;

        tnorm = obs_data.tobs;
        y = log(obs_data.Ct_obs./obs_data.Ct_obs(1));
        final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
        k = (y(end) - y(1)) ./ tnorm(end);
        terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time) / k;
        n = terminalTime ./ obs_data.tobs(end);
        t = 0:tnorm(end) * n;

%         figure;
%         subplot(211)
%         scatter(obs_data.tobs,obs_data.Ct_obs,'DisplayName','TC'); hold on
%         plot(t,obs_data.Ct_obs(1).*exp(k.*t))
%         scatter(obs_data.tobs_Co,obs_data.Co_obs,'DisplayName','AromaticC'); hold on
%         xlabel("time"); ylabel('gC/bag')
% 
%         subplot(212)
%         plot(1-obs_data.Ct_obs([1,idtemp])./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1),'-o'); hold on
%         xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
%         xlim([0,1]);ylim([0,inf])
%         title([num2str(ii),'-',temp{ii}])
        CO_N= max(obs_data.Co_obs./obs_data.Co_obs(1));

        if(CO_N>1.8)
            mo_f=3;
            disp("mo_f=2")
        else
            mo_f = 1;
        end
        param.lb = [0.0001, param.mo*mo_f, 0.1];
        param.ub = [0.04, param.mo*mo_f, 400];

        param.emax = emax_fun(init_CN);
        param.CO_0 = obs_data.Co_obs(1);
        param.CT_0 = obs_data.Ct_obs(1);
        fig = figure(2);
        init_guess = [-k, 0.5, 1]; % [vh_max, mo,ro];
        [par, sol, rmse, rsquare] = find_parameter(obs_data, param, ...
            init_guess, g, fig, @makeplot_state_space, @ysim_state_space, n);
        param.rmse = rmse;
        param.rsquare = rsquare;
        makeplot_state_space(sol, '', param, obs_data, g, fig)
        exportgraphics(fig, fig_path_all+"Preston09_" +sites{i}+"_" +temp{j}+".png", 'Resolution', 100)
        LC0 = AIS_C(1) / Ct_obs(1);
        par_ = [par_; i, j, par, LC0, init_CN, LC0 * init_CN, rmse, rsquare];
        disp("current simulation is ii =" +j)
    end
end
save_fname = excel_path + "fitted_par_Preston09_b.txt";
save(save_fname, "par_", '-ascii', '-double', '-tabs')

%%
% massLoss = readtable('data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'M2:X9');
% massLoss = removevars(massLoss, {'Cdc'}); % because for cdc lignin at t=2 N/A
% massLoss = removevars(massLoss, {'Gfh'}); % excluding grass
%
% close all
% spcode = massLoss.Properties.VariableNames(2:end);
% for i = 1:length(spcode)
%     sample_Preston2009_NMR
% end
%
% close all
% fig = figure;
% par_Preston09_NMR = [];
% for i = 1:length(spcode)
%     sample_Preston2009_NMR
%     param.emax = emax_fun(init_CN);
%     param.CO_0 = obs_data.Co_obs(1);
%     param.CT_0 = obs_data.Ct_obs(1);
%
%     init_guess = [0.001, 0.5, 100]; % [vh_max, mo,ro];
%     [par, sol, rmse, rsquare] = find_parameter(obs_data, param, ...
%         init_guess, g, fig, @makeplot_state_space, @ysim_state_space, n);
%     param.rmse = rmse;
%     param.rsquare = rsquare;
%     makeplot_state_space(sol, '', param, obs_data, g, fig)
%     exportgraphics(fig, fig_path_all+"Preston09_NMR_" +temp{i}+".png", 'Resolution', 100)
%     LC0 = lignin_C(1) / Ct_obs(1);
%     par_Preston09_NMR = [par_Preston09_NMR; 0, i, par, LC0, init_CN, LC0 * init_CN, rmse, rsquare]; % 0=MAR site
%     disp("current simulation is i =" +i)
%
% end
% save_fname = excel_path + "fitted_par_Preston09_NMR.txt";
% save(save_fname, "par_Preston09_NMR", '-ascii', '-double', '-tabs')
