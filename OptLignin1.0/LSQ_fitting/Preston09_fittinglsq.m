clear all
clc
close all

%%

set_up
close all
set_up
study ="Preston09";
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\";
fig_path= "fig\";
% param.lb = [0.0005,0.0001, param.mo, 0.01]; param.ub = [0.04,0.04,
% param.mo, 400];

%% Read_data MAR site data
PA_initial = readtable('..\data/preston_2009.xlsx', 'Sheet', 'PA_initial');
PAdata = readtable('..\data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'A1:J60');
PAnmr = readtable('..\data/preston_2009.xlsx', 'Sheet', 'NMR75', 'Range', 'A2:K35');
sites = {'MAR'};
PA = PAdata(strcmpi(PAdata.site, sites{1}), :);

massLoss = readtable('..\data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'M2:X9');
massLoss = removevars(massLoss, {'Cdc'}); % because for cdc lignin at t=2 N/A
massLoss = removevars(massLoss, {'Gfh'}); % excluding grass

spcode = massLoss.Properties.VariableNames(2:end);

%% Preston09 fitting
close all
fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];
par_ = [];
for i = 1:length(spcode)
    sample_Preston2009
    param.emax = emax_fun(init_CN);
    init_guess = [0.003	0.003 0.2	50];
    CO_N= max(obs_data.Co_obs./obs_data.Co_obs(1));

    if(CO_N>1.8)
        mo_f=3;
        disp("mo_f=2")
    else
        mo_f = 1;
    end
        param.lb = [0.0001,0.0001, param.mo*mo_f, 0.1];
        param.ub = [0.04,0.04, param.mo*mo_f, 400];

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

    LC0 = AIS_C(1) / Ct_obs(1);
    par_ = [par_; par, LC0, init_CN, LC0 * init_CN,...
        r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau]; % 0=MAR site
    disp("current simulation is i =" +i)
    %     if(r2_co<0)
    %         error("rsquare negative")
    %     end
end

final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);

writetable(final_table,excel_path+"data summary_LSQ.xlsx")

set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)

%% 'CBR','NH1','TER'
study = "Preston09_b";
fig = figure;
fig.Position = [100, 100, 800, 800];
tiledlayout('flow');
tau = [];
Ctau=[];
sites = {'CBR', 'NH1', 'TER'};
massLoss = readtable('..\data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'AC2:AF23');

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
        param.emax = emax_fun(init_CN);

        CO_N= max(obs_data.Co_obs./obs_data.Co_obs(1));

        if(CO_N>1.8)
            mo_f=3;
            disp("mo_f=2")
        else
            mo_f = 1;
        end
        param.lb = [0.0001,0.0001, param.mo*mo_f, 0.1];
        param.ub = [0.04,0.04, param.mo*mo_f, 400];

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

        LC0 = AIS_C(1) / Ct_obs(1);
        par_ = [par_; par, LC0, init_CN, LC0 * init_CN,...
            r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau]; %
        disp("current simulation is i =" +i)
    end
end
final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
id = strcmp(final_table.study_name, study);
final_table(id, 9:26) = array2table(par_);

writetable(final_table,excel_path+"data summary_LSQ.xlsx")

set(gcf,'color','w')
exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)

%% NMR
% 
% study = "Preston09NMR";
% 
% massLoss = readtable('..\data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'M2:X9');
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
% fig.Position = [100, 100, 800, 800];
% tiledlayout('flow');
% tau = [];
% Ctau=[];
% 
% par_ = [];
% for i = 1:length(spcode)
%     sample_Preston2009_NMR
%     param.emax = emax_fun(init_CN);
% 
%     init_guess = [0.003	0.003 0.2	50];
%     [par, sol] = find_parameter_lsq(obs_data, param, ...
%         init_guess, @ysim_state_space_lsq)  ;
%     [r2, rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT]= est_r2_rmse(obs_data,sol);
% 
%     nexttile
%     fit_plotting(sol, obs_data);lh = legend('show');lh.Box = "off";
%     title("i=" +i)
%     CT=sol.y(1,:);Co=sol.y(2,:);
% 
%     dcodt = diff(sol.y(2,:))./diff(sol.x);
%     neg_idx  = find(dcodt<0, 1);
%     tau=sol.x(neg_idx);
%     CTN=CT./CT(1); Ctau=CTN(neg_idx);
% 
%     LC0 = lignin_C(1) / Ct_obs(1);
%     par_ = [par_; par, LC0, init_CN, LC0 * init_CN,...
%         r2,rmse, r2_co, rmse_co, r2_CT,rmse_CT,AIC,AIC_Co,AIC_CT,tau,Ctau]; %
%     disp("current simulation is i =" +i)
% 
% end
% 
% final_table = readtable(excel_path+ "data summary_LSQ.xlsx");
% id = strcmp(final_table.study_name, study);
% final_table(id, 9:26) = array2table(par_);
% writetable(final_table,excel_path+"data summary_LSQ.xlsx")
% set(gcf,'color','w')
% exportgraphics(fig, fig_path+study+".png", 'Resolution', 300)

