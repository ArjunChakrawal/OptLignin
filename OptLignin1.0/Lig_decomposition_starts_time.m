fid = fopen('simulation_details.txt', 'a');
fprintf(fid, 'Running Lig_decomposition_starts_time...................');
time = datestr(clock, 'YYYY/mm/dd HH:MM:SS:FFF');
fprintf(fid, '%23s\n', time);

%%
clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Berg";
%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
par_est = load("est_params\" +scenario_name+"\fitted_par_Berg.txt");
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

berg_data = readtable('data\Berg and McClaugherty 1989 - Suppl Mat.xlsx');
datacode = unique(berg_data.DatasetCode);

ix = 0;
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
        par = par_est(ix, 3:5);
        [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
        [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
        stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

        vo = sol.NumericalResults.Control;
        time = sol.NumericalResults.Independent;
        CT = sol.NumericalResults.State(1, :);

        voN = vo ./ max(vo);
        id = find(voN > vo_thres);
        Ligdecomposition_Starts(ix) = time(id(1));
        C_remain_lig_dec_start(ix) = CT(id(1)) / CT(1);
        avgvo(ix) = mean(vo(time < obs_data.tobs(end)));
        init_vo(ix) = vo(1);
        max_vo(ix) = max(vo);
        vo_at_Ligdecomposition_Starts(ix) = vo(id(1));

        %         figure(fig_vo);
        %         nexttile
        %         plot(time, voN, 'linewidth', 2); hold on
        %         plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
        %         xlabel('day'); ylabel('vo/max(vo)')
        %         title("i=" +i+", j=" +j+"(" +decdata2{1, 1}+" " +decdata2{1, 2}+")");
    end
end
%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)

temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'j', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', ...
    'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%%
clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Hirobe2004";
%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
par_est = load("est_params\" +scenario_name+"\fitted_par_Hirobe.txt");
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

Hirobe04_data = readtable('data/Hirobe04.xls', 'Sheet', 'C and Mineral');
organic_hirobe = readtable('data/Hirobe04.xls', 'Sheet', 'Organic');
datacode = unique(organic_hirobe.Species);

for i = 1:length(datacode)
    sample_hirobe
    CN0 = Ct_obs(1) / Nt_obs(1);
    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    param.emax = emax_fun(CN0);
    par = par_est(i, 2:4);
    [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];
    %     makeplot_state_space(sol, '',param,obs_data,g,fig)

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);

    voN = vo ./ max(vo);
    id = find(voN > vo_thres);
    Ligdecomposition_Starts(i) = time(id(1));
    C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
    avgvo(i) = mean(vo(time < obs_data.tobs(end)));
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i) = vo(id(1));

    %     figure(fig_vo);
    %     nexttile
    %     plot(time, voN, 'linewidth', 2); hold on
    %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
    %     xlabel('day'); ylabel('vo/max(vo)')
    %     title("i=" +i+"(" +Hirobe04_data.SpeciesName{i}+")")
end
%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)

temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Osono 2004 Upper
clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Osono2004Upper";
par_est = load("est_params\" +scenario_name+"\fitted_par_Osono04U.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;


treename = {'FC', 'QC', 'AM', 'AR', 'SC', 'BG', 'MJ', 'PH', 'CC', 'CL', 'MO', 'PR', 'AT', 'CJ'};
osono_t = table2array(readtable('data\Osono04.xls', "Sheet", 'Sampling_Dates', "Range", 'C2:C12'));
g_bag = table2array(readtable('data\Osono04.xls', "Sheet", '1)Content', "Range", 'B4:O14'));
mgAUR_gdLitter = table2array(readtable('data\Osono04.xls', "Sheet", '1)Content', "Range", 'B32:O42'));
mgN_gdLitter_U = table2array(readtable('data\Osono04.xls', "Sheet", '1)Content', "Range", 'B144:O154'));

amount_Lig_C0 = (0.001 * mgAUR_gdLitter(1, :))'; % initial fraction of lignin C in litter C
amount_C = g_bag .* initial_fraction_of_C; %gC/ bag
amount_N0 = mgN_gdLitter_U(1, :) .* g_bag(1, :) * 0.001; % g N/ bag
CN0_osono = (initial_fraction_of_C * g_bag(1, :)) ./ (amount_N0);
amount_AUR_C_osono04 = fraction_of_C_in_AUR .* mgAUR_gdLitter .* g_bag * 0.001; % g lignin / bag
sz = size(amount_C);

for i = 1:sz(2)
    sample_Osono04
    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    param.emax = emax_fun(CN0_osono(i));
    par = par_est(i, 2:4);
    [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    voN = vo ./ max(vo);
    id = find(voN > vo_thres);
    Ligdecomposition_Starts(i) = time(id(1));
    C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
    avgvo(i) = mean(vo(time < obs_data.tobs(end)));
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i) = vo(id(1));
    %     figure(fig_vo);
    %     nexttile
    %     plot(time, voN, 'linewidth', 2); hold on
    %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
    %     xlabel('day'); ylabel('vo/max(vo)')
    %     title("i=" +i)

end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)

temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Osono 2004 Lower
clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Osono2004Lower";
par_est = load("est_params\" +scenario_name+"\fitted_par_Osono04L.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

osono_t = table2array(readtable('data\Osono04.xls', "Sheet", 'Sampling_Dates', "Range", 'C2:C12'));
g_bag = table2array(readtable('data\Osono04.xls', "Sheet", '1)Content', "Range", 'B17:O27'));
mgAUR_gdLitter = table2array(readtable('data\Osono04.xls', "Sheet", '1)Content', "Range", 'B45:O55'));
mgN_gdLitter_U = table2array(readtable('data\Osono04.xls', "Sheet", '1)Content', "Range", 'B157:O167'));
amount_Lig_C0 = (0.001 * mgAUR_gdLitter(1, :))'; % initial fraction of lignin C in litter C

amount_Lig_C0 = (0.001 * mgAUR_gdLitter(1, :))'; % initial fraction of lignin C in litter C
amount_C = g_bag .* initial_fraction_of_C; %gC/ bag
amount_N0 = mgN_gdLitter_U(1, :) .* g_bag(1, :) * 0.001; % g N/ bag
CN0_osono = (initial_fraction_of_C * g_bag(1, :)) ./ (amount_N0);
amount_AUR_C_osono04 = fraction_of_C_in_AUR .* mgAUR_gdLitter .* g_bag * 0.001; % g lignin / bag

sz = size(amount_C);
for i = 1:sz(2)
    sample_Osono04
    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    param.emax = emax_fun(CN0_osono(i));
    par = par_est(i, 2:4);
    [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    voN = vo ./ max(vo);
    id = find(voN > vo_thres);
    Ligdecomposition_Starts(i) = time(id(1));
    C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
    avgvo(i) = mean(vo(time < obs_data.tobs(end)));
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i) = vo(id(1));
    %     figure(fig_vo);
    %     nexttile
    %     plot(time, voN, 'linewidth', 2); hold on
    %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
    %     xlabel('day'); ylabel('vo/max(vo)')
    %     title("i=" +i)

end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Osono2017
clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Osono17";
par_est = load("est_params\" +scenario_name+"\fitted_par_Osono17.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

t = readtable('data\Osono17.xlsx', "Sheet", 'Sampling_Dates', "Range", 'B2:B12');
Osono2017 = readtable('data\Osono17.xlsx', 'Sheet', 'Origianal_data');
iddata = find(Osono2017.Collection == 0);
day = [0, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36]' * 30;
treename_Osono2017 = unique(Osono2017.Tree);

numSpecies = length(treename_Osono2017);
for i = 1:numSpecies
    sample_Osono17
    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    CN0 = amount_C(tsrt) / amount_N(tsrt);
    param.emax = emax_fun(CN0);
    par = par_est(i, 2:4);
    [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];
    %     makeplot_state_space(sol, '',param,obs_data,g,figure)
    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    voN = vo ./ max(vo);
    id = find(voN > 0.005);
    Ligdecomposition_Starts(i) = time(id(1));
    C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
    avgvo(i) = mean(vo(time < obs_data.tobs(end)));
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i) = vo(id(1));
    %     figure(fig_vo);
    %     nexttile
    %     plot(time, voN, 'linewidth', 2); hold on
    %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
    %     xlabel('day'); ylabel('vo/max(vo)')
    %     title("i=" +i)
    %     pause
end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Preston09_fitting

clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Preston09";
par_est = load("est_params\" +scenario_name+"\fitted_par_Preston09.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

PA_initial = readtable('data/preston_2009.xlsx', 'Sheet', 'PA_initial');
PAdata = readtable('data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'A1:J60');
PAnmr = readtable('data/preston_2009.xlsx', 'Sheet', 'NMR75', 'Range', 'A2:K35');

sites = {'MAR'};
PA = PAdata(strcmpi(PAdata.site, sites{1}), :);
massLoss = readtable('data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'M2:X9');
massLoss = removevars(massLoss, {'Cdc'}); % because for cdc lignin at t=2 N/A
massLoss = removevars(massLoss, {'Gfh'}); % excluding grass
spcode = massLoss.Properties.VariableNames(2:end);

for i = 1:length(spcode)
    sample_Preston2009
    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    param.emax = emax_fun(init_CN);
    par = par_est(i, 3:5);
    [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    voN = vo ./ max(vo);
    id = find(voN > 0.005);
    Ligdecomposition_Starts(i) = time(id(1));
    C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
    avgvo(i) = mean(vo(time < obs_data.tobs(end)));
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i) = vo(id(1));
    %     figure(fig_vo);
    %     nexttile
    %     plot(time, voN, 'linewidth', 2); hold on
    %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
    %     xlabel('day'); ylabel('vo/max(vo)')
    %     title("i=" +i)
end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'j', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

% %% Preston09_fittingNMR data
% clearvars;
% close all;
% set_up;param.numiter=50;
% stat_par = [];
% excel_path = "est_params\" +scenario_name + "\";
% %fig_pathCombined "fig\" +scenario_name + "\combined\";
% study = "Preston09NMR";
% par_est = load("est_params\" +scenario_name+"\fitted_par_Preston09_NMR.txt");
%
% PA_initial = readtable('data/preston_2009.xlsx', 'Sheet', 'PA_initial');
% PAdata = readtable('data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'A1:J60');
% PAnmr = readtable('data/preston_2009.xlsx', 'Sheet', 'NMR75', 'Range', 'A2:K35');
%
% sites = {'MAR'};
% PA = PAdata(strcmpi(PAdata.site, sites{1}), :);
% massLoss = readtable('data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'M2:X9');
% massLoss = removevars(massLoss, {'Cdc'}); % because for cdc lignin at t=2 N/A
% massLoss = removevars(massLoss, {'Gfh'}); % excluding grass
% spcode = massLoss.Properties.VariableNames(2:end);
%
% %fig_vo = figure;
% %fig_vo.Position = [100, 100, 1200, 800];
% %tiledlayout('flow');
% Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
% C_remain_lig_dec_start = Ligdecomposition_Starts;
% avgvo = Ligdecomposition_Starts;
% init_vo = avgvo;
% max_vo = init_vo;
% vo_at_Ligdecomposition_Starts = init_vo;
%
%
% for i = 1:length(spcode)
%     sample_Preston2009_NMR
%     param.CO_0 = obs_data.Co_obs(1);
%     param.CT_0 = obs_data.Ct_obs(1);
%     param.emax = emax_fun(init_CN);
%     par = par_est(i, 3:5);
%     [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
%     [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
%     stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];
%
%     vo = sol.NumericalResults.Control;
%     time = sol.NumericalResults.Independent;
%     CT = sol.NumericalResults.State(1, :);
%     voN = vo ./ max(vo);
%     id = find(voN > 0.005);
%     Ligdecomposition_Starts(i) = time(id(1));
%     C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
%     avgvo(i) = mean(vo(time < obs_data.tobs(end)));
%     init_vo(i) = vo(1);
%     max_vo(i) = max(vo);
%     vo_at_Ligdecomposition_Starts(i) = vo(id(1));
%     %     figure(fig_vo);
%     %     nexttile
%     %     plot(time, voN, 'linewidth', 2); hold on
%     %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
%     %     xlabel('day'); ylabel('vo/max(vo)')
%     %     title("i=" +i)
% end
%
% %exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
% temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
%     vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
% par_sav = array2table(temp_par, ...
%     'VariableNames', {'i', 'j', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
%     'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
%     'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
% writetable(par_sav, excel_path+study+"_final.xlsx")

%% Preston 2009 other site
clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Preston09_b";
par_est = load("est_params\" +scenario_name+"\fitted_par_Preston09_b.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

PA_initial = readtable('data/preston_2009.xlsx', 'Sheet', 'PA_initial');
PAdata = readtable('data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'A1:J60');

sites = {'CBR', 'NH1', 'TER'};
massLoss = readtable('data/preston_2009.xlsx', 'Sheet', 'PA', 'Range', 'AC2:AF23');

Initial_litter_mass = 1; % grams
tobs = massLoss.t_year_(1:7);
ix = 0;
for i = 1:3
    PA = PAdata(strcmpi(PAdata.site, sites{i}), :);
    mltemp = massLoss(strcmpi(massLoss.site, sites{i}), :);
    temp = massLoss.Properties.VariableNames(3:end);
    for j = 1:2
        ix = ix + 1;
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

        param.CO_0 = obs_data.Co_obs(1);
        param.CT_0 = obs_data.Ct_obs(1);
        param.emax = emax_fun(init_CN);
        par = par_est(ix, 3:5);
        [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
        [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
        stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

        vo = sol.NumericalResults.Control;
        time = sol.NumericalResults.Independent;
        CT = sol.NumericalResults.State(1, :);

        voN = vo ./ max(vo);
        id = find(voN > vo_thres);
        Ligdecomposition_Starts(ix) = time(id(1));
        C_remain_lig_dec_start(ix) = CT(id(1)) / CT(1);
        avgvo(ix) = mean(vo(time < obs_data.tobs(end)));
        init_vo(ix) = vo(1);
        max_vo(ix) = max(vo);
        vo_at_Ligdecomposition_Starts(ix) = vo(id(1));

        %         figure(fig_vo);
        %         nexttile
        %         plot(time, voN, 'linewidth', 2); hold on
        %         plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
        %         xlabel('day'); ylabel('vo/max(vo)')
        %         title("i=" +i+", j=" +j);
    end
end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'j', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% He_2019
clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "He2019";
par_est = load("est_params\" +scenario_name+"\fitted_par_He2019.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

day = table2array(readtable('data/He et al 2019.xlsx', 'Range', 'H3:H9')); %gC litter
C_remaining = table2array(readtable('data/He et al 2019.xlsx', 'Range', 'M3:P9'));
lignin_C_remaining = table2array(readtable('data/He et al 2019.xlsx', 'Range', 'M12:P18'));
init_CN = [37.10, 37.1, 32.7, 32.7];
SP = readtable('data/He et al 2019.xlsx', 'Range', 'I2:L2').Properties.VariableNames;

for i = 1:length(SP)
    param.emax = emax_fun(init_CN(i));
    obs_data = [];
    obs_data.tobs = day;
    obs_data.Ct_obs = C_remaining(:, i);
    obs_data.Co_obs = lignin_C_remaining(:, i);

    tnorm = obs_data.tobs;
    f = fit(obs_data.tobs, obs_data.Ct_obs, 'exp1');
    final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
    k = f.b;
    terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time) / k;
    n = terminalTime ./ obs_data.tobs(end);
    tt = 0:tnorm(end) * n;

    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);

    par = par_est(i, 2:4);
    [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    voN = vo ./ max(vo);
    id = find(voN > 0.005);
    Ligdecomposition_Starts(i) = time(id(1));
    C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
    avgvo(i) = mean(vo(time < obs_data.tobs(end)));
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i) = vo(id(1));
    %     figure(fig_vo);
    %     nexttile
    %     plot(time, voN, 'linewidth', 2); hold on
    %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
    %     xlabel('day'); ylabel('vo/max(vo)')
    %     title("i=" +i)
end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Ruzek2021

clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Ruzek2021";
par_est = load("est_params\" +scenario_name+"\fitted_par_Ruzek_Rooibos.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;
Cday = table2array(readtable('data/Ruzek et al 2021.xlsx', 'sheet', ...
    'processed_data', 'Range', 'B6:B10'));
C_remaining = table2array(readtable('data/Ruzek et al 2021.xlsx', 'sheet', ...
    'processed_data', 'Range', 'G6:J10'));
Lday = table2array(readtable('data/Ruzek et al 2021.xlsx', 'sheet', ...
    'processed_data', 'Range', 'L6:L9'));
lignin_C_remaining = table2array(readtable('data/Ruzek et al 2021.xlsx', 'sheet', ...
    'processed_data', 'Range', 'Q6:T9'));
SP = readtable('data/Ruzek et al 2021.xlsx', 'sheet', ...
    'processed_data', 'Range', 'M11:P11').Properties.VariableNames;
init_CN = 53.2;


for i = 1:size(C_remaining, 2)
    param.emax = emax_fun(init_CN);
    obs_data = [];
    obs_data.tobs = Cday;
    obs_data.Ct_obs = C_remaining(:, i);
    obs_data.tobs_Co = Lday;
    obs_data.Co_obs = lignin_C_remaining(:, i);
    idtemp = [];
    for j = 1:length(Lday)
        idtemp(j) = find(Cday == Lday(j));
    end
    tnorm = obs_data.tobs;
    f = fit(obs_data.tobs, obs_data.Ct_obs, 'exp1');
    final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
    k = f.b;
    terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time) / k;
    n = terminalTime ./ obs_data.tobs(end);
    tt = 0:tnorm(end) * n;

    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);

    par = par_est(i, 2:4);
    [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    voN = vo ./ max(vo);
    id = find(voN > 0.005);
    Ligdecomposition_Starts(i) = time(id(1));
    C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
    avgvo(i) = mean(vo(time < obs_data.tobs(end)));
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i) = vo(id(1));
    %     figure(fig_vo);
    %     nexttile
    %     plot(time, voN, 'linewidth', 2); hold on
    %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
    %     xlabel('day'); ylabel('vo/max(vo)')
    %     title("i=" +i)
end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Herzog

% clearvars;
% close all;
% set_up;stat_par=[];
% excel_path = "est_params\" +scenario_name + "\";
% %fig_pathCombined "fig\" +scenario_name + "\combined\";
% study = "Herzog2019";
% par_est = load("est_params\" +scenario_name+"\fitted_par_Herzog.txt");
%
% %fig_vo = figure;
% %fig_vo.Position = [100, 100, 1200, 800];
% %tiledlayout('flow');
% Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
% C_remain_lig_dec_start = Ligdecomposition_Starts;
% avgvo = Ligdecomposition_Starts;
% init_vo = avgvo;
% max_vo = init_vo;
% vo_at_Ligdecomposition_Starts = init_vo;
%
% day = table2array(readtable('data/Herzog et al 2019.xlsx', 'Range', 'B36:B41')); % day
% C_remaining = table2array(readtable('data/Herzog et al 2019.xlsx', 'Range', 'K36:L41')); % day
% lignin_C_remaining = table2array(readtable('data/Herzog et al 2019.xlsx', 'Range', 'M36:N41')); % day
% init_CN = [57.7, 69.5]; % [Control  irrigated]
% SP = {'Control', 'irrigated'};
%
% for i = 1:2
%     param.emax = emax_fun(init_CN(i));
%     obs_data = [];
%     obs_data.tobs = day;
%     obs_data.Ct_obs = C_remaining(:, i);
%     obs_data.Co_obs = lignin_C_remaining(:, i); % lignin C as VSC
%
%     tnorm = obs_data.tobs;
%     f = fit(obs_data.tobs, obs_data.Ct_obs, 'exp1');
%     final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
%     k = f.b;
%     terminalTime = log(final_C*0.75) / k;
%     n = terminalTime ./ obs_data.tobs(end);
%     tt = 0:tnorm(end) * n;
%
%     param.CO_0 = obs_data.Co_obs(1);
%     param.CT_0 = obs_data.Ct_obs(1);
%
%     par = par_est(i, 2:4);
%     [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);[r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT,AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);stat_par = [stat_par;[ r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT,AIC,AIC_Co,AIC_CT]];
%
%     vo = sol.NumericalResults.Control;
%     time = sol.NumericalResults.Independent;
%     CT = sol.NumericalResults.State(1, :);
%     voN = vo ./ max(vo);
%     id = find(voN > 0.005);
%     Ligdecomposition_Starts(i) = time(id(1));
%     C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
%     avgvo(i) = mean(vo(time<obs_data.tobs(end)));
%     init_vo(i) = vo(1);
%     max_vo(i) = max(vo);
%     vo_at_Ligdecomposition_Starts(i) = vo(id(1));
%     figure(fig_vo);
%     nexttile
%     plot(time, voN, 'linewidth', 2); hold on
%     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
%     xlabel('day'); ylabel('vo/max(vo)')
%     title("i=" +i)
% end
%
% %exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
% temp_par = [par_est(:,1:end-2),stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
%     vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
% par_sav = array2table(temp_par, ...
%     'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
%     'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT','AIC', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
%     'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
% writetable(par_sav, excel_path+study+"_final.xlsx")

%% Zhou2017

clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Zhou2017";
par_est = load("est_params\" +scenario_name+"\fitted_par_Zhou2017.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

data = readtable('data/Zhou etal 2017.xlsx', 'Range', 'A4:G56');
CN0 = 51.68;
AISC0 = 0.219;
treatementName = unique(data.treatment);

for i = 1:length(treatementName)
    temp = data(data.treatment == treatementName(i), :);
    obs_data = [];
    obs_data.tobs = temp.day;
    obs_data.Ct_obs = temp.CG;
    obs_data.Co_obs = aromatic_fraction_inAIS(temp.LigninCG);
    tnorm = obs_data.tobs;
    y = log(obs_data.Ct_obs./obs_data.Ct_obs(1));
    final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
    k = (y(end) - y(1)) ./ tnorm(end);
    terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time) / k;
    n = terminalTime ./ obs_data.tobs(end);
    t = 0:tnorm(end) * n;

    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    param.emax = emax_fun(CN0);

    par = par_est(i, 2:4);
    [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    voN = vo ./ max(vo);
    id = find(voN > 0.005);
    Ligdecomposition_Starts(i) = time(id(1));
    C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
    avgvo(i) = mean(vo(time < obs_data.tobs(end)));
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i) = vo(id(1));
    %     figure(fig_vo);
    %     nexttile
    %     plot(time, voN, 'linewidth', 2); hold on
    %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
    %     xlabel('day'); ylabel('vo/max(vo)')
    %     title("i=" +i)
end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Huang2021

clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Huang2021";
par_est = load("est_params\" +scenario_name+"\fitted_par_Huang2021.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

data = readtable('data/Huang2021.xlsx', 'Range', 'A10:Q38');
CN0 = 29;
AISC0 = 0.219;
treatementName = unique(data.treatment);

for i = 1:length(treatementName)
    temp = data(data.treatment == treatementName(i), :);
    obs_data = [];
    obs_data.tobs = temp.day;
    obs_data.Ct_obs = temp.CInG;
    obs_data.Co_obs = aromatic_fraction_inAIS(temp.ligninCInG);
    tnorm = obs_data.tobs;
    y = log(obs_data.Ct_obs./obs_data.Ct_obs(1));
    final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
    k = (y(end) - y(1)) ./ tnorm(end);
    terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time) / k;
    n = terminalTime ./ obs_data.tobs(end);
    t = 0:tnorm(end) * n;

    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    param.emax = emax_fun(CN0);

    par = par_est(i, 2:4);
    [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    voN = vo ./ max(vo);
    id = find(voN > 0.005);
    Ligdecomposition_Starts(i) = time(id(1));
    C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
    avgvo(i) = mean(vo(time < obs_data.tobs(end)));
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i) = vo(id(1));
    %     figure(fig_vo);
    %     nexttile
    %     plot(time, voN, 'linewidth', 2); hold on
    %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
    %     xlabel('day'); ylabel('vo/max(vo)')
    %     title("i=" +i)
end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Kou2015

clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Kou2015";
par_est = load("est_params\" +scenario_name+"\fitted_par_Kou2015.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

data = readtable('data/Kou2015.xlsx', 'Range', 'A4:M58');
Littertype = unique(data.type);
treatementName = unique(data.treatment);

ix = 0;
for i = 1:length(Littertype)
    temp = data(data.type == Littertype(i), :);
    for j = 1:length(treatementName)
        ix = ix + 1;
        temp2 = temp(temp.treatment == treatementName(j), :);
        obs_data = [];
        obs_data.tobs = temp2.day;
        obs_data.Ct_obs = temp2.CInG;
        obs_data.Co_obs = aromatic_fraction_inAIS(temp2.LigninCInG);
        CN0 = temp2.CN0(1);
        tnorm = obs_data.tobs;
        y = log(obs_data.Ct_obs./obs_data.Ct_obs(1));
        final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
        k = (y(end) - y(1)) ./ tnorm(end);
        terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time) / k;
        n = terminalTime ./ obs_data.tobs(end);
        t = 0:tnorm(end) * n;

        param.CO_0 = obs_data.Co_obs(1);
        param.CT_0 = obs_data.Ct_obs(1);
        param.emax = emax_fun(CN0);

        par = par_est(ix, 3:5);
        [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
        [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
        stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

        vo = sol.NumericalResults.Control;
        time = sol.NumericalResults.Independent;
        CT = sol.NumericalResults.State(1, :);

        voN = vo ./ max(vo);
        id = find(voN > vo_thres);
        Ligdecomposition_Starts(ix) = time(id(1));
        C_remain_lig_dec_start(ix) = CT(id(1)) / CT(1);
        avgvo(ix) = mean(vo(time < obs_data.tobs(end)));
        init_vo(ix) = vo(1);
        max_vo(ix) = max(vo);
        vo_at_Ligdecomposition_Starts(ix) = vo(id(1));

        %         figure(fig_vo);
        %         nexttile
        %         plot(time, voN, 'linewidth', 2); hold on
        %         plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
        %         xlabel('day'); ylabel('vo/max(vo)')
        %         title("i=" +i+", j=" +j);
    end
end
%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'j', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Yue et al 2016

clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Yue2016";
par_est = load("est_params\" +scenario_name+"\fitted_par_Yue.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

day = table2array(readtable('data/Yue et al 2016.xlsx', 'Range', 'A30:A35')); %gC litter
C_remaining = table2array(readtable('data/Yue et al 2016.xlsx', 'Range', 'J30:M35')); %gC litter
lignin_C_remaining = table2array(readtable('data/Yue et al 2016.xlsx', 'Range', 'N30:Q35')); %gC lignin
init_C = table2array(readtable('data/Yue et al 2016.xlsx', 'Range', 'B23:E23'));
init_CN = table2array(readtable('data/Yue et al 2016.xlsx', 'Range', 'B25:E25'));
SP = readtable('data/Yue et al 2016.xlsx', 'Range', 'B22:E22').Properties.VariableNames;


for i = 1:length(SP)
    param.emax = emax_fun(init_CN(i));
    obs_data = [];
    obs_data.tobs = day;
    obs_data.Ct_obs = C_remaining(:, i);
    obs_data.Co_obs = aromatic_fraction_inAIS(lignin_C_remaining(:, i));

    tnorm = obs_data.tobs;
    f = fit(obs_data.tobs, obs_data.Ct_obs, 'exp1');
    final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
    k = f.b;
    terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time) / k;
    n = terminalTime ./ obs_data.tobs(end);
    tt = 0:tnorm(end) * n;

    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    param.emax = emax_fun(init_CN(i));
    if (i == 1)
        n = 1;
    end
    par = par_est(i, 2:4);
    [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
    [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
    stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);

    voN = vo ./ max(vo);
    id = find(voN > vo_thres);
    Ligdecomposition_Starts(i) = time(id(1));
    C_remain_lig_dec_start(i) = CT(id(1)) / CT(1);
    avgvo(i) = mean(vo(time < obs_data.tobs(end)));
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i) = vo(id(1));
    %     figure(fig_vo);
    %     nexttile
    %     plot(time, voN, 'linewidth', 2); hold on
    %     plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
    %     xlabel('day'); ylabel('vo/max(vo)')
    %     title("i=" +i)
end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Tu2011
clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Tu2011";
par_est = load("est_params\" +scenario_name+"\fitted_par_Tu2011.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

SP = readtable('data/Tu et al 2011.xlsx', 'Range', 'A10:A33');
obsday = table2array(readtable('data/Tu et al 2011.xlsx', 'Range', 'B10:B33'));
C_remaining = readtable('data/Tu et al 2011.xlsx', 'Range', 'P10:S33'); %gC litter
lignin_C_remaining = readtable('data/Tu et al 2011.xlsx', 'Range', 'T10:W33');
init_chem = readtable('data/Tu et al 2011.xlsx', 'Range', 'A3:J6'); %gC litter
treatementName = C_remaining.Properties.VariableNames;
SPname = unique(SP.Type);
ix = 0;
for i = 1:3
    idd = find(strcmpi(SP.Type, SPname{i}));
    for j = 1:4
        ix = ix + 1;
        init_CN = init_chem.initC_N(strcmpi(init_chem.Type, SPname{i}));
        param.emax = emax_fun(init_CN);
        obs_data = [];
        obs_data.tobs = obsday(idd);
        obs_data.Ct_obs = C_remaining(idd, j).Variables;
        obs_data.Co_obs = aromatic_fraction_inAIS(lignin_C_remaining(idd, j).Variables);

        tnorm = obs_data.tobs;
        f = fit(obs_data.tobs, obs_data.Ct_obs, 'exp1');
        final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
        k = f.b;
        terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time) / k;
        n = terminalTime ./ obs_data.tobs(end);
        tt = 0:tnorm(end) * n;

        param.CO_0 = obs_data.Co_obs(1);
        param.CT_0 = obs_data.Ct_obs(1);
        par = par_est(ix, 3:5);
        [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
        [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
        stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

        vo = sol.NumericalResults.Control;
        time = sol.NumericalResults.Independent;
        CT = sol.NumericalResults.State(1, :);

        voN = vo ./ max(vo);
        id = find(voN > vo_thres);
        Ligdecomposition_Starts(ix) = time(id(1));
        C_remain_lig_dec_start(ix) = CT(id(1)) / CT(1);
        avgvo(ix) = mean(vo(time < obs_data.tobs(end)));
        init_vo(ix) = vo(1);
        max_vo(ix) = max(vo);
        vo_at_Ligdecomposition_Starts(ix) = vo(id(1));

        %         figure(fig_vo);
        %         nexttile
        %         plot(time, voN, 'linewidth', 2); hold on
        %         plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
        %         xlabel('day'); ylabel('vo/max(vo)')
        %         title("i=" +i+", j=" +j);
    end
end
%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'j', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Tu2014

clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "Tu2014";
par_est = load("est_params\" +scenario_name+"\fitted_par_Tu2014.txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1200, 800];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

data = readtable('data/Tu_et_al2014.xlsx', 'Range', 'B16:O360');
CN0 = table2array(readtable('data/Tu_et_al2014.xlsx', 'Range', 'K2:K11'));
AISC0 = table2array(readtable('data/Tu_et_al2014.xlsx', 'Range', 'N2:N11'));
SP = readtable('data/Tu_et_al2014.xlsx', 'Range', 'A1:A11');

treatementName = unique(data.treatment);
SPname = unique(data.SpeciesCode);
ix = 0;
for i = 1:length(SPname)
    tempdata1 = data(data.SpeciesCode == SPname(i), :);
    for j = 1:length(treatementName)
        ix = ix + 1;
        tempdata2 = tempdata1(tempdata1.treatment == treatementName(j), :);
        C_2 = tempdata2.C_C_0_(2);
        N_2 = tempdata2.N_N_0_(2);
        
        if (N_2 < 0.9 && N_2 < C_2 || C_2 < 0.7)
            tsrt = 2;
            temp = 1;
        else
            tsrt = 1;
            temp = 0;
        end

        obs_data=[];
        obs_data.tobs  = tempdata2.days(tsrt:end)-tempdata2.days(tsrt);
        obs_data.Ct_obs = tempdata2.CInG(tsrt:end);
        obs_data.Co_obs  = aromatic_fraction_inAIS(tempdata2.ligninCInG(tsrt:end));
        CT_N = obs_data.Ct_obs./obs_data.Ct_obs(1);
        iddiff = [false;diff(CT_N)./CT_N(1:end-1) >0.1];
        obs_data.tobs(iddiff)=[];obs_data.Ct_obs(iddiff)=[];obs_data.Co_obs(iddiff)=[];

        tnorm = obs_data.tobs;
        f = fit(obs_data.tobs, obs_data.Ct_obs, 'exp1');
        final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
        k = f.b;
        terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time) / k;
        n = terminalTime ./ obs_data.tobs(end);
        tt = 0:tnorm(end) * n;

        param.CO_0 = obs_data.Co_obs(1);
        param.CT_0 = obs_data.Ct_obs(1);
        param.emax = emax_fun(CN0(i));
        par = par_est(ix, 3:5);
        [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
        [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
        stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];

        vo = sol.NumericalResults.Control;
        time = sol.NumericalResults.Independent;
        CT = sol.NumericalResults.State(1, :);

        voN = vo ./ max(vo);
        id = find(voN > vo_thres);
        Ligdecomposition_Starts(ix) = time(id(1));
        C_remain_lig_dec_start(ix) = CT(id(1)) / CT(1);
        avgvo(ix) = mean(vo(time < obs_data.tobs(end)));
        init_vo(ix) = vo(1);
        max_vo(ix) = max(vo);
        vo_at_Ligdecomposition_Starts(ix) = vo(id(1));

        %         figure(fig_vo);
        %         nexttile
        %         plot(time, voN, 'linewidth', 2); hold on
        %         plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
        %         xlabel('day'); ylabel('vo/max(vo)')
        %         title("i=" +i+", j=" +j);
    end
end
%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)
temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'j', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% MagilandAber1998
clearvars;
close all;
set_up;param.numiter=50;
stat_par = [];
fig_path_all = "fig\" +scenario_name + "\";

excel_path = "est_params\" +scenario_name + "\";
%fig_pathCombined "fig\" +scenario_name + "\combined\";
study = "MagilandAber1998";
par_est = load(excel_path+"fitted_par_" +study+".txt");

%fig_vo = figure;
%fig_vo.Position = [100, 100, 1000, 600];
%tiledlayout('flow');
Ligdecomposition_Starts = zeros(size(par_est, 1), 1);
C_remain_lig_dec_start = Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo = avgvo;
max_vo = init_vo;
vo_at_Ligdecomposition_Starts = init_vo;

data = readtable('data/Magil and Aber 1998.xlsx', 'Range', 'A8:Q242');

ix = 0;
for i = 1:2
    sitedata = data(data.SiteCode == i, :);
    spcode = unique(sitedata.SpeciesCode);
    for j = 1:length(spcode)
        spdata = sitedata(sitedata.SpeciesCode == spcode(j), :);
        treatmentcode = unique(spdata.treatmentCtrl_0_LN_1_HN_2);
        for k = 1:length(treatmentcode)
            ix = ix + 1;
            litterdata = spdata(spdata.treatmentCtrl_0_LN_1_HN_2 == treatmentcode(k), :);
            C_2 = 0.01 * litterdata.M_M_0_(2);
            if (C_2 < 0.7)
                tsrt = 2;
                temp = 1;
            else
                tsrt = 1;
                temp = 0;
            end
            AISC0 = litterdata.ligninCInG(tsrt)./litterdata.CInG(tsrt);
            CN0 = litterdata.CN0(tsrt);
            obs_data=[];
            obs_data.tobs  = (litterdata.month(tsrt:end) -litterdata.month(tsrt)) *30;
            obs_data.Ct_obs = litterdata.CInG(tsrt:end);
            obs_data.Co_obs  = aromatic_fraction_inAIS(litterdata.ligninCInG(tsrt:end));

            CT_N = obs_data.Ct_obs./obs_data.Ct_obs(1);
            iddiff = [false;diff(CT_N)./CT_N(1:end-1) >0.1];
            obs_data.tobs(iddiff)=[];obs_data.Ct_obs(iddiff)=[];obs_data.Co_obs(iddiff)=[];

            tnorm = obs_data.tobs;
            f = fit(obs_data.tobs, obs_data.Ct_obs./obs_data.Ct_obs(1), 'exp1');
            final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
            decayk = f.b;
            terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time) / decayk;
            n = terminalTime ./ obs_data.tobs(end);
            tt = 0:tnorm(end) * n;
            param.CO_0 = obs_data.Co_obs(1);
            param.CT_0 = obs_data.Ct_obs(1);
            param.emax = emax_fun(CN0);

            par = par_est(ix, 3:5);
            [~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
            [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT] = est_r2_rmse(obs_data, sol);
            stat_par = [stat_par; [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT, AIC,AIC_Co,AIC_CT]];
            param.rsquare = r2;
            param.rmse = rmse;
            fig = figure(2);
            clf;
            makeplot_state_space(sol, '', param, obs_data, g, fig)
            exportgraphics(fig, fig_path_all+study+"_" +sitedata.Site(1)+"_" + ...
                spdata.Var3(1)+"_k" +k+".png", 'Resolution', 100)

            vo = sol.NumericalResults.Control;
            time = sol.NumericalResults.Independent;
            CT = sol.NumericalResults.State(1, :);

            voN = vo ./ max(vo);
            id = find(voN > vo_thres);
            Ligdecomposition_Starts(ix) = time(id(1));
            C_remain_lig_dec_start(ix) = CT(id(1)) / CT(1);
            avgvo(ix) = mean(vo(time < obs_data.tobs(end)));
            init_vo(ix) = vo(1);
            max_vo(ix) = max(vo);
            vo_at_Ligdecomposition_Starts(ix) = vo(id(1));

            %             figure(fig_vo);
            %             nexttile
            %             plot(time, voN, 'linewidth', 2); hold on
            %             plot([1, 1]*time(id(1)), [0, 1], '--k', 'linewidth', 2)
            %             xlabel('day'); ylabel('vo/max(vo)')
            %             title("i" +i+" j" +j+" k" +k)

        end
    end
end

%exportgraphics(figure(fig_vo), fig_pathCombined+study+"_vo.png", 'Resolution', 300)

temp_par = [par_est(:, 1:end-2), stat_par, Ligdecomposition_Starts, avgvo, init_vo, max_vo, ...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];

par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'j', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'r2', 'rmse', 'r2_co', 'rmse_co', 'r2_CT', 'rmse_CT', 'AIC','AIC_Co','AIC_CT', 'LigDecStartDay', 'avg_vo', 'init_vo', 'max_vo', ...
    'vo_at_Ligdecomposition_Starts', 'C_remain_lig_dec_start'});
writetable(par_sav, excel_path+study+"_final.xlsx")

%% Snajdr 2011

% Snajdr_fitting;

%% McKee2016

% McKee2016_fitting

%% He_2016

% He_2016_fitting

%% Fioretto 2005

% Fioretto_2005_fitting

%%
fid = fopen('simulation_details.txt', 'a');
fprintf('Finished Lig_decomposition_starts_time...................\n');
time = datestr(clock, 'YYYY/mm/dd HH:MM:SS:FFF');
fprintf(fid, '%23s\n', time);

%%
fid = fopen('simulation_details.txt', 'a');
fprintf(fid, 'Running feed_excel ...................');
time = datestr(clock, 'YYYY/mm/dd HH:MM:SS:FFF');
fprintf(fid, '%23s\n', time);
feed_excel
fid = fopen('simulation_details.txt', 'a');
fprintf(fid, 'Finished feed_excel...................');
time = datestr(clock, 'YYYY/mm/dd HH:MM:SS:FFF');
fprintf(fid, '%23s\n', time);
