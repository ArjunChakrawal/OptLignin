clear all
clc
close all
%%
set_up
disp("this is the g-function with "+ scenario_name+" scenario")
fig_path_all = "fig\"+scenario_name+"\";
fig_pathCombined = "fig\"+scenario_name+"\combined\";
excel_path = "est_params\"+scenario_name+"\";

%% Read_data
Read_data

%% Berg
parFname = "est_params\"+scenario_name+"\fitted_par_Berg.txt";
par_Berg = load(parFname);
close all
fig1 = figure;
tiledlayout(8, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig1.Position = [552, 50, 859, 1066];
fig2 = figure;
tiledlayout(8, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig3 = figure;
tiledlayout(8, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig4 = figure;
tiledlayout(8, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig5 = figure;
tiledlayout(8, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig6 = figure;
tiledlayout(8, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact');

fig2.Position = fig1.Position;
fig3.Position = fig1.Position;
fig4.Position = fig1.Position;
fig5.Position = fig1.Position;
fig6.Position = fig1.Position;
fig7=figure;
xfig=figure;

ix = 0;
Ligdecomposition_Starts = ones(length(par_Berg), 1);
avgvo = Ligdecomposition_Starts;
min_vo= avgvo;max_vo=min_vo;
C_loss_lig_dec_start=avgvo;

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
        par = par_Berg(ix, 3:5);
        [~,sol] =  opt_con(param,g,par,obs_data.tobs(end).*n);
        %         IC=[param.CT_0, param.CO_0];
        %         [~, sol] = ysim_state_space(par, ocp, obs_data, IC);
        makeplot_state_space(sol, '',param,obs_data,g,xfig)
%         exportgraphics(xfig, fig_path_all+decdata2.Dataset{1}+"_j"+j+".png",'Resolution',100)

        vo = sol.NumericalResults.Control;
        time = sol.NumericalResults.Independent;
        CT = sol.NumericalResults.State(1, :);
        Co = sol.NumericalResults.State(2, :);
        L=Co./CT; AIS = Co*4;L_AIS=AIS./CT;
        e = param.emax - par(3) * vo;
        gL=g(L,a, b);
        vhmax=sol.NumericalResults.Parameter(1)  ;
        Dh = vhmax.*gL.*(CT-Co);
        Do = vo.*Co;

        voN = vo ./ max(vo);
        id = find(voN > 0.01);
        Ligdecomposition_Starts(ix) = time(id(1)) ;
        C_loss_lig_dec_start(ix) = CT(id(1))/CT(1);

        avgvo(ix) = mean(vo);
        title_name="i=" +i+", j=" +j+"(" +decdata2{1, 1}+" " +decdata2{1, 2}+")";
%         post_plots
% 
%         figure(fig6);
%         nexttile
%         plot(L, Do./(Dh+Do),'o','MarkerFaceColor',LC(i,:));
%         ylabel('Do./(Dh+Do)')
%         xlabel('LCI=Co/CT')
%         title(title_name)
% 
%         figure(fig7);
%         scatter(L_AIS, Do./(Dh+Do),'o','MarkerFaceColor',LC(i,:));hold on
%         ylabel('Do./(Dh+Do)')
%         xlabel('LCI=AISC/CT');
    end
end

% figureHandles = findobj('Type', 'figure');
% set(figureHandles,'color','w')
% exportgraphics(figure(fig1), fig_pathCombined+'Berg_state_space.png', 'Resolution', 300)
% exportgraphics(figure(fig2), fig_pathCombined+'Berg_CUE.png', 'Resolution', 300)
% exportgraphics(figure(fig3), fig_pathCombined+'Berg_vo.png', 'Resolution', 300)
% exportgraphics(figure(fig4), fig_pathCombined+'Berg_Lig_Startday.png', 'Resolution', 300)
% exportgraphics(figure(fig5), fig_pathCombined+'Berg_L_vs_vo.png', 'Resolution', 300)
% exportgraphics(figure(fig6), fig_pathCombined+'Berg_L_vs_relDo.png', 'Resolution', 300)
% exportgraphics(figure(fig7), fig_pathCombined+'Berg_L_vs_relDo_Scatter.png', 'Resolution', 300)

temp_par = [par_Berg, Ligdecomposition_Starts, avgvo,C_loss_lig_dec_start];
berg_par_sav = array2table(temp_par, ...
    'VariableNames', {'i', 'j', 'vhamx', 'mo', 'ro', 'AISC0', 'CN0', 'LN0', ...
    'rmse', 'rsquare', 'LigDecStartDay', 'avg_vo','C_loss_lig_dec_start'});
writetable(berg_par_sav, excel_path+"berg_final_table.xlsx")

%% HIROBE
close all
parFname = "est_params\"+scenario_name+"\fitted_par_Hirobe.txt";
par_hirobe = load(parFname);

xfig=figure;

% fig1 = figure;
% tiledlayout(5, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact')
% fig1.Position = [327, 185, 726, 829];
% fig2 = figure;
% tiledlayout(5, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% fig3 = figure;
% tiledlayout(5, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% fig4 = figure;
% tiledlayout(5, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% fig5 = figure;
% 
% fig2.Position = fig1.Position;
% fig3.Position = fig1.Position;
% fig4.Position = fig1.Position;

datacode = unique(organic_hirobe.Species);
Initial_litter_mass = 5; % grams

Ligdecomposition_Starts = ones(length(datacode), 1);
avgvo = Ligdecomposition_Starts;C_loss_lig_dec_start=Ligdecomposition_Starts;
for i = 1:length(datacode)
    sample_hirobe
    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    param.emax = emax_fun(CN0);
    par = par_hirobe(i, 2:4);
    [~,sol] =  opt_con(param,g,par,obs_data.tobs(end).*n);
    makeplot_state_space(sol, '',param,obs_data,g,xfig)
%     exportgraphics(xfig, fig_path_all+"Hirobe_"+i+"_"+Hirobe04_data.SpeciesName{i}+".png",'Resolution',100)

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    Co = sol.NumericalResults.State(2, :);
    e = param.emax - par(3) * vo;
    L=Co./CT; AIS = Co*4;L_AIS=AIS./CT;
    gL=g(L,a, b);
    vhmax=sol.NumericalResults.Parameter(1)  ;
    Dh = vhmax.*gL.*(CT-Co);
    Do = vo.*Co;

    voN = vo ./ max(vo);
    id = find(voN > 0.01);
    Ligdecomposition_Starts(i) = time(id(1));
    avgvo(i) = mean(vo); C_loss_lig_dec_start(i) = CT(id(1))/CT(1);
%     title_name = "i=" +i+"(" +Hirobe04_data.SpeciesName{i}+")";
%     post_plots

end
temp_par = [par_hirobe, Ligdecomposition_Starts, avgvo,C_loss_lig_dec_start];
hirobe_data = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx(d-1)', 'mo', 'ro(d)', 'AISC0', 'CN0', 'LN0', ...
    'rmse', 'rsquare', 'LigDecStartDay', 'avg_vo','C_loss_lig_dec_start'});
writetable(hirobe_data, excel_path+'hirobe_final_table.xlsx')
% figureHandles = findobj('Type', 'figure');
% set(figureHandles,'color','w')
% exportgraphics(figure(fig1), fig_pathCombined+'hirobe_state_space.png', 'Resolution', 300)
% exportgraphics(figure(fig2), fig_pathCombined+'hirobe_CUE.png', 'Resolution', 300)
% exportgraphics(figure(fig3), fig_pathCombined+'hirobe_vo.png', 'Resolution', 300)
% exportgraphics(figure(fig4), fig_pathCombined+'hirobe_Lig_Startday.png', 'Resolution', 300)
% exportgraphics(figure(fig5), fig_pathCombined+'hirobe_L_vs_relDo_Scatter.png', 'Resolution', 300)

%% OSONO04
close all
par_Osono = load("est_params\"+scenario_name+"\fitted_par_Osono04U.txt");

xfig = figure;
fig1 = figure;
tiledlayout(5, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact')
fig1.Position = [327, 185, 726, 829];
fig2 = figure;
tiledlayout(5, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig3 = figure;
tiledlayout(5, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig4 = figure;
tiledlayout(5, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig5 = figure;

fig2.Position = fig1.Position;
fig3.Position = fig1.Position;
fig4.Position = fig1.Position;

sz = size(amount_C);
Ligdecomposition_Starts = ones(sz(2), 1);
avgvo = Ligdecomposition_Starts;C_loss_lig_dec_start=Ligdecomposition_Starts;
for i = 1:sz(2)
    sample_Osono04U
    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    par = par_Osono(i, 2:4);
    [~,sol] =  opt_con(param,g,par,obs_data.tobs(end).*n);
    makeplot_state_space(sol, '',param,obs_data,g,xfig)
    exportgraphics(xfig, fig_path_all+"Osono04U_"+i+".png",'Resolution',100)
    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    Co = sol.NumericalResults.State(2, :);
    e = param.emax - par(3) * vo;
    L=Co./CT; AIS = Co*4;L_AIS=AIS./CT;
    gL=g(L,a, b);
    vhmax=sol.NumericalResults.Parameter(1)  ;
    Dh = vhmax.*gL.*(CT-Co);
    Do = vo.*Co;
    voN = vo ./ max(vo);
    id = find(voN > 0.01);
    Ligdecomposition_Starts(i) = time(id(1)) ;
    avgvo(i) = mean(vo);C_loss_lig_dec_start(i) = CT(id(1))/CT(1);
    title_name="i=" +i+"(" +treename{i}+")";
    post_plots

end
temp_par = [par_Osono, Ligdecomposition_Starts, avgvo,C_loss_lig_dec_start];
Osono_data = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx(d-1)', 'mo', 'ro(d)', 'AISC0', 'CN0', 'LN0', ...
    'rmse', 'rsquare', 'LigDecStartDay', 'avg_vo','C_loss_lig_dec_start'});
writetable(Osono_data, excel_path+'Osono_final_table.xlsx')
figureHandles = findobj('Type', 'figure');
set(figureHandles,'color','w')

exportgraphics(figure(fig4), fig_pathCombined+'Osono04_Lig_Startday.png', 'Resolution', 300)
exportgraphics(figure(fig1), fig_pathCombined+'Osono04_state_space.png', 'Resolution', 300)
exportgraphics(figure(fig2), fig_pathCombined+'Osono04_CUE.png', 'Resolution', 300)
exportgraphics(figure(fig3), fig_pathCombined+'Osono04_vo.png', 'Resolution', 300)
exportgraphics(figure(fig5), fig_pathCombined+'Osono04_L_vs_relDo_Scatter.png', 'Resolution', 300)

%% OSONO17
close all
parFname = "est_params\"+scenario_name+"\fitted_par_Osono17.txt";
par_Osono17 = load(parFname);

xfig = figure;

fig1 = figure;
tiledlayout(4, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact')
fig1.Position = [327, 185, 726, 829];
fig2 = figure;
tiledlayout(4, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig3 = figure;
tiledlayout(4, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig4 = figure;
tiledlayout(4, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
fig5 = figure;

fig2.Position = fig1.Position;
fig3.Position = fig1.Position;
fig4.Position = fig1.Position;

numSpecies = 12;
Ligdecomposition_Starts = ones(numSpecies, 1);
avgvo = Ligdecomposition_Starts;
C_loss_lig_dec_start=Ligdecomposition_Starts;
for i = 1:numSpecies
    sample_Osono17
    param.CO_0 = obs_data.Co_obs(1);
    param.CT_0 = obs_data.Ct_obs(1);
    par = par_Osono17(i, 2:4);
    [~,sol] =  opt_con(param,g,par,obs_data.tobs(end).*n);
    makeplot_state_space(sol, '',param,obs_data,g,xfig)
    exportgraphics(xfig, fig_path_all+"Osono17_"+i+".png",'Resolution',100)
    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    Co = sol.NumericalResults.State(2, :);
    e = param.emax - par(3) * vo;
    L=Co./CT; AIS = Co*4;L_AIS=AIS./CT;
    gL=g(L,a, b);
    vhmax=sol.NumericalResults.Parameter(1)  ;
    Dh = vhmax.*gL.*(CT-Co);
    Do = vo.*Co;
    voN = vo ./ max(vo);
    id = find(voN > 0.01);
    Ligdecomposition_Starts(i) = time(id(1)) ;
    C_loss_lig_dec_start(i) = CT(id(1))/CT(1);
    avgvo(i) = mean(vo);
    post_plots

end
temp_par = [par_Osono17, Ligdecomposition_Starts, avgvo,C_loss_lig_dec_start];
Osono17_data = array2table(temp_par, ...
    'VariableNames', {'i', 'vhamx(d-1)', 'mo', 'ro(d)', 'AISC0', 'CN0', 'LN0', ...
    'rmse', 'rsquare', 'LigDecStartDay', 'avg_vo','C_loss_lig_dec_start'});
writetable(Osono17_data, excel_path+'Osono17_final_table.xlsx')
figureHandles = findobj('Type', 'figure');
set(figureHandles,'color','w')

exportgraphics(figure(fig1), fig_pathCombined+'Osono17_state_space.png', 'Resolution', 300)
exportgraphics(figure(fig2), fig_pathCombined+'Osono17_CUE.png', 'Resolution', 300)
exportgraphics(figure(fig3), fig_pathCombined+'Osono17_vo.png', 'Resolution', 300)
exportgraphics(figure(fig4), fig_pathCombined+'Osono17_Lig_Startday.png', 'Resolution', 300)
exportgraphics(figure(fig5), fig_pathCombined+'Osono17_L_vs_relDo_Scatter.png', 'Resolution', 300)
%%
for i =1:10 
    fprintf("\n")
end
disp(upper("finished Post Processing"))


