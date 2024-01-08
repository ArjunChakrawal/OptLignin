clear all
clc
close all

set_up
excel_path = "est_params\" +scenario_name + "\";

fig=figure;fig.Position=[100   150   900   220];
fig.Color='w' ;
tiledlayout('flow',TileSpacing='tight',Padding='compact')
for i=1:3
    ax(i)=nexttile;
end
LC=[44,127,184; ...
    mean([44,127,184; 255, 255, 255], 1); ...
    mean([226, 50, 1; 255, 255, 255], 1); ...
    226, 50, 1] / 255;
Rsquare=[];RMSE=[];ARC0=[];L0=[];climate=[];
%% boreal
berg_data = readtable('data\Berg and McClaugherty 1989 - Suppl Mat.xlsx');
par_est = load("est_params\" +scenario_name+"\fitted_par_Berg.txt");
par = par_est(1, 3:5);

i = 1;j=1;
id = find(berg_data.DatasetCode == i);
dataset = berg_data(id(1):id(end), :);
spcode = unique(dataset.SpeciesCode);
sample_berg
param.emax = emax_fun(CN0);
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);
% init_guess = par; % [vh_max, mo,ro];
% param.lb=[0.0005, param.mo,0.5];
% param.ub=[0.04, param.mo,400];
% [par,sol,~,~] = find_parameter(obs_data,...
%     param,init_guess,g,figure,@makeplot_state_space,@ysim_state_space,n);
[~,sol] =  opt_con(param,g,par,obs_data.tobs(end)*n);
plot_figure5(sol,obs_data, ax,LC(1,:),param)
[r2, rmse] = est_r2_rmse(obs_data, sol);
Rsquare =[Rsquare,r2];
RMSE =[RMSE,rmse];
ARC0=[ARC0;param.CO_0];
L0=[L0;param.CO_0/param.CT_0];
climate = [climate;"Boreal"];
%% warm temperate
par_est = load("est_params\" +scenario_name+"\fitted_par_Tu2014.txt");
data = readtable('data/Tu_et_al2014.xlsx', 'Range', 'B16:O360');
CN0 = table2array(readtable('data/Tu_et_al2014.xlsx', 'Range', 'K2:K11'));
AISC0 = table2array(readtable('data/Tu_et_al2014.xlsx', 'Range', 'N2:N11'));
SP = readtable('data/Tu_et_al2014.xlsx', 'Range', 'A1:A11');

treatementName = unique(data.treatment);
SPname = unique(data.SpeciesCode);
i=8;j=1;
tempdata1 = data(data.SpeciesCode == SPname(i), :);
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
obs_data = [];
obs_data.tobs = tempdata2.days(tsrt:end) - tempdata2.days(tsrt);
obs_data.Ct_obs = tempdata2.CInG(tsrt:end);
obs_data.Co_obs = aromatic_fraction_inAIS(tempdata2.ligninCInG(tsrt:end));

tnorm = obs_data.tobs;
f = fit(obs_data.tobs, obs_data.Ct_obs, 'exp1');
final_C = obs_data.Ct_obs(end) ./ obs_data.Ct_obs(1);
k = f.b;
terminalTime = log(final_C*0.5) / k;
n = terminalTime ./ obs_data.tobs(end);
tt = 0:tnorm(end) * n;
param.CO_0 = obs_data.Co_obs(1);
param.CT_0 = obs_data.Ct_obs(1);
param.emax = emax_fun(CN0(i));
% par = [0.0021    0.1000    4.1516];
par = par_est(29, 3:5);
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);

plot_figure5(sol,obs_data, ax,LC(2,:),param)
[r2, rmse] = est_r2_rmse(obs_data, sol);
Rsquare =[Rsquare,r2];
RMSE =[RMSE,rmse];ARC0=[ARC0;param.CO_0];
L0=[L0;param.CO_0/param.CT_0];climate = [climate,"Warm temp."];
%% cold temperate
par_est = load("est_params\" +scenario_name+"\fitted_par_Osono04U.txt");
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
i=14;sample_Osono04
param.CO_0 = obs_data.Co_obs(1);
param.CT_0 = obs_data.Ct_obs(1);
param.emax = emax_fun(CN0_osono(i));
par = par_est(i, 2:4);
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
plot_figure5(sol,obs_data, ax,LC(3,:),param)
[r2, rmse] = est_r2_rmse(obs_data, sol);
Rsquare =[Rsquare,r2];
RMSE =[RMSE,rmse];ARC0=[ARC0;param.CO_0];
L0=[L0;param.CO_0/param.CT_0];
climate = [climate,"Cold temp."];

%% Tropical
Hirobe04_data = readtable('data/Hirobe04.xls', 'Sheet', 'C and Mineral');
organic_hirobe = readtable('data/Hirobe04.xls', 'Sheet', 'Organic');
datacode = unique(organic_hirobe.Species);
par_est = load("est_params\" +scenario_name+"\fitted_par_Hirobe.txt");

i=3;sample_hirobe;
CN0 = Ct_obs(1) / Nt_obs(1);
param.CO_0 = obs_data.Co_obs(1);
param.CT_0 = obs_data.Ct_obs(1);
param.emax = emax_fun(CN0);
par = par_est(i, 2:4);
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);
plot_figure5(sol,obs_data, ax,LC(4,:),param)
[r2, rmse] = est_r2_rmse(obs_data, sol);
Rsquare =[Rsquare,r2];
RMSE =[RMSE,rmse];ARC0=[ARC0;param.CO_0];
L0=[L0;param.CO_0/param.CT_0];climate = [climate,"Tropical"];
%%
for i=1:3
xlabel(ax(i),'C loss (% of initial)')
end
lstr= ""+num2str(L0,1);
p=[];str=[];

for i=1:4
    p(i)=plot(ax(1),nan, nan,'LineWidth',2,'Color',LC(i,:));
    hold(ax(1),'on')
    str = [str,lstr(i)+" ("+climate(i)+")"];
%     str=[str,climate(i)+" ({\it rmse} = "+ num2str(RMSE(i),2)+", {\itr^2}="+num2str(Rsquare(i),2)+")"];
end

str=["\it P. sylvestris","\it E. grandis","\it C. japonica","\it S. beccariana"];
lh=legend(ax(3),p,str);
set(ax, 'Fontsize', 10, 'LineWidth', 0.45,'box','on', ...
    'Xcolor', [1, 1, 1]*0, 'Ycolor', [1, 1, 1]*0,...
    'Xgrid','off','Ygrid','off')
lh.Box='off';lh.Location='bestoutside';lh.FontSize=11;lh.NumColumns=1;
% title(lh,{'{\it ARC}_0','[gC g^{-1} totalC]'})
title(lh,{''})
str='abc';
for j=1:3
        text(ax(j),0.87,0.9,"("+str(j)+")",'Units','normalized')
end
exportgraphics(fig, 'results/Figure4.png', Resolution = 300)

%%
function plot_figure5(sol,obs_data,ax,LC,param)

ro = sol.NumericalResults.Parameter(3); %Y
% time = sol.NumericalResults.Independent ;
CT = sol.NumericalResults.State(1, :);
Co = sol.NumericalResults.State(2, :);
vo = sol.NumericalResults.Control(1, :); %per Y
voN = vo ./ max(vo);
id = find(voN > param.vo_thres);
x=(1-CT./CT(1))*100;
tau=x(id(1));

% x=time;
plot(ax(1),(1-CT./CT(1))*100, Co./Co(1), 'Linewidth', 2, 'Color',LC);hold(ax(1),'on')

try
    [~,idtemp] = intersect(obs_data.tobs,obs_data.tobs_Co,'stable');
    Ct_obs = obs_data.Ct_obs(idtemp);
    scatter(ax(1),(1-Ct_obs./Ct_obs(1))*100, ...
        obs_data.Co_obs./obs_data.Co_obs(1), ...
        40, 'DisplayName', "CTobs", 'MarkerFaceColor', LC,...
        'MarkerEdgeColor','k');
catch
    scatter(ax(1),(1-obs_data.Ct_obs./obs_data.Ct_obs(1))*100, ...
        obs_data.Co_obs./obs_data.Co_obs(1), ...
        40, 'DisplayName', "CTobs",  'MarkerFaceColor', LC,...
        'MarkerEdgeColor','k');
end
ylabel(ax(1),'Aromatic C [-]')

plot(ax(2),x, vo,'-','LineWidth',2,'Color',LC);hold(ax(2),'on')
xline(ax(2),tau,'k:','LineWidth',1.5,'Color',LC)
text(ax(2),tau+2,4e-3,"\tau="+num2str(tau,'%1.0f'),'FontSize',11, ...
    'Color',LC,'FontWeight','normal', Rotation=90)
xlim(ax(2),[0 inf])
ylabel(ax(2),'{\it v_O} [d^{-1}]','Interpreter','tex')

plot(ax(3),x, param.emax - ro .* vo,'LineWidth',2,'Color',LC);hold(ax(3),'on')
ylabel(ax(3),'{\it CUE} [-]','Interpreter','tex')
ylim(ax(3),[0 inf])
% [r2, rmse, r2_co, rmse_co, r2_CT, rmse_CT,AIC] = est_r2_rmse(obs_data, sol);
% text(ax(1),10,5e-5,'\tau','FontSize',12,'Color',LC,'FontWeight','bold')

end










