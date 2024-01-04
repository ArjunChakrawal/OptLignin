clearvars;
close all;
set_up;param.numiter=50;
par_est = load("est_params\" +scenario_name+"\fitted_par_Huang2021.txt");
data = readtable('data/Huang2021.xlsx', 'Range', 'A10:W38');


CN0 = 29;
AISC0 = 0.219;
treatementName = unique(data.treatment);
LC=linspecer(length(treatementName));
fig=figure; set(gcf, Color='w')
fig.Position(3:4)=[675,528];
voNorm=[];PeroxidaseNorm=[];
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

    vo = sol.NumericalResults.Control;
    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    voN = interp1(time, vo ./ max(vo),obs_data.tobs);

    voNorm=[voNorm;voN];
    PeroxidaseNormalized0_1 = (temp.PeroxidaseMl_g_day- min(temp.PeroxidaseMl_g_day))...
    ./(max(temp.PeroxidaseMl_g_day) -min(temp.PeroxidaseMl_g_day));
    PeroxidaseNorm=[PeroxidaseNorm;PeroxidaseNormalized0_1];
    subplot(221)
    plot(time, vo ./ max(vo), 'linewidth', 2,Color=LC(i,:)); hold on
    scatter(obs_data.tobs,PeroxidaseNormalized0_1,'linewidth', 2, MarkerEdgeColor=LC(i,:))
    xlim([0, obs_data.tobs(end)])
    xlabel('Time [d]')
    ylabel('Normalized enzyme activity')
    subplot(223)
    scatter(voN,PeroxidaseNormalized0_1,'linewidth', 2, MarkerEdgeColor=LC(i,:)); hold on 
    xlabel('Simulated')
    ylabel('Observed')
end

subplot(223);plot([0,1],[0,1], '--k', LineWidth=1.5)


% res=voNorm-PeroxidaseNorm;
% dev= mean(PeroxidaseNorm) -PeroxidaseNorm;
% SST = sum(dev.^2);
% SSR = sum(res.^2);
% r2 = 1 - SSR / SST;
[r,pval]=corr(voNorm,PeroxidaseNorm);
subplot(223);text(0.1,0.8, "[corr, pval] =["+num2str(r,2)+"," +num2str(pval,2)+"]",'Unit','Normalized')


%%
clearvars;
% close all;
LC=linspecer(2);

set_up;param.numiter=50;
init_CN = 49;
t= table2array(readtable('data/Snajdr et al 2011.xlsx','Range','P4:P12')); %months
C_remaining=table2array(readtable('data/Snajdr et al 2011.xlsx','Range','M4:M12')); %gC litter
lignin_C_remaining=table2array(readtable('data/Snajdr et al 2011.xlsx','Range','F22:F26')); % g C lignin
tobs_Co = table2array(readtable('data/Snajdr et al 2011.xlsx','Range','L22:L26'));  %months
aromaticC = aromatic_fraction_inAIS(lignin_C_remaining);
param.emax = emax_fun(init_CN);

perOx= table2array(readtable('data/Snajdr et al 2011.xlsx','Range','N4:N12'));
lacc= table2array(readtable('data/Snajdr et al 2011.xlsx','Range','R4:R12'));
ligox= perOx+lacc;

obs_data=[];
obs_data.tobs  =  t;
obs_data.Ct_obs = C_remaining;
obs_data.tobs_Co  = tobs_Co;
obs_data.Co_obs  =aromaticC;

f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=f.b ;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);

param.emax = emax_fun(init_CN);
param.CO_0=obs_data.Co_obs(1);
param.CT_0=obs_data.Ct_obs(1);

par_est = load("est_params\" +scenario_name+"\fitted_par_Snajdr.txt");
par = par_est(1, 2:4);
[~, sol] = opt_con(param, g, par, obs_data.tobs(end).*n);

vo = sol.NumericalResults.Control;
time = sol.NumericalResults.Independent;
CT = sol.NumericalResults.State(1, :);
voN =vo ./ max(vo);
subplot(222)
plot(time, voN, 'linewidth', 2); hold on
xlim([0, obs_data.tobs(end)])
scatter(obs_data.tobs, ligox./max(ligox),'linewidth', 2, MarkerEdgeColor=LC(1,:)); hold on
xlabel('Time [d]')
ylabel('Normalized enzyme activity')

subplot(224)
scatter(interp1(time, vo ./ max(vo),obs_data.tobs),ligox./max(ligox),'linewidth', 2, ...
    MarkerEdgeColor=LC(1,:))
[r,pval]=corr(ligox,interp1(time, vo ./ max(vo),obs_data.tobs));

text(0.1,0.8, "[corr, pval] =["+num2str(r,2)+"," +num2str(pval,2)+"]",'Unit','Normalized')
xlabel('Simulated')
ylabel('Observed')
hold on;plot([0,1],[0,1], '--k', LineWidth=1.5)

% legend('box','off')
ax=findall(gcf, 'type','axes');
set(ax,'box','on')
text(subplot(221),0.05,0.9,'(a)','Unit','Normalized')
text(subplot(222),0.05,0.9,'(b)','Unit','Normalized')
text(subplot(223),0.05,0.9,'(c)','Unit','Normalized')
text(subplot(224),0.05,0.9,'(d)','Unit','Normalized')

title(subplot(221),"Huang et al. 2021")
title(subplot(222),"Šnajdr et al. 2011")
fig.Position=[172   346   864   688];
exportgraphics(gcf,"results/FigureS6.png",Resolution=400)


