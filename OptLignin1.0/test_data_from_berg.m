datacode = unique(berg_data.DatasetCode);
i = 1;
id = find(berg_data.DatasetCode == i);
dataset = berg_data(id(1):id(end), :);
spcode = unique(dataset.SpeciesCode);
j = 1;
idsp = find(dataset.SpeciesCode == spcode(j));
decdata2 = dataset(idsp(1):idsp(end), :);
id = decdata2.NMg_g < 0;
decdata2(id, :) = [];
total_C = (1-decdata2.MassLoss_*0.01)*initial_fraction_of_C; % gC/ bag
amount_AUR_C = (decdata2.AISMg_gInitialLitter*0.001*0.5).*(1-decdata2.MassLoss_(1)*0.01);  % gligC/ bag
aromaticC = aromatic_fraction_inAIS(amount_AUR_C); % true lignin


C_2 = decdata2.gC_gCInitial(2);
N_2 = decdata2.gN_gNInitial(2);
if (N_2 < 0.9 && N_2 < C_2 || C_2 < 0.7)
    tsrt = 2;
    temp = 1;
else
    tsrt = 1;
    temp = 0;
end

obs_data = [];
Time_days_= decdata2.Time_days_(tsrt:end) - decdata2.Time_days_(tsrt);
obs_data.tobs = Time_days_;
obs_data.Ct_obs = total_C(tsrt:end);
obs_data.Co_obs = aromaticC(tsrt:end);

tnorm = obs_data.tobs;
y=log(obs_data.Ct_obs./obs_data.Ct_obs(1));
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=(y(end)-y(1))./tnorm(end);
t=0:tnorm(end)*2;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);


figure;
tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
nexttile
plot(obs_data.tobs,obs_data.Ct_obs,'-o','DisplayName','total C'); hold on
plot(t,obs_data.Ct_obs(1).*exp(k.*t),'DisplayName','exp fit')
scatter(terminalTime,obs_data.Ct_obs(1).*exp(k.*terminalTime),100,'Marker','*', ...
    'MarkerEdgeColor','red','DisplayName','simulation period')
plot(obs_data.tobs,obs_data.Co_obs,'-o','DisplayName','AromaticC'); hold on
xlabel("time"); ylabel('gC/bag')

nexttile
plot(obs_data.Ct_obs,obs_data.Co_obs,'-o'); hold on
CN0=decdata2.gCInitial_gNInitial(tsrt);

%{
figure;
tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
nexttile
scatter(obs_data.tobs,obs_data.Ct_obs./obs_data.Ct_obs(1),100,'o', ...
    'DisplayName','C remaining','LineWidth',1); hold on
scatter(obs_data.tobs,amount_AUR_C(2:end)./amount_AUR_C(2),100,'d', ...
    'DisplayName','AUR remaining','LineWidth',1); hold on

set(gcf,'color','w')
set(gca,'FontSize',20)
ylim([0, inf])
lh=legend('show');lh.Box='off';lh.FontSize=24;

axis tight


% scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1)); hold on
% xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')

%}

