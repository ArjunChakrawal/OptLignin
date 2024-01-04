
Initial_litter_mass = 1; % grams

temp=massLoss.Properties.VariableNames(2:end);

initial_fraction_of_C = 0.001*PA_initial.TotalC_g_kg_1_(strcmpi(PA_initial.CIDETCode,temp{i}));
init_CN = PA_initial.C_N_atomic_(strcmpi(PA_initial.CIDETCode,temp{i}));

id= strcmpi(PAnmr.SPCode,temp{i});
lignin_C_conc= PAnmr.AROM(id)+PAnmr.PHEN(id); % g/kg litter or mg/g litter

ml= table2array(massLoss(:,temp{1}));
Ct_obs = Initial_litter_mass.*0.01 .*ml.*initial_fraction_of_C; % g litter/ bag * fraction of litter*fraction of C = (gC/g bag)
tobs  =  massLoss.t_year_;
tobs_Co  = PAnmr.year(id);

for j = 1:length(tobs_Co)
    idtemp(j) = find(tobs == tobs_Co(j));
end
lignin_C =Ct_obs(idtemp) .*lignin_C_conc.* 0.01;


obs_data=[];
obs_data.tobs  =  massLoss.t_year_*365;
obs_data.Ct_obs = Ct_obs;
obs_data.tobs_Co  = tobs_Co*365;
obs_data.Co_obs  =lignin_C;

tnorm = obs_data.tobs;
y=log(obs_data.Ct_obs./obs_data.Ct_obs(1));
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=(y(end)-y(1))./tnorm(end);
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);
t=0:tnorm(end)*n;


% figure
% nexttile
% scatter(obs_data.tobs,obs_data.Ct_obs,'DisplayName','TC'); hold on
% plot(t,obs_data.Ct_obs(1).*exp(k.*t))
% scatter(obs_data.tobs_Co,obs_data.Co_obs,'DisplayName','AromaticC'); hold on
% xlabel("time"); ylabel('gC/bag')
% 
% nexttile
% plot(1-obs_data.Ct_obs(idtemp)./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1),'-o'); hold on
% xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
% xlim([0,1]);ylim([0,inf])
% title([num2str(i),'-',temp{i}])
% 
