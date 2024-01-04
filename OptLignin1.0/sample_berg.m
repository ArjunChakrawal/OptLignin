
idsp = find(dataset.SpeciesCode == spcode(j));
decdata2 = dataset(idsp(1):idsp(end), :);
id = decdata2.NMg_g < 0;
decdata2(id, :) = [];
total_C = (1-decdata2.MassLoss_*0.01)*initial_fraction_of_C; % gC/g litter
amount_AUR_C = (decdata2.AISMg_gInitialLitter*0.001*0.5).*(1-decdata2.MassLoss_(1)*0.01);  % gligC/g litter
aromaticC = aromatic_fraction_inAIS(amount_AUR_C); % true lignin

C_2 = decdata2.gC_gCInitial(2);
% N_2 = decdata2.gN_gNInitial(2);
if (C_2 < 0.7)
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
% fun = @(b,x) exp(b*x);
% k = nlinfit(obs_data.tobs ,obs_data.Ct_obs./obs_data.Ct_obs(1),fun,0.001);
f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
k=f.b ;

final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);
tt=0:tnorm(end)*n;

AISC0=decdata2.gLigninC_gC(tsrt);
CN0=1./decdata2.gN_gC(tsrt);
LN0=decdata2.gligC_gN(tsrt);