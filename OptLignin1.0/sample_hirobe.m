Initial_litter_mass = 5; % grams
id = find(Hirobe04_data.Species == datacode(i));
ddata = Hirobe04_data(id, :);
obs_data=[];
obs_data.tobs  = ddata.Time*30;
Nt_obs = Initial_litter_mass.*(0.01 .* ddata.weightRemain_ )...
    .* (ddata.N_ .* 0.01);
Ct_obs = Initial_litter_mass.*(0.01 .* ddata.weightRemain_ )...
    .* (ddata.C_ .* 0.01); % g litter/ bag * fraction of litter*fraction of C = (gC/g bag)
obs_data.Ct_obs = Ct_obs;

datacode = unique(organic_hirobe.Species);
idL = find(organic_hirobe.Species == datacode(i));
data_Lig_Carb = organic_hirobe(idL, :);
obs_data.tobs_Co = data_Lig_Carb.Period.*30;

for j = 1:length(data_Lig_Carb.Period)
    idtemp(j) = find(ddata.Time == data_Lig_Carb.Period(j));
end
% g litter/ bag * fraction of litter*fraction of lignin C (g lignin C/ g litter)= g lignin C/ bag
AIS_C = Initial_litter_mass.*0.01 * ddata.weightRemain_(idtemp) .* ...
    data_Lig_Carb.Lignin_mg_g_ .* 0.001*fraction_of_C_in_AUR;
aromaticC = aromatic_fraction_inAIS(AIS_C); % true lignin
obs_data.Co_obs=aromaticC;

tnorm = obs_data.tobs;
y=log(obs_data.Ct_obs./obs_data.Ct_obs(1));
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=(y(end)-y(1))./tnorm(end);
t=0:tnorm(end)*2;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);