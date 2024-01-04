ddata=Osono2017(strcmpi(Osono2017.Tree,treename_Osono2017{i}),:);

id_t = ddata.Collection;
Time_days_ = day(1:length(id_t));
amount_C = ddata.MassG .* ddata.C_ .* 0.01; % g C/ bag
amount_N = ddata.MassG .* ddata.N_ .* 0.01; % g N/ bag
amount_AUR_C = fraction_of_C_in_AUR.*ddata.MassG .* ddata.AUR_ .* 0.01; % g lignin / bag

gC_gCInitial = amount_C./amount_C(1);
C_2 = gC_gCInitial(2);
if ( C_2 < 0.7)
    tsrt = 2;
    temp = 1;
else
    tsrt = 1;
    temp = 0;
end
tt=Time_days_(tsrt:end)-Time_days_(tsrt);
obs_data = [];
obs_data.tobs = tt;
obs_data.Ct_obs = amount_C(tsrt:end); %Remaining mass of litter (gC / bag)
obs_data.Co_obs = aromatic_fraction_inAIS(amount_AUR_C(tsrt:end));


tnorm = obs_data.tobs;
y=log(obs_data.Ct_obs./obs_data.Ct_obs(1));
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=(y(end)-y(1))./tnorm(end);
t=0:tnorm(end)*2;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);


% figure;
% subplot(211)
% scatter(obs_data.tobs,obs_data.Ct_obs,'DisplayName','TC'); hold on
% plot(t,obs_data.Ct_obs(1).*exp(k.*t))
% scatter(obs_data.tobs,obs_data.Co_obs,'DisplayName','AromaticC'); hold on
% xlabel("time"); ylabel('gC/bag')
% 
% subplot(212)
% plot(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1),'-o'); hold on
% xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
% xlim([0,1]);ylim([0,inf])
% title([num2str(i),'-',treename_Osono2017{i}])