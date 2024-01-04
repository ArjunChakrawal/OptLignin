gC_gCInitial = amount_C(:, i)./amount_C(1, i);
C_2 = gC_gCInitial(2);
if ( C_2 < 0.7)
    tsrt = 2;
    temp = 1;
else
    tsrt = 1;
    temp = 0;
end
tt=osono_t(tsrt:end)-osono_t(tsrt);
obs_data = [];
obs_data.tobs = tt;
obs_data.Ct_obs = amount_C(tsrt:end, i); %Remaining mass of litter (gC / bag)
obs_data.Co_obs = aromatic_fraction_inAIS(amount_AUR_C_osono04(tsrt:end, i));

tnorm = obs_data.tobs;
y=log(obs_data.Ct_obs./obs_data.Ct_obs(1));
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=(y(end)-y(1))./tnorm(end);
t=0:tnorm(end)*2;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);

AISC0=  amount_AUR_C_osono04(tsrt, i)/obs_data.Ct_obs(1);