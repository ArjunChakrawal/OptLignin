
obs_data=[];
obs_data.tobs  =  t(1:tend);
obs_data.Ct_obs = table2array(gCt_obs(1:tend,i));
obs_data.Co_obs  =aromatic_fraction_inAIS(table2array(gCo_obs(1:tend,i)));

tnorm = obs_data.tobs;
f = fit(obs_data.tobs ,obs_data.Ct_obs,'exp1');
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);
k=f.b ;
tt=0:tnorm(end)*2;
terminalTime = log(final_C*fraction_of_final_C_at_Terminal_time)/k;
n=terminalTime./obs_data.tobs(end);

figure;
subplot(121)
plot(obs_data.tobs,obs_data.Ct_obs./obs_data.Ct_obs(1),'o-','DisplayName','TC'); hold on
plot(tt,1*exp(k.*tt))
scatter(terminalTime,exp(k.*terminalTime),100,'Marker','*', ...
    'MarkerEdgeColor','red')
plot(obs_data.tobs,obs_data.Co_obs./obs_data.Co_obs(1),'o-','DisplayName','ligninC'); hold on
xlabel("time"); ylabel('gC/bag')

subplot(122)
scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1),obs_data.Co_obs./obs_data.Co_obs(1));
xlabel("mass loss (fraction of initial)"); ylabel('Aromatic (normalized to initial value)')
