function double_exp(obs_data)



f = fit(obs_data.tobs ,obs_data.Ct_obs./obs_data.Ct_obs(1),'exp2');
final_C = obs_data.Ct_obs(end)./obs_data.Ct_obs(1);

figure;
scatter(obs_data.tobs,obs_data.Ct_obs./obs_data.Ct_obs(1));hold on
plot(1:100:10000,f(1:100:10000))

fun=@(x) f(x)-final_C*fraction_of_final_C_at_Terminal_time;
terminalTime = fsolve(fun,1000);

end