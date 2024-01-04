function fit_plotting(sol, obs_data)
LC=linspecer(2);
CT=sol.y(1,:);Co=sol.y(2,:);
plot(sol.x,CT./CT(1),'DisplayName',"total C");hold on
plot(sol.x,Co./Co(1),'DisplayName',"Aromaic C");hold on
scatter(obs_data.tobs, obs_data.Ct_obs./obs_data.Ct_obs(1), ...
   30, 'DisplayName', "CTobs", 'MarkerFaceColor', LC(1, :));
try
   scatter(obs_data.tobs_Co, obs_data.Co_obs./obs_data.Co_obs(1), ...
      30, 'DisplayName', "COobs", 'MarkerFaceColor', LC(2, :)); hold on
catch
   scatter(obs_data.tobs, obs_data.Co_obs./obs_data.Co_obs(1), ...
      30, 'DisplayName', "COobs", 'MarkerFaceColor', LC(2, :)); hold on
end

xlabel('day'); ylabel('C [g])')
end