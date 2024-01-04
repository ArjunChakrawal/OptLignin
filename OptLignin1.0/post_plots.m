% state space
figure(fig1);
nexttile
plot(1-CT./CT(1), Co./Co(1), 'linewidth', 2, 'DisplayName', 'model'); hold on

try
    for j = 1:length(obs_data.tobs_Co)
        idtemp(j) = find(obs_data.tobs == obs_data.tobs_Co(j));
    end
    Ct_obs = obs_data.Ct_obs(idtemp);
    scatter(1-Ct_obs./Ct_obs(1), obs_data.Co_obs./obs_data.Co_obs(1), 20, 'DisplayName',"data"); hold on
catch
    scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1), obs_data.Co_obs./obs_data.Co_obs(1), 20, 'DisplayName',"data"); hold on
end
xlabel("mass loss"); ylabel('AIS/AIS_0')
title(title_name)

% cue
figure(fig2);
nexttile
plot(time, e, 'linewidth', 2)
xlabel('time (day)')
ylabel('CUE')
ylim([0, param.emax])
title(title_name)

% vo vs time
figure(fig3);
nexttile
stairs(time, vo, 'linewidth', 2, 'DisplayName', 'vO');
ylabel('vo')
xlabel('time (day)')
ylim([0, inf])
title(title_name)

% lig start time
figure(fig4);
nexttile
plot(time, voN); hold on
plot([1, 1]*time(id(1)), [0, 1])
xlabel('day'); ylabel('vo/max(vo)')
title(title_name)

% LAISC  vs relative lig decompostion rate
figure(fig5);
scatter(L_AIS, Do./(Dh+Do),'o');hold on
ylabel('Do./(Dh+Do)')
xlabel('LCI=AISC/CT');
