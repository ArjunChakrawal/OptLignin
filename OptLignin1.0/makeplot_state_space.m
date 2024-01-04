function makeplot_state_space(sol, ~, param, obs_data, g, fig)
LC = linspecer(4);


f = param.f;

emax = param.emax;
a = param.a;
b = param.b;

vhmax = sol.NumericalResults.Parameter(1);
mo = sol.NumericalResults.Parameter(2);
ro = sol.NumericalResults.Parameter(3);
T = sol.NumericalResults.Parameter(4);

time = sol.NumericalResults.Independent * T ./ 365;

CT = sol.NumericalResults.State(1, :);
Co = sol.NumericalResults.State(2, :);
vo = sol.NumericalResults.Control(1, :);


L = Co ./ CT;
gL = g(L, a, b);
Dh = vhmax .* gL .* (CT - Co);
Do = vo .* Co;
e = emax - ro * vo;
G = e .* (Dh + f .* Do);

figure(fig);
clf;
fig.Position = [160, 140, 850, 650];
tiledlayout('flow', 'TileSpacing', 'Compact', 'Padding', 'Compact');
set(gcf, 'Color', 'w')

nexttile
plot(time, CT, 'linewidth', 2, 'DisplayName', 'hydrolyzable'); hold on
scatter(obs_data.tobs./365, obs_data.Ct_obs, ...
    30, 'DisplayName', "CTobs", 'MarkerFaceColor', LC(1, :));
plot(time, Co, 'linewidth', 2, 'DisplayName', 'Unhydrolyzable');
try
    scatter(obs_data.tobs_Co./365, obs_data.Co_obs, ...
        30, 'DisplayName', "COobs", 'MarkerFaceColor', LC(2, :)); hold on
catch
    scatter(obs_data.tobs./365, obs_data.Co_obs, ...
        30, 'DisplayName', "COobs", 'MarkerFaceColor', LC(2, :)); hold on
end
% lh = legend('show');
% lh.Box = "off";
xlabel('time (Y)')
grid on
try
    title({"vh_{max}=" +vhmax+", ro=" +ro+", mo =" +mo, "rmse="+param.rmse+" r2=" + param.rsquare })
catch
    title("vh_{max}=" +vhmax+"ro=" +ro+", mo =" +mo)
end
nexttile
plot(time, CT./CT(1), 'linewidth', 2, 'color', LC(1, :), 'DisplayName', 'CT'); hold on
scatter(obs_data.tobs./365, obs_data.Ct_obs./obs_data.Ct_obs(1), ...
    30, 'DisplayName', "CTobs", 'MarkerFaceColor', LC(1, :)); hold on
plot(time, Co./Co(1), 'linewidth', 2, 'color', LC(2, :), 'DisplayName', 'Co./Co(1)');
try
    scatter(obs_data.tobs_Co./365, obs_data.Co_obs./obs_data.Co_obs(1), ...
        30, 'DisplayName', "COobs", 'MarkerFaceColor', LC(2, :)); hold on
catch
    scatter(obs_data.tobs./365, obs_data.Co_obs./obs_data.Co_obs(1), ...
        30, 'DisplayName', "COobs", 'MarkerFaceColor', LC(2, :)); hold on
end
% legend('show')
xlabel('time (Y)')
ylabel('normalized to initial')
grid on
p1 = plot(nan, nan, 'linewidth', 2, 'color', LC(1, :));
p2 = plot(nan, nan, 'linewidth', 2, 'color', LC(2, :));
lh = legend([p1, p2], {'totalC', 'AromaticC'});
lh.Box = "off";
% ylim([0, 1.5])
nexttile
plot(1-CT./CT(1), Co./Co(1), 'linewidth', 2, 'DisplayName', 'model'); hold on
try
    [~,idtemp] = intersect(obs_data.tobs,obs_data.tobs_Co,'stable');
    Ct_obs = obs_data.Ct_obs(idtemp);
    scatter(1-Ct_obs./Ct_obs(1), obs_data.Co_obs./obs_data.Co_obs(1), 20, 'DisplayName', "data"); hold on
catch
    scatter(1-obs_data.Ct_obs./obs_data.Ct_obs(1), obs_data.Co_obs./obs_data.Co_obs(1), 20, 'DisplayName', "data"); hold on
end
% legend('show')
xlabel("mass loss"); ylabel('AIS (normalized to initial value)')
grid on
% ylim([0, 1.5])

nexttile
% stairs(time, vh, 'linewidth', 2,'DisplayName','vH');
plot(time, vo, 'linewidth', 2, 'DisplayName', 'vO');
ylabel('vo')
xlabel('time (Y)')
ylim([0, inf])
grid on

% nexttile
% plot(time, Do, 'linewidth', 2,'DisplayName','Do = f vo Co');hold on
% plot(time, Dh,'.-','DisplayName','Dh = vh_{max}Ch' )
% xlabel('time (Y)')
% legend('show')
% grid on

nexttile
plot(time, e, 'linewidth', 2)
xlabel('time (Y)')
ylabel('CUE')
grid on

nexttile
plot(time, cumtrapz(time, G), 'linewidth', 2)
xlabel('time (Y)')
ylabel('\int Growth rate')
grid on

nexttile
scatter(time, L, 'linewidth', 2); hold on
try
    scatter(obs_data.tobs(idtemp)./365,obs_data.Co_obs./Ct_obs)
catch
    scatter(obs_data.tobs./365,obs_data.Co_obs./obs_data.Ct_obs)
end
xlabel('time (Y)')
ylabel('L')
grid on

% nexttile
% plot(L, Do./(Dh+Do), 'linewidth', 2)
% xlabel('L')
% ylabel('Dh./(Dh+Do)')
% grid on

% nexttile
% plot(L, Do, 'linewidth', 2,'DisplayName','Do = f vo Co');hold on
% plot(L, Dh,'.-','DisplayName','Dh = vh_{max}Ch' )
% xlabel('L')
% legend('show')
% grid on

nexttile
plot(time, gL, 'linewidth', 2)
xlabel('time (Y)')
ylabel('g')
ylim([0, 1])
grid on

nexttile
plot(L, gL, 'linewidth', 2)
xlabel('L')
ylabel('g')
ylim([0, 1])
grid on

nexttile
try
    plot(CT,Co, 'linewidth', 2); hold on
    scatter(obs_data.Ct_obs,obs_data.Co_obs); hold on
catch
end

h = findobj('Type','axes');
for i = 1:numel(h)
    if strcmp(h(i).XLabel.String, 'time (Y)')
        h(i).XLim = [0 obs_data.tobs(end)/365]; % Set x-axis limits
    end
end

%
% nexttile
% scatter(vo,gL); hold on
% xlabel('vo');ylabel('g');
% grid on
%
% nexttile
% scatter(time, Dh);hold on
% xlabel('time (Y)');ylabel('Dh = vh_{max}(1-g)Ch');grid on
% set(gcf, 'Color','w')

% nexttile
% scatter(vo, G);hold on
% xlabel('vo');ylabel('G');grid on
% set(gcf, 'Color','w')
% exportgraphics(gcf,fname, 'Resolution',300)
end