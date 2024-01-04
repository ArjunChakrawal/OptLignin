function makeplot(sol,param,obs_data,g,a,b)
% close all
emax= param.emax;
f=param.f;

time = sol.NumericalResults.Independent./365;
Ch = sol.NumericalResults.State(1,:);
Co = sol.NumericalResults.State(2,:);
vo= sol.NumericalResults.Control(1,:);
L = Co ./ (Ch + Co);
gL=1-g(L,a, b);

vh=sol.NumericalResults.Parameter(1) .* gL;
mo = sol.NumericalResults.Parameter(2);
ro = sol.NumericalResults.Parameter(3);

Dh = vh.*Ch;
Do = vo.*Co;
e= emax - ro * vo;
G = e .* (Dh + f*Do);
fig=figure;
fig.Position=[ 160 140 850 650];
tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
set(gcf, 'Color','w')

nexttile
plot(time, Ch, 'linewidth', 2,'DisplayName','labile');hold on
plot(time, Co, 'linewidth', 2,'DisplayName','oxidizable');
legend('show')
xlabel('time (Y)')
grid on
title("ro="+ro + ", mo ="+mo+", f ="+param.f)

CT=Ch+Co;
nexttile
plot(time, CT./CT(1), 'linewidth', 2,'DisplayName','Ch+Co');hold on
scatter(obs_data.tobs./365, obs_data.Ct_obs./obs_data.Ct_obs(1), 20,'DisplayName',"CTobs"); hold on
plot(time, Co./Co(1), 'linewidth', 2,'DisplayName','Co');
try
    scatter(obs_data.tobs_Co./365, obs_data.Co_obs./obs_data.Co_obs(1), 20, 'DisplayName',"COobs"); hold on
catch
    scatter(obs_data.tobs./365, obs_data.Co_obs./obs_data.Co_obs(1), 20, 'DisplayName',"COobs"); hold on
end
% legend('show')
xlabel('time (Y)')
grid on
ylim([0,inf])

nexttile
% stairs(time, vh, 'linewidth', 2,'DisplayName','vH');
stairs(time, vo, 'linewidth', 2,'DisplayName','vO');
ylabel('vo')
xlabel('time (Y)')
ylim([0, inf])
grid on

nexttile
plot(time, Do, 'linewidth', 2,'DisplayName','Do = f vo Co');hold on
plot(time, Dh,'.-','DisplayName','Dh = vh_{max}Ch' )
xlabel('time (Y)')
legend('show')
grid on

nexttile
plot(time, e, 'linewidth', 2)
xlabel('time (Y)')
ylabel('CUE')
grid on

nexttile
plot(time, G, 'linewidth', 2)
xlabel('time (Y)')
ylabel('Growth rate')
grid on

nexttile
plot(time, L, 'linewidth', 2)
xlabel('time (Y)')
ylabel('L')
grid on

nexttile
plot(L, Dh./(Dh+Do), 'linewidth', 2)
xlabel('L')
ylabel('Dh./(Dh+Do)')
grid on

nexttile
plot(L, Do, 'linewidth', 2,'DisplayName','Do = f vo Co');hold on
plot(L, Dh,'.-','DisplayName','Dh = vh_{max}Ch' )
xlabel('L')
legend('show')
grid on

nexttile
plot(time, gL, 'linewidth', 2)
xlabel('time')
ylabel('g');ylim([0, 1])
grid on

nexttile
plot(L, gL, 'linewidth', 2)
xlabel('L')
ylabel('g')
ylim([0, 1])
grid on

dummy=0:0.01:1;
nexttile
plot(dummy, 1-g(dummy,a, b), 'linewidth', 2)
xlabel('L')
ylabel('g')
ylim([0, 1])
grid on

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