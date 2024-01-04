function par_est_vo(g,a,b,obs_data,fCO_0, lw, col,Topt_end,param,ha)
% clc
% parameters
vo_max=0.0009;
vh_max=0.002;

fCH_0 = 1 - fCO_0; % initial fraction of CH in litter[-]
fun = @(est_par, t) output_Var_(est_par, param, fCH_0, fCO_0,obs_data,g,a,b);
options1 = optimoptions('lsqcurvefit', 'Display', 'final', ...
    ...%      'PlotFcn', {@optimplotx, @optimplotfirstorderopt, @optimplotresnorm}, ...
    'MaxFunEvals', 2000, 'MaxIter', 1000, 'TolFun', 1e-12, 'TolX', 1e-12,...
    'FiniteDifferenceStepSize', 1e-9, 'FiniteDifferenceType', 'forward');
lb = [0.0001, 0.000];
ub = [0.05, 0.05];

init_guess = [vh_max, vo_max];
[LS_est_par, ~, ~ , ~, ~, ~, ~] = lsqcurvefit(fun, init_guess, ...
    obs_data.tobs, [obs_data.Ct_obs; obs_data.Co_obs], lb, ub, options1);
disp(LS_est_par)

%%
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
myode=@(t,y)mass_balance_ode(y, LS_est_par, param,g,a,b);
[~,State] = ode23s(myode,obs_data.tobs,[fCH_0; fCO_0],options);
Ch = State(:,1);
Co = State(:,2);
[rmse_ch,rmse_co, rsquare_ch,rsquare_co]=  mystat(Ch+Co, Co, obs_data );


[time,State] = ode23s(myode,[0 Topt_end],[fCH_0; fCO_0],options);
Ch = State(:,1);
Co = State(:,2);
vo = LS_est_par(1);
vh_max = LS_est_par(2);
L = Co ./ (Ch + Co);
Dh = vh_max .* g(L,a, b).*Ch;
Do = vo.*Co;
e= param.emax - param.ro .* vo;
G = e .* (Dh + Do);
axes(ha(1))
plot(time, Ch+Co, 'linewidth', lw,'Color',col(3,:),'DisplayName','Ch+Co');hold on
scatter(obs_data.tobs, obs_data.Ct_obs, 50,'Marker','diamond',...
    'MarkerEdgeColor','k','MarkerFaceColor','none','DisplayName',"CTobs");
plot(time, Co./Co(1),':', 'linewidth', lw,'Color',col(3,:),'DisplayName','Co./Co(1)');
scatter(obs_data.tobs, obs_data.Co_obs, 50, 'DisplayName',"COobs",...
    'MarkerEdgeColor','k','MarkerFaceColor','none');
% legend('show')
grid('off')
% p1=plot(nan, nan, 'linewidth', lw,'Color',col(3,:));
% p2=plot(nan, nan, ':','linewidth', lw,'Color',col(3,:));
% lh=legend([p1,p2],["total litter C ","fraction of structural C"]);
% lh.Box='off';lh.FontSize=14;lh.Location='southeast';
set(gca, 'FontSize',12,'Color','w')
title(['rmse Ct =',num2str(rmse_ch, '%1.3f'),' rmse Co =',num2str(rmse_co, '%1.3f'),...
    ',\newline r2 Ct =',num2str(rsquare_ch, '%1.3f'),', r2 Co=',num2str(rsquare_co, '%1.3f'),...
    ',\newline vo max =',num2str(LS_est_par(1), '%1.2E'), ' d^{-1}',', vh=',num2str(LS_est_par(2), '%1.2E'), ' d^{-1}'] ...
    , 'FontSize',11,'FontWeight','normal');
% exportgraphics(gca, 'CtCo_LS.png','Resolution',150)

axes(ha(6))
plot(time, cumtrapz(time,G), '-','linewidth', lw,'Color',col(3,:),...
    'DisplayName',"Least-square Optimization");hold on
xlabel('time (days)')
ylabel('Cumulative growth rate (\int_0^T G dt)')
grid off
lh=legend("show");lh.Box='off';lh.FontSize=12;lh.Location="southeast";
set(gca, 'FontSize',12, 'box','on')
% exportgraphics(gca, 'int_G_int.png','Resolution',150)

figure(7)
plot(L, g(L,a, b),'-','Color',col(3,:), ...
    'LineWidth',0.5,'DisplayName',"Least-square Optimization"); hold on 
title("g(L) = 1-L^2")   

% title("$g_1(L) = 1-L^2$,  "+ ...
%     "$g_2(L,a,b) =1-\frac{1}{1+exp \left(-\frac{L-0.6}{0.08}\right)}$", ...
%     'interpreter','latex')
legend('interpreter','tex','FontSize',12)
ylabel("Factor reducing carbohydrate uptake")
lh=legend("show");lh.Box='off';lh.FontSize=12;lh.Location="best"  ;
set(gca,'FontSize',12,'Color','w')
xlim([0,1])
% exportgraphics(gca, 'g_func.png','Resolution',150)

%%
    function ysim= output_Var_(est_par, param, fCH_0, fCO_0, obs_data,g,a,b)
        [~,yf] = ode23s(@(t,y)mass_balance_ode(y, est_par, param,g,a,b), ...
            obs_data.tobs ,[fCH_0; fCO_0], ...
            odeset('RelTol',1e-8,'AbsTol',1e-10));
        Chsim = yf(:,1);
        Cosim = yf(:,2);
        ysim = [Chsim + Cosim; Cosim ./ Cosim(1)];
    end

    function dx = mass_balance_ode(y, est_par,param,g,a,b)
        Chf = y(1);
        Cof = y(2);
        vof = est_par(1);
        vh = est_par(2);
        mu = 1;
%         mof = est_par(3);
        Lf = Cof / (Chf + Cof);
        Dhf = vh * g(Lf,a, b)*Chf;
        Dof = vof*Cof;
        ef= param.emax - param.ro * vof;
        Gf = ef * (Dhf + Dof);
        dx = [-Dhf + (1 - param.mo) *mu* Gf; ...
            -Dof + param.mo *mu* Gf];
    end

    

end
