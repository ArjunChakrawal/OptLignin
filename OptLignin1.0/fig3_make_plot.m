function fig3_make_plot(sol, i, LC,lw,axf, linstl,alpha)

% emax = param.emax;
% a = param.a;
% b = param.b;
% 
% vhmax = sol.NumericalResults.Parameter(1);
% mo = sol.NumericalResults.Parameter(2);
% ro = sol.NumericalResults.Parameter(3);
T = sol.NumericalResults.Parameter(4);
time = sol.NumericalResults.Independent * T;
CT = sol.NumericalResults.State(1, :);
Co = sol.NumericalResults.State(2, :);
vo = sol.NumericalResults.Control(1, :);

% L = Co ./ CT;
% gL = g(L, a, b);
% Dh = vhmax .* gL .* (CT - Co);
% Do = vo .* Co;
% e = emax - ro * vo;
% G = e .* (Dh +  Do);
% voN = vo ./ max(vo);

plot(axf(2), time, vo, linstl{i}, 'LineWidth', 1.5, 'Color', [[0,0,0]./255,alpha]);
hold(axf(2), 'on')
% [[0,0,0]./255,alpha]
plot(axf(1), time, CT./CT(1), linstl{i}, 'LineWidth', 1.5, 'Color', [LC(1,:),alpha]);hold(axf(1), 'on')
plot(axf(1), time, Co./Co(1), linstl{i}, 'LineWidth', 1.5, 'Color', [LC(2,:),alpha]);

end
