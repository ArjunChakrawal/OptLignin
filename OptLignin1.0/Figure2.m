clear all
clc
close all


set_up
param.a=a;param.b=b;param.g=g;
param.vo_thres=vo_thres;
disp("this is the g-function with " +scenario_name+" scenario")
excel_path = "est_params\" +scenario_name + "\";
fig_path_all = "fig\" +scenario_name + "\";


fig=figure;
fig.Position=[100,100,810,340];
tiledlayout(2, 4, TileSpacing="compact", Padding="compact");
set(gcf, 'Color', 'w')
for i = 1:2
    for j = 1:4
        ax(i, j) = nexttile;
    end
end


CN0 = 50;
param.CT_0 = 0.5;
param.TT = 2000;
linestyle = {':', '-'};

ARC0_loop = 0.05;
ro_loop = 25;
vh_loop = 0.0015;
pname={'vhmax','ro','ARC0','mo'};
tr=[0.5,1];

% ARC0 = [ARC0_loop/2,ARC0_loop,ARC0_loop*2];
% vh=[vh_loop/2,vh_loop,vh_loop*2]; ro =[ro_loop/2,ro_loop,ro_loop*2];
ARC0 = [ARC0_loop,ARC0_loop*2];
vh=[vh_loop,vh_loop*3];
ro =[ro_loop,ro_loop*2];
mo = [0,param.mo];
tstep=1./vh';
tstep(tstep./20>1)=100;

LC=flipud(lines(length(vh)));

% LC = [255, 102, 102; % red
%     0, 153, 153]/255; % blue

tic
for i =1:length(CN0)
    param.emax = emax_fun(CN0(i));
    for j=1:4
        plotting_function(param,ARC0_loop,ro_loop,vh_loop, pname{j}, ...
            vh, ro, ARC0,mo,ax(:,j),linestyle,LC,1.5)
    end
end

%% figure 2 
title(ax(1,1),"varying \it v_H", 'FontWeight', 'normal')
title(ax(1,2),"varying \it r_O", 'FontWeight', 'normal')
title(ax(1,3),"varying \it ARC_0", 'FontWeight', 'normal')
title(ax(1,4),"varying \it \mu_O", 'FontWeight', 'normal')

set(ax(1,1:4),'XtickLabel',[])
xlabel(ax(2,1:4),'time [d]');

ylabel(ax(1,1),'C remaining [-]');
ylabel(ax(2,1), '{\it v_O} [d^{-1}]', 'Interpreter', 'tex')


leg=["{\it v_H} = ","{\it r_O} = ","{\it ARC_0} = ","{\it \mu_O = }"];
pvar={vh, ro, ARC0,mo};

for j=1:4
    tt=pvar{j};
    for i=1:2
        p(i)=plot(ax(2,j),nan,nan,linestyle{i},'Color',[0,0,0,tr(i)],...
            'LineWidth',2);

        str2(i)=leg(j)+tt(i);
    end
    lh=legend(ax(2,j),p,str2);
    lh.Box = 'off';lh.Location="northwest";lh.FontSize=9;
    lh.Position(2)= lh.Position(2)+0.025;lh.Position(1)= lh.Position(1)-0.01;
%     lh.Position=[0.7470    0.7992    0.1283    0.1010];
end

set(ax(1,1:4),'YLim' ,[0,1.75])
% set(ax(2,1:4),'YLim' ,[0,2])
set(ax(2,1:4),'YLim' ,[0,0.004])

tL.TileSpacing='compact';
str='abcdefghijkl';
ix=0;
for i =1:2
    for j=1:4
        ix=ix+1;
        text(ax(i,j),0.85,0.9,"("+str(ix)+")",'Units','normalized')
    end
end
set(ax, 'Fontsize', 10, 'LineWidth', 0.45,'box','on', ...
    'Xcolor', [1, 1, 1]*0.25, 'Ycolor', [1, 1, 1]*0.25,...
    'Xgrid','off','Ygrid','off')

p1=plot(ax(1,4),nan,nan,'Color',LC(1,:),'LineWidth',2);
p2=plot(ax(1,4),nan,nan,'Color',LC(2,:),'LineWidth',2);
lh=legend(ax(1,4),[p1,p2],["Total C","Aromatic C"]);
lh.Box = 'off';lh.Location="northwest";lh.FontSize=10;
lh.Position(2)= lh.Position(2)+0.035;lh.Position(1)= lh.Position(1)-0.01;

set(ax(1:2,2:4),'YtickLabel',[])

exportgraphics(fig, 'results/Figure2.png', 'Resolution','300')


%%

function plotting_function(param,ARC0_loop,ro_loop,...
    vh_loop,pnamef,vh,ro,ARC0,mo,axf1,linestyle,LC,lw)
g =param.g;
param.CO_0 = param.CT_0 .*ARC0_loop;
[ocp, sol] = opt_con(param, g, [vh_loop, param.mo, ro_loop], param.TT);

switch pnamef
    case 'vhmax'
        nloop=length(vh);
    case 'ro'
        nloop=length(ro);
    case 'ARC0'
        nloop=length(ARC0);
    case 'mo'
        nloop=length(mo);
end
%     tau_L = zeros(nloop, 1);

for i = 1:nloop
    switch pnamef
        case 'vhmax'
            par=[vh(i), param.mo, ro_loop];
            IC = [param.CT_0, param.CT_0 .* ARC0_loop];
            ocp.Parameter(1:3).setBound(par)
            ocp.State(1:2).setInitial(IC)
            ocp.parameterize();
            sol = ocp.optimize( ...
                'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
                );
        case 'ro'
            par=[vh_loop, param.mo, ro(i)];
            IC = [param.CT_0, param.CT_0 .* ARC0_loop];
            ocp.Parameter(1:3).setBound(par)
            ocp.State(1:2).setInitial(IC)
            ocp.parameterize();
            sol = ocp.optimize( ...
                'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
                );

        case 'ARC0'
            par=[vh_loop, param.mo, ro_loop];
            IC = [param.CT_0, param.CT_0 .* ARC0(i)];
            ocp.Parameter(1:3).setBound(par)
            ocp.State(1:2).setInitial(IC)
            ocp.parameterize();
            sol = ocp.optimize( ...
                'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
                );
        case 'mo'
            par=[vh_loop, mo(i), ro_loop];
            IC = [param.CT_0, param.CT_0 .* ARC0_loop];
            ocp.Parameter(1:3).setBound(par)
            ocp.State(1:2).setInitial(IC)
            ocp.parameterize();
            sol = ocp.optimize( ...
                'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
                );
    end

    time = sol.NumericalResults.Independent;
    vo = sol.NumericalResults.Control(1, :);
    CT = sol.NumericalResults.State(1, :);
    Co = sol.NumericalResults.State(2, :);
    vhmax_sol=sol.NumericalResults.Parameter(1);
    ro_sol = sol.NumericalResults.Parameter(3);
    e = param.emax - ro_sol * vo;
    L = Co ./ CT;
    gL = g(L, param.a, param.b);

    Dh = vhmax_sol .* gL .* (CT - Co);
    Do = vo .* Co;
    DT=(1-e).*(Dh+Do);

    plot(axf1(1), time, CT./CT(1), linestyle{i}, 'LineWidth', lw, ...
        'Color', LC(1,:));
    hold(axf1(1), 'on')
    plot(axf1(1), time, Co./Co(1), linestyle{i}, 'LineWidth', lw, ...
        'Color', LC(2,:));

    plot(axf1(2), time, vo, linestyle{i}, 'LineWidth', lw, ...
        'Color', 'k');
    hold(axf1(2), 'on')


end
end









