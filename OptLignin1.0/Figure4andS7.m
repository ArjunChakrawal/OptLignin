clear all
clc
close all

%% figure 4
set_up
param.a=a;param.b=b;param.g=g;
param.vo_thres=vo_thres;
disp("this is the g-function with " +scenario_name+" scenario")
excel_path = "est_params\" +scenario_name + "\";
fig_path_all = "fig\" +scenario_name + "\";
fname= excel_path+"fit_results_table.xlsx";
fit_results_table=readtable(fname);

CN0 = [20,100];
param.CT_0 = 0.5;
param.TT = 2000;
linestyle = {'-', ':'};

ARC0_loop = 0.05;
ro_loop = 25;
vh_loop = 0.0015;
param.emax = emax_fun(CN0(2));
data=readtable('results/Figure3_datatable.xlsx');

LC=flipud(lines(2));
varstr={'vh','ro','ARC0','mo'};
yvarstr={'tau','avgvo','obj'};

fig = figure;
clf;
fig.Position(3:4)=[908,496];

t = tiledlayout(3, 4);
set(gcf, 'Color', 'w')
for i = 1:t.GridSize(1)
    for j = 1:t.GridSize(2)
        ax(i, j) = nexttile;
%         axis square
    end
end

for i=1:length(yvarstr)
    for j=1:length(varstr)
        temp=data(strcmp(data.var,varstr{j}),:);
        plot(ax(i,j),temp.value(temp.CN==20),temp.(yvarstr{i})(temp.CN==20),':', ...
            'LineWidth',1.5,'Color',LC(1,:)); hold(ax(i,j), 'on');
%         xlabel(ax(i,j),varstr{j});ylabel(ax(i,j),yvarstr{i})
        plot(ax(i,j),temp.value(temp.CN==100),temp.(yvarstr{i})(temp.CN==100), '-',...
            'LineWidth',1.5,'Color',LC(2,:))
%          xlabel(ax(i,j),varstr{j});ylabel(ax(i,j),yvarstr{i})
%         axis(ax(i,j),'tight');
    end
end
for j=1:t.GridSize(1)
    ax(j,1).XAxis.Scale='log';
end


title(ax(1,1),"varying \it v_H", 'FontWeight', 'normal')
title(ax(1,2),"varying \it r_O", 'FontWeight', 'normal')
title(ax(1,3),"varying \it ARC_0", 'FontWeight', 'normal')
title(ax(1,4),"varying \it \mu_O", 'FontWeight', 'normal')

% ylabel(ax(1,1),'\int_0^T {\it G }dt [gC d^{-1}]')
ylabel(ax(3,1),["Biomass =","\int_0^T {\it G }dt [gC]"])
ylabel(ax(1, 1), '\tau [d]')
ylabel(ax(2, 1), '$\bar{v}_O$ [d$^{-1}$]', Interpreter='latex')

% ylabel(ax(2, 1), 'max $({v}_O)$ [d$^{-1}$]', Interpreter='latex')

xlabel(ax(3,1),'{\it v_H} [d^{-1}]');
xlabel(ax(3,2),'{\it r_O} [d]');
xlabel(ax(3,3), ' {\it ARC_0} [gC g^{-1} totalC]');
xlabel(ax(3, 4), '{\it \mu_O}')



str='abcdefghijkl';
ix=0;
for i =1:3
    for j=1:4
        ix=ix+1;
        text(ax(i,j),0.05,0.9,"("+str(ix)+")",'Units','normalized')
    end
end
set(ax(3,1:4),'Ylim',[0,0.35])
set(ax(2,1:4),'Ylim',[0,2e-3])
set(ax(1:3,1),'Xlim',[0,0.1], 'Xtick',[1e-4,1e-3,1e-2,1e-1])
set(ax(1:3,4),'Xlim',[0,0.5],'Xtick',[0 0.25 0.5])

p1=plot(ax(2,4),nan,nan,':','LineWidth',1.5,'Color',LC(1,:));
p2=plot(ax(2,4),nan,nan,'-','LineWidth',1.5,'Color',LC(2,:));
lh=legend(ax(2,4),[p1,p2],["Low CN","High CN"]);
lh.Box='off';lh.Location='best';
set(ax, 'Fontsize', 10, 'LineWidth', 0.45,'box','on', ...
    'Xcolor', [1, 1, 1]*0.25, 'Ycolor', [1, 1, 1]*0.25,...
    'Xgrid','off','Ygrid','off')
t.TileSpacing="compact";
% t.Padding="loose";

yL = ylim(ax(1,1));
x = [(min(fit_results_table.vhmax)) (max(fit_results_table.vhmax)) (max(fit_results_table.vhmax)) (min(fit_results_table.vhmax))];
y = [0 0 yL(2) yL(2)];
patch(ax(1,1),x,y,'black','FaceAlpha',0.1,'EdgeColor','none')

yL = ylim(ax(2,1));
y = [0 0 yL(2) yL(2)];
patch(ax(2,1),x,y,'black','FaceAlpha',0.1,'EdgeColor','none')

yL = ylim(ax(3,1));
y = [0 0 yL(2) yL(2)];
patch(ax(3,1),x,y,'black','FaceAlpha',0.1,'EdgeColor','none')


yL = ylim(ax(1,2));
x = [(min(fit_results_table.ro)) (max(fit_results_table.ro)) (max(fit_results_table.ro)) (min(fit_results_table.ro))];
y = [0 0 yL(2) yL(2)];
patch(ax(1,2),x,y,'black','FaceAlpha',0.1,'EdgeColor','none')

yL = ylim(ax(2,2));
y = [0 0 yL(2) yL(2)];
patch(ax(2,2),x,y,'black','FaceAlpha',0.1,'EdgeColor','none')

yL = ylim(ax(3,2));
y = [0 0 yL(2) yL(2)];
patch(ax(3,2),x,y,'black','FaceAlpha',0.1,'EdgeColor','none')

exportgraphics(fig, 'results/Figure4.png', Resolution = 300)
%% figureS7
ARC0 = linspace(0,0.2,21);
ro = [0,0.1, 1, 10, 100];
LC=copper(length(ro));
figure;
for i =1:length(ro)
[tau,C_tau,avgvo,init_vo,max_vo,vo_at_tau,obj] = est_tau(param,ARC0_loop, ro(i),...
            vh_loop, 'ARC0', vh_loop, ro, ARC0,param.mo,100);
subplot(2,2,1)
plot(ARC0, tau/365, LineWidth=2,Color=LC(i,:), DisplayName="r_O="+ro(i));
legend('show')
xlabel('Initial fraction of aroamtic C [gC g^{-1} total C]');ylabel('Lag time [Y]');hold on
subplot(2,2,2)
plot(ARC0, avgvo, LineWidth=2,Color=LC(i,:), DisplayName="ro="+ro(i)); 
xlabel('Initial fraction of aroamtic C [gC g^{-1} total C]');ylabel('avg vo [d^{-1}]');hold on
subplot(2,2,3)
plot(ARC0, max_vo, LineWidth=2,Color=LC(i,:), DisplayName="ro="+ro(i)); 
xlabel('Initial fraction of aroamtic C [gC g^{-1} total C]');ylabel('max vo [d^{-1}]');hold on
subplot(2,2,4)
plot(ARC0, vo_at_tau, LineWidth=2,Color=LC(i,:), DisplayName="ro="+ro(i));
xlabel('Initial fraction of aroamtic C [gC g^{-1} total C]');ylabel('vo(t=\tau) [d^{-1}]');hold on

end
ax=gca;
exportgraphics(gcf, 'results/figureS7.png', Resolution = 300)



