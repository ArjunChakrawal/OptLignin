clear all
clc
% close all

%%
set_up
param.a=a;param.b=b;param.g=g;
param.vo_thres=vo_thres;
disp("this is the g-function with " +scenario_name+" scenario")
excel_path = "est_params\" +scenario_name + "\";
fig_path_all = "fig\" +scenario_name + "\";

%%
% close all
CN0 = [20,100];
param.CT_0 = 0.5;
param.TT = 3000;
linestyle = {'-', ':'};

ARC0_loop = 0.05;
ro_loop = 25;
vh_loop = 0.0015;
param.emax = emax_fun(CN0(1));

param.CO_0 = param.CT_0 .*ARC0_loop*3;
% vh=20*10^-4;
vh=0.007666822074546;
tstep=1./vh';
tstep(tstep./20>1)=100;
tstep(tstep./20<1)=round(linspace(150,500,sum((tstep./20<1))));

[ocp, sol] = opt_con(param, param.g, [vh_loop, param.mo, ro_loop*2], param.TT);
% ocp.solve('controlIntervals', tstep, ...
%                 'collocationPoints', 'legendre', ... %
%                 'polynomialDegree', 3, ...
%                 'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
%                 );
time = sol.NumericalResults.Independent;
vo = sol.NumericalResults.Control(1, :);
CT = sol.NumericalResults.State(1, :);
Co = sol.NumericalResults.State(2, :);

dcodt = diff(Co)./diff(time);
neg_idx  = find(dcodt<0, 1);
voN = vo ./ max(vo);
if(isempty(neg_idx ))
    id=length(time);
else
    id = find(voN > vo_thres);
end
tau_L = time(id(1))
figure; tiledlayout('flow');
nexttile;plot(time, CT); hold on; plot(time, Co); plot(time, CT-Co);
nexttile; plot(time, vo);hold on;plot(vo*0+tau_L,vo); title("tau="+tau_L+"  vo(tau)="+vo(id(1)))


% 
% % vh= 0.05;
% % LCol={copper(length(vh)),winter(length(vh))};
% % 
% % ro =200;
% % pnamef='ro';
% % tic
% % [tau_L, ~, avgvo, ~,obj,sol] = est_tau(param,ARC0_loop,ro_loop,...
% %     vh_loop, pnamef, vh, ro, ARC0_loop,param.mo,200);
% % toc
% % figure;plot(sol.NumericalResults.Independent,sol.NumericalResults.State(1,:));hold on
% % plot(sol.NumericalResults.Independent,sol.NumericalResults.State(2,:))
% % vo=sol.NumericalResults.Control ;
% % figure;plot(sol.NumericalResults.Independent,vo./max(vo))
% % dcodt = diff(sol.NumericalResults.State(2,:))./diff(sol.NumericalResults.Independent);
% % figure;plot(sol.NumericalResults.Independent(2:end),dcodt)
% 

%%
% close all

mo=linspace(0.01,0.5,100);

ARC0 = linspace(0.01,0.2,100);
vh= 10.^([linspace(-4,-3.5,15),linspace(-3.5,-2,40),...
    linspace(-2,-1,20),linspace(-1,0,10)]);
vh=unique(vh);
tstep=1./vh';
tstep(tstep./20>1)=100;
tstep(tstep./20<1)=round(linspace(150,500,sum((tstep./20<1))));

ro = linspace(1,200,100);
pvar={vh, ro, ARC0,mo};
pname={'vh','ro','ARC0','mo'};

nn = length(vh)+length(ro)+length(ARC0)+length(mo);
CN = [repmat(CN0(1), nn, 1)];
var2 = {repmat('vh', length(vh), 1); repmat('ro', length(ro), 1); repmat('ARC0', length(ARC0), 1)};
var=repmat(' ' , nn, 1);
value = zeros(nn, 1);
tau = zeros(nn, 1);
avgvo = zeros(nn, 1);
obj = zeros(nn, 1);
avgvo_vh = zeros(nn, 1);
vo_at_tau = zeros(nn, 1);
init_vo = zeros(nn, 1);
max_vo = zeros(nn, 1);

T = table(CN, var, value, tau,obj, avgvo,vo_at_tau,init_vo,max_vo,avgvo_vh, ...
    'VariableNames', {'CN', 'var', 'value', ...
    'tau','obj', 'avgvo','vo_at_tau','init_vo','max_vo','avgvo_vh'});

T.var = repmat(' ', nn, length('ARC0'));
T.var(1:length(vh), 1:length('vh')) = repmat('vh', length(vh), 1);
T.var(length(vh)+1:length(vh)+length(ro), 1:length('rO')) = repmat('ro', length(ro), 1);
T.var(length(vh)+length(ro)+1:length(vh)+length(ro)+length(ARC0), 1:length('ARC0')) = repmat('ARC0', length(ARC0), 1);
T.var(length(vh)+length(ro)+length(ARC0)+1:end, 1:length('mo')) = repmat('mo', length(mo), 1);

T.value=[vh,ro,ARC0,mo]';
T2=T;
T2.CN=T2.CN*0+CN0(2);

T=[T;T2];
var=[];
for i=1:length(T.var)
    var=  [var;{strrep(T.var(i,:), ' ', '')}];
end
T.var=var;

%%
data=[];
tic
for i =1:length(CN0)
    param.emax = emax_fun(CN0(i));

    for j=1:4
%         nn = length(pvar{j});
%         tau = zeros(nn, 1)+999*i;
%         avgvo = zeros(nn, 1)+888*i;
%         obj = zeros(nn, 1)+777*i;
%         avgvo_vh = zeros(nn, 1)+666*i;
%         vo_at_tau = zeros(nn, 1)+555*i;
%         vo_init = zeros(nn, 1)+444*i;
%         vo_init = zeros(nn, 1)+444*i;
%         max_vo = zeros(nn, 1)+333*i;
        [tau,C_tau,avgvo,init_vo,max_vo,vo_at_tau,obj] = est_tau(param,ARC0_loop,ro_loop,...
            vh_loop, pname{j}, vh, ro, ARC0,mo,tstep);
        if(j==1)
            avgvo_vh=avgvo./pvar{j}';
        else
            avgvo_vh=avgvo./vh_loop;
        end
        temp=[tau,obj, avgvo,vo_at_tau,init_vo,max_vo,avgvo_vh];
        temp2=T(strcmp(T.var,pname{j}),:);
        temp2=temp2(temp2.CN==CN0(i),:);
        temp2(:,4:end)=array2table(temp);
        data=[data;temp2];
    end
end

% T.tau=temp(:,1); T.avgvo=temp(:,2); T.obj=temp(:,3);
% T.avgvo_vh=temp(:,4);

toc

%%
writetable(data,'results/Figure3_datatable.xlsx')




