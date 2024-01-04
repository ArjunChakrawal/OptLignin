clear all
clc
close all
%%
set_up
disp("this is the g-function with "+ scenario_name+" scenario")
excel_path = "est_params\"+scenario_name+"\";
fig_path_all = "fig\"+scenario_name+"\";
vomax =0.1;
%%
num=20000;num2=500;
txt={'vhmax', 'mo', 'ro' ,'CN0','CT0','L0'};

X = lhsdesign(num,6);
lb = [1/3650, 0.05,0.1,10,0.5,0.01];
ub=[1/60,0.2,100,200,0.5,0.25];
rescaled_x=X;

for i =1:6
    x=X(:,i);
    rescaled_x(:,i) = (x - min(x)) / (max(x) - min(x)) * (ub(i) - lb(i)) +  lb(i);
end
e = emax_fun(rescaled_x(:,5)) - rescaled_x(:,1).*rescaled_x(:,3) ;

id=e>0;
rescaled_x2=rescaled_x(id,:);
figure(1);tiledlayout('flow'); set(gcf,'Color','w')
for i=1:6
    nexttile
    plot(X(:,i),'.');title(txt(i))
    nexttile
    plot(rescaled_x(:,i),'.');title(txt(i)); hold on
    plot(rescaled_x2(:,i),'.');title(txt(i))
    xlim([0, length(rescaled_x2)*2])
end
nexttile;plot(e,'.');hold on
plot(e(id),'.');
ylim([0,0.5])
% xlim([0, length(rescaled_x2)*2])

nexttile;plot(rescaled_x2(:,1),rescaled_x2(:,3),'.')
id2 = randsample(length(rescaled_x2),num2);
rescaled_x2=rescaled_x2(id2,:);
e = emax_fun(rescaled_x2(:,4)) - rescaled_x2(:,1).*rescaled_x2(:,3) ;
df=array2table([rescaled_x2,e]);
df.Properties.VariableNames={'vhmax', 'mo', 'ro' ,'CN0','init_C0','aromaticC0','e'};
% figure;corrplot(df)
AR0=rescaled_x2(:,5).*rescaled_x2(:,6);
figure; plot(rescaled_x2(:,5),AR0,'.')
figure; plot(rescaled_x2(:,6),AR0,'.')
%%
close all
i=1;
init_guess = rescaled_x2(i,:);
% init_guess = [0.00096,0.2,13,131.5,0.5,0.0223];
CT_0= init_guess(5);
CO_0= CT_0*init_guess(6);

f=param.f;
a = param.a;
b = param.b;
TT=365*5;

% vomax=param.vomax;
% Create the Yop system
sys = YopSystem( ...
    'states', 2, ...
    'controls', 1, ...
    'parameters', 4);
% Symbolic variables
time = sys.t;
x = sys.x;
u = sys.u;
vh_max = sys.p(1);
mo=sys.p(2);
ro=sys.p(3);
CN0=sys.p(4);
emax=emax_fun(CN0);

CT = x(1);
Co = x(2);
vo = u(1);

L = Co /CT;
Dh = vh_max * g(L, a, b) * (CT-Co);
Do = vo * Co;
e = emax - ro * vo;
G = e * (Dh + Do);
% G = e * (Dh + f*Do)/( vo);

dx = [-Dh-Do +  G; ...
    -Do + mo * G];

sys.set('ode', dx);

% Formulate optimal control problem
ocp = YopOcp();
ocp.max({timeIntegral(G)});
ocp.st( ...
    'systems', sys, ...
    ... % state bounds
    {0, '<=', CT}, ...
    {0, '<=', Co}, ...
    ... %     {0, '<=', Co/(Co+Ch), '<=', L_max}, ...
    ... % control bounds
    {0, '<=', vo, '<=',vomax}, ...
    ... % parameter bounds
    {init_guess(1), '==', vh_max}, ...
    {init_guess(2), '==', mo}, ...
    {init_guess(3), '==', ro}, ...
    {init_guess(4), '==', CN0}, ...
    ... % Initial conditions
    {0, '==', t_0(time)}, ...
    {CT_0, '==', t_0(CT)}, ...
    {CO_0, '==', t_0(Co)}, ...
    ... % Terminal conditions
    {TT, '==', t_f(time)} ...
    ...%     {fCH_0 * ft, '==', t_f(Ch)} ...
    );
% Solving the OCP
tic
sol = ocp.solve( ... %     'initialGuess', initialGuess, ...
    'controlIntervals', 30, ...
    'collocationPoints', 'legendre', ... %
    'polynomialDegree', 3, ...
    'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
    );
toc

time = sol.NumericalResults.Independent;
CT = sol.NumericalResults.State(1, :);
Co = sol.NumericalResults.State(2, :);
vo = sol.NumericalResults.Control(1, :);

emax_fun(init_guess(4)) - init_guess(3).*vo(1)

figure;subplot(211);
plot(time, CT./CT(1), 'linewidth', 2, 'DisplayName', 'CT'); hold on
plot(time, Co./Co(1), 'linewidth', 2, 'DisplayName', 'Co');legend('show')
subplot(212);
plot(time, vo, 'linewidth', 2, 'DisplayName', 'vo'); hold on

%%
close all

% fig=figure;
Ligdecomposition_Starts = zeros(num2, 1);
C_remain_lig_dec_start=Ligdecomposition_Starts;
avgvo = Ligdecomposition_Starts;
init_vo= avgvo;
max_vo=init_vo;
vo_at_Ligdecomposition_Starts=init_vo;
tic
for i =1:num2

    par =rescaled_x2(i,1:4);
    ocp.Parameter(1:4).setBound(par)
    CT_0= init_guess(5);
    CO_0= CT_0*init_guess(6);
    ocp.State(1:2).setInitial([CT_0,CO_0]);
    ocp.parameterize();
    sol = ocp.optimize( ... %     'initialGuess', initialGuess, ...
        'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
        );

    time = sol.NumericalResults.Independent;
    CT = sol.NumericalResults.State(1, :);
    Co = sol.NumericalResults.State(2, :);
    vo = sol.NumericalResults.Control(1, :);
%     figure(fig);clf;subplot(211);
%     plot(time, CT./CT(1), 'linewidth', 2, 'DisplayName', 'CT'); hold on
%     plot(time, Co./Co(1), 'linewidth', 2, 'DisplayName', 'Co');legend('show')
%     subplot(212);
%     plot(time, vo, 'linewidth', 2, 'DisplayName', 'vo'); hold on


    voN = vo ./ max(vo);
    if(isempty(isnan(voN)))
        id=1;
    else
        id = find(voN > vo_thres);
    end
    Ligdecomposition_Starts(i) = time(id(1)) ;
    C_remain_lig_dec_start(i) = CT(id(1))/CT(1);
    avgvo(i) = mean(vo);
    init_vo(i) = vo(1);
    max_vo(i) = max(vo);
    vo_at_Ligdecomposition_Starts(i)=vo(id(1));
    i
    %     pause
end
toc
%%
temp_par = [rescaled_x2, Ligdecomposition_Starts, avgvo,init_vo,max_vo,...
    vo_at_Ligdecomposition_Starts, C_remain_lig_dec_start];
par_sa = array2table(temp_par, ...
    'VariableNames', [txt,{'tauL', 'avg_vo','initvo','max_vo', ...
    'votau','Ctau'}]);
writetable(par_sa, excel_path+"sensitivity.xlsx")

%%
par_sa=readtable(excel_path+"sensitivity.xlsx");
par_sa.AR0=par_sa.CT0.*par_sa.L0;
colname={'vhmax', 'mo', 'ro' ,'CN0','tauL','Ctau','AR0'};% 
figure; 
[SpearmanR, Pvalue]=corrplot(par_sa(:,colname),Type="Pearson",TestR="on", Alpha=0.05);

x=categorical({'\it vh_{max}', '\it m_o', '\it r_o' , ...
    '\it CN_0','\it \tau_L','\it C_\tau','\it ARC_0'});
figure;b=bar(x,[R.("Ctau "),R.("tauL ")]);
b(1).EdgeColor='none';b(2).EdgeColor='none';
ylabel('Pearson correlation');
lh=legend({'\it C_\tau','\it \tau_L'});lh.Location='northwest';lh.Box='on';
set(gca, 'Fontsize',12, 'LineWidth', 0.45, ...
    'Xcolor', [1, 1, 1]*0, 'Ycolor', [1, 1, 1]*0)
set(gcf,'color','w')
grid on 
exportgraphics(gcf, 'results\sensitivity analysis.png','Resolution',300)

% figure; tiledlayout('flow')
% for i=1:size(par_sa,2)
% nexttile;boxplot(par_sa(:,i).Variables)
% xlabel(par_sa(:,i).Properties.VariableNames)
% end

