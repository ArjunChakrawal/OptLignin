clearvars;
close all;
clc

% parameters
mO = 0.2;  
Kr = 0.01;      
L_max = 1;
alpha = 2.0;
emax = 0.5;
mu=0.25;
ro =2;
vh=0.1;
fo = 0.1;
param = [mO, L_max,alpha, emax,mu,vh,ro,fo];
% microbial operating parameters
control=0.1; %vO
% simulation parameters
fCH_0 = 0.5;       % initial fraction of CH in litter[-]
fCO_0 = 0.25;       % initial fraction of CO in litter[-]
t_initial = 0.0; % initial time [day]
t_final = 400.0;  % final time
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
myode=@(time,state)mass_bal(time, state, control, param);

[time,State] = ode23s(myode,[t_initial t_final],[fCH_0; fCO_0],options);
Ch = State(:,1);
Co = State(:,2);
vo= control(1);


L= Co./(Ch+Co);
go = (L/L_max).^alpha;
Dh = vh.*Ch.*(1-go);
Do = vo.*Co*fo;
e = (emax- ro*vo).*ones(length(time),1);

G = e.*(Dh+Do);

int_G = cumtrapz(time,G);

% time=time./max(time);
% Plot the results from dynamic simulation
figure
clf
t = tiledlayout(2,4,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(time, State(:,1), 'linewidth', 2,'DisplayName','CH'); hold on 
plot(time,  State(:,2), 'linewidth', 2,'DisplayName','CO');
legend('show')
xlabel('Time/Tmax')
grid on

nexttile
plot(time, G, 'linewidth', 2)
xlabel('Time')
ylabel('Growth rate')
grid on

nexttile
plot(time, e, 'linewidth', 2)
xlabel('Time')
ylabel('CUE')
grid on

nexttile
plot(L, go, 'linewidth', 2,'DisplayName','gO');hold on 
plot(L, 1-go, 'linewidth', 2,'DisplayName','gH')
xlabel('L')
legend('show')
grid on

nexttile
plot(time, L, 'linewidth', 2); hold on
plot(time,1- L, 'linewidth', 2)
xlabel('Time')
ylabel('L')
grid on

nexttile
plot(time, int_G, 'linewidth', 2)
xlabel('Time')
ylabel('\intgrowth rate')
grid on

nexttile
plot(time, Dh, 'linewidth', 2);hold on 
plot(time, Do, 'linewidth', 2)
xlabel('Time')
ylabel('D')
legend('DH','DO')
grid on


%%

function dstate_dt = mass_bal(time, state, control, param)
% global param
% param = [mO,M,Kr, L_max,alpha, emax];

mO = param(1);
L_max = param(2);
alpha = param(3);
emax = param(4);
mu = param(5);
vh = param(6); 
ro= param(7); 
fo= param(8); 

Ch = state(1);
Co = state(2);
vo= control(1);

L= Co/(Ch+Co);
gL = (L/L_max)^alpha;
Dh = vh*Ch*(1-gL);
Do = vo*Co*fo;
e = emax- ro*vo;
dstate_dt = [-Dh + (1-mO)*mu*e*Dh; ...
    -Do + mO*mu*e*Dh];

end