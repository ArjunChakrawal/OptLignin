clearvars
close all
clc
% format short

%% Solve the optimal control problem using Pontryagin's maximum principle
% necessary conditions for optimality
% dhdu = -dlambdadt
% dHdu=0


% State equations
% state costate and control
syms C  p  u
% system parameters
syms b r e g mu;

% Cost function inside the integral
G = e * b * (u - r) / (u + b);

% System dynamics
dC = (-u - g * C + mu * G);

% Hamiltonian
H = G + p * dC;


% Costate or adjoint
dHdC = diff(H, C); % eq 3
% state
dHdu = diff(H, u);
dHdp = diff(H, p); % gives back the dynamic constraint (not needed)

% solve for control u
u_opt = solve(dHdu == 0, u);
u_opt = u_opt(2);
u_fun = matlabFunction(u_opt);

% substitute the value of optimal u in dxdt and dpdt
dCdt = subs(dC, u, u_opt);
dpdt = subs(-dHdC, u, u_opt);
H_T = subs(H, u, u_opt);

syms C(t)  p(t) T
dCdt = subs(dCdt, [C, p], [C(t), p(t)]);
dpdt = subs(dpdt, [C, p], [C(t), p(t)]);

% following lines convert above system of odes to anonymous function
eqs = [diff(C(t), t) == T* dCdt, ...
    diff(p(t), t) == T * dpdt];

vars = [C(t), p(t)];
[M, F] = massMatrixForm(eqs, vars);
f = M \ F;
odefun = odeFunction(f, vars, T, b, r, e, g, mu);

% get the Hamiltonian as a fucntion of final terminal state and costate
syms xtf pf
H_T = subs(H_T, [C, p], [xtf, pf]);
H_Tfun = matlabFunction(H_T);

%% Set parameter values
b = 0.5;g = 0.5;r = 0.1;e = 0.5;mu = 0.5;
C0=1;
CT=0;
sol= bvp4c_solver(b, r, e, g, mu,odefun,H_Tfun,C0,CT);
[time,G,u]=make_plots(sol,u_fun,b, r, e, g, mu);

disp("J is " + trapz(time,G))
disp("H_T is "+H_Tfun(b,e,g,mu,sol.y(2, end),r,sol.y(1, end)))

% b = 0.5;g = 0.1;r = 0.1;e = 0.5;mu = 0.5;
% sol= bvp4c_solver(b, r, e, g, mu,odefun,H_Tfun,C0,CT);
% make_plots(sol,u_fun,b, r, e, g, mu)
% 
% b = 0.5;g = 0.1;r = 0.1;e = 0.5;mu = 0.5;
% sol= bvp4c_solver(b, r, e, g, mu,odefun,H_Tfun,C0,CT);
% make_plots(sol,u_fun,b, r, e, g, mu)
%%
function sol= bvp4c_solver(b, r, e, g, mu,odefun,H_Tfun,C0,CT)
solinit = bvpinit(linspace(0, 1, 100), ...
    [0.1, 0.1],10);

options = bvpset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Stats', 'on');
BVP_bc2 = @(ya,yb, b, g, r, e) [ya(1) - C0; ... % C(0)=1
    yb(1)-CT; ... %C(Tf)=0
    H_Tfun(b, e, g, mu, yb(2), r, yb(1))...%H(Tf)=0
    ];

BVP_bc = @(ya, yb,T) BVP_bc2(ya,yb, b, g, r, e);
BVP_ode = @(t, y,T)odefun(t, y, T,b, r, e, g, mu);

sol = bvp4c(BVP_ode, BVP_bc, solinit, options);
end

function [time,gg, ut]=make_plots(sol,u_fun,b, r, e, g, mu)
y = sol.y;
time = sol.parameters.*sol.x;
p=sol.y(2,:);
ut = u_fun(b,e,mu,p,r);
G = @(u) e .* b .* (u - r) ./ (u + b);

figure(1)
set(gcf,'color','w')
subplot(2, 2, 1)
plot(time, sol.y(1, :), '-', 'LineWidth', 2, ...
    'DisplayName',"\beta="+b+" \gamma="+g+" \mu ="+mu + " r="+r); hold on
ylabel('C(t)')
xlabel('time');
grid on
hl=legend('show');
hl.Box ='off';
subplot(2, 2, 2)
plot(time, sol.y(2, :), '-', 'LineWidth', 2); hold on
ylabel('\lambda(t)')
xlabel('time');
grid on
subplot(2, 2, 3)
plot(time, ut, '-', 'LineWidth', 2); hold on
xlabel('time'); ylabel('u_{opt}')
grid on
subplot(2, 2, 4)
plot(time, G(ut), '-', 'LineWidth', 2); hold on
xlabel('time'); ylabel('G_{opt}')
grid on
saveas(gcf, 'indirect_bvp4c_Fig1.png')

figure(2)
set(gcf,'color','w')
subplot(3,1, 1)
plot(sol.y(1, :), G(ut), '-', 'LineWidth', 2); hold on
xlabel('C');
ylabel('G_{opt}');
grid on
subplot(3,1, 2)
plot(sol.y(1, :), ut, '-', 'LineWidth', 2); hold on
xlabel('C');
ylabel('u_{opt}');
grid on
subplot(3,1, 3)
plot(ut, G(ut), '-', 'LineWidth', 2); hold on
xlabel('u_{opt}');
ylabel('G_{opt}');
grid on
saveas(gcf, 'indirect_bvp4c_Fig2.png')
gg =G(ut);



end