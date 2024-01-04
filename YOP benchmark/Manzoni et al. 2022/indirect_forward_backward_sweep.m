%-------------------------------------------------------------------------%
% Author: Arjun Chakrawal
% Copyright (cc-by 4.0).
%-------------------------------------------------------------------------%

clearvars
close all
clc
format short

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
syms H;
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

% Costate dynamics
dpdt = matlabFunction(-dHdC);

% get the Hamiltonian as a fucntion state and costate
H_Tfun = matlabFunction(H);

%%
b = 0.5;
e = 0.5;
r = 0.1*0 ;
g = 0.5;
mu = 0.5 ;
C0 =1;
CT=0;
u_ub = 1;
u_lb = 0;
% G as a function of u
G = @(u) e .* b .* (u - r) ./ (u + b);

% dCdt as a function of C and u
dCdt = @(C,u) (-u - g * C + mu * G(u)); %eq1
% use root finding algorithm to multishoot for final lamda and terminal
% time
myfun = @ (x)Ham_zero_fun(x,b,e,r,g,mu,C0,u_ub,u_lb,dCdt,dpdt,u_fun,H_Tfun,CT);
T = fzero(myfun,[1,2]);
%%
myfun = @ (x)fun(x,b,e,r,g,mu,C0,u_ub,u_lb,dCdt,dpdt,u_fun,T,CT);
theta = fzero(myfun, 0.5);


%%
y1 = sweep(theta,b,e,r,g,mu,C0,u_ub,u_lb,dCdt,dpdt,u_fun,T);
time= y1(1, :);
Copt = y1(2,:);
lambda_ = y1(4,:);
uopt = y1(3,:);

figure
set(gcf, 'color', 'w')
subplot(2, 2, 1)
plot(time, Copt, '-', 'LineWidth', 2); hold on
ylabel('C(t)')
xlabel('time');
grid on
subplot(2, 2, 2)
plot(time, lambda_, '-', 'LineWidth', 2); hold on
ylabel('\lambda(t)')
xlabel('time');
grid on
subplot(2, 2, 3)
plot(time, uopt, '-', 'LineWidth', 2); hold on
xlabel('time'); ylabel('u_{opt}')
grid on
subplot(2, 2, 4)
plot(time, G(uopt), '-', 'LineWidth', 2); hold on
xlabel('time'); ylabel('G_{opt}')
grid on
saveas(gcf, 'indirect_forward_backward_sweep_Fig1.png')

figure
set(gcf, 'color', 'w')
subplot(3, 1, 1)
plot(Copt, G(uopt), '-', 'LineWidth', 2); hold on
xlabel('C');
ylabel('G_{opt}');
grid on
subplot(3, 1, 2)
plot(Copt, uopt, '-', 'LineWidth', 2); hold on
xlabel('C');
ylabel('u_{opt}');
grid on
subplot(3, 1, 3)
plot(uopt, G(uopt), '-', 'LineWidth', 2); hold on
xlabel('u_{opt}');
ylabel('G_{opt}');
grid on
saveas(gcf, 'indirect_forward_backward_sweep_Fig2.png')

%%

function y = sweep(theta,b,e,r,g,mu,C0,u_ub,u_lb,dCdt,dpdt,u_fun,T)
% forward backward sweep
test = -1;
delta = 0.0001;
N = 1000;
t = linspace(0,T,N+1);
h = T/N;
h2 = h/2;

C = zeros(1,N+1);
C(1) = C0;


lambda = zeros(1,N+1);
lambda(N+1) = theta;
u = zeros(1,N+1);

while(test < 0)
    
    oldu = u;
    oldC = C;
    oldlambda = lambda;
    
    for i=1:N
        m11 =  dCdt(C(i),u(i));
        m21 = dCdt(C(i)+m11*h2,u(i));
        m31 = dCdt(C(i)+m21*h2,u(i));
        m41 = dCdt(C(i)+m31*h,u(i));
        C(i+1) = C(i) + (h/6)*(m11 + 2*m21 + 2*m31 + m41);
    end
    
    for i = 1:N
        j = N + 2 - i;
        m11 = dpdt(g,lambda(j));
        m21 = dpdt(g,lambda(j)+m11*h2);
        m31 = dpdt(g,lambda(j)+m21*h2);
        m41 = dpdt(g,lambda(j)+m31*h);
        lambda(j-1) = lambda(j) - (h/6)*(m11 + 2*m21 + 2*m31 + m41);
    end
        
    temp = u_fun(b,e,mu,lambda,r);
    u1 = min(u_ub,max(u_lb,temp));
    u = 0.5*(u1 + oldu);
    
    temp1 = delta*sum(abs(u)) - sum(abs(oldu - u));
    temp2 = delta*sum(abs(C)) - sum(abs(oldC - C));
    temp3 = delta*sum(abs(lambda)) - sum(abs(oldlambda - lambda));
    test = min(temp1, min(temp2, min(temp3)));
end

y(1,:) = t;
y(2,:) = C;
y(3,:) = u;
y(4,:) = lambda;
end



function res = fun(theta,b,e,r,g,mu,C0,u_ub,u_lb,dCdt,dpdt,u_fun,T,B)
% find root for terminal lambda
z = sweep(theta,b,e,r,g,mu,C0,u_ub,u_lb,dCdt,dpdt,u_fun,T);
res = z(2, 1001) - B;
end


function res = Ham_zero_fun(T,b,e,r,g,mu,C0,u_ub,u_lb,dCdt,dpdt,u_fun,H_Tfun,CT)
% find root for terminal time

myfun = @ (x)fun(x,b,e,r,g,mu,C0,u_ub,u_lb,dCdt,dpdt,u_fun,T,CT);
theta = fzero(myfun, 0.5);

z = sweep(theta,b,e,r,g,mu,C0,u_ub,u_lb,dCdt,dpdt,u_fun,T);
uT=z(3,1001);
res =H_Tfun(CT,b,e,g,mu,theta,r,uT);

end














