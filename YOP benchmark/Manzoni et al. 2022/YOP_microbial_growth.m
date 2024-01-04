%-------------------------------------------------------------------------%
% Author: Arjun Chakrawal
% Copyright (cc-by 4.0).
%-------------------------------------------------------------------------%
clearvars;
% close all;
clc

% parameters
tic
beta = 0.5;
gamma = 0.5;
r = 0.1;
emax = 0.5;
mu = 0.5;
k=0.0;
G = @(u,emax,r,k) beta .* emax .* ((u - r) ./ (u + beta))-k*u;
% figure(3);
% u=0:0.01:2.5;
% plot(u,G(u,emax,0.1,0.1));hold on 
% plot(u,G(u,emax,0.1,0.2));
% plot(u,G(u,emax,0.1,0.3));
% grid on

% test for a value of gamma
param = [beta, emax, r, mu, gamma,k];

sol = yop_opt(param);
time = sol.Variables.Independent;

% Plot the results
% figure
% sol.plot(sol.NlpSolution.x,sol.NlpSolution.lam_x, 'linewidth', 2)

figure
% clf
subplot(311); hold on
sol.plot(time, sol.Variables.State, 'linewidth', 2)
xlabel('Time')
ylabel('Cs')
grid on
subplot(312); hold on
sol.stairs(time, sol.Variables.Control, 'linewidth', 2)
xlabel('Time')
ylabel('uptake rate (Control)')
grid on
subplot(313); hold on
sol.stairs(time, sol.Variables.Control, 'linewidth', 2)
xlabel('Time')
ylabel('uptake rate (Control)')
grid on

growth = G(sol.NumericalResults.Control,emax,r,k);
trapz(sol.NumericalResults.Independent,growth)
figure;
% clf;
hold on
plot(sol.NumericalResults.Independent, growth, 'o')
xlim([0, inf])
xlabel('Time')
ylabel('G opt')
grid on
