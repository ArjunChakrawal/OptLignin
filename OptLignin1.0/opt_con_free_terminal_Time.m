function [ocp,sol] = opt_con_free_terminal_Time(param,g,init_guess)
ft=0.1;
emax= param.emax;
f=param.f;
a = param.a;
b = param.b;
CO_0= param.CO_0;
CT_0= param.CT_0;
Tmax=param.Tmax;
vomax=param.vomax;
% Create the Yop system
sys = YopSystem( ...
    'states', 2, ...
    'controls', 1, ...
    'parameters', 4);
% Symbolic variables
time = sys.t;
vh_max = sys.p(1);
mo=sys.p(2);
ro=sys.p(3);
Terminal_time=sys.p(4);

CT = sys.x(1);
Co = sys.x(2);
vo = sys.u(1);

L = Co /CT;
Dh = vh_max * g(L, a, b) * (CT-Co);
Do = vo * Co;
e = emax - ro * vo;
G = e * (Dh + f*Do);
dx = Terminal_time.*[-Dh-Do +  G; ...
    -Do + mo * G];
y.CT = CT;
y.Co = Co;
y.vo = vo;

sys.set('ode', dx);
sys.set('y', y);

% Formulate optimal control problem
ocp = YopOcp();
ocp.max({timeIntegral(Terminal_time*G)});
ocp.st( ...
    'systems', sys, ...
    ... % state bounds
    {0, '<=', CT}, ...
    {0, '<=', Co}, ...
    ... %     {0, '<=', Co/(Co+Ch), '<=', L_max}, ...
    ... % control bounds
    {0, '<=', vo, '<=',vomax}, ...
    {0, '<=', Terminal_time, '<=',Tmax}, ...
    ... % parameter bounds
    {init_guess(1), '==', vh_max}, ...
    {init_guess(2), '==', mo}, ...
    {init_guess(3), '==', ro}, ... 
    ... % Initial conditions
    {0, '==', t_0(time)}, ...
    {CT_0, '==', t_0(CT)}, ...
    {CO_0, '==', t_0(Co)}, ...
    ... % Terminal conditions
    {1, '==', t_f(time)}, ...
    {CT_0 * ft, '==', t_f(CT)} ...
    );
% Solving the OCP
tic
sol = ocp.solve( ... %     'initialGuess', initialGuess, ...
    'controlIntervals', 50, ...
    'collocationPoints', 'legendre', ... %
    'polynomialDegree', 5, ...
    'ipopt', struct('max_iter', 2000, "print_level", 1, "max_cpu_time", 50) ...
    );
toc
% disp("The final value of objective function is " +-sol.NumericalResults.Objective)
% disp("The final time is " +sol.NumericalResults.Independent(end))

end