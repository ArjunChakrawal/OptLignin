function sol = yop_function(vh_max, param, init_cond, vo_max, T,mo)

ro = param.ro;
L_max = param.L_max;
alpha = param.alpha;
emax = param.emax;
ft= param.ft;
mu = param.mu;
fCH_0 = init_cond(1); % initial fraction of CH in litter[-]
fCO_0 = init_cond(2); % initial fraction of CO in litter[-]
fCH_T=fCH_0*ft;

% Create the Yop system
sys = YopSystem( ...
    'states', 2, ...
    'controls', 1);
% Symbolic variables
time = sys.t;
x = sys.x;
u = sys.u;

Ch = x(1);
Co = x(2);
vo = u(1);
L = Co / (Ch + Co);
gL = (L / L_max)^alpha;
Dh = vh_max * (1- gL)*Ch;
Do = vo*Co;    
e= emax - ro * vo;
G = e * (Dh + Do);
dx = [-Dh + (1 - mo) * mu* G; ...
    -Do + mo * mu * G];
y.Ch = Ch;
y.Co = Co;
y.vo = vo;

sys.set('ode', dx);
sys.set('y', y);

% Formulate optimal control problem
ocp = YopOcp();
ocp.max({timeIntegral(G)});
ocp.st( ...
    'systems', sys, ...
    ... % state bounds
    {0, '<=', Ch}, ...
    {0, '<=', Co}, ...
    {0, '<=', Co/(Co+Ch), '<=', L_max}, ...
    ... % control bounds
    {0, '<=', vo, '<=', vo_max}, ...
    ... % Initial conditions
    {0, '==', t_0(time)}, ...
    {fCH_0, '==', t_0(Ch)}, ...
    {fCO_0, '==', t_0(Co)}, ...
    ... % Terminal conditions
    {0, '<=', t_f(time), '<=', T}, ...
    {fCH_T, '==', t_f(Ch)} ...
    ... %     {0.01, '==', t_f(Co)} ...
    );

% Solving the OCP
sol = ocp.solve( ... %     'initialGuess', initialGuess, ...
    'controlIntervals', 200, ...
    'collocationPoints', 'legendre', ... %
    'polynomialDegree', 3, ...
    'ipopt', struct('max_iter', 2000,"print_level",1,"max_cpu_time",50) ...
    );

% sol = ocp.solve('controlIntervals', 50, 'collocationPoints', 'radau', ...
%     'ipopt', struct('max_iter', 5000));
% sol = ocp.solve('initialGuess', sol);
% sol = ocp.solve('controlIntervals', 500);
   

end