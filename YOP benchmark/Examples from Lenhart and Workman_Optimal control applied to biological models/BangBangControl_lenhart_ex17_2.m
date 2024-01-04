%% A Linear Problem With Bang Bang Control
% Create the Yop system
sys = YopSystem('states', 1, 'controls', 1);
% Symbolic variables
t = sys.t;
x = sys.x;
u = sys.u;

% Model
xdot = x + u;
sys.set('ode',xdot)

% Formulate optimal control problem
ocp = YopOcp();
ocp.max({ timeIntegral( 2*x-3*u ) });
ocp.st(...
     'systems', sys, ...
     ... % Initial conditions
    {  0  '==' t_0( t )   }, ...
    {  5  '==' t_0( x )   }, ...
    ... % Constraints
    {-inf '<=' x '<=' inf }, ...
    { 0 '<='  u   '<=' 2   }, ...
    ... % Terminal conditions
    {  2  '==' t_f( t ) } ...
    );

% Solving the OCP
sol = ocp.solve('controlIntervals', 30);

%% Plot the results
figure
subplot(121); hold on
sol.stairs(sys.t, sys.u)
xlabel('Time')
ylabel('Control')

subplot(122); hold on
sol.plot(sys.t, sys.x(1))
xlabel('Time')
ylabel('state')


