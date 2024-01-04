%% A Linear Problem With Bang Bang Control
% Create the Yop system
sys = YopSystem('states', 1, 'controls', 1);
% Symbolic variables
t = sys.t;
x = sys.x;
u = sys.u;

% Model
xdot = u;
sys.set('ode',xdot)

% Formulate optimal control problem
ocp = YopOcp();
ocp.min({ timeIntegral( (x-t^2)^2) });
ocp.st(...
     'systems', sys, ...
     ... % Initial conditions
    {  0  '==' t_0( t )   }, ...
    {  1  '==' t_0( x )   }, ...
    ... % Constraints
    {-inf '<=' x '<=' inf }, ...
    { 0 '<='  u   '<=' 4   }, ...
    ... % Terminal conditions
    {  2  '==' t_f( t ) } ...
    );

% Solving the OCP
sol = ocp.solve('controlIntervals', 60);

%% Plot the results
figure
subplot(121); hold on
sol.plot(sys.t, sys.u,'.')
ylim([0,5])
grid on 
xlabel('Time')
ylabel('Control')

subplot(122); hold on
sol.plot(sys.t, sys.x(1))
ylim([0,5])
xlabel('Time')
ylabel('state')
grid on 


