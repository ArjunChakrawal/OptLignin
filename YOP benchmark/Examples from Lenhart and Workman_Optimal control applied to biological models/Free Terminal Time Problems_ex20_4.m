clearvars; close all;clc
x0 = 10;
% Create the Yop system
sys = YopSystem(...
    'states', 2, ...
    'controls', 1, ...
    'model', @trolleyModel ...
    );

time = sys.t;
y1 = sys.y;

ocp = YopOcp();
ocp.min({ timeIntegral( 1) });
ocp.st(...
    'systems', sys, ...
    {-1,'<=',y1.acceleration,'<=',1},...
    ... % Initial conditions
    { 0  '==' t_0(time)            }, ...
    {  x0  '==' t_0( y1.position ) }, ...
    {  0  '==' t_0( y1.speed    ) }, ...
    ... % Terminal conditions
    {  0 '<=' t_f( time) '<=' inf }, ...
    {  0  '==' t_f( y1.position ) },...
    {  0  '<=' t_f( y1.speed ) }...
    );

% Solving the OCP
tic
sol = ocp.solve('controlIntervals', 50);
toc
disp("analytical terminal time is " + num2str(2*x0^0.5))
%% Plot the results
figure(1)
subplot(311); hold on
sol.plot(time, y1.position)
xlabel('Time')
ylabel('Position')
grid on
subplot(312); hold on
sol.plot(time, y1.speed)
xlabel('Time')
ylabel('Velocity')
grid on
subplot(313)
sol.stairs(time, y1.acceleration)
xlabel('Time')
ylabel('Acceleration (Control)')
grid on
%%
function [dx, y] = trolleyModel(time, state, control)

position = state(1);
speed = state(2);
acceleration = control;
dx = [speed; acceleration];

y.position = position;
y.speed = speed;
y.acceleration = acceleration;

end