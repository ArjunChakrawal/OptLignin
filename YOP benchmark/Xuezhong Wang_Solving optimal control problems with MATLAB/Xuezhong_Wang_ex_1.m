clearvars; close all;clc

% Create the Yop system
bdSystem = YopSystem(...
    'states', 2, ...
    'controls', 1, ...
    'model', @trolleyModel ...
    );

time = bdSystem.t;
trolley = bdSystem.y;

ocp = YopOcp();
ocp.min({ timeIntegral( 1/2*trolley.acceleration^2 ) });
ocp.st(...
    'systems', bdSystem, ...
    ... % Initial conditions
    {  0  '==' t_0( trolley.position ) }, ...
    {  0  '==' t_0( trolley.speed    ) }, ...
    ... % Terminal conditions
    {  2  '==' t_f( time ) }, ...
    {  5  '==' t_f( trolley.position ) }, ...
    { 2  '==' t_f( trolley.speed    ) }...
    );

% Solving the OCP
tic
sol = ocp.solve('controlIntervals', 25);
toc
% Plot the results
figure(1)
subplot(211); hold on
sol.plot(time, trolley.position)
xlabel('Time')
ylabel('Position')
grid on
subplot(212); hold on
sol.plot(time, trolley.speed)
xlabel('Time')
ylabel('Velocity')
grid on
figure(2); hold on
sol.stairs(time, trolley.acceleration)
xlabel('Time')
ylabel('Acceleration (Control)')
grid on

function [dx, y] = trolleyModel(time, state, control)

position = state(1);
speed = state(2);
acceleration = control;
dx = [speed; -speed+acceleration];

y.position = position;
y.speed = speed;
y.acceleration = acceleration;

end