%% Starting Drop
% Clear
clc; clear; close all;

% Constants
g = -9.81; % [m/s^2]

% Define X and Z
xDrop = linspace(-63.28,-28.28,1000);
zDrop = -xDrop-28.28;

% Plot Drop
plot3(xDrop, ones(1,length(xDrop)), zDrop)
grid on

% Connecting Arc
% Clear
% clc; clear; close all;

% Define X and Z
xArc = linspace(-28.28,0,1000);
zArc = -sqrt(400-((xArc+14.14).^2))+14.14;

% Plot arc
hold on 
plot3(xArc, ones(1,length(xArc)), zArc)
grid on

% Calculate Arc Length
lengthArc = 2*pi*20*(90/360)

% Zero G Parabola
% Define Constants
g = -9.81; % [m/s^2]
t0 = 0; % [s]
tf = 5; % [s]
t = linspace(t0,tf,10000); % [s]
V0 = 25; % [m/s]
theta = 45; % [degrees]
x0 = 0; % [m]
z0 = 0; % [m]

% Calculate X and Z
xParabola = (V0.*t*cosd(theta)) + x0;
zParabola = z0 + (V0.*t*sind(theta)) + (0.5*g.*(t.^2));

% Plot Parabola
hold on
plot3(xParabola, ones(1,length(t)), zParabola)
grid on;
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

% Calculate Arc Length
f = @(t) sqrt((V0*cosd(theta)).^2 + (V0*sind(theta)+(g.*t)).^2);
lengthParabola = integral(f, t0, tf)

% Calculate G's

