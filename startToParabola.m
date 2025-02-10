%% Starting Drop
% Clear
clc; clear; close all;

% Constants
g = -9.81; % [m/s^2]

% Define X and Z
xDrop = linspace(-63.28,-28.28,1000);
zDrop = -xDrop-28.28;

% Plot Drop
plot3(xDrop, ones(1,length(xDrop)), zDrop, 'b')
grid on

% Calculate length
lengthDrop = (63.28-28.28)*sqrt(2);

% Connecting Arc
% Clear
% clc; clear; close all;

% Define X and Z
xArc = linspace(-28.28,0,1000);
zArc = -sqrt(400-((xArc+14.14).^2))+14.14;

% Plot arc
hold on 
plot3(xArc, ones(1,length(xArc)), zArc, 'b')
grid on

% Calculate Arc Length
lengthArc = 2*pi*20*(90/360);

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
plot3(xParabola, ones(1,length(t)), zParabola, 'b')
grid on;
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

% Calculate Arc Length
f = @(t) sqrt((V0*cosd(theta)).^2 + (V0*sind(theta)+(g.*t)).^2);
lengthParabola = integral(f, t0, tf);

% Transtition to banked turn
% Arc of a circle with radius 50, starting at an angle of -60.6 degrees
% with respect to the verticle
xTurn_Trans = linspace(max(xParabola), 131.9488, 1000);
zTurn_Trans = -sqrt((50^2) - (xTurn_Trans - 131.9448).^2) - 9.2367;
hold on 
plot3(xTurn_Trans, ones(1, length(xTurn_Trans)), zTurn_Trans, 'b')
length_Turn_Trans = (60.6*(pi/180))*50;

% Calculate total length
lengthTot = lengthDrop + lengthArc + lengthParabola + length_Turn_Trans

