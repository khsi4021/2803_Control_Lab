%% Starting Drop
% Clear
clc; clear; close all;

% Constants
g = -9.81; % [m/s^2]

% Define X and Z
xDrop = linspace(-63.28,-28.28,1000);
zDrop = -xDrop-28.28+90;

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
zArc = -sqrt(400-((xArc+14.14).^2))+14.14+90;

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
zParabola = z0 + (V0.*t*sind(theta)) + (0.5*g.*(t.^2))+90;

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
zTurn_Trans = -sqrt((50^2) - (xTurn_Trans - 131.9448).^2) - 9.2367+90;
hold on 
plot3(xTurn_Trans, ones(1, length(xTurn_Trans)), zTurn_Trans, 'b')
length_Turn_Trans = (60.6*(pi/180))*50;

% clear
% clc
% close all

%banked turn
rBanked = 50; %radius of turn
g = 9.81;
vBanked = 43; %velocity
thetaBanked = 70; %degrees from upright of roller cart

%placeholder values, delete later
m = 10;

F_r = (m*vBanked.^2)/rBanked; %force in radial direction
F_g = m*g; %force of gravity
F_N = F_r * sind(thetaBanked) - F_g * cosd(thetaBanked); %force normal to the seat of the cart
F_L = F_r * cosd(thetaBanked) + F_g * sind(thetaBanked); %force perpendicular to the normal

gs_banked_up = F_N / (m*g); %gs through the seat
    %limit is 6g
gs_banked_lateral = (F_L)/(m*g); %gs through sidebar
    %limit is 3g

S_banked = pi*rBanked; %path length

theta_circ = 0:pi/50:pi; %in radians
yBanked = (rBanked * cos(theta_circ))-49;
xBanked = (rBanked * sin(theta_circ))+max(xTurn_Trans);
zBanked = zeros(1,51)+min(zTurn_Trans);

hold on
plot3(xBanked,yBanked,zBanked, 'b');

% Loop
rLoop = 40;
thetaLoop = linspace(0, 2*pi, 1000);
xLoop = (rLoop .* cos(thetaLoop)) + min(xBanked);
yLoop = min(yBanked) * ones(1,length(thetaLoop));
zLoop = (rLoop .* sin(thetaLoop))+min(yBanked)+2*rLoop+90;
lengthLoop = 2*pi*rLoop;

hold on
plot3(xLoop, yLoop, zLoop, 'b')

% Reform Axes
xlim([-100 200])
ylim([-150 150])
zlim([-150 150])

% Calculate Loop G's
sLoop = linspace(0,2*pi*rLoop);
gLoop = ((2.*(125-(min(zBanked)+rLoop.*(1-cos(sLoop./rLoop)))))/rLoop)-sin((sLoop./rLoop)+(3*pi/2));
figure
plot(sLoop, gLoop)
grid on
xlabel('Arc Length')
ylabel("G's")
xlim([0 2*pi*rLoop])
ylim([-1 6])

% Calculate total length
lengthTot = lengthDrop + lengthArc + lengthParabola + length_Turn_Trans + S_banked+lengthLoop

