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

% Transition from loop to brake section
% Arc
xStrTrans = linspace(min(xBanked), min(xBanked - 50), 1000);
hold on
plot3(xStrTrans, min(yBanked) * ones(1, length(xStrTrans)), min(zBanked) * ones(1, length(xStrTrans)), 'b')

rTrans = ((125-min(zBanked))/3)+5;
thetaTrans = linspace(pi/2, 2*pi/3, 1000);
xTransLoop = (rTrans.*cos(thetaTrans))+min(xStrTrans);
zTransLoop = (rTrans.*sin(thetaTrans))-rTrans+min(zBanked);
plot3(xTransLoop, min(yBanked)*ones(1,length(xTransLoop)), zTransLoop, 'b')
lengthTransLoop = rTrans*pi/6;

% % G's along arc
 sTrans = linspace(0, rTrans*pi/6, 1000);
 GArcTrans = ((2*(125-min(zBanked)+rTrans*(1-sin((pi/2)-(sTrans/rTrans)))))/rTrans)-sin((pi/2)-(sTrans/rTrans));
% figure
% plot(sTrans, GArcTrans)

% Straight descent
xStrDec = linspace(min(xTransLoop), sqrt(3)*(50*(1-cos(pi/6))-((-1/sqrt(3))*min(xTransLoop)+min(zTransLoop))));
zStrDec = (1/sqrt(3))*xStrDec-(1/sqrt(3))*min(xTransLoop)+min(zTransLoop);
plot3(xStrDec, min(yBanked)*ones(1,length(xStrDec)), zStrDec, 'b')
lengthStrDec = sqrt(((max(xStrDec)-min(xStrDec))^2+(max(zStrDec)-min(zStrDec))^2));

% Arc to z=0
rArcFinal = 50;
thetaArcFinal = linspace(-pi/3, -pi/2, 1000);
xArcFinal = rArcFinal.*cos(thetaArcFinal)+min(xStrDec)-50*sin(pi/6);
zArcFinal = rArcFinal.*sin(thetaArcFinal)+50;
plot3(xArcFinal, min(yBanked)*ones(1,length(xArcFinal)), zArcFinal, 'b')
lengthArcFinal = rArcFinal*pi/6;

% Calculate total length
lengthTot = lengthDrop + lengthArc + lengthParabola+ length_Turn_Trans + S_banked+lengthLoop + lengthTransLoop + lengthStrDec + lengthArcFinal


%% Calculate Loop G's
 sLoop = linspace(0,2*pi*rLoop);
 gLoop = ((2.*(125-(min(zBanked)+rLoop.*(1-cos(sLoop./rLoop)))))/rLoop)-sin((sLoop./rLoop)+(3*pi/2));
 figure
 plot(sLoop, gLoop)
 grid on
 xlabel('Arc Length')
 ylabel("G's")
 xlim([0 2*pi*rLoop])
 ylim([-1 6])
 title("G's vs arc length, loop section")

%% plots for gs during banked turn
pathlength_banked = linspace(0,S_banked,100);
gs_banked_lateral_vec = linspace(gs_banked_lateral,gs_banked_lateral,100);
gs_banked_up_vec = linspace(gs_banked_up,gs_banked_up,100);

figure();
subplot(2,1,1)
hold on
plot(pathlength_banked, gs_banked_lateral_vec);
xlabel('pathlength');
ylabel('gs lateral')
title('gs lateral for banked turn');
hold off

subplot(2,1,2);
hold on
plot(pathlength_banked,gs_banked_up_vec);
xlabel('pathlength');
ylabel('gs up through the seat');
title('gs up through the seat for the banked turn');
hold off

%% gs calculations

%theta for initial section section
z_delta = zDrop(1,1) - zDrop(1,1000);
x_delta = xDrop(1,1) - xDrop(1,1000);
theta_b = abs(atand(z_delta / x_delta));

%velocity function for second section
h0_red = 125 - zDrop(1,1); %z-position where red section starts
v = sqrt(2*9.81 * (125 - h0_red - zArc)); %should it be -zArc or plus ???
r_red = 20; %radius of part of loop leading in into the parabola

%first section
gs_up_blue = cosd(theta_b);
gs_forward_blue = sind(theta_b);

%second section
for i=1:999
    x_change(i) = xArc(i+1)-xArc(i);
    z_change(i) = zArc(i+1)-zArc(i);
end

theta_red = atand((abs(x_change)) ./ (abs(z_change)));

%gs second section
gs_up_red = (v(1:999).^2)/(r_red*9.81) + sind(theta_red);
gs_forward_red = cosd(theta_red(1:499));
gs_back_red = cosd(theta_red(501:999));

%second section max
gs_up_red_max = max(gs_up_red);
gs_forward_red_max = max(gs_forward_red);
gs_back_red_max = max(gs_back_red);

% path length
%first section
l_b = sqrt( xDrop(1:1)^2 + zDrop(1:1)^2);
length_blue = linspace(0,l_b,1000);

gs_up_blue_vec = zeros(1,length(length_blue));
gs_up_blue_vec(:) = gs_up_blue;
gs_forward_blue_vec = zeros(1,length(length_blue));
gs_forward_blue_vec(:) = gs_forward_blue;

for i=1:999
    l_r(i) = sqrt( (xArc(i+1)-xArc(i))^2 + (zArc(i+1)-zArc(i))^2 );
end
l_r = sum(l_r,'all');
length_red = linspace(0,l_r,length(gs_up_red));

%% total figures

pathlength_para = linspace(0,lengthParabola,100);
gs_up_para = linspace(0,0,100);

figure();
%gs up through riders seat
hold on
plot(length_blue, gs_up_blue_vec);
plot(length_red + l_b, gs_up_red);
plot(l_b + l_r + pathlength_para,gs_up_para);
%add in transition from parabola to banked turn
plot(l_b + l_r + lengthParabola + pathlength_banked, gs_banked_up_vec);
plot(l_b + l_r + lengthParabola + pathlength_banked(1,100) + sLoop, gLoop);
plot(l_b + l_r + lengthParabola + pathlength_banked(1,100) + sLoop(1,100) + sTrans,GArcTrans);

xline([l_b l_b+l_r l_b+l_r+lengthParabola],'--');
ylabel('gs up');
xlabel('pathlength');
title('gs up through seat');

figure();
%gs from seat restraint(blue and red section)
hold on
plot(length_blue, abs(gs_forward_blue_vec));
plot(length_red(1:499) + l_b, gs_forward_red);
plot([(.5*l_r)+l_b ; l_r+l_b], [0;0]);
ylabel('gs_forward');
xlabel('pathlength');
title('gs back from seat restraint');

figure();
%gs from back of seat pushing
hold on
plot([0;(l_b+(l_r*.5))], [0;0]);
plot(length_red(501:999) + l_b, gs_back_red);
ylabel('gs back of seat');
xlabel('pathlength');
title('gs through back of seat');

function v = calcVelocity(h0, h)
v = sqrt(2*9.81*(h0-h))
end
