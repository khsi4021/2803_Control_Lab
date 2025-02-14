%% Starting Drop
% Clear
clc; clear; close all;
% Constants
g = -9.81; % [m/s^2]
% Define X and Z
xDrop = linspace(-63.28,-28.28,1000);
yDrop = zeros(1, 1000);
zDrop = -xDrop-28.28+90;
vDrop = sqrt(2 * -g * (125 - zDrop));
X = [xDrop; xDrop];
Y = [yDrop; yDrop];
Z = [zDrop; zDrop];
C = [vDrop; vDrop];
% Plot Drop
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet); 
grid on
% Calculate length
lengthDrop = (63.28-28.28)*sqrt(2);
% Connecting Arc
% Clear
% clc; clear; close all;
% Define X and Z
xArc = linspace(-28.28,0,1000);
yArc = zeros(1, 1000);
zArc = -sqrt(400-((xArc+14.14).^2))+14.14+90;
vArc = sqrt(2 * -g .* (125 - zArc));
X = [xArc; xArc];
Y = [yArc; yArc];
Z = [zArc; zArc];
C = [vArc; vArc];
% Plot arc
hold on 
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet); 
grid on
% Calculate Arc Length
lengthArc = 2*pi*20*(90/360);
% Zero G Parabola
% Define Constants
g = -9.81; % [m/s^2]
t0 = 0; % [s]
tf = 5.35; % [s]
t = linspace(t0,tf,10000); % [s]
V0 = calcVelocity(125, max(zArc)); % [m/s]
theta = 51.9; % [degrees]
x0 = 0; % [m]
z0 = 0; % [m]
% Calculate X and Z
xParabola = (V0.*t*cosd(theta)) + x0;
yParabola = zeros(1, 10000);
zParabola = z0 + (V0.*t*sind(theta)) + (0.5*g.*(t.^2))+90;
vParabola = sqrt(2 * -g .* (125 - zParabola));
X = [xParabola; xParabola];
Y = [yParabola; yParabola];
Z = [zParabola; zParabola];
C = [vParabola; vParabola];
% Plot Parabola
hold on
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet);
grid on;
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
% Calculate Arc Length
f = @(t) sqrt((V0*cosd(theta)).^2 + (V0*sind(theta)+(g.*t)).^2);
lengthParabola = integral(f, t0, tf);
%% Arc Trans. Arc of a circle with radius 50, starting at an angle of -60.6 degrees
% with respect to the verticle
r_Arc_Trans = 50;
xTurn_Trans = linspace(max(xParabola), 131.9488, 1000);
yTurn_Trans = zeros(1, 1000);
zTurn_Trans = -sqrt((50^2) - (xTurn_Trans - 131.9448).^2) - 9.2367+90;
vTurn_Trans = sqrt(2 * -g .* (125 - zTurn_Trans));
X = [xTurn_Trans; xTurn_Trans];
Y = [yTurn_Trans; yTurn_Trans];
Z = [zTurn_Trans; zTurn_Trans];
C = [vTurn_Trans; vTurn_Trans];
hold on 
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet);
length_Turn_Trans = (60.6*(pi/180))*50;
theta_Arc_Trans = atand(xTurn_Trans ./ zTurn_Trans);
v_Arc_Trans = calcVelocity(125, zTurn_Trans);
length_Turn_Trans_vec = linspace(0,length_Turn_Trans,1000);
for i=1:999
    atan_Arc_Trans(i) = ((v_Arc_Trans(1,i+1))^2 - (v_Arc_Trans(1,i))^2) / (2*length_Turn_Trans_vec(i));%calculating tangent acceleration
end
gs_up_Arc_Trans = cosd(theta_Arc_Trans) + (v_Arc_Trans.^2)/(r_Arc_Trans*9.81);
gs_back__Arc_Trans = sind(theta_Arc_Trans(1:999)) + atan_Arc_Trans./9.81;
%% banked turn
rBanked = 50; %radius of turn
g = 9.81;
thetaBanked = 70; %degrees from upright of roller cart
%placeholder values, delete later
m = 10;
theta_circ = 0:pi/50:pi; %in radians
yBanked = (rBanked * cos(theta_circ))-50;
xBanked = (rBanked * sin(theta_circ))+max(xTurn_Trans);
zBanked = zeros(1,51)+min(zTurn_Trans);
vBanked = calcVelocity(125,zBanked(1,1)); %velocity
F_r = (m*vBanked.^2)/rBanked; %force in radial direction
F_g = m*g; %force of gravity
F_N = F_r * sind(thetaBanked) - F_g * cosd(thetaBanked); %force normal to the seat of the cart
F_L = F_r * cosd(thetaBanked) + F_g * sind(thetaBanked); %force perpendicular to the normal
gs_banked_up = F_N / (m*g); %gs through the seat
    %limit is 6g
gs_banked_lateral = (F_L)/(m*g); %gs through sidebar
    %limit is 3g
S_banked = pi*rBanked; %path length
vBanked(1, 1:51) = vBanked;
X = [xBanked; xBanked];
Y = [yBanked; yBanked];
Z = [zBanked; zBanked];
C = [vBanked; vBanked];
hold on
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet);
yBanked = (rBanked * cos(theta_circ))-49;
xBanked = (rBanked * sin(theta_circ))+max(xTurn_Trans);
zBanked = zeros(1,51)+min(zTurn_Trans);
vBanked = calcVelocity(125,zBanked(1,1)); %velocity
%% Loop
rLoop = 40;
thetaLoop = linspace(0, 2*pi, 100);
xLoop = (rLoop .* cos(thetaLoop)) + min(xBanked);
yLoop = min(yBanked) * ones(1,length(thetaLoop));
zLoop = (rLoop .* sin(thetaLoop)) + zBanked(51) + rLoop;
lengthLoop = 2*pi*rLoop;
% velocity for loop section
v_loop = calcVelocity(125,zLoop);

% Duplicate for 2-row surface
X = [xLoop; xLoop];
Y = [yLoop; yLoop];
Z = [zLoop; zLoop];
C = [v_loop; v_loop];
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet); 
% Reform Axes
xlim([-100 200])
ylim([-150 150])
zlim([-150 150])
%%  Arc after loop. Transition from loop to brake section
% Arc
xStrTrans = linspace(min(xBanked), min(xBanked - 50), 1000);
yStrTrans = min(yBanked) * ones(1, length(xStrTrans));
zStrTrans = min(zBanked) * ones(1, length(xStrTrans));
vStrTrans = sqrt(2 * g .* (125 - zStrTrans));
X = [xStrTrans; xStrTrans];
Y = [yStrTrans; yStrTrans];
Z = [zStrTrans; zStrTrans];
C = [vStrTrans; vStrTrans];
 
hold on
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet);
rTrans = ((125-min(zBanked))/3)+5;
thetaTrans = linspace(pi/2, 2*pi/3, 1000);
xTransLoop = (rTrans.*cos(thetaTrans))+min(xStrTrans);
yTransLoop = min(yBanked) * ones(1, length(xTransLoop));
zTransLoop = (rTrans.*sin(thetaTrans))-rTrans+min(zBanked);
vTransLoop = sqrt(2 * g * (125 - zTransLoop));
X = [xTransLoop; xTransLoop];
Y = [yTransLoop; yTransLoop];
Z = [zTransLoop; zTransLoop];
C = [vTransLoop; vTransLoop];
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet);
lengthTransLoop = rTrans*pi/6;
%G's along arc
 sTrans = linspace(0, rTrans*pi/6, 1000);
 GArcTrans = ((2*(125-min(zBanked)+rTrans*(1-sin((pi/2)-(sTrans/rTrans)))))/rTrans)-sin((pi/2)-(sTrans/rTrans));
%% Straight descent
xStrDec = linspace(min(xTransLoop), sqrt(3)*(50*(1-cos(pi/6))-((-1/sqrt(3))*min(xTransLoop)+min(zTransLoop))));
yStrDec = min(yBanked) * ones(1, length(xStrDec));
zStrDec = (1/sqrt(3))*xStrDec-(1/sqrt(3))*min(xTransLoop)+min(zTransLoop);
vStrDec = sqrt(2 * g .* (125 - zStrDec));
X = [xStrDec; xStrDec];
Y = [yStrDec; yStrDec];
Z = [zStrDec; zStrDec];
C = [vStrDec; vStrDec];
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet);
lengthStrDec = sqrt(((max(xStrDec)-min(xStrDec))^2+(max(zStrDec)-min(zStrDec))^2));
StrDec_deltaz = zStrDec(1,1) - zStrDec(1,100);
StrDec_deltax = xStrDec(1,1) - xStrDec(1,100);
theta_StrDec = atand(StrDec_deltaz / StrDec_deltax);
gs_StrDec_up = cosd(theta_StrDec);
gs_StrDec_back = sin(theta_StrDec);
%% Arc to z=0
rArcFinal = 50;
thetaArcFinal = linspace(-pi/3, -pi/2, 1000);
xArcFinal = rArcFinal.*cos(thetaArcFinal)+min(xStrDec)-50*sin(pi/6);
yArcFinal = min(yBanked)*ones(1,length(xArcFinal));
zArcFinal = rArcFinal.*sin(thetaArcFinal)+50;
vArcFinal = sqrt(2 * g .* (125 - zArcFinal));
X = [xArcFinal; xArcFinal];
Y = [yArcFinal; yArcFinal];
Z = [zArcFinal; zArcFinal];
C = [vArcFinal; vArcFinal];
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet);
lengthArcFinal = rArcFinal*pi/6;
theta_Arcto0 = atand(xArcFinal ./ zArcFinal);
v_Arcto0 = calcVelocity(125, zArcFinal);
lengthArcFinal_vec = linspace(0, lengthArcFinal, 1000);
for i=1:999
    atan_Arcto0(i) = ((v_Arcto0(1,i+1))^2 - (v_Arcto0(1,i))^2) / (2*lengthArcFinal_vec(i)); %calculating tangent acceleration
end
gs_Arcto0_up = cosd(theta_Arcto0) + (v_Arcto0.^2)/(rArcFinal*9.81);
gs_Arcto0_back = sind(theta_Arcto0(1:999)) + atan_Arcto0./9.81;
%% break section
a_brake=-9.81;
function v_brake = V_perSec(t_brake,a_brake)
    acceleration = a_brake;
    v_brake = sqrt(2*9.81*125) +acceleration*(t_brake);
end
s_brake = (2*9.81*125)/(2*(-a_brake));
t_brake = linspace(0,6,100); %can change the 5 to a 6 and then add in a cutoff to get v=0, and change 20 to 50
g_brake = a_brake/9.81;
v_brake = V_perSec(t_brake,a_brake);
for i=1:length(v_brake)
    if v_brake(i) <= 0
        v_brake(i) = 0;
    end
end
I = find(v_brake==0,1);
dist_brake = v_brake(1,1)*t_brake(1,I) + .5*a_brake* (t_brake(1,I))^2;
xBrake = (xArcFinal(1,1000) - linspace(0,dist_brake,100));
yBrake = min(yBanked)*ones(1,length(xBrake));
zBrake = zeros(1,100);
X = [xBrake; xBrake];
Y = [yBrake; yBrake];
Z = [zBrake; zBrake];
C = [v_brake; v_brake];
surf(X, Y, Z, C,'EdgeColor','interp','FaceColor','none','LineWidth',2);
colormap(jet);
c = colorbar();
c.Label.String = 'Speed (m/s)';
%gs for breaking section
gs_braking_back = linspace(1,1,100);
gs_braking_up = linspace(1,1,100);
dist_braking_vec = linspace(0,dist_brake,100);
%% Calculate total length
lengthTot = lengthDrop + lengthArc + lengthParabola+ length_Turn_Trans + S_banked + lengthLoop + lengthTransLoop + lengthStrDec + lengthArcFinal + dist_brake; 
%% Calculate and plot Loop G's
sLoop = linspace(0,2*pi*rLoop);
 % G's on the loop 
G_loop_forward = sin(thetaLoop); 
G_loop_back = -G_loop_forward; 
G_loop_up = (v_loop.^2) / (rLoop*g) + cos(thetaLoop); 
G_loop_down = -G_loop_up;
% fix entries 
G_loop_forward(G_loop_forward < 0) = 0; % G force exerted by the back of the seat
G_loop_back(G_loop_back < 0) = 0; % G force pushing rider back
G_loop_up(G_loop_up < 0) = 0; % G force exerted by the bottom of the seat
G_loop_down(G_loop_down < 0) = 0; % G force pushing the rider down
figure();
subplot(2,1,1)
plot(thetaLoop, G_loop_forward);
title("Arc Length vs Forward Loop G-Force");
xlabel("Arc Length");
ylabel("G-Force");
subplot(2,1,2)
plot(thetaLoop, G_loop_back);
title("Arc Length vs Backward Loop G-Force");
xlabel("Arc Length");
ylabel("G-Force");
figure();
subplot(2,1,1)
plot(thetaLoop, G_loop_up);
title("Arc Length vs Upward Loop G-Force");
xlabel("Arc Length");
ylabel("G-Force");
subplot(2,1,2)
plot(thetaLoop, G_loop_down);
title("Arc Length vs Downward Loop G-Force");
xlabel("Arc Length");
ylabel("G-Force");
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
%% plots for gs during braking section
figure();
subplot(2,1,1);
plot(dist_braking_vec,gs_braking_up);
ylabel('gs up');
xlabel('pathlength');
title('gs up through seat for breaking section');
subplot(2,1,2);
plot(dist_braking_vec,gs_braking_back);
ylabel('gs back');
xlabel('pathlength');
title('gs back through restraint for breaking section');
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
l_b = lengthDrop;
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
lengthStrDec_vec = linspace(0,lengthStrDec,100);
gs_StrDec_up_vec = linspace(gs_StrDec_up, gs_StrDec_up,100);
gs_StrDec_back_vec = linspace(gs_StrDec_back,gs_StrDec_back,100);
figure();
%gs up through riders seat
hold on
plot(length_blue, gs_up_blue_vec,'b');
plot(length_red + l_b, gs_up_red,'b');
plot(l_b + l_r + pathlength_para,gs_up_para,'b');
plot(l_b + l_r + lengthParabola + length_Turn_Trans_vec, gs_up_Arc_Trans,'b');
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked, gs_banked_up_vec,'b');
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked + sLoop(1,100), G_loop_up,'b');
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked(1,100) + sLoop(1,100) + sTrans,GArcTrans,'b');
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked(1,100) + sLoop(1,100) + sTrans(1,1000) + lengthStrDec_vec, gs_StrDec_up_vec,'b');
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked(1,100) + sLoop(1,100) + sTrans(1,1000) + lengthStrDec + lengthArcFinal_vec, gs_Arcto0_up,'b');
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked(1,100) + sLoop(1,100) + sTrans(1,1000) + lengthStrDec + lengthArcFinal + dist_braking_vec,gs_braking_up,'b');
%xlines to devide each section of the track
xline([l_b l_b+l_r l_b+l_r+lengthParabola l_b+l_r+lengthParabola+length_Turn_Trans l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100) l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)],'b--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal+dist_brake],'--');
yline(6,'--','G Limit');
ylabel('G-Force (N)');
xlabel('Pathlength (m)');
title('G-Force Pushing Up On Rider By Seat');
%%
figure();
%gs from seat restraint(blue and red section)
hold on
plot(length_blue, abs(gs_forward_blue_vec));
plot(length_red(1:499) + l_b, gs_forward_red);
plot([(.5*l_r)+l_b ; l_r+l_b], [0;0]);
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked + sLoop(1, 100), G_loop_back,'b');
%xlines to devide each section of the track
xline([l_b l_b+l_r l_b+l_r+lengthParabola l_b+l_r+lengthParabola+length_Turn_Trans l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100) l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)],'b--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal+dist_brake],'--');
yline(4,'--','G Limit');
ylabel('G-Force (N)');
xlabel('Path-Length (m)');
title('G-Force Pushing Back On Rider By Lap Bar');
%%
% gs from seat restraint pushing down on rider
figure();
hold on; 
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked + sLoop(1, 100), G_loop_down,'b');
%xlines to devide each section of the track
xline([l_b l_b+l_r l_b+l_r+lengthParabola l_b+l_r+lengthParabola+length_Turn_Trans l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100) l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)],'b--');
xline(l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000),'--');
xline(l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec,'--');
xline(l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal,'--');
xline(l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal+dist_brake,'--');
yline(1,'--','G Limit');
ylabel('G-Force (N)');
xlabel('Pathlength (m)');
title('G-Force Pushing Down On Rider By Lap Bar');
hold off;
%% ???
%We are still missing some of the gs for this part, also I might be
%confusing the gs from the back of the seat with the gs from the seat
%restraint. Even so, the shape of a lot of the g plots dont make sense. We
%probably need to fix/change this.
figure();
%gs from back of seat pushing
hold on
plot([0;(l_b+(l_r*.5))], [0;0]);
plot(l_b + length_red(501:999), gs_back_red);
%plot(l_b + l_r + pathlength_para,,'b');
plot(l_b + l_r + lengthParabola + length_Turn_Trans_vec(3:1000), gs_back__Arc_Trans(2:999),'b');
%plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked, ,'b');
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked(1,100) + sLoop, G_loop_forward ,'b');
%plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked(1,100) + sLoop(1,100) + sTrans,,'b');
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked(1,100) + sLoop(1,100) + sTrans(1,1000) + lengthStrDec_vec, gs_StrDec_back_vec);
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked(1,100) + sLoop(1,100) + sTrans(1,1000) + lengthStrDec + lengthArcFinal_vec(2:1000), gs_Arcto0_back);
plot(l_b + l_r + lengthParabola + length_Turn_Trans + pathlength_banked(1,100) + sLoop(1,100) + sTrans(1,1000) + lengthStrDec + lengthArcFinal + dist_braking_vec, gs_braking_back);
%xlines to devide each section of the track
xline([l_b l_b+l_r l_b+l_r+lengthParabola l_b+l_r+lengthParabola+length_Turn_Trans l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100) l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)],'b--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal+dist_brake],'--');
yline(4,'--','G Limit');
ylabel('G-Force (N)');
xlabel('Pathlength (m)');
title('G-Force Pushing Forward On Rider By Seat');
%% plot for lateral gs for whole coaster
length_1 = l_b+l_r+lengthParabola+length_Turn_Trans; %pathlength up until banked turn
length_1_vec = linspace(0,length_1,100);
length_2 = sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal+dist_brake; %pathlength after banked turn
length_2_vec = linspace(0,length_2,100);
gs_lat_0 = zeros(1,100);
figure();
hold on
plot(length_1_vec, gs_lat_0);
plot(length_1+pathlength_banked, gs_banked_lateral_vec);
plot(length_1 + S_banked + length_2_vec, gs_lat_0);
xline([l_b l_b+l_r l_b+l_r+lengthParabola l_b+l_r+lengthParabola+length_Turn_Trans l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100) l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)],'b--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal],'--');
xline([l_b+l_r+lengthParabola+length_Turn_Trans+pathlength_banked(1,100)+sLoop(1,100)+sTrans(1,1000)+lengthStrDec+lengthArcFinal+dist_brake],'--');
yline(3,'--','G Limit');
xline(length_1,'--');
xline(length_1+S_banked,'--');
yline(3,'--','g limit');
ylabel('G-Force (N)');
xlabel('Pathlength (m)');
title('gs lateral');
%% Velocity function
function v = calcVelocity(h0, h)
v = sqrt(2*9.81*(h0-h));
end