clear;
clc;
close all;
% Dynamics and Control Lab
% Nathan Stahl

% Code for loop
% V = sqrt(v_0^2 + 2g(h_0 - h))

% Initial Conditions
r = 10; % Radius of loop (m)
v_0 = 20; % Velocity of roller coaster going into the loop (m/s)
g = 9.81; % Acceleration due to gravity (m/s^2)
minv = sqrt(2 * g * 2 * r); % Minimum velocity to complete loop with current radius

% length of track 
d_loop = 2 * pi * r; % length of track in meters

% graph for loop
theta = linspace(0, 2*pi, 200);
x_loop = r .* sin(theta); % x distance
y_loop = zeros(1, 200); % y change = 0
z_loop = r .* sin(theta - (pi/2)) + r; % change in height over time

% velocity for loop section
vel = sqrt(v_0^2 - 2 * g .* z_loop);

figure();
plot(theta, vel);
title("Distance vs Velocity");
xlabel("X-Distance");
ylabel("Velocity");



% G's on the track 
G_loop_back = sin(theta); % G force exerted by the back of the seat
G_loop_seat = (vel.^2) / (r*g) + cos(theta); % G force exerted by the bottom of the seat

figure();
plot(theta, G_loop_back);
title("Distance vs Back of Seat G's");
xlabel("X-Distance");
ylabel("G-Force");

figure();
plot(theta, G_loop_seat);
title("Distance vs Bottom of Seat G's");
xlabel("X-Distance");
ylabel("G-Force");



%plot
figure();
plot3(x_loop, y_loop, z_loop);

%% 7) 3D Plot with Speed-Colored Line
% "Surface Trick": replicate each array in 2 rows so that 'surf' draws 
% a thin strip whose "edges" are color-interpolated by velocity.

% Convert to row vectors (in case they're column vectors):
X = x_loop(:).';
Y = y_loop(:).';
Z = z_loop(:).';
C = vel(:).';   % color array (speed) must match dimension

% Duplicate for 2-row surface
X2 = [X; X];
Y2 = [Y; Y];
Z2 = [Z; Z];
C2 = [C; C];

figure('Name','3D Loop, Colored by Speed');
surf(X2, Y2, Z2, C2,'EdgeColor','interp', ...   % interpolate color along edges
    'FaceColor','none',  ...    % no colored faces, just edges
    'LineWidth',2);
axis equal; grid on;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Vertical Loop in 3D (Line Colored by Speed)');
colormap(jet); 
clim([min(C) max(C)]);  % color scale matches velocity range
colorbar('Ticks',linspace(min(C),max(C),5), ...
         'Label','Speed (m/s)');
view(45,20);

