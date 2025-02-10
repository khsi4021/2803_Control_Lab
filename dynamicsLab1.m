clear;
clc;
close all;
% Dynamics and Control Lab
% Nathan Stahl

% Code for loop
% V = sqrt(v_0^2 + 2g(h_0 - h))

h_0 = 125;
h_bottom = 100;
str_angle = 10;
r = 10;
v_0 = 0;



% velocity for loop section
height = [(h_0 - h_bottom):0.1:(h_0 - h_bottom + r), (h_0 - h_bottom + r):0.1:(h_0-h_bottom)];
v_a = sqrt(v_0^2 + 2 * 9.81 * (h_0 - h_bottom));
v_b = sqrt(v_0^2 + 2 * 9.81 * (h_0 - (h_bottom + r)));
v_c = sqrt(v_0^2 + 2 * 9.81 * (h_0 - (h_bottom + 2 * r)));
v_d = sqrt(v_0^2 + 2 * 9.81 * (h_0 - (h_bottom + r)));

% Velocity for straight section
distance = 0:0.1:((h_0 - h_bottom)/sind(str_angle));
v_straight = sqrt(v_0^2 / 9.81 .* distance); 
g_str_tot = cosd(str_angle);
g_str_back = g_str_tot * cosd(str_angle);
g_str_up = g_str_tot * sind(str_angle);


% length of track 
d_str = (h_0 - h_bottom)/sind(str_angle);
d_loop = 2 * pi * r;

d_tot = d_str + d_loop;

% G's on the track 
G_a_back = (v_a^2 / (9.81 * r)) + 1;  
G_b_back = (v_b^2 / (9.81 * r));
G_c_back = (v_c^2 / (9.81 * r)) - 1;
G_d_back = (v_d^2 / (9.81 * r));


%graph for straight section
x_str = linspace(0, (h_0 - h_bottom) / tand(str_angle), 100);
y_str = zeros(1, 100);
z_str = linspace(0, (h_0-h_bottom), 100);

%graph for loop
theta_cir = linspace(0, 2*pi, 200);
x_loop = r .* cos(theta_cir);
y_loop = zeros(1, 200);
z_loop =  r + r .* sin(theta_cir);

%plot
plot3(x_str, y_str, z_str, x_loop, y_loop, z_loop);





% 
% figure();
% hold on;
% r = 10;
% xc = 0;
% yc = r;
% 
% theta = linspace(0,2*pi);
% x = r*cos(theta) + xc;
% y = r*sin(theta) + yc;
% plot(x,y)
% 
% 
% axis equal
