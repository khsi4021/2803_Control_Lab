clear
clc
close all

%banked turn
r = 50; %radius of turn
g = 9.81;
v = 50; %velocity
theta = 70; %degrees from upright of roller cart

%placeholder values, delete later
m = 10;

F_r = (m*v.^2)/r; %force in radial direction
F_g = m*g; %force of gravity
F_N = F_r * sind(theta) - F_g * cosd(theta); %force normal to the seat of the cart
F_L = F_r * cosd(theta) + F_g * sind(theta); %force perpendicular to the normal

gs_banked_up = F_N / (m*g); %gs through the seat
    %limit is 6g
gs_banked_lateral = (F_L)/(m*g); %gs through sidebar
    %limit is 3g

S_banked = pi*r; %path length

theta_circ = 0:pi/50:pi; %in radians
x = r * cos(theta_circ);
y = r * sin(theta_circ);
z = zeros(1,51);

plot3(x,y,z);
