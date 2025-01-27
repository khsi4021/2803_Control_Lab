clear
clc
close all

%banked turn
r = 10;
g = 9.81;
v = 20;

gs_banked = v^2 /(g*r);

S_banked = pi*r;

theta = 0:pi/50:pi;
x = r * cos(theta);
y = r * sin(theta);
z = zeros(1,51);

plot3(x,y,z);