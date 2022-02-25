% Analytical series solution

% Initialization

clear all;
close all;
clc;
set(0,'DefaultFigureWindowStyle','docked');

% Solving V=V0 @ x=0 and V=0 @ x=L in region LxW
% Implement funtion 'pbaspect' to fix Z aspect ratio

L = 90;
W = 2/3 * L;
V0 = 1;

nx = 100;
ny = 100;

La = linspace(-L/2,L/2,nx);
Wa = linspace(0, W, ny);

V = zeros(nx,ny);
Vs = zeros(nx,ny);
n = 1;
f = 0;

for n=1:2:210
    for i=1:nx
        for j = 1:ny
            V(i,j) = (4*V0/pi)*((1/n)*(cosh(n*pi*La(i)/(W))/cosh(n*pi*(L/2)/(W)))*sin(n*pi*Wa(j)/(W)));              
        end
    end
    Vs = Vs + V;
    surf(Vs)
    pbaspect([1 1 0.5])
    pause(0.01)
end

