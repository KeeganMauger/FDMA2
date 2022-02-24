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

nx = L;
ny = W;

V = zeros(nx,ny);
n = 1;
f = 0;



for i=1:nx
    for j = 1:ny
        for k=1:10
            f = f + ((1/n)*(cosh(n*pi*i/nx)/cosh(n*pi*ny/nx))*sin(n*pi*j/nx));
            n = n+2;
        end
        V(i,j) = f * (4*V0/pi);
    end
end

surf(V)