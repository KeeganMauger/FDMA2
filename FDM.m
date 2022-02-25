%% Elec 4700 Assignment 2
%%% Finite Difference Method for Solving Laplace
% Keegan Mauger
% 101042551

%% Section 1: Electrostatic Potential in a Rectangular Region
%
% Using Laplace's equation by the finite difference method, an
% electrostatic potential problem was to be solved. The problem was
% modelled by an orthogonal resistor network in a region of W by L, choosen
% to be 90 by 60 for this problem. By using a mesh of resistors, boundary
% conditions and intrusions became easier to model.
%
% Using the matrix form GV=F, the electrostatic potential in the
% rectangular region was solved using $del^2*V = 0$. 
%
%%% Section 1a: $V = V_0 @ x = 0$ and $V = 0 @ x = L$
% For this problem, the boundary conditions were set such that $V = V_0 = 1$
% at the left region boundary, and $V = 0$ at the right region boundary.
% The y-axis boundary conditions were not set. This case was then solved,
% with the code and results shown below. Note that the region dimensions
% are taken to be unitless.

clear all;
close all;
clc;
%set(0,'DefaultFigureWindowStyle','docked');

% Solving V=V0 @ x=0 and V=0 @ x=L in region LxW
% Implement funtion 'pbaspect' to fix Z aspect ratio

L = 90;
W = 2/3 * L;
V0 = 1;

fMesh = 1;                  % Mesh factor
nx = fMesh*L;
ny = fMesh*W;
G = sparse(nx*ny);
%V = sparse(nx,ny);
F = sparse(1,nx*ny);

La = linspace(0,L,nx);
Wa = linspace(0,W,ny);

for i = 1:nx                %Iteration through length
    for j = 1:ny            %Iteration through width
        n = j + (i-1)*ny;

        if i == 1          % x=0 BCs
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 1;
        elseif i == nx     % x=1 BCs
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 0;       %F(n)=0 sets z at final width to 0
              
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            G(n,n) = -(4);
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end
            
    end
end
% figure(1)
% spy(G)

V = G\F';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        Vmap(i,j) = V(n);
    end
end

figure(1)
surf(Vmap)
pbaspect([1 1 0.5])
view(90,270)
title('2D plot of Electrostatic Potential over a Rectangular Region')
xlabel('Region Width')    %Axis swapped due to view change, named correctly
ylabel('Region Length')
zlabel('Voltage (V)')
saveas(gcf,'Figure1')

%%% Results
% On the plot, the colourmap represents the voltage of the region. As
% expected, the voltage begins in the region at 1V and linearly decreses
% over the length to zero volts.