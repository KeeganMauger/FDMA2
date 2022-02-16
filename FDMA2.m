% ELEC4700 Assignment 2: Finite Difference Method

% Initialization

clear all;
close all;
clc;
set(0,'DefaultFigureWindowStyle','docked');

% Solving V=V0 @ x=0 and V=0 @ x=L in region LxW

L = 90;
W = 2/3 * L;
V0 = 1;

nx = L;
ny = W;
G = sparse(nx*ny);
V = sparse(nx,ny);
B = sparse(1,nx*ny);

for i = 1:nx                %Iteration through length
    for j = 1:ny            %Iteration through width
        n = j + (i-1)*ny;

        if i == 1;
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 1;
        elseif i == nx
            G(n,:) = 0;
            G(n,n) = 1;
%         elseif j == 1
%             nxm = j + (i-2)*ny;
%             nxp = j + (i)*ny;
%             nyp = j+1 + (i-1)*ny;
%             
%             G(n,n) = 1;
%         elseif j == ny
%             nxm = j + (i-2)*ny;
%             nxp = j + (i)*ny;
%             nym = j-1 + (i-1)*ny;
%             
%             G(n,n) = 1;
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
figure(1)
spy(G)