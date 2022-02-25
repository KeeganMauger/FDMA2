% ELEC4700 Assignment 2: Finite Difference Method

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

fMesh = 0.5:0.5:3;
for k = 1:length(fMesh)
    
    nx = round(fMesh(k)*L);
    ny = round(fMesh(k)*W);
    G = sparse(nx*ny);
    %V = sparse(nx,ny);
    F = sparse(1,nx*ny);
    
    
    Acond = 1;              % background conductivity of region, low resistance
    Bcond = 1e-2;           % Conductivity of boxes, highly resistive
    BN = 0:1:10;                % Changing bottleneck
    cMap = zeros(nx,ny);
    Lb = 20;
    Wb = 20;

    for u = 1:nx
        for v = 1:ny
            if (u >= 35 && u <= 55)
                if v >= 0 && v <= 20
                    cMap(u,v) = Bcond;
                elseif v >= 40 && v <= 60
                    cMap(u,v) = Bcond;
                else
                    cMap(u,v) = Acond;
                end
            else
                cMap(u,v) = Acond;
            end
        end
    end





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

                % COMMENT BELOW FOR 1a
                %F(n) = 1;       %F(n)=1 sets z at final width to 1

            elseif j == 1                 % y=0 BCs
                nxm = j + (i-2)*ny;
                nxp = j + (i)*ny;
                nyp = j+1 + (i-1)*ny;

                rxm = (cMap(i,j) + cMap(i-1,j))/2;
                rxp = (cMap(i,j) + cMap(i+1,j))/2;
                ryp = (cMap(i,j) + cMap(i,j+1))/2;

                G(n,n) = -(rxm+rxp+ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nyp) = ryp;

            elseif j == ny                % y=1 BCs
                nxm = j + (i-2)*ny;
                nxp = j + (i)*ny;
                nym = j-1 + (i-1)*ny;

                rxm = (cMap(i,j) + cMap(i-1,j))/2;
                rxp = (cMap(i,j) + cMap(i+1,j))/2;
                rym = (cMap(i,j) + cMap(i,j-1))/2;

                G(n,n) = -(rxm+rxp+rym);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;


                % COMMENT ABOVE FOR 1a

            else
                nxm = j + (i-2)*ny;
                nxp = j + (i)*ny;
                nym = j-1 + (i-1)*ny;
                nyp = j+1 + (i-1)*ny;

                rxm = (cMap(i,j) + cMap(i-1,j))/2;
                rxp = (cMap(i,j) + cMap(i+1,j))/2;
                rym = (cMap(i,j) + cMap(i,j-1))/2;
                ryp = (cMap(i,j) + cMap(i,j+1))/2;

                G(n,n) = -(rxm+rxp+rym+ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;
                G(n,nyp) = ryp;
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


    for i = 1:nx
        for j = 1:ny
            if i == 1
                Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
            elseif i == nx
                Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
            else
                Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
            end
            if j == 1
                Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
            elseif j == ny
                Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
            else
                Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
            end
        end
    end

    Ex = -Ex;
    Ey = -Ey;

    eFlowx = cMap .* Ex;        %Jx
    eFlowy = cMap .* Ey;        %Jy


    C0 = sum(eFlowx(1, :));
    Cnx = sum(eFlowx(nx, :));
    Curr(k) = (C0 + Cnx) * 0.5;
end

plot(fMesh,Curr)
xlabel('Mesh factor multiplier, where default mesh is 90x60')
ylabel('Current')