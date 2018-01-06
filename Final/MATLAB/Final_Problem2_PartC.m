clear;
clc;
close all;

% domain
xL = 0;
xR = 12;
yB = 0;
yT = 3;

% diffusion coefficient
D = 0.2;

% Velocities
Vx = -.8;
Vy = -.4;

% boundary conditions
c_bc = @(t,x,y) 0;
g = @(t,x,y) 0;

% initial conditions
c_start = @(x,y) 0;

% time interval
t_start = 0;
t_final = 10;
dt = .1;
t_all = zeros(1, (t_final - t_start)/dt);
t_count = 1;

% exact solution
c_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% space discretization
Nx = 160;
Ny = 40;

num_splits = 3;
x0 = 4;
y0 = 0;

% time-step

subplot_num = 1;

figure('rend','painters','pos',[100 100 1400 800])

x = linspace(xL, xR, Nx);
y = linspace(yB, yT, Ny);
dx = (xR-xL)/(Nx-1);
dy = (yT-yB)/(Ny-1);

beach_1 = zeros(1, (t_final - t_start)/dt);
beach_2 = zeros(1, (t_final - t_start)/dt);
beach_3 = zeros(1, (t_final - t_start)/dt);

c_old = zeros(Ny,Nx);
c_new = zeros(Ny,Nx);
c_exa = zeros(Ny,Nx);

for i = 1:Nx
    for j =1:Ny
        c_old(j,i) = c_start(x(i), y(j));
    end
end

t = t_start;

% Create sparse matrix and allocate memory for right-hand side
RHS = zeros(Nx*Ny,1);

% Calculate the matrix before the while-loop to save time
% internal points
A = make_A(1,Nx,Ny,dx,dy,dt,D,Vx,Vy);


while t < t_final
    
    if (round(t,1) == 1) || (round(t,1) == 4) || (round(t,1) == 7)
        subplot(3,3,subplot_num)
        contourf(x,y,c_old, 100, 'LineColor','none');
        title(['Water at t = ' num2str(t)]);
        colorbar;
        caxis([0,.03]);
        hold on
        subplot_num = subplot_num + 1;
    end
    
    if t + dt > t_final
        dt = t_final-t;
        
        A = make_A(1,Nx,Ny,dx,dy,dt,D,Vx,Vy);
    end
    
    % internal points
    for i = 1:Nx
        for j = 1:Ny
            p = (j-1)*Nx+i;
            if (i==1) || (i==Nx) || (j==Ny)
                RHS(p) = c_bc(t+dt,x(i),y(j));
            elseif (j==1)
                RHS(p) = c_old(j,i) + dt*f(t+dt,x(i),y(j)) - Vx*dt*((c_old(j,i+1) - c_old(j,i))/(dx)) - Vy*dt*((c_old(j+1,i) - c_old(j,i))/(dy)) - (2*dt/dy)*g(t+dt,x(i),y(j));
            else
                RHS(p) = c_old(j,i) + dt*f(t+dt,x(i),y(j)) - Vx*dt*((c_old(j,i+1) - c_old(j,i))/(dx)) - Vy*dt*((c_old(j+1,i) - c_old(j,i))/(dy));
            end
        end
    end
    
    % solve system of equations
    c_new = reshape(A\RHS,Nx,Ny)';
    
    beach_1(t_count) = c_old(1,54);
    beach_2(t_count) = c_old(1,80);
    beach_3(t_count) = c_old(1,107);
    
    
    c_old = c_new;
    t_all(t_count) = t;
    t = t+dt;
    t_count = t_count + 1;
    
end

t_all(t_count) = t;
beach_1(t_count) = c_old(1,54);
beach_2(t_count) = c_old(1,80);
beach_3(t_count) = c_old(1,107);

subplot(3,3,[4,9]);
plot(t_all, beach_1, '-', t_all, beach_2, '-', t_all, beach_3, '-')
title('Concentration vs Time at Beaches');
xlabel('Time');
ylabel('Concentration');
legend('Beach at (4,0)', 'Beach at (6,0)', 'Beach at (8,0)', 'Location', 'southeast');

function [val] = f(t,x,y)

x_s = 10;
r_s = .1;
epsilon = .1;

if (t > .5)
    val = 0;
else
    val = .5*(1-tanh((sqrt((x-x_s)^2 + y^2) - r_s)/epsilon));
end


end