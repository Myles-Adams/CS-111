clear;
clc;
close all;

% Domain
xL = -2;
xR = 2;
yB = -2.5;
yT = 2.5;

% Diffusion coefficient
lambda = 0.0015;

% boundary conditions
T_bc = @(t,x,y) min(20 + 80*t/60, 100);

% Source term
f = @(t,x, y) 0;

% Initial conditions
T_start = @(x,y) 20;

% Time
t_start = 0;
t_final = 1500;
dt = 5;
t_all = zeros(1, (t_final - t_start)/dt);
T_center_all = zeros(1, (t_final - t_start)/dt);
t_count = 1;

subplot_num = 1;

% Space discretization
Nx = 80;
Ny = 100;

% Discretize
x = linspace(xL, xR, Nx);
y = linspace(yB, yT, Ny);
dx = (xR-xL)/(Nx-1);
dy = (yT-yB)/(Ny-1);

% Pre-allocate arrays
T_old = zeros(Ny,Nx);
T_new = zeros(Ny,Nx);
T_exa = zeros(Ny,Nx);

% Populate T_old with start values
for i = 1:Nx
    for j =1:Ny
        T_old(j,i) = T_start(x(i), y(j));
    end
end

t = t_start;

% Create sparse matrix and allocate memory for right-hand side
RHS = zeros(Nx*Ny,1);

% Calculate A matrix before the while-loop to save time
A = make_A(0,Nx,Ny,dx,dy,dt,lambda);

% Create figure
figure('rend','painters','pos',[400 100 800 1000])

% Loop through all values of t
while t < t_final
    
    % Plot potato temp at different times
    if (t == 0) || (t == 200) || (t == 400) || (t == 600)
        subplot(4,2,subplot_num)
        contourf(x,y,T_old, 100, 'LineColor','none');
        title(['      Temperature of potato at t = ' num2str(t)]);
        colorbar;
        caxis([20,100]);
        hold on
        subplot_num = subplot_num + 1;
    end
    
    % Check for overshoot
    if t + dt > t_final
        dt = t_final-t;
        
        % Recalculate A
        A = make_A(0,Nx,Ny,dx,dy,dt,lambda);
    end
    
    % Calculate RHS array
    for i = 1:Nx
        for j = 1:Ny
            p = (j-1)*Nx+i;
            if (i==1) || (i==Nx) || (j==1) || (j==Ny)
                RHS(p) = T_bc(t+dt,x(i),y(j));
            else
                RHS(p) = T_old(j,i) + dt*f(t+dt,x(i),y(j));
            end
        end
    end
    
    % Solve system of equations
    T_new = reshape(A\RHS,Nx,Ny)';
    
    % Store current time and center temperature
    t_all(t_count) = t;
    T_center_all(t_count) = T_old(Ny/2, Nx/2);
    
    % Move new to old
    T_old = T_new;
    
    % Iterate time
    t = t+dt;
    t_count = t_count + 1;
    
end

% Plot temp vs time
subplot(4,2,[5,8]);
plot(t_all, T_center_all, '-')
title('Temperature vs Time at Center of Potato');
text(890,65,'\leftarrow Reaches 65^{\circ}C at t=890', 'fontsize', 14)
xlabel('Time (s)');
ylabel('Temperature (^{\circ}C)');