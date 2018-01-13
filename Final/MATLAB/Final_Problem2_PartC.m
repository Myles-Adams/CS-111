clear;
clc;
close all;

% domain
xL = 0;
xR = 12;
yB = 0;
yT = 3;

% Diffusion coefficient
D = 0.2;

% Velocities
Vx = -.8;
Vy = -.4;

% Boundary conditions
c_bc = @(t,x,y) 0;
g = @(t,x,y) 0;

% Initial conditions
c_start = @(x,y) 0;

% Time
t_start = 0;
t_final = 10;
dt = .1;
t_all = zeros(1, (t_final - t_start)/dt);
t_count = 1;

% Exact solution
c_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% Space discretization
Nx = 160;
Ny = 40;

% Number of discretizations 
num_splits = 3;

subplot_num = 1;

% Create Figure
figure('rend','painters','pos',[100 100 1400 800])

% Discretize area
x = linspace(xL, xR, Nx);
y = linspace(yB, yT, Ny);
dx = (xR-xL)/(Nx-1);
dy = (yT-yB)/(Ny-1);

% Preallocate arrays for beach values
beach_1 = zeros(1, (t_final - t_start)/dt);
beach_2 = zeros(1, (t_final - t_start)/dt);
beach_3 = zeros(1, (t_final - t_start)/dt);

% preallocate arrays for concentration values
c_old = zeros(Ny,Nx);
c_new = zeros(Ny,Nx);
c_exa = zeros(Ny,Nx);

% Populate c_old with start values
for i = 1:Nx
    for j =1:Ny
        c_old(j,i) = c_start(x(i), y(j));
    end
end

t = t_start;

% Create sparse matrix and allocate memory for right-hand side
RHS = zeros(Nx*Ny,1);

% Calculate A matrix before the while-loop to save time
A = make_A(1,Nx,Ny,dx,dy,dt,D,Vx,Vy);

% Iterate through all values of t
while t < t_final
    
    % Plot concentrations at certain times
    if (round(t,1) == 1) || (round(t,1) == 4) || (round(t,1) == 7)
        subplot(3,3,subplot_num)
        contourf(x,y,c_old, 100, 'LineColor','none');
        title(['Water at t = ' num2str(t)]);
        colorbar;
        caxis([0,.03]);
        hold on
        subplot_num = subplot_num + 1;
    end
    
    % Check for overshoot
    if t + dt > t_final
        dt = t_final-t;
        
        % Recalculate A
        A = make_A(1,Nx,Ny,dx,dy,dt,D,Vx,Vy);
    end
    
    % Calculate RHS Array
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
    
    % Solve system of equations
    c_new = reshape(A\RHS,Nx,Ny)';
    
    % Store concentration values at beaches
    beach_1(t_count) = c_old(1,54);
    beach_2(t_count) = c_old(1,80);
    beach_3(t_count) = c_old(1,107);
    
    % Store values and prepare for next iteration
    c_old = c_new;
    t_all(t_count) = t;
    t = t+dt;
    t_count = t_count + 1;
    
end

% Plot beach concentrations over time
subplot(3,3,[4,9]);
plot(t_all, beach_1, '-', t_all, beach_2, '-', t_all, beach_3, '-')
title('Concentration vs Time at Beaches');
text(1.3,.006273,'\leftarrow Reaches .006273', 'fontsize', 14)
text(1.52,.0055,'at t = 1.3', 'fontsize', 14)
text(3.3,.006549,'\leftarrow Reaches .006549', 'fontsize', 14)
text(3.52,.0058,'at t = 3.3', 'fontsize', 14)
text(5.2,.00597,'\leftarrow Reaches .00597', 'fontsize', 14)
text(5.42,.0052,'at t = 5.2', 'fontsize', 14)
xlabel('Time (days)');
ylabel('Concentration');
legend('Beach at (4,0)', 'Beach at (6,0)', 'Beach at (8,0)', 'Location', 'northeast');

% Function for source term ??
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