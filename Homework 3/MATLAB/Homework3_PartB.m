% Homework 3 Part b
clear;
clc;
close all;

% Boundaries 
xL = -1;
xR = 1;

yB = -1;
yT = 1;

% Constants
r = .2;
r_c = .2;
eps = .05;
x0 = -.7;
y0 = -.2;
alpha = .2;
U = .1;

% Initial condition function
c0 = @(x,y) .5*(1 - tanh((sqrt((x-x0)^2 + (y-y0)^2) - r)/eps));

% Exact solution function
c_exact = @(t,x,y) .5*(1-tanh((sqrt((x-R*cos(t))^2 + (y-R*sin(t))^2) - r)/eps));

% Boundary condition function
c_bc = @(t,x,y) 0;

% Times
t_start = 0;
t_final = 15;

% space discretization
Nx = 100;
Ny = 100;

% Courant number
C = 0.75;
    
% Create x and y arrays
x = linspace(xL, xR, Nx);
y = linspace(yB, yT, Ny);

% Define dx and dy
dx = x(2)-x(1);
dy = y(2)-y(1);

% High number to start as minimum value
dt_min = 100000;

% Find smallest dt that satisfies equation for various dx/dy and vx/vy
for i=1:Nx
    for j=1:Ny
        [vx,vy] = v_xy(0,x(i),y(j));
        dt = C/(vx/dx + vy/dy);
        if (dt > 0) && (dt < dt_min)
            dt_min = dt;
        end
    end
end

dt = dt_min;

% Impose initial conditions
for i=1:Nx
    for j=1:Ny
        c_old(i,j) = c0(x(i),y(j));
    end
end
c_new = c_old;

t = t_start;

% Loop through all values of t from t_start to t_final by dt
while t < t_final
    
    % Check for overshoot
    if t+dt > t_final
        dt = t_final - t;
    end
    
    % Upwind scheme
    for i = 1:Nx
        for j = 1:Ny
            
            % Calculate velocities
            [vel_x, vel_y] = v_xy(t,x(i),y(j));
            
            % Check if boundary condition is needed based on vel_x, if not
            % calculate dc/dx
            if vel_x >=0
                if i > 1
                    dcdx = (c_old(i,j)-c_old(i-1,j))/dx;
                else
                    c_new(i,j) = c_bc(t,x(i),y(j));
                    continue
                end
            else
                if i < Nx
                    dcdx = (c_old(i+1,j)-c_old(i,j))/dx;
                else
                    c_new(i,j) = c_bc(t,x(i),y(j));
                    continue
                end
            end
            
            % Check if boundary condition is needed based on vel_y, if not
            % calculate dc/dy
            if vel_y >= 0
                if j > 1
                    dcdy = (c_old(i,j)-c_old(i,j-1))/dy;
                else
                    c_new(i,j) = c_bc(t,x(i),y(j));
                    continue
                end
            else
                if j < Ny
                    dcdy = (c_old(i,j+1)-c_old(i,j))/dy;
                else
                    c_new(i,j) = c_bc(t,x(i),y(j));
                    continue
                end
            end
            
            % If boundary condition is not needed, calculate c_new
            c_new(i,j) = c_old(i,j) - dt*vel_x*dcdx - dt*vel_y*dcdy;
        end
    end
    
    % Step to next time t
    t = t+dt;
    c_old = c_new;
    
end

% Graph at t_final
contourf(x,y,c_old','LevelList', linspace(0.0,1.0,50), 'LineColor', 'none');
colorbar;
hold on

nx_quiver = 20; % coarse grid for plotting velocity field
ny_quiver = 20;

x_quiver = linspace(xL, xR, nx_quiver);
y_quiver = linspace(yB, yT, ny_quiver);

for i = 1:nx_quiver
    for j = 1:ny_quiver
        [vx_quiver(i,j),vy_quiver(i,j)] = v_xy(0,x_quiver(i),y_quiver(j));
    end
end
quiver(x_quiver, y_quiver, vx_quiver', vy_quiver','Color','w');
% -----------------------------------------------------------------------
draw_disk(0,0,.2);
hold off
axis([xL xR yB yT]);
axis equal;