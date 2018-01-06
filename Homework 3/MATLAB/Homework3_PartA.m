% two dimensional advection
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
R = .5;
eps = .1;

% Initial condition function
c0 = @(x,y) .5*(1 - tanh((sqrt((x-R)^2 + y^2) - r)/eps));

% Velocity functions
v_x = @(t,x,y) -y;
v_y = @(t,x,y) x;

% Exact solution function
c_exact = @(t,x,y) .5*(1-tanh((sqrt((x-R*cos(t))^2 + (y-R*sin(t))^2) - r)/eps));

% Boundary condition function
c_bc = @(t,x,y) 0;

% Times
t_start = 0;
t_final = .2;

% Space discretization
Nx = 50;
Ny = 50;

% Courant number
C = 0.75;

% How many different space discretizations we will be calculating for
num_splits = 4;

for k=1:num_splits
    
    % Create array for x and y
    x = linspace(xL, xR, Nx);
    y = linspace(yB, yT, Ny);
    
    % Define dx and dy
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    
    % Calculate dt for maximum vx and vy
    vx_max = 1;
    vy_max = 1;
    dt = C/((vx_max/dx) + (vy_max/dy));
    
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
                vel_x = v_x(t,x(i),y(j));
                vel_y = v_y(t,x(i),y(j));
                
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
    
    % Calculate exact solution
    for i=1:Nx
        for j=1:Ny
            c_exa(i,j) = c_exact(t,x(i),y(j));
        end
    end
    
    % Calculate max error
    error(k) = max(max(abs(c_exa - c_new)));
    NxNy_all(k) = Nx;
    
    % Step to next values for Nx and Ny
    Nx = Nx*2;
    Ny = Ny*2;
end

% Calculate order of accuracy
for i=2:num_splits
    order(i) = log(error(i-1)/error(i))/log(2);
end

% Print out data
fprintf('Resolution \t Error \t Order\n');
for i = 1:num_splits
    fprintf('%g \t %g \t %g\n', NxNy_all(i), error(i), order(i));
end