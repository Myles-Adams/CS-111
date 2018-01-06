clear;
clc;
close all;

% Domain
xL = -1;
xR = 3;
yB = -1.5;
yT = 1.5;

% Diffusion coefficient
D = 0.7;

% Velocities
Vx = -.8;
Vy = -.4;

% Boundary conditions
c_bc = @(t,x,y) sin(x)*cos(y)*exp(-t);
g = @(t,x,y) -D*sin(x)*sin(y)*exp(-t) - Vy*sin(x)*cos(y)*exp(-t);

% Source term
f = @(t,x, y) exp(-t)*(-sin(x)*cos(y)+Vx*cos(x)*cos(y)-Vy*sin(x)*sin(y)+2*D*sin(x)*cos(y));

% Initial conditions
c_start = @(x,y) sin(x)*cos(y);

% Time interval
t_start = 0;
t_final = 1;

% Exact solution
c_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% Space discretization
Nx = 20;
Ny = 15;

% Number of discretizations
num_splits = 4;

% Iterate through discretizations
for k = 1:num_splits
    
    % Discretize area
    x = linspace(xL, xR, Nx);
    y = linspace(yB, yT, Ny);
    dx = (xR-xL)/(Nx-1);
    dy = (yT-yB)/(Ny-1);
    
    % Time-step
    dt = .5*dx;
    
    % Pre-allocate arrays
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
        
        % Check for overshoot
        if t + dt > t_final
            dt = t_final-t;
            
            % Recalculate A
            A = make_A(1,Nx,Ny,dx,dy,dt,D,Vx,Vy);
        end
        
        % Calculate RHS array
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
        
        % Prepare for next iteration and iterate t
        c_old = c_new;
        t = t+dt;
        
    end
    
    % Calculate c_exact
    for i = 1:Nx
        for j = 1:Ny
            c_exa(j,i) = c_exact(t,x(i),y(j));
        end
    end
    
    % Calculate error and store
    error(k) = max(max(abs(c_exa - c_new)));
    
    % Store discretization sizes
    Nx_all(k) = Nx;
    Ny_all(k) = Ny;
    
    % Double discretization size
    Nx = Nx*2;
    Ny = Ny*2;
end

% Calculate order of accuracy
for i=2:num_splits
    order(i) = log(error(i-1)/error(i))/log(2);
end

% Print out data
fprintf('Resolution\tError\t\tOrder\n');
for i = 1:num_splits
    fprintf('%gx%g\t\t%g\t%g\n', Nx_all(i), Ny_all(i), error(i), order(i));
end