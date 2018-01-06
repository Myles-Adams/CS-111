clear;
clc;
close all;

% Domain
xL = -1;
xR = 1;
yB = -.5;
yT = 1.7;

% Diffusion coefficient
lambda = 0.75;

% Boundary conditions
T_bc = @(t,x,y) sin(x)*cos(y)*exp(-t);

% Source term
f = @(t,x, y) (2*lambda - 1)*sin(x)*cos(y)*exp(-t);

% initial conditions
T_start = @(x,y) sin(x)*cos(y);

% Time interval
t_start = 0;
t_final = 1;

% Exact solution
T_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% Space discretization
Nx = 25;
Ny = 30;

num_splits = 3;

% Run problem for each discretization
for k = 1:num_splits
    
    % Discretize the potato
    x = linspace(xL, xR, Nx);
    y = linspace(yB, yT, Ny);
    dx = (xR-xL)/(Nx-1);
    dy = (yT-yB)/(Ny-1);
    
    % Time-step
    dt = .5*dx;
    
    % Allocate arrays
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
    
    % Loop through all values of t
    while t < t_final
        
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
        
        % solve system of equations
        T_new = reshape(A\RHS,Nx,Ny)';
        
        % Prepare for next iteration
        T_old = T_new;
        
        % Incremenet t
        t = t+dt;
        
    end
    
    % Calculate t_exact
    for i = 1:Nx
        for j = 1:Ny
            T_exa(j,i) = T_exact(t,x(i),y(j));
        end
    end
    
    % Solve error and store
    error(k) = max(max(abs(T_exa - T_new)));

    % Store discretization sizes
    Nx_all(k) = Nx;
    Ny_all(k) = Ny;
    
    % Change discretization sizes
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