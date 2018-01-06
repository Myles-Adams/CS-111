clear;
clc;
close all;

% domain
xL = -1;
xR = 1;
yB = -.5;
yT = 1.7;

% diffusion coefficient
lambda = 0.75;

% boundary conditions
T_bc = @(t,x,y) sin(x)*cos(y)*exp(-t);

% source term
f = @(t,x, y) (2*lambda - 1)*sin(x)*cos(y)*exp(-t);

% initial conditions
T_start = @(x,y) sin(x)*cos(y);

% time interval
t_start = 0;
t_final = 1;

% exact solution
T_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% space discretization
Nx = 25;
Ny = 30;

num_splits = 3;

for k = 1:num_splits
    
    x = linspace(xL, xR, Nx);
    y = linspace(yB, yT, Ny);
    dx = (xR-xL)/(Nx-1);
    dy = (yT-yB)/(Ny-1);
    
    % time-step
    dt = .5*dx;
        
    T_old = zeros(Ny,Nx);
    T_new = zeros(Ny,Nx);
    T_exa = zeros(Ny,Nx);
    
    for i = 1:Nx
        for j =1:Ny
            T_old(j,i) = T_start(x(i), y(j));
        end
    end
    
    t = t_start;
    
    % Create sparse matrix and allocate memory for right-hand side
    RHS = zeros(Nx*Ny,1);
    
    % Calculate the matrix before the while-loop to save time
    % internal points
    A = make_A(0,Nx,Ny,dx,dy,dt,lambda);
    
    while t < t_final
        
        if t + dt > t_final
            dt = t_final-t;
            
            A = make_A(0,Nx,Ny,dx,dy,dt,lambda);
        end
        
        % Make RHS
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
        
        T_old = T_new;
        
        for i = 1:Nx
            for j = 1:Ny
                T_exa(j,i) = T_exact(t+dt,x(i),y(j));
            end
        end
        
        
        
        t = t+dt;
        
    end
    
    error(k) = max(max(abs(T_exa - T_new)));
    Nx_all(k) = Nx;
    Ny_all(k) = Ny;
    
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