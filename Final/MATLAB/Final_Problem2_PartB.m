clear;
clc;
close all;

% domain
xL = -1;
xR = 3;
yB = -1.5;
yT = 1.5;

% diffusion coefficient
D = 0.7;

% Velocities
Vx = -.8;
Vy = -.4;

% boundary conditions
c_bc = @(t,x,y) sin(x)*cos(y)*exp(-t);
g = @(t,x,y) -D*sin(x)*sin(y)*exp(-t) - Vy*sin(x)*cos(y)*exp(-t);

% source term
f = @(t,x, y) exp(-t)*(-sin(x)*cos(y)+Vx*cos(x)*cos(y)-Vy*sin(x)*sin(y)+2*D*sin(x)*cos(y));

% initial conditions
c_start = @(x,y) sin(x)*cos(y);

% time interval
t_start = 0;
t_final = 1;

% exact solution
c_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% space discretization
Nx = 20;
Ny = 15;

num_splits = 4;

for k = 1:num_splits
    
    x = linspace(xL, xR, Nx);
    y = linspace(yB, yT, Ny);
    dx = (xR-xL)/(Nx-1);
    dy = (yT-yB)/(Ny-1);
    
    % time-step
    dt = .5*dx;
    
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
        
        c_old = c_new;
        t = t+dt;
        
        for i = 1:Nx
            for j = 1:Ny
                c_exa(j,i) = c_exact(t,x(i),y(j));
            end
        end
        error_test = abs(c_exa - c_new);
        
    end
    

    
    error(k) = max(max(abs(c_exa - c_new)));
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