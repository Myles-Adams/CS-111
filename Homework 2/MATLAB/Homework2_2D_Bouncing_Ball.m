% Outline of the code simulating a ball in a closed container using the
% Euler method (it is not a complete code and does not run)
clear;
clc;
close all;
 
% boundaries of the container
a = 0;
b = 1;
c = 0;
d = 1;
 
% damping coefficients
alpha = 0.8;
beta = 0.9;
 
% free-fall acceleration
g = 9.81;
 
% radius of the ball
r = 0.05;

% Functions for calculating values
f_x = @(dt, vx) dt*vx;
f_y = @(dt, vy) dt*vy;
f_vy = @(dt, vy_old, g) vy_old - g*dt;
 
% Initial conditions
t_start = 0;
 
x_start = 0.1;
y_start = 0.7;
 
vx_start = 3;
vy_start = 1;
 
% Final time and time-step
dt = 0.01;
t_final = .931;
 
% Impose initial conditions
t = t_start;
 
x_old = x_start;
y_old = y_start;
 
vx_old = vx_start;
vy_old = vy_start;
 
% Exact solution at the final time (doesn't work for more than three bounces)
[x_exact, y_exact, vx_exact, vy_exact] = bouncing_ball_exact(t,a,b,c,d,r,g,alpha,beta,x_start,y_start,vx_start,vy_start);
 
% Auxiliary arrays to plot the trajectory
x_all = x_start;
y_all = y_start;

% Max error values
x_maxerror = 0;
y_maxerror = 0;
vx_maxerror = 0;
vy_maxerror = 0;
 
while t < t_final
   
    % Reset dt_now value
    dt_now = dt;
   
    % Check for overshoot of t_final
    if t + dt > t_final
        dt_now = t_final-t;
    end
   
    % Trapezoidal step
    % Calculate new x value
    x_new = x_old + f_x(dt_now, vx_old);
    
    % Calculate new vy and y values
    vy_temp = f_vy(dt_now, vy_old, g);
    y_new = y_old + .5*(f_y(dt_now, vy_old)+f_y(dt_now, vy_temp));
   
    % Reset the vx and vy values for collision checks
    vx_new = vx_old;
    vy_new = vy_old;
 
   
    % check for a collision with the right wall
    if x_new > b-r
       
        % Calculate dt between old time and wall collision
        dt_now = dt*((b-r)-x_old)/(x_new-x_old);
        
        % Calculate new vy right before collision
        vy_new = f_vy(dt_now, vy_old, g);
           
        % Calculate the new x and y values
        x_new = x_old + .5*(f_x(dt_now, vx_old)+f_x(dt_now, vx_old));
        y_new = y_old + .5*(f_y(dt_now, vy_old)+f_y(dt_now, vy_new));
       
        % Calculate vx and vy right after collision
        vx_new = -alpha*vx_new;
        vy_new =  beta*vy_new;
        
        % Check for a collision with the left wall
    elseif x_new < r + a
        
        % Calculate dt between old time and wall collision
        dt_now = dt*(x_old - (r + a))/(x_old - x_new);
        
        % Calculate new vy right before collision
        vy_new = f_vy(dt_now, vy_old, g);
           
        % Calculate the new x and y values
        x_new = x_old + .5*(f_x(dt_now, vx_old)+f_x(dt_now, vx_old));
        y_new = y_old + .5*(f_y(dt_now, vy_old)+f_y(dt_now, vy_new));
        
        % Calculate vx and vy right after collision
        vx_new = -alpha*vx_new;
        vy_new =  beta*vy_new;
        
       % Check for a collision with the top wall
    elseif y_new > d - r
       
        % Calculate dt between old time and wall collision
        dt_now = dt_now*((d-r) - y_old)/(y_new - y_old);
        
        % Calculate new x and y values
        x_new = x_old + .5*(f_x(dt_now, vx_old)+f_x(dt_now, vx_old));
        y_new = (d-r);
        
        % Calculate new vy value right before collision
        vy_new = f_vy(dt_now, vy_old, g);
        
        % Calculate new vx and vy values right after collision
        vx_new = beta*vx_new;
        vy_new = -alpha*(vy_new);
       
        % Check for a collision with the bottom wall
    elseif y_new < r + c
       
        % Calculate dt between old time and wall collision
        dt_now = dt_now*(y_old-(r+c))/(y_old-y_new);
        
        % Calculate new x and y values
        x_new = x_old + .5*(f_x(dt_now, vx_old)+f_x(dt_now, vx_old));
        y_new = (r+c);
        
        % Calculate new vy value right before collision
        vy_new = f_vy(dt_now, vy_old, g);
        
        % Calculate new vx and vy values right after collision
        vx_new = beta*vx_new;
        vy_new = -alpha*(vy_new);
        
       % No collision 
    else
        
        % Calculate the new vy
        vy_new = vy_temp;
        
    end
    
    % Save trajectory points to plot
    x_all = [x_all x_new];
    y_all = [y_all y_new];
     
    % Advance time
    t = t+dt_now;
   
    % Plot the ball and its trajectory
    draw_disk(x_new, y_new, r); axis([a b, c d]); axis square;
    [x_exact, y_exact, vx_exact, vy_exact] = bouncing_ball_exact(t,a,b,c,d,r,g,alpha,beta,x_start,y_start,vx_start,vy_start);
    
    hold on
    plot(x_all, y_all,'.-');
    hold off
    pause(dt);
    
    % Print values to compare numerical and exact solutions
    fprintf("Time: %f\n", t);
    fprintf("x_exact = %f x_new = %f \n", x_exact, x_new);
    fprintf("y_exact = %f y_new = %f \n", y_exact, y_new);
    fprintf("vx_exact = %f vx_new = %f \n", vx_exact, vx_new);
    fprintf("vy_exact = %f vy_new = %f \n", vy_exact, vy_new);
   
    % Prepare values for next iteration
    x_old = x_new;
    y_old = y_new;
   
    vx_old = vx_new;
    vy_old = vy_new;
   
end

% Print the error values
fprintf("X Error: %g\n", abs(x_exact - x_new));
fprintf("y Error: %g\n", abs(y_exact - y_new));
fprintf("vX Error: %g\n", abs(vx_exact - vx_new));
fprintf("vy Error: %g\n", abs(vy_exact - vy_new));