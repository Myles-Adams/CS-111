clear;
clc;
close all;

% Constants
g = 9.81;
m = .437;
c_drag = 0.0057;
c_lift = 0.0061;
l_goal = 7.32;
x_goal_start = -.5*l_goal;
x_goal_end = .5*l_goal;
h_goal = 2.44;
h_wall = 2;

% Time
t_start = 0;
t_final = 1;
dt = .1;

% Initial values for example
y_initial = [0,1,0,0,2,1];

% Function for example system
f_example = @(t, y) f_example_sys(t,y);

% Solve example system using RK4
[t_num, sol_num] = solve_sys_rk4(f_example, y_initial, t_start, t_final, dt);

% Solve exact solution of example system
[t_exact, sol_exact] = solve_sys_exact(y_initial, t_start, t_final, dt);

% Create array for error values
error_sum = zeros(1, length(sol_num));

% Find error values
for n = 1:length(sol_num)
    for i = 1:6
        error_sum(n) = error_sum(n) + (sol_exact(i,n) - sol_num(i,n))^2;
    end
    error_sum(n) = sqrt(error_sum(n));
end

% Find max error
e_max = max(error_sum);
fprintf("Max Error: %g\n\n", e_max);

% Change time values
t_start = 0;
t_final = 2;
dt = .02;

% Initial values and constants for each case

% Case 1
x0 = -15;
y0 = 23;
z0 = 0;

vx0 = 18;
vy0 = -23;
vz0 = 8.5;

sx = .1;
sy = .1;
sz = -.99;

x_wall_start = -11.2;
y_wall_start = 14.7;
x_wall_end = -8.6;
y_wall_end = 16;

% Case 2
% x0 = 25.5;
% y0 = 15;
% z0 = 0;
% 
% vx0 = -28;
% vy0 = -11;
% vz0 = 7.5;
% 
% sx = 0;
% sy = 0;
% sz = 1;
% 
% x_wall_start = 18.1;
% y_wall_start = 9.7;
% x_wall_end = 17;
% y_wall_end = 11.3;

% Case 3
% x0 = -4;
% y0 = 35;
% z0 = 0;
% 
% vx0 = -8;
% vy0 = -36;
% vz0 = 5;
% 
% sx = -.32;
% sy = 0;
% sz = .95;
% 
% x_wall_start = -4.3;
% y_wall_start = 25.9;
% x_wall_end = -.8;
% y_wall_end = 25.9;

% Case 4
% x0 = 16;
% y0 = 28;
% z0 = 0;
% 
% vx0 = -25;
% vy0 = -20;
% vz0 = 8;
% 
% sx = .1;
% sy = .15;
% sz = .98;
% 
% x_wall_start = 12.6;
% y_wall_start = 19.6;
% x_wall_end = 9.9;
% y_wall_end = 20.8;

% Arrays for initial conditions and spin constants
init_cond = [x0, y0, z0, vx0, vy0, vz0];
s = [sx; sy; sz];

% Function for solving free kick system of equations
f_free_kick = @(t, v) f_soccer(t, v, s, c_drag, c_lift, m, g);

% Solve free kick system of equations
[t, output_soccer] = solve_sys_rk4(f_free_kick, init_cond, t_start, t_final, dt);

% Plot the free kick
plot_trajectory(output_soccer(1,:), output_soccer(2,:), output_soccer(3,:), x_wall_start, y_wall_start, x_wall_end, y_wall_end, h_wall)

% Function for checking if ball is passing wall
f_wall = @(x) y_wall_start + (x - x_wall_start)*((y_wall_end - y_wall_start)/(x_wall_end - x_wall_start));

hitWall = 0;

% Checking if ball is hitting wall and going in the goal at each time step
for i=1:(length(output_soccer) -1)
    
    % Check if hits wall
    if (output_soccer(2,i) > f_wall(output_soccer(1,i))) && (output_soccer(2,i+1) < f_wall(output_soccer(1,i+1)))
        
        % Find position of the ball as it passes the wall
        pos_star_wall = f_pos_star_wall(output_soccer, i, x_wall_start, y_wall_start, x_wall_end, y_wall_end);
        x_star_wall = pos_star_wall(1);
        y_star_wall = pos_star_wall(2);
        z_star_wall = pos_star_wall(3);
        
        % Check if it is in bounds of the wall
        if ((x_star_wall >= x_wall_start) && (x_star_wall <= x_wall_end)) || ((x_star_wall <= x_wall_start) && (x_star_wall >= x_wall_end))
            if ((y_star_wall >= y_wall_start) && (y_star_wall <= y_wall_end)) || ((y_star_wall <= y_wall_start) && (y_star_wall >= y_wall_end))
                if ((z_star_wall >= 0) && (z_star_wall <= h_wall))
                    fprintf("Hit Wall\n");
                    hitWall = 1;
                end
            end
        end
    
    % Check if going in goal
    elseif (output_soccer(2,i) > 0) && (output_soccer(2,i+1) < 0)
        
        % Find position of ball as it crosses goal line
        pos_star_goal = f_pos_star_goal(output_soccer, t, i);
        x_star_goal = pos_star_goal(1);
        z_star_goal = pos_star_goal(3);
        
        % Check if its in bounds of the goal
        if ((x_star_goal > x_goal_start) && (x_star_goal < x_goal_end)) && ((z_star_goal > 0) && (z_star_goal < h_goal))
            fprintf("Goal Scored\n");
            break;
        else
            fprintf("Missed Goal\n");
            break;
        end
    end
end

if hitWall == 0
    fprintf("Missed Wall\n");
end