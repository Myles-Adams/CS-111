% Falling Object - Homework #1

clear;

% Set parameters for this problem
g = 9.81; % acceleration of free-fall
m = 75; % mass
cd = 0.25; % drag coefficient

% Function for free-fall
f_parachute = @(t,v) g - (cd/m)*v^2;

% Set final and initial time, start velocity, and time-step
t_start = 0;
t_final = 15;
v_start = 0;
dt = 3.3;

% Find data using all three methods
[t,v_euler_dt] = solve_ode_euler(f_parachute, v_start, t_start, t_final, dt);
[t,v_trap_dt] = solve_ode_trapezoid(f_parachute, v_start, t_start, t_final, dt);
[t,v_rk4_dt] = solve_ode_rk4(f_parachute, v_start, t_start, t_final, dt);

% Create arrays for exact solution and errors for the three methods and
% value for size of arrays
n = ceil((t_final - t_start)/dt);
v_exact_dt = zeros(1, n);
error_euler_dt = zeros(1, n);
error_trap_dt = zeros(1, n);
error_rk4_dt = zeros(1, n);

% Calculate the exact solution
for i = 1:length(t)
    v_exact_dt(i) = sqrt(g*m/cd)*tanh(t(i)*sqrt(g*cd/m));
    error_euler_dt(i) = abs(v_euler_dt(i) - v_exact_dt(i));
    error_trap_dt(i) = abs(v_trap_dt(i) - v_exact_dt(i));
    error_rk4_dt(i) = abs(v_rk4_dt(i) - v_exact_dt(i));
end

% Plot euler, trapezoidal, RK4, and exact on same graph
plot(t, v_euler_dt, 'o', t, v_trap_dt, 'o', t, v_rk4_dt, 'o', t, v_exact_dt)
xlabel('Time');
ylabel('Velocity');
legend('Euler', 'Trapezoidal', 'RK4', 'Exact', 'Location', 'southeast');

% Put values for dt into an array
dt_values = [.1, .05, .025, .0125];

fprintf('Euler''s Method\n______________\n');
fprintf('Time-Step\tMax Error\tOrder\n');
for i = 1:4
    
    % Create arrays for exact values and errors for dt
    n = ceil((t_final - t_start)/dt_values(i));
    v_exact_dt = zeros(1, n);
    error_euler_dt = zeros(1, n);
    
    % Create arrays for exact values and errors for dt/2
    x = ceil((t_final - t_start)/(dt_values(i)/2));
    v_exact_dt2 = zeros(1,x);
    error_euler_dt2 = zeros(1, x);
    
    % Calculate numerical values for Euler's Method using dt and dt/2
    [t_dt,v_euler_dt] = solve_ode_euler(f_parachute, v_start, t_start, t_final, dt_values(i));
    [t_dt2,v_euler_dt2] = solve_ode_euler(f_parachute, v_start, t_start, t_final, dt_values(i)/2);
    
    % Find the exact value and error at every point for dt
    for j = 1:length(t_dt)
        v_exact_dt(j) = sqrt(g*m/cd)*tanh(t_dt(j)*sqrt(g*cd/m));
        error_euler_dt(j) = abs(v_euler_dt(j) - v_exact_dt(j));
    end
    
    % Find the max error for dt
    emax_euler_dt = max(error_euler_dt);
    
    % Find the exact value and error at every point for dt/2
    for j = 1:length(t_dt2)
        v_exact_dt2(j) = sqrt(g*m/cd)*tanh(t_dt2(j)*sqrt(g*cd/m));
        error_euler_dt2(j) = abs(v_euler_dt2(j) - v_exact_dt2(j));
    end
    
    % Calculate max error for (dt/2)
    emax_euler_dt2 = max(error_euler_dt2);
    
    % Calculate order of accuracy for Euler's method
    order_euler = log(emax_euler_dt/emax_euler_dt2)/log(2);
    
    % Print out the information
    fprintf('%f\t%g\t%f\n', dt_values(i), emax_euler_dt, order_euler);
end

fprintf('\n');
fprintf('Trapezoidal Method\n______________\n');
fprintf('Time-Step\tMax Error\tOrder\n');

for i = 1:4
    
    % Create arrays for exact values and errors for dt
    n = ceil((t_final - t_start)/dt_values(i));
    v_exact_dt = zeros(1, n);
    error_trap_dt = zeros(1, n);
    
    % Create arrays for exact values and errors for dt/2
    x = ceil((t_final - t_start)/(dt_values(i)/2));
    v_exact_dt2 = zeros(1, x);
    error_trap_dt2 = zeros(1, x);
    
    % Calculate numerical values for Trapezoidal Method using dt and dt/2
    [t_dt,v_trap_dt] = solve_ode_trapezoid(f_parachute, v_start, t_start, t_final, dt_values(i));
    [t_dt2,v_trap_dt2] = solve_ode_trapezoid(f_parachute, v_start, t_start, t_final, dt_values(i)/2);
    
    % Find the exact value and error at every point for dt
    for j = 1:length(t_dt)
        v_exact_dt(j) = sqrt(g*m/cd)*tanh(t_dt(j)*sqrt(g*cd/m));
        error_trap_dt(j) = abs(v_trap_dt(j) - v_exact_dt(j));
    end
    
    % Find the max error for dt
    emax_trap_dt = max(error_trap_dt);
    
    % Find the exact value and error at every point for dt/2
    for j = 1:length(t_dt2)
        v_exact_dt2(j) = sqrt(g*m/cd)*tanh(t_dt2(j)*sqrt(g*cd/m));
        error_trap_dt2(j) = abs(v_trap_dt2(j) - v_exact_dt2(j));
    end
    
    % Calculate max error for (dt/2)
    emax_trap_dt2 = max(error_trap_dt2);
    
    % Calculate order of accuracy for the trapezoidal method
    order_trap = log(emax_trap_dt/emax_trap_dt2)/log(2);
    
    % Print out the information
    fprintf('%f\t%g\t%f\n', dt_values(i), emax_trap_dt, order_trap);
end

fprintf('\n');
fprintf('RK4 Method\n___________\n');
fprintf('Time-Step\tMax Error\tOrder\n');

for i = 1:4
    
    % Create arrays for exact values and errors for dt
    n = ceil((t_final - t_start)/dt_values(i));
    v_exact_dt = zeros(1, n);
    error_rk4_dt = zeros(1, n);
    
    % Create arrays for exact values and errors for dt/2
    x = ceil((t_final - t_start)/(dt_values(i)/2));
    v_exact_dt2 = zeros(1, x);
    error_rk4_dt2 = zeros(1, x);
    
    % Calculate numerical values for RK4 method using dt and dt/2
    [t_dt,v_rk4_dt] = solve_ode_rk4(f_parachute, v_start, t_start, t_final, dt_values(i));
    [t_dt2,v_rk4_dt2] = solve_ode_rk4(f_parachute, v_start, t_start, t_final, dt_values(i)/2);
    
    % Find the exact value and error at every point for dt
    for j = 1:length(t_dt)
        v_exact_dt(j) = sqrt(g*m/cd)*tanh(t_dt(j)*sqrt(g*cd/m));
        error_rk4_dt(j) = abs(v_rk4_dt(j) - v_exact_dt(j));
    end
    
    % Find the max error for dt
    emax_rk4_dt = max(error_rk4_dt);
    
    % Find the exact value and error at every point for dt/2
    for j = 1:length(t_dt2)
        v_exact_dt2(j) = sqrt(g*m/cd)*tanh(t_dt2(j)*sqrt(g*cd/m));
        error_rk4_dt2(j) = abs(v_rk4_dt2(j) - v_exact_dt2(j));
    end
    
    % Calculate max error for (dt/2)
    emax_rk4_dt2 = max(error_rk4_dt2);
    
    % Calculate order of accuracy for the RK4 method
    order_rk4 = log(emax_rk4_dt/emax_rk4_dt2)/log(2);
    
    % Print out the information
    fprintf('%f\t%g\t%f\n', dt_values(i), emax_rk4_dt, order_rk4);
end