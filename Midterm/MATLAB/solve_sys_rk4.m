function [t, output] = solve_sys_rk4(f, init_cond, t_start, t_final, dt)

% Find number of steps
n_total = ceil((t_final - t_start)/dt) + 1;

% Find number of initial conditions
m = length(init_cond);

% Preallocated array
t = zeros(1, n_total);
output = zeros(m, n_total);

% Put initial conditions in output array
output(:,1) = init_cond;

n = 1;

% Loop through all values of t
while t(n) < t_final
    
    % Check for overshoot
    if (t(n) + dt) > t_final
        dt = t_final - t(n);
    end
    
    % Calculate new time
    t(n+1) = t(n) + dt;
    
    % Calculate all values for RK4 formula
    k1_v = f(t(n), output(:, n));
    k2_v = f(t(n) + (1/2)*(dt), output(:, n) + (1/2)*(dt)*(k1_v));
    k3_v = f(t(n) + (1/2)*(dt), output(:, n) + (1/2)*(dt)*(k2_v));
    k4_v = f(t(n) + dt, output(:, n) + (dt)*(k3_v));
    
    % Calculate output using RK4 formula
    output(:, n+1) = output(:, n) + (dt)*((k1_v+(2*k2_v)+(2*k3_v)+k4_v)/6);
    
    % Increase total count of iterations
    n = n + 1;
end
end