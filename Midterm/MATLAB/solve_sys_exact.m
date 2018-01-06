function [t, y] = solve_sys_exact(y0, t_start, t_final, dt)

% Find number of steps
n_total = ceil((t_final - t_start)/dt) + 1;

% Preallocated arrays
t = zeros(1, n_total);
y = zeros(6, n_total);

% Put initial values into y
y(:,1) = y0;

n = 1;

% Loop through all values of t
while t(n) < t_final
    
    % Check for overshoot
    if (t(n) + dt) > t_final
        dt = t_final - t(n);
    end
    
    % Calculate new time
    t(n+1) = t(n) + dt;
    
    % Calculate values using exact formulas
    y(1, n+1) = sin(t(n+1));
    y(2, n+1) = cos(t(n+1));
    y(3, n+1) = log(1 + t(n+1));
    y(4, n+1) = (.5*(t(n+1))^2 + t(n+1))*sin(t(n+1));
    y(5, n+1) = 1 + exp(-1*t(n+1));
    y(6, n+1) = 1/(1 + (t(n+1))^3);
    
    
    % Increase total count of iterations
    n = n + 1;
     
end
end