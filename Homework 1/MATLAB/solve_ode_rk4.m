function [t,v] = solve_ode_rk4(f,v_start, t_start, t_final, dt)
%solve_ode_rk4

% Find size of arrays
n_total = ceil((t_final - t_start)/dt);

% Create preallocated arrays
v = zeros(1, n_total);
t = zeros(1, n_total);

% Put initial values into arrays
v(1) = v_start;
t(1) = t_start;

% Initialize auxiliary variable to keep track of the number of iterations
n = 1;

% Use RK4 method iteratively until we reach given t_final
while t(n) < t_final
    
    % Check for overshooting
    if t(n)+dt > t_final
        % If new time t(n+1)=t(n)+dt is greater than t_final we reduce dt
        % so that that new time t(n+1) would be exactly t_final
        dt = t_final-t(n);
    end
    
    % Calculate new time
    t(n+1) = t(n) + dt;
    
    % Calculate all values for RK4 formula
    k1 = f(t(n), v(n));
    k2 = f(t(n) + (1/2)*(dt), v(n) + (1/2)*(dt)*(k1));
    k3 = f(t(n) + (1/2)*(dt), v(n) + (1/2)*(dt)*(k2));
    k4 = f(t(n) + dt, v(n) + (dt)*(k3));
    
    % Calculate velocity using RK4 formula
    v(n+1) = v(n) + (dt)*((k1+(2*k2)+(2*k3)+k4)/6);
    
    % Increase total count of iterations
    n = n + 1;
end

end

