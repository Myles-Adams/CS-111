function [t,v] = solve_ode_trapezoid(f,v_start, t_start, t_final, dt)
%solve_ode_trapezoid

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

% Use trapezoidal method iteratively until we reach given t_final
while t(n) < t_final
    
    % Check for overshooting
    if t(n)+dt > t_final
        % If new time t(n+1)=t(n)+dt is greater than t_final we reduce dt
        % so that that new time t(n+1) would be exactly t_final
        dt = t_final-t(n);
    end
    
    % Calculate new time
    t(n+1) = t(n) + dt;
    
    % Calculate slope at time t(n)
    slope_fn = f(t(n), v(n));
    
    % Approximate v at n+1
    v(n+1) = v(n) + dt*slope_fn;
    
    % Calculate slope at time t(n+1)
    slope_fn1 = f(t(n+1), v(n+1));
    
    % Calculate v at n+1
    v(n+1) = v(n) + dt*(1/2)*(slope_fn + slope_fn1);
    
    % Increase total count of iterations
    n = n + 1;
end


end

