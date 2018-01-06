function [sol] = f_example_sys(t, y)

% Preallocate array
sol = zeros(6,1);

% Use t and y to find solutions to system of equations
sol(1,1) = y(2,1);
sol(2,1) = sin(1 - exp(y(3,1)));
sol(3,1) = 1/(1 - log(y(5,1) - 1));
sol(4,1) = (.5*t^2+t)*y(2,1) + exp(y(3,1))*y(1,1);
sol(5,1) = 1 - y(5,1);
sol(6,1) = -3*((exp(y(3,1)) - 1)/(1+t^3))^2;

end