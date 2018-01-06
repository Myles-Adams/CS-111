function [output] = f_soccer(t, input, s, c_drag, c_lift, m, g)

% Preallocated output array
output = zeros(6,1);

% Create variables for ease of reading
v_x = input(4,1);
v_y = input(5,1);
v_z = input(6,1);
v_mag = sqrt((v_x)^2 + (v_y)^2 + (v_z)^2);
s_x = s(1,1);
s_y = s(2,1);
s_z = s(3,1);

% Set values for output array according to system of equations.
output(1,1) = v_x;
output(2,1) = v_y;
output(3,1) = v_z;
output(4,1) = (-c_drag*v_mag*v_x + c_lift*v_mag*(s_y*v_z - s_z*v_y))/m;
output(5,1) = (-c_drag*v_mag*v_y + c_lift*v_mag*(s_z*v_x - s_x*v_z))/m;
output(6,1) = (-g*m - c_drag*v_mag*v_z + c_lift*v_mag*(s_x*v_y - s_y*v_x))/m;

end