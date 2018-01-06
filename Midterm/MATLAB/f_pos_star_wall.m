function [pos_star] = f_pos_star_wall(pos, ndx, x_wall_start, y_wall_start, x_wall_end, y_wall_end)

% Create variables for ease of reading
x_n = pos(1,ndx);
x_n1 = pos(1,ndx+1);
y_n = pos(2,ndx);
y_n1 = pos(2,ndx+1);
z_n = pos(3,ndx);
z_n1 = pos(3, ndx+1);

% Solve small part of equations to simplify equations for x_star, y_star,
% and z_star
A = (y_n1 - y_n)/(x_n1 - x_n);
B = (y_wall_end - y_wall_start)/(x_wall_end - x_wall_start);

% Calculate values for x_star, y_star, and z_star 
x_star = (y_wall_start - y_n + x_n*A - x_wall_start*B)/(A-B);
y_star = y_n + (x_star - x_n)*((y_n1 - y_n)/(x_n1 - x_n));
z_star = z_n + (x_star - x_n)*((z_n1 - z_n)/(x_n1 - x_n));

% Put values into array for returning
pos_star = [x_star;y_star;z_star];

end