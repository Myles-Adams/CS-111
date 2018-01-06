function [vx,vy] = v_xy(t,x,y)

% Constants
r_c = .2;
alpha = .1;
U = .1;

% Check if within bounds of circle, if not calculate value of vx and vy
if sqrt((x^2 + y^2)) < r_c
    vx = 0;
    vy = 0;
else
    vx = U*(1 + ((r_c)^2)*((y^2 - x^2)/((x^2 + y^2)^2))) + 2*alpha*U*r_c*((y)/((x^2 + y^2)^(3/2)));
    vy = -2*U*((r_c)^2)*((x*y)/((x^2 + y^2)^2)) - 2*alpha*U*r_c*((x)/((x^2 + y^2)^(3/2)));
end

