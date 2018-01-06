function [A] = make_A(type,Nx,Ny,dx,dy,dt,lambda,Vx,Vy)

A = sparse(Nx*Ny,Nx*Ny);

if (type == 0)
    
    a_c = 1 + 2*((lambda*dt)/(dx^2)) + 2*((lambda*dt)/(dy^2));
    a_l = -(lambda*dt)/(dx^2);
    a_r = a_l;
    a_b = -(lambda*dt)/(dy^2);
    a_t = a_b;
    
    for i = 1:Nx
        for j = 1:Ny
            p = (j-1)*Nx+i;
            if (i==1) || (i==Nx) || (j==1) || (j==Ny)
                A(p,p) = 1;
            else
                A(p,p) = a_c;
                A(p,p-1) = a_l;
                A(p,p+1) = a_r;
                A(p,p-Nx) = a_b;
                A(p,p+Nx) = a_t;
            end
        end
    end
    
elseif (type == 1)
    
    a_c = 1 + 2*((lambda*dt)/(dx^2)) + 2*((lambda*dt)/(dy^2));
    a_l = -(lambda*dt)/(dx^2);
    a_r = a_l;
    a_b = -(lambda*dt)/(dy^2);
    a_t = a_b;
    
    a_c_r = a_c + (2*Vy*dt)/(dy);
    a_l_r = -(lambda*dt)/(dx^2);
    a_r_r = a_l_r;
    a_t_r = a_t + a_b;
    
    for i = 1:Nx
        for j = 1:Ny
            p = (j-1)*Nx+i;
            if (i==1) || (i==Nx) || (j==Ny)
                A(p,p) = 1;
            elseif (j==1)
                A(p,p) = a_c_r;
                A(p,p+1) = a_r_r;
                A(p,p-1) = a_l_r;
                A(p,p+Nx) = a_t_r;
            else
                A(p,p) = a_c;
                A(p,p-1) = a_l;
                A(p,p+1) = a_r;
                A(p,p-Nx) = a_b;
                A(p,p+Nx) = a_t;
            end
            
        end
    end
    
    
end
end

