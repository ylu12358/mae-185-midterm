function [U] = prim2cons(rho,u,v,T,cv)
%PRIM2CONS Convert primitive vars to conservative.
%   [U] = prim2cons(rho,u,v,T,cv)

    % Get values of conservative variables
    rhou = rho.*u;
    rhov = rho.*v;
    e = cv*T;
    Et = rho.*(e+(u.^2+v.^2)/2);

    % Assign data to U
    [nx,ny] = size(u);
    U = zeros(4,nx,ny);
    U(1,:,:) = rho;
    U(2,:,:) = rhou;
    U(3,:,:) = rhov;
    U(4,:,:) = Et;
end