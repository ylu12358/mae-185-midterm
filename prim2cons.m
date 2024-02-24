% Converts primitive variables into conservative variables
function U = prim2cons(rho,u,v,T,cv)
    % Get values of conservative variables
    rhou = rho.*u;
    rhov = rho.*v;
    e = cv.*T;
    Et = rho.*(e+(u.^2+v.^2)./2);
    % Initialize arrays/values for looping through data
    U = zeros([4, size(rho)]);
    nx = size(rho,1);
    ny = size(rho,2);
    for i = 1:nx
        for j = 1:ny
            % Input data into U for return
            U(:,i,j) = [rho(i,j); rhou(i,j); rhov(i,j); Et(i,j)];
        end
    end
end