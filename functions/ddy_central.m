function dfdy = ddy_central(f,dy)
%DDY_CENTRAL Computes second-order central first difference function in y
%   dfdy = ddy_central(f,dy)
    
    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny);
    
    % central difference
    for i=1:nx
        for j=2:ny-1
            dfdy(i,j) = (f(i,j+1)-f(i,j-1))/2/dy;
        end
    end
    
    % forward difference for first point
    j = 1;
    for i=1:nx
        dfdy(i,j) = (-3*f(i,j)+4*f(i,j+1)-f(i,j+2))/2/dy;
    end
    
    % backward difference for last point
    j = ny;
    for i=1:nx
        dfdy(i,j) = (3*f(i,j)-4*f(i,j-1)+f(i,j-2))/2/dy;
    end
    
end