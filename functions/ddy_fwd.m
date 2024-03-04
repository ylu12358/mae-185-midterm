function dfdy = ddy_fwd(f,dy,bc)
%DDY_FWD Computes first-order forward difference function in y
%   dfdy = ddy_fwd(f,dy,bc)

    % set default value for 'bc'
    if nargin<3, bc = 'one-sided'; end

    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny);

    % forward difference
    for i=1:nx
        for j=1:ny-2
            dfdy(i,j) = (-3*f(i,j)+4*f(i,j+1)-f(i,j+2))/2/dy;
        end
    end

    switch bc
        case 'periodic'

            % assuming periodicity (top boundary)
            j = ny;
            for i=1:nx
                dfdy(i,j) = (f(i,j)-f(i,ny))/dy;
            end

        otherwise

            % central difference for second to last point
            j = ny-1;
            for i=1:nx
                dfdy(i,j) = (f(i,j+1)-f(i,j-1))/2/dy;
            end

            % backward difference for last point
            j = ny;
            for i=1:nx
                dfdy(i,j) = (3*f(i,j)-4*f(i,j-1)+f(i,j-2))/2/dy;
            end
    end
end