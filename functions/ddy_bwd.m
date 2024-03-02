function dfdy = ddy_bwd(f,dy,bc)
%DDY_BWD Computes first-order backward difference function in y
%   dfdy = ddy_bwd(f,dy,bc)

    % set default value for 'bc'
    if nargin<3, bc = 'one-sided'; end

    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdy        = zeros(nx,ny);

    % backward difference
    for i=1:nx
        for j=2:ny
            dfdy(i,j) = (f(i,j)-f(i,j-1))/dy;
        end
    end

    switch bc
        case 'periodic'

            % assuming periodicity (bottom boundary)
            j = 1;
            for i=1:nx
                dfdy(i,j) = (f(i,j)-f(i,ny))/dy;
            end

        otherwise

            % forward difference for first point
            j = 1;
            for i=1:nx
                dfdy(i,j) = (f(i,j+1)-f(i,j))/dy;
            end
    end
end