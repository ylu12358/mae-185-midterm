function dfdx = ddx_bwd(f,dx,bc)
%DDX_BWD Computes first-order backward difference function in x
%   dfdx = ddx_bwd(f,dx,bc)

    % set default value for 'bc'
    if nargin<3, bc = 'one-sided'; end

    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdx        = zeros(nx,ny);

    % backward difference
    for i=2:nx
        for j=1:ny
            dfdx(i,j) = (f(i,j)-f(i-1,j))/dx;
        end
    end

    switch bc
        case 'periodic'

            % assuming periodicity (left boundary)
            i = 1;
            for j=1:ny
                dfdx(i,j) = (f(i,j)-f(end,j))/dx;
            end

        otherwise

            % forward difference for first point
            i = 1;
            for j=1:ny
                dfdx(i,j) = (f(i+1,j)-f(i,j))/dx;
            end      
    end
end