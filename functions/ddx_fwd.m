function dfdx = ddx_fwd(f,dx,bc)
%DDX_FWD Computes first-order forward difference function in x
%   dfdx = ddx_fwd(f,dx,bc)

    % set default value for 'bc'
    if nargin<3, bc = 'one-sided'; end

    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdx        = zeros(nx,ny);

    % forward difference
    for i=1:nx-2
        for j=1:ny
            dfdx(i,j) = (-3*f(i,j)+4*f(i+1,j)-f(i+2,j))/2/dx;
        end
    end

    switch bc
        case 'periodic'

            % assuming periodicity (right boundary)
            i = nx;
            for j=1:ny
                dfdx(i,j) = (f(1,j)-f(i,j))/dx;
            end

        otherwise

            % central difference for second to last point
            i = nx-1;
            for j=1:ny
                dfdx(i,j) = (f(i+1,j)-f(i-1,j))/2/dx;
            end    

            % backward difference for last point
            i = nx;
            for j=1:ny
                dfdx(i,j) = (3*f(i,j)-4*f(i-1,j)+f(i-2,j))/2/dx;
            end      
    end
end