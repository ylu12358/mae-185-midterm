% First-order forward difference difference function in x
function dfdx = ddx_fwd(f,dx,bc)

    % set default value for 'bc'
    if nargin<3, bc = 'one-sided'; end

    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    dfdx        = zeros(nx,ny);

    % forward difference
    for i=1:nx-1
        for j=1:ny
            dfdx(i,j) = (f(i+1,j)-f(i,j))/dx;
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

            % backward difference for last point
            i = nx;
            for j=1:ny
                dfdx(i,j) = (f(i,j)-f(i-1,j))/dx;
            end      
    end
end