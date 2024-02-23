% Second-order second central difference function in x
function d2fdx2 = d2dx2(f,dx,bc)
    
    % set default value for 'bc'
    if nargin<3, bc = 'one-sided'; end
    
    % determine field size
    [nx,ny]     = size(f);

    % allocate return field
    d2fdx2      = zeros(nx,ny);
    
    dx2         = dx^2;
    
    % central difference
    for i=2:nx-1
        for j=1:ny
            d2fdx2(i,j) = (f(i+1,j)-2*f(i,j)+f(i-1,j))/dx2;
        end
    end
    
    switch bc
        case 'periodic'
            
            i = 1;
            for j=1:ny
                d2fdx2(i,j) = (f(i+1,j)-2*f(i,j)+f(nx,j))/dx2;
            end
            
            i = nx;
            for j=1:ny
                d2fdx2(i,j) = (f(1,j)-2*f(i,j)+f(i-1,j))/dx2;
            end
            
        otherwise
            % forward difference for first point
            i = 1;
            for j=1:ny
                d2fdx2(i,j) = (2*f(i,j)-5*f(i+1,j)+4*f(i+2,j)-f(i+3,j))/dx2;
            end
            
            % backward difference for last point
            i = nx;
            for j=1:ny
                d2fdx2(i,j) = (2*f(i,j)-5*f(i-1,j)+4*f(i-2,j)-f(i-3,j))/dx2;
            end
    end
    
end