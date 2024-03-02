function d2fdy2 = d2dy2(f,dy,bc)
%D2DY2 Computes second-order second central difference function in y
%   d2fdy2 = d2dy2(f,dy,bc)
    
    % set default value for 'bc'
    if nargin<3, bc = 'one-sided'; end
    
    % determine field size
    [nx,ny]     = size(f);    
    
    % allocate return field
    d2fdy2      = zeros(nx,ny);
    
    dy2         = dy^2;
    
    % central difference
    for i=1:nx
        for j=2:ny-1
            d2fdy2(i,j) = (f(i,j+1)-2*f(i,j)+f(i,j-1))/dy2;
        end
    end
    
    switch bc
        case 'periodic'
            
            j = 1;
            for i=1:nx
                d2fdy2(i,j) = (f(i,j+1)-2*f(i,j)+f(i,ny))/dy2;
            end
            
            j = ny;
            for i=1:nx
                d2fdy2(i,j) = (f(i,1)-2*f(i,j)+f(i,j-1))/dy2;
            end
            
        otherwise
            % forward difference for first point (bottom)
            j = 1;
            for i=1:nx
                d2fdy2(i,j) = (2*f(i,j)-5*f(i,j+1)+4*f(i,j+2)-f(i,j+3))/dy2;
            end
            
            % backward difference for last point (top)
            j = ny;
            for i=1:nx
                d2fdy2(i,j) = (2*f(i,j)-5*f(i,j-1)+4*f(i,j-2)-f(i,j-3))/dy2;
            end
    end
    
end