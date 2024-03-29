function [u,v,p,T,rho] = enforceBC(u,v,p,T,rho,pinf,Tinf,uinf,R,bc)
%enforceBC Enforces boundary conditions.
%   [u,v,p,T,rho] = enforceBC(u,v,p,T,rho,pinf,Tinf,uinf,R,bc)
    
    % Enforce BCs on primitive variables

    % Left boundary (inlet)
    u(1,:) = uinf;
    v(1,:) = 0;
    p(1,:) = pinf;
    T(1,2:end) = Tinf;
    rho(1,:) = p(1,:)./(R.*T(1,:));
    
    % Right boundary (outlet)
    u(end,:) = 2*u(end-1,:)-u(end-2,:);
    v(end,:) = 2*v(end-1,:)-v(end-2,:);
    p(end,:) = 2*p(end-1,:)-p(end-2,:);
    T(end,2:end) = 2*T(end-1,2:end)-T(end-2,2:end);
    rho(end,:) = p(end,:)./(R.*T(end,:));
    
    % Top boundary (far-field)
    u(:,end) = uinf;
    v(:,end) = 0;
    p(:,end) = pinf;
    T(:,end) = Tinf;
    rho(:,end) = p(:,end)./(R.*T(:,end));
    
    % Bottom boundary (wall)
    u(:,1) = 0;
    v(:,1) = 0;
    p(:,1) = 2*p(:,2)-p(:,3);
    rho(:,1) = p(:,1)./(R.*T(:,1));

    % If statement for temperature at wall BC
    if bc == "isothermal"
        T(:,1) = Tinf;
    end
    
    % Bottom left corner point (leading edge)
    u(1,1) = 0;
    v(1,1) = 0;
    p(1,1) = pinf;
    T(1,1) = Tinf;
end