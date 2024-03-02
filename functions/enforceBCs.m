function [u,v,p,T,rho] = enforceBCs(u,v,p,T,rho,pinf,Tinf,uinf,R)
    
    % Enforce BCs on primitive variables
    % Left boundary
    u(1,:) = uinf;
    v(1,:) = 0;
    p(1,:) = pinf;
    T(1,:) = Tinf;
    rho(1,:) = p(1,:)./(R.*T(1,:));
    
    % Right boundary
    u(end,:) = 2*u(end-1,:)-u(end-2,:);
    v(end,:) = 2*v(end-1,:)-v(end-2,:);
    p(end,:) = 2*p(end-1,:)-p(end-2,:);
    T(end,:) = 2*T(end-1,:)-T(end-2,:);
    rho(end,:) = p(end,:)./(R.*T(end,:));
    
    % Top boundary: overwrites top right corner point
    u(:,end) = uinf;
    v(:,end) = 0;
    p(:,end) = pinf;
    T(:,end) = Tinf;
    rho(:,end) = p(:,end)./(R.*T(:,end));
    
    % Bottom boundary: overwrites bottom right corner point
    u(:,1) = 0;
    v(:,1) = 0;
    p(:,1) = 2*p(:,2)-p(:,3);
    T(:,1) = Tinf;
    rho(:,1) = p(:,1)./(R.*T(:,1));
    
    % Bottom left corner point
    u(1,1) = 0;
    v(1,1) = 0;
    p(1,1) = pinf;
    T(1,1) = Tinf;
end