function [Ubar, Ebar, Fbar] = predictor(U, E, F, R, cv, Pr, dx, dy, dt,...
    uinf, pinf, Tinf)

    %% Setup
    % Load functions
    addpath('functions');
    
    % Extract primitive variables from conservative
    [rho,u,v,T,p,~,Et] = cons2prim(U,R,cv);
    
    % Update all necessary physical parameters
    mu = sutherland(T);
    k = mu.*(R+cv)./Pr;

    % Preallocate space for Ubar
    Ubar = U;
    
    %% Compute and assemble flux arrays
    % Compute partial derivatives of primitive variables needed to assemble flux array E

    % Compute the 2D stresses and heat flux
    tau_xx = 2.*mu.*(ddx_bwd(u,dx) - (ddx_bwd(u,dx) + ddy_central(v,dy))./3);
    tau_xy = mu.*(ddy_central(u,dy) + ddx_bwd(v,dx));
    qdot_x = -k.*ddx_bwd(T,dx);
    
    % Update E
    E(1,:,:) = squeeze(U(2,:,:));
    E(2,:,:) = rho.*u.^2 + p - tau_xx;
    E(3,:,:) = rho.*u.*v - tau_xy;
    E(4,:,:) = (Et+p).*u - u.*tau_xx - v.*tau_xy + qdot_x;
    
    % Compute partial derivatives of primitive variables needed to assemble flux array F
    
    % Compute the 2D stresses and heat flux
    tau_yy = 2.*mu.*(ddy_bwd(v,dy) - (ddx_central(u,dx) + ddy_bwd(v,dy))./3);
    tau_xy = mu.*(ddy_bwd(u,dy) + ddx_central(v,dx));
    qdot_y = -k.*ddy_bwd(T,dy);
    
    % Update F
    F(1,:,:) = squeeze(U(3,:,:));
    F(2,:,:) = squeeze(E(3,:,:));
    F(3,:,:) = rho.*v.^2 + p - tau_yy;
    F(4,:,:) = (Et+p).*v - v.*tau_yy - u.*tau_xy + qdot_y;
    
    %% Compute Ubar 
    % Compute Ubar (using forward differences for derivatives of  E and F)
    for ind = 1:4
        Ubar(ind,:,:) = U(ind,:,:) - dt.*... 
            (ddx_fwd(squeeze(E(ind,:,:)),dx) - ...
            ddy_fwd(squeeze(F(ind,:,:)),dy));
    end
    
    % Update primitive variables
    [rho,u,v,T,p,~,~] = cons2prim(Ubar,R,cv);
    
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

    % Update Ubar, update Ebar and Fbar for export
    Ubar = prim2cons(rho,u,v,T,cv);
    Ebar = E;
    Fbar = F;
end