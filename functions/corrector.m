function [U] = corrector(U, Ubar, Ebar, Fbar, Pr, dx, dy, dt, R, cv, cp,...
    uinf, pinf, Tinf)
    
    % Extract primitive variables at predictor level from Ubar
    [rho,u,v,T,p,~,Et] = cons2prim(Ubar,R,cv);

    %% Update all necessary physical parameters
    mu = sutherland(T);
    k = (cp/Pr)*mu;

    %% Compute partial derivatives of primitive variables needed to assemble Ebar and Fbar
    
    % Reminder: corrector applied bwd FD in x on Ebar, so use fwd FD in x 
    % and central FD in y on internal derivatives

    % Compute the 2D stresses and heat flux for Ebar 
    tau_xx = 2.*mu.*(ddx_fwd(u,dx) - (ddx_fwd(u,dx) + ddy_central(v,dy))./3);
    tau_xy = mu.*(ddy_central(u,dy) + ddx_fwd(v,dx));
    qdot_x = -k.*ddx_fwd(T,dx);
    
    % Update Ebar
    Ebar(1,:,:) = squeeze(Ubar(2,:,:));
    Ebar(2,:,:) = rho.*u.^2 + p - tau_xx;
    Ebar(3,:,:) = rho.*u.*v - tau_xy;
    Ebar(4,:,:) = (Et+p).*u - u.*tau_xx - v.*tau_xy + qdot_x;

    % Reminder: corrector applied bwd FD in y on Fbar, so use fwd FD in y 
    % and central FD in x on internal derivatives

    % Compute the 2D stresses and heat flux for Fbar 
    tau_yy = 2.*mu.*(ddy_fwd(v,dy) - (ddx_central(u,dx) + ddy_fwd(v,dy))./3);
    tau_xy = mu.*(ddy_fwd(u,dy) + ddx_central(v,dx));
    qdot_y = -k.*ddy_fwd(T,dy);
    
    % Update Fbar
    Fbar(1,:,:) = squeeze(Ubar(3,:,:));
    Fbar(2,:,:) = squeeze(Ebar(3,:,:));
    Fbar(3,:,:) = rho.*v.^2 + p - tau_yy;
    Fbar(4,:,:) = (Et+p).*v - v.*tau_yy - u.*tau_xy + qdot_y;
    
    %% Compute U using backward FDs for derivatives of Ebar and Fbar
    U(1,:,:) = 0.5.*(U(1,:,:) + Ubar(1,:,:) ...
        - dt.*ddx_bwd(squeeze(Ebar(1,:,:)),dx) ...
        - dt.*ddy_bwd(squeeze(Fbar(1,:,:)),dy)); % 1st conservative var (rho)
   
    U(2,:,:) = 0.5.*(U(2,:,:) + Ubar(2,:,:) ...
        - dt.*ddx_bwd(squeeze(Ebar(2,:,:)),dx) ...
        - dt.*ddy_bwd(squeeze(Fbar(2,:,:)),dy)); % 2nd conservative var (rho*u)
    
    U(3,:,:) = 0.5.*(U(3,:,:) + Ubar(3,:,:) ...
        - dt.*ddx_bwd(squeeze(Ebar(3,:,:)),dx) ...
        - dt.*ddy_bwd(squeeze(Fbar(3,:,:)),dy)); % 3rd conservative var (rho*v)

    U(4,:,:) = 0.5.*(U(4,:,:) + Ubar(4,:,:) ...
        - dt.*ddx_bwd(squeeze(Ebar(4,:,:)),dx) ...
        - dt.*ddy_bwd(squeeze(Fbar(4,:,:)),dy)); % 4th conservative var (Et)

    %% Update primitive variables
    [rho,u,v,T,p,~,~] = cons2prim(U,R,cv);

    %% Enforce BCs on primitive variables 
    % Also enforce BCs on rho due to dependence on p and T 

    % Inlet (left)
    u(1,:) = uinf;
    v(1,:) = 0;
    T(1,:) = Tinf;
    p(1,:) = pinf;
    rho(1,:) = p(1,:)./(R.*T(1,:));

    % Outlet (right)
    u(end,:) = 2*u(end-1,:)-u(end-2,:);
    v(end,:) = 2*v(end-1,:)-v(end-2,:);
    T(end,:) = 2*T(end-1,:)-T(end-2,:);
    p(end,:) = 2*p(end-1,:)-p(end-2,:);
    rho(end,:) = p(end,:)./(R.*T(end,:));
    
    % Wall (bottom)
    u(:,1) = 0;
    v(:,1) = 0;
    T(:,1) = Tinf;
    p(:,1) = 2*p(:,2)-p(:,3); % 2nd order extrapolation
    rho(:,1) = p(:,1)./(R.*T(:,1));
    
    % Far-field (top)
    u(:,end) = uinf;
    v(:,end) = 0;
    T(:,end) = Tinf;
    p(:,end) = pinf;
    rho(:,end) = p(:,end)./(R.*T(:,end));

    % Leading edge (bottom left point) 
    u(1,1) = 0;

    %% Update U

    U = prim2cons(rho,u,v,T,cv);

end