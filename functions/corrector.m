function [U] = corrector(U, Ubar, Ebar, Fbar, Pr, dx, dy, dt, R, cv, cp,...
    uinf, pinf, Tinf)

    %% Setup

    % Extract primitive variables at predictor level from Ubar
    [rho,u,v,T,p,~,Et] = cons2prim(Ubar,R,cv);

    % Update all necessary physical parameters
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
    U(1,:,:) = 0.5.*(squeeze(U(1,:,:)) + squeeze(Ubar(1,:,:)) ...
        - dt.*ddx_bwd(squeeze(Ebar(1,:,:)),dx) ...
        - dt.*ddy_bwd(squeeze(Fbar(1,:,:)),dy)); % 1st conservative var (rho)
   
    U(2,:,:) = 0.5.*(squeeze(U(2,:,:)) + squeeze(Ubar(2,:,:)) ...
        - dt.*ddx_bwd(squeeze(Ebar(2,:,:)),dx) ...
        - dt.*ddy_bwd(squeeze(Fbar(2,:,:)),dy)); % 2nd conservative var (rho*u)
    
    U(3,:,:) = 0.5.*(squeeze(U(3,:,:)) + squeeze(Ubar(3,:,:)) ...
        - dt.*ddx_bwd(squeeze(Ebar(3,:,:)),dx) ...
        - dt.*ddy_bwd(squeeze(Fbar(3,:,:)),dy)); % 3rd conservative var (rho*v)

    U(4,:,:) = 0.5.*(squeeze(U(4,:,:)) + squeeze(Ubar(4,:,:)) ...
        - dt.*ddx_bwd(squeeze(Ebar(4,:,:)),dx) ...
        - dt.*ddy_bwd(squeeze(Fbar(4,:,:)),dy)); % 4th conservative var (Et)

    %% Update primitive variables
    [rho,u,v,T,p,~,~] = cons2prim(U,R,cv);

    %% Enforce BCs on primitive variables 
    [u,v,~,T,rho] = enforceBCs(u,v,p,T,rho,pinf,Tinf,uinf,R);


    %% Update U

    U = prim2cons(rho,u,v,T,cv);

end