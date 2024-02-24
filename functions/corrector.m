function [U] = corrector(U, Ubar, Ebar, Fbar, Pr, dx, dy, dt, R, cv, cp)
    
    %% Setup
    % Extract primitive variables at predictor level from Ubar
    [rho,u,v,T,p,e,Et] = cons2prim(Ubar,R,cv);

    % Update all necessary physical parameters
    mu = sutherland(T);
    k = (cp/Pr)*mu;

    %% Compute and assemble flux arrays
    
    U = 0.5*(U+Ubar - dt*ddx_bwd(Ebar,dx) - dt*ddy_bwd(Fbar,dy));

    % Compute partial derivatives of primitive variables needed to assemble
    % flux array Ebar

    % Compute the 2D stresses and heat flux
    tau_xx = 2.*mu.*(ddx_fwd(u,dx) - (ddx_fwd(u,dx) + ddy_central(v,dy))./3);
    tau_xy = mu.*(ddy_central(u,dy) + ddx_fwd(v,dx));
    qdot_x = k.*ddx_fwd(T,dx);
    
    % Update Ebar
    Ebar(1,:,:) = squeeze(Ubar(2,:,:));
    Ebar(2,:,:) = rho.*u.^2 + p - tau_xx;
    Ebar(3,:,:) = rho.*u.*v - tau_xy;
    Ebar(4,:,:) = (Et+p).*u - u.*tau_xx - v.*tau_xy + qdot_x;
    
    % Compute partial derivatives of primitive variables needed to assemble
    % flux array Fbar
    
    % Compute the 2D stresses and heat flux
    tau_yy = 2.*mu.*(ddy_fwd(v,dy) - (ddx_central(u,dx) + ddy_fwd(v,dy))./3);
    tau_xy = mu.*(ddy_fwd(u,dy) + ddx_central(v,dx));
    qdot_y = k.*ddy_fwd(T,dy);
    
    % Update Fbar
    Fbar(1,:,:) = squeeze(U(3,:,:));
    Fbar(2,:,:) = squeeze(Ebar(3,:,:));
    Fbar(3,:,:) = rho.*v.^2 + p - tau_yy;
    Fbar(4,:,:) = (Et+p).*v - v.*tau_yy - u.*tau_xy + qdot_y;

    


    % enforce BCs on primitive variables 

end