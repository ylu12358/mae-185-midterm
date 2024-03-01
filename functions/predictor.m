function [Ubar, Ebar, Fbar] = predictor(U, E, F, R, cv, Pr, dx, dy, dt,...
    uinf, pinf, Tinf)

    %% Setup
    
    % Extract primitive variables from conservative
    [rho,u,v,T,p,~,Et] = cons2prim(U,R,cv);
    
    % Update all necessary physical parameters
    mu = sutherland(T);
    k = mu.*(R+cv)/Pr;

    % Preallocate space for Ubar
    Ubar = U;
    
    %% Compute and assemble flux arrays
    % Compute partial derivatives of primitive variables needed to assemble flux array E

    % Compute the 2D stresses and heat flux
    tau_xx = 2*mu.*(ddx_bwd(u,dx) - (ddx_bwd(u,dx) + ddy_central(v,dy))./3);
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
        Ubar(ind,:,:) = squeeze(U(ind,:,:)) - dt*... 
            (ddx_fwd(squeeze(E(ind,:,:)),dx) - ...
            ddy_fwd(squeeze(F(ind,:,:)),dy));
    end
    
    % Update primitive variables
    [rho,u,v,T,p,~,~] = cons2prim(Ubar,R,cv);
    
    % Enforce BCs on primitive variables
    [u,v,~,T,rho] = enforceBCs(u,v,p,T,rho,pinf,Tinf,uinf,R);


    % Update Ubar, update Ebar and Fbar for export
    Ubar = prim2cons(rho,u,v,T,cv);
    Ebar = E;
    Fbar = F;
end