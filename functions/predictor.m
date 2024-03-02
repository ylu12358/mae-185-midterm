function [Ubar, Ebar, Fbar] = predictor(U,E,F,R,cp,cv,Pr,dx,dy,dt,uinf,pinf,Tinf,bc)
%PREDICTOR Executes predictor step of MacCormack method.
%   [Ubar, Ebar, Fbar] = predictor(U,E,F,R,cp,cv,Pr,dx,dy,dt,uinf,pinf,Tinf)

    %% Setup
    
    % Extract primitive variables from conservative
    [rho,u,v,T,p,~,Et] = cons2prim(U,R,cv);
    
    % Update all necessary physical parameters
    mu = sutherland(T);
    k = (cp/Pr)*mu;

    % Preallocate space for Ubar
    Ubar = U;
    
    % Create temperature gradient array 
    dTdy = ddy_bwd(T,dy);

    
    %% Compute partial derivatives of primitive variables needed to assemble E and F
    
    if bc == "adiabatic" 
        
        % Enforce adiabatic wall BC on temperature gradient array 
        dTdy(:,1) = 0;
    
    end

    % Compute the 2D stresses and heat flux for E
    tau_xx = 2*mu.*(ddx_bwd(u,dx) - (ddx_bwd(u,dx) + ddy_central(v,dy))/3);
    tau_xy = mu.*(ddy_central(u,dy) + ddx_bwd(v,dx));
    qdot_x = -k.*ddx_bwd(T,dx);
    
    % Update E
    E(1,:,:) = squeeze(U(2,:,:));
    E(2,:,:) = rho.*u.^2 + p - tau_xx;
    E(3,:,:) = rho.*u.*v - tau_xy;
    E(4,:,:) = (Et+p).*u - u.*tau_xx - v.*tau_xy + qdot_x;
    
    % Compute the 2D stresses and heat flux for F
    tau_yy = 2*mu.*(ddy_bwd(v,dy) - (ddx_central(u,dx) + ddy_bwd(v,dy))/3);
    tau_xy = mu.*(ddy_bwd(u,dy) + ddx_central(v,dx));
    qdot_y = -k.*dTdy;
    
    % Update F
    F(1,:,:) = squeeze(U(3,:,:));
    F(2,:,:) = squeeze(E(3,:,:));
    F(3,:,:) = rho.*v.^2 + p - tau_yy;
    F(4,:,:) = (Et+p).*v - v.*tau_yy - u.*tau_xy + qdot_y;
    
    %% Compute Ubar using forward FDs for derivatives of E and F
    for ind = 1:4
        Ubar(ind,:,:) = squeeze(U(ind,:,:)) ...
            - dt*ddx_fwd(squeeze(E(ind,:,:)),dx) ...
            - dt*ddy_fwd(squeeze(F(ind,:,:)),dy);
    end
    
    %% Update primitive variables
    [rho,u,v,T,p,~,~] = cons2prim(Ubar,R,cv);
    
    %% Enforce BCs on primitive variables
    [u,v,~,T,rho] = enforceBC(u,v,p,T,rho,pinf,Tinf,uinf,R,bc);

    %% Update Ubar, update Ebar and Fbar for export
    Ubar = prim2cons(rho,u,v,T,cv);
    Ebar = E;
    Fbar = F;
    
end