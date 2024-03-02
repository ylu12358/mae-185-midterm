function [U] = corrector(U,Ubar,Ebar,Fbar,R,cp,cv,Pr,dx,dy,dt,uinf,pinf,Tinf,bc)
%CORRECTOR Executes corrector step of MacCormack method.
%   [U] = corrector(U,Ubar,Ebar,Fbar,R,cp,cv,Pr,dx,dy,dt,uinf,pinf,Tinf)

    %% Setup

    % Extract primitive variables at predictor level from Ubar
    [rho,u,v,T,p,~,Et] = cons2prim(Ubar,R,cv);

    % Update all necessary physical parameters
    mu = sutherland(T);
    k = (cp/Pr)*mu;

    % Create temperature gradient array 
    dTdy = ddy_fwd(T,dy);
    %% Compute partial derivatives of primitive variables needed to assemble Ebar and Fbar

    % set default value for 'bc'
    if nargin<3
        
        bc = 'isothermal'; 
   
    elseif bc == "adiabatic" 
        
        % Enforce adiabatic wall BC on temperature gradient array 
        dTdy(:,1) = 0;
    
    end

    % Compute the 2D stresses and heat flux for Ebar
    tau_xx = 2*mu.*(ddx_fwd(u,dx) - (ddx_fwd(u,dx) + ddy_central(v,dy))/3);
    tau_xy = mu.*(ddy_central(u,dy) + ddx_fwd(v,dx));
    qdot_x = -k.*ddx_fwd(T,dx);
    
    % Update Ebar
    Ebar(1,:,:) = squeeze(Ubar(2,:,:));
    Ebar(2,:,:) = rho.*u.^2 + p - tau_xx;
    Ebar(3,:,:) = rho.*u.*v - tau_xy;
    Ebar(4,:,:) = (Et+p).*u - u.*tau_xx - v.*tau_xy + qdot_x;

    % Compute the 2D stresses and heat flux for Fbar 
    tau_yy = 2*mu.*(ddy_fwd(v,dy) - (ddx_central(u,dx) + ddy_fwd(v,dy))/3);
    tau_xy = mu.*(ddy_fwd(u,dy) + ddx_central(v,dx));
    qdot_y = -k.*dTdy;
    
    % Update Fbar
    Fbar(1,:,:) = squeeze(Ubar(3,:,:));
    Fbar(2,:,:) = squeeze(Ebar(3,:,:));
    Fbar(3,:,:) = rho.*v.^2 + p - tau_yy;
    Fbar(4,:,:) = (Et+p).*v - v.*tau_yy - u.*tau_xy + qdot_y;
    
    %% Compute U using backward FDs for derivatives of Ebar and Fbar
    for ind = 1:4
        U(ind,:,:) = 0.5*( ...
            squeeze(U(ind,:,:)) ...
            + squeeze(Ubar(ind,:,:)) ...
            - dt.*ddx_bwd(squeeze(Ebar(ind,:,:)),dx) ...
            - dt.*ddy_bwd(squeeze(Fbar(ind,:,:)),dy) ...
            );
    end
   
    %% Update primitive variables
    [rho,u,v,T,p,~,~] = cons2prim(U,R,cv);

    %% Enforce BCs on primitive variables 
    [u,v,~,T,rho] = enforceBC(u,v,p,T,rho,pinf,Tinf,uinf,R);


    %% Update U
    U = prim2cons(rho,u,v,T,cv);

end