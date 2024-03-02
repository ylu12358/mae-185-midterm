function S = schlieren(rho,dx,dy)
%SCHLIEREN Generates schlieren-like image
%   S = schlieren(rho,dx,dy)

    %% Parameters
    beta = 0.8;
    kappa = 10;
    
    %% Calculated gradient
    rhograd = abs(sqrt((ddx_central(rho,dx)).^2 + (ddy_central(rho,dy)).^2));
    S = beta*exp(-kappa*rhograd/max(rhograd,[],"all"));
end