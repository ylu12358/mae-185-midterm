function [mu] = sutherland(T)
%SUTHERLAND Uses Sutherland's law to compute the dynamic viscosity of air.
%   [mu] = sutherland(T)

    % Initialize reference values (using room temp in K for To)
    T0 = 288.15;
    mu0 = 1.735e-5;
    S1 = 110.4;
    mu = mu0.*((T./T0).^1.5).*(T0+S1)./(T+S1);
end