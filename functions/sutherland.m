% Uses Sutherland's law to compute the dynamic viscosity of air
function [mu] = sutherland(T)
    % Initialize reference values (using room temp in K for To)
    To = 293.15;
    muo = 1.735E-5;
    fprintf('Assuming %d K for reference temperature and %f N•s/m^2 for reference dynamic viscosity\n', To, muo);
    S1 = 110.4;
    mu = muo.*((T./To).^1.5).*(To+S1)./(T+S1);
end