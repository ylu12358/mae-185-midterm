% Converts conservative variables into primitive variables
function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)
    % Retrieve whole values of conservative variables from U
    rho = squeeze(U(1,:,:));
    Et = squeeze(U(4,:,:));
    % Arithmetically compute remaining values from given data
    u = squeeze(U(2,:,:))./rho;
    v = squeeze(U(3,:,:))./rho;
    e = Et./rho-(u.^2+v.^2)/2;
    T = e./cv;
    p = rho.*R.*T;
end