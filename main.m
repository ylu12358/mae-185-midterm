%% clean slate
clc;
clear all;
close all;

%% load functions
addpath('functions')

%% physical and simulation parameters
cp = 1005;
cv = 718;
R = 287;
gamma = 1.4;
mu0 = 1.735*10^(-5);
S1 = 110.4;
Pr = 0.71;

L = 10^(-5);
H = 8*10^(-6);

nx = 75;
ny = 80;
nt = 1500;
t = 0;

%% initialize grids using ndgrid
[xx,yy] = ndgrid(linspace(0,L,nx),linspace(0,H,ny));

%% initialize primitive variables
emptyM = zeros(nx,ny);
identityM = ones(nx,ny);

rho0 = 1.225;
rho = rho0*identityM;

u = emptyM;
v = emptyM;

T0 = 288.15;
T = T0*identityM;

p0 = 101300;
p = p0*identityM;

M = 4;

%% preallocate other arrays
emptyM3D = zeros(4,nx,ny);
E = emptyM3D;
F = emptyM3D;
Ubar = emptyM3D;

%% initial physical parameters
a = sqrt(gamma*R*T0);
uinf = M*a;

%mu = sutherland()

%% dx, dy, dt
dx = diff(xx);
dx = dx(1);
dy = diff(yy');
dy = dy(1);
dt = 2.35*10^(-11);

%% time loop
for j = 1:nt
    % U = corrector(predictor(U,E,F,...));

    % enforce BCs on correctors
    u(:,[1 end]) = 0;
    v(:,[1 end]) = 0;    

    % draw graph
    figure(1);
    f1.Position = [100,100,3400,600];

    cblabels = {'\rho [kg/m^3]',
        'u [m/s]',
        'v [m/s]',
        'e [J/kg = m^3/s^3]',
        'p [Pa = kg/m/s^2]',
        'T [K]'};
    
    subtitles = {'density',
        'velocity - x',
        'velocity - y',
        'specific internal energy',
        'pressure',
        'temperature'};

    sgtitle(['MacCormack for Compressible Navier-Stokes: t = ' num2str(t)]);
    for j = 1:size(varsplot,1)
        subplot(2,3,j);
        pcolor(xx,yy,squeeze(varsplot(j,:,:)));
        shading interp;
        cb(j) = colorbar;
        ylabel(cb(j),cblabels(j));
        hold on;
        axis equal tight;
        title(subtitles(j));
        xlabel('x'); ylabel('y');
        set(gca,'FontSize',14);
    end

    % update figure
    drawnow;
    
    % increment t
    t = t + dt;
end