%% clean slate
clc;
clear all;
close all;

%% begin timer
tic

%% load functions
addpath('functions')

%% define physical and simulation parameters
% air at standard conditions
pinf = 101300;
Tinf = 288.15;
rho0 = 1.225;
R = 287;
cp = 1005;
cv = 718;
gamma = 1.4;
Pr = 0.71;

% physical parameters
L = 10^(-5);
H = 8*10^(-6);

M = 4;

% simulation parameters
nx = 75;
ny = 80;
nt = 1500;
t = 0;

%% initialize grid
[xx,yy] = ndgrid(linspace(0,L,nx),linspace(0,H,ny));

%% initialize primitive variables w/ ICs
% compute u_infinity
a = sqrt(gamma*R*Tinf);
uinf = M*a;

% initialize variables
rho = rho0*ones(nx,ny);

u = uinf*ones(nx,ny);
v = zeros(nx,ny);

T = Tinf*ones(nx,ny);

p = pinf*ones(nx,ny);

%% initialize conservative variables
U = prim2cons(rho,u,v,T,cv);

%% preallocate other arrays
E = zeros(4,nx,ny);
F = zeros(4,nx,ny);
Ubar = zeros(4,nx,ny);

%% compute initial physical parameters
mu = sutherland(T);

%% dx, dy, dt
dx = diff(xx);
dx = dx(1);
dy = diff(yy');
dy = dy(1);
dt = 2.35*10^(-11);

%% boundary conditions
% wall
u(:,1) = 0;
v(:,1) = 0;
T(:,1) = Tinf;
% inlet and far field
u(1,:) = uinf;
u(:,end) = uinf;
v(1,:) = 0;
v(:,end) = 0;
p(1,:) = pinf;
p(:,end) = pinf;
T(1,:) = Tinf;
T(:,end) = Tinf;
% leading edge
u(1,1) = 0;
v(1,1) = 0;
p(1,1) = pinf;
T(1,1) = Tinf;

%% time loop
for i = 1:nt
    % increment t
    t = t + dt;

    % predictor and corrector step
    [Ubar, Ebar, Fbar] = predictor(U,E,F,R,cv,Pr,dx,dy,dt,uinf,pinf,Tinf);
    U = corrector(U,Ubar,Ebar,Fbar,Pr,dx,dy,dt,R,cv,cp,uinf,pinf,Tinf);
    
    % extract conservative variables for plot
    [rho,u,v,T,p,e,~] = cons2prim(U,R,cv);

    varsplot = zeros(6,nx,ny);
    varsplot(1,:,:) = rho;
    varsplot(2,:,:) = u;
    varsplot(3,:,:) = v;
    varsplot(4,:,:) = e;
    varsplot(5,:,:) = p;
    varsplot(6,:,:) = T;

    % draw graph
    figure(1);
    f1.Position = [100,100,3400,600];

    cblabels = {'\rho [kg/m^3]',...
        'u [m/s]',...
        'v [m/s]',...
        'e [J/kg = m^3/s^3]',...
        'p [Pa = kg/m/s^2]',...
        'T [K]'};
    
    subtitles = {'density',...
        'velocity - x',...
        'velocity - y',...
        'specific internal energy',...
        'pressure',...
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

    % update figure every 20 iterations
    if mod(i,20) == 0
        drawnow;
    end
end

%% end timer
toc