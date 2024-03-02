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
gamma = cp/cv;
Pr = 0.71;

% physical parameters
L = 10^(-5);
H = 8*10^(-6);

M = 4;

% simulation parameters
nx = 75;
ny = 80;
nt = 1500;
figureskipped = 50; % only plots iterations at increments of this value
t = 0;

plotschlieren = false;
plotadiabaticwall = false;
bc = "isothermal";

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

%% preallocate other arrays
E = zeros(4,nx,ny);
F = zeros(4,nx,ny);
Ubar = zeros(4,nx,ny);

%% compute initial physical parameters
mu = sutherland(T);

%% dx, dy, dt, safety factor sf
dx = diff(xx);
dx = dx(1);
dy = diff(yy');
dy = dy(1);
sf = 1.2;
dt = 2.35*10^(-11)/sf;

%% initialize conservative variables
U = prim2cons(rho,u,v,T,cv);

%% initialize graph window and parameters
f1 = figure(1);
f1.Position = [100,100,1400,1000];

varsplot = zeros(6,nx,ny);

cblabels = {'\rho [kg/m^3]',...
    'u [m/s]',...
    'v [m/s]',...
    'e [J/kg = m^3/s^3]',...
    'p [Pa = kg/m/s^2]',...
    'T [K]'};

subtitles = {'density',...
    'velocity - x',...
    'velocity - y',...
    'spec. int. energy',...
    'pressure',...
    'temperature'};

if plotschlieren == true
    cblabels{1} = ' ';
    subtitles{1} = 'schlieren';
end

if plotadiabaticwall == true
    
    cblabels{5,6} = ' ';
    subtitles{5} = 'pressure; adiabatic wall';
    subtitles{6} = 'temperature; adiabatic wall';
end

convergencevar = rho(40,20);
convergencet = t;

%% time loop
for iter = 1:nt
    % increment t
    t = t + dt;

    % predictor and corrector step
    [Ubar, Ebar, Fbar] = predictor(U,E,F,R,cp,cv,Pr,dx,dy,dt,uinf,pinf,Tinf);
    U = corrector(U,Ubar,Ebar,Fbar,R,cp,cv,Pr,dx,dy,dt,uinf,pinf,Tinf);

    % update figure every few iterations
    if mod(iter,figureskipped) == 0 || iter == nt
            
        % extract conservative variables for plot
        [rho,u,v,T,p,e,~] = cons2prim(U,R,cv);
    
        varsplot(1,:,:) = rho;
        if plotschlieren == true
            varsplot(1,:,:) = schlieren(rho,dx,dy);
        end
        varsplot(2,:,:) = u;
        varsplot(3,:,:) = v;
        varsplot(4,:,:) = e;
        varsplot(5,:,:) = p;
        varsplot(6,:,:) = T;
        if plotadiabaticwall == true
            varsplot(5,:,:) = p./pinf;
            varsplot(6,:,:) = T./Tinf;
        end
    
        % overall title
        sgtitle(['MacCormack for Compressible Navier-Stokes: t = ' num2str(t) ', n = ' num2str(iter)]);
    
        % pcolor
        for j = 1:size(varsplot,1)
            ax = subplot(3,3,j);
            pcolor(xx,yy,squeeze(varsplot(j,:,:)));
            shading interp;
            cb(j) = colorbar;
            ylabel(cb(j),cblabels(j));
            if j == 1 && plotschlieren == true % schlieren
                colormap(ax,"gray");
                clim([0 1]);
            elseif (j == 5 || j == 6) && plotadiabaticwall == true % adiabatic wall
                
            else
                colormap(ax,"turbo");
            end
            axis equal tight;
            title(subtitles(j));
            xlabel('x'); ylabel('y');
            set(gca,'FontSize',12);
        end

        % convergence
        convergencevar = [convergencevar rho(40,20)];
        convergencet = [convergencet  t];
        subplot(3,3,[7 8 9]);
        plot(convergencet,convergencevar,'r-','LineWidth',2);
        hold on; box on; grid on;
        xlim([0 t]);
        title(['convergence of density at x = ' num2str(40*dx) ...
            ', y = ' num2str(20*dy)]);
        xlabel('t'); ylabel('\rho');
        set(gca,'FontSize',14);

        drawnow;
    end
end

%% end timer
toc