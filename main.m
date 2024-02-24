%% clean slate
clc;
clear all;
close all;

%% load functions
addpath('functions')

%% constants
cp = 1005;
cv = 718;
R = 287;
gamma = 1.4;
mu0 = 1.735*10^(-5);
S1 = 110.4;
Pr = 0.71;

L = 10^(-5);
H = 8*10^(-6);

%% initialize
M = 4;

nx = 75;
ny = 80;
xx = ndgrid(1:nx,1:ny);
empty = zeros(nx,ny);
u = empty;
v = empty;

E = empty;
F = empty;
Ubar = empty;

p0 = 101300;
p = p0*ones(nx,ny);

T0 = 288.15;
T = T0*ones(nx,ny);

rho0 = 1.225;
rho = rho0*ones(nx,ny);

%% calorically perfect gas equations
a = sqrt(gamma*R*T0);
uinf = M*a;

%% dx, dy, dt
dx = L/nx;
dy = H/ny;
dt = 2.35*10^(-11);
