%% clean slate
clc;
clear all;
close all;

%% load data
adiabatic = load("Tadiabatic.mat");
isothermal = load("Tisothermal.mat");

%% plot temperature at wall
f1 = figure(1);
f1.Position = [100,100,1200,800];

plot(adiabatic.xx(:,1),adiabatic.T(:,1),'LineWidth',2);
hold on; grid on; box on;
plot(isothermal.xx(:,1),isothermal.T(:,1),'LineWidth',2);

title('Wall temperature for adiabatic and isothermal boundary conditions');
xlabel('x'); ylabel('T [K]');
legend('adiabatic','isothermal','Location','east');
set(gca,'FontSize',14);