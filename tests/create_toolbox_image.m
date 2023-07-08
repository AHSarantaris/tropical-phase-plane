clc
clear
close all
a = 0.25;
F = [1 a-1 -1 0;-1 -1 0 0; -1 -1 0 2];
G = [-1 0 0 0; 1 0 1 -1; 1 0 1 1];
ulim = [-2 1];
vlim = [-1.5 1.5];
% Trajectory parameters
w0 = [ulim(1); 0.25];
eps = 0.08;
tspan = [0 exp(1/eps)];
% Plot phase plane with trajectory
f = figure("Visible","on","Units","normalized","OuterPosition",[0.2 0.2 .4 .6]);
TropicalPhasePlane(F,G,ulim,vlim)
% xticks([])
% yticks([])
hold on
[T,W] = ode15s(@tropicalODE,tspan,w0,[],F,G,eps);
plot(W(:,1),W(:,2),'Color',[0 0.5 0],'LineWidth',1.5)

exportgraphics(f,'toolbox_image.png','Resolution',300) 