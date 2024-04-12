clc
clearvars
close all

interpreter = "tex";
% interpreter = "latex";

save_figures = true;

% Load solutions
s1 = load('../ctcs/recent_solution.mat');

s2 = load('../node-only-cnstr/recent_solution_cvx.mat');

prb = s1.prb;

% Position
fig = figure('Position',[214,19,407,407]);

plot(linspace(0,4000),tand(prb.gam_gs)*linspace(0,4000),'-','LineWidth',3.5,'Color',[1,0.5,0.5])
hold on
plot(s1.x(1,:),s1.x(3,:),'-k');
plot(s1.xbar(1,:),s1.xbar(3,:),'.k');
plot(s2.x_sim(1,:),s2.x_sim(3,:),'--','Color',[0.7,0.7,0.7]);
plot(s2.x(1,:),s2.x(3,:),'.','Color',[0.7,0.7,0.7]);
ylim([0,2000]);
xlim([0,4000]);
ylabel('Altitude [m]');
xlabel('Downrange [m]');
title('Position');

ax = gca;
ax.FontSize = 20;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];

%%% Magnified inset
annotation(fig,"rectangle",[0.6151, 0.3661, 0.0998, 0.0584],'LineWidth',1);
annotation(fig,"arrow",'Position',[0.7199,0.43,0.0541,0.0934],...
           'LineWidth',1,'HeadStyle','plain','HeadLength',7,'HeadWidth',4);

axes('Position',[0.6552,0.4976,0.25,0.25])
plot(linspace(0,4000),tand(prb.gam_gs)*linspace(0,4000),'-','LineWidth',3.5,'Color',[1,0.5,0.5])
hold on
plot(s1.x(1,:),s1.x(3,:),'-k');
plot(s1.xbar(1,:),s1.xbar(3,:),'.k');
plot(s2.x_sim(1,:),s2.x_sim(3,:),'--','Color',[0.7,0.7,0.7]);
plot(s2.x(1,:),s2.x(3,:),'.','Color',[0.7,0.7,0.7]);
ylim([157,457]);
xlim([2420,2920]);
xticks([2500,2800]);
yticks([200,400]);

ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];
ax.FontSize = 16;
set(ax,'LineWidth',1);
%%%

if save_figures 
    exportgraphics(fig,'traj.pdf','ContentType','vector');
    savefig(fig,'traj.fig');
end

mass_bar1 = exp(s1.xbar(7,:));
mass_bar2 = exp(s2.x(7,:));
mass_sim = exp(s2.x_sim(7,:));

nrm_thrust1 = misc.compute_vec_norm(s1.u(1:3,:)) .* s1.mass;
nrm_thrust2 = misc.compute_vec_norm(s2.u_sim(1:3,:)) .* mass_sim;

nrm_thrust_bar1 = misc.compute_vec_norm(s1.ubar(1:3,:)) .* mass_bar1;
nrm_thrust_bar2 = misc.compute_vec_norm(s2.u(1:3,:)) .* mass_bar2;

sig1 = s1.u(4,:) .* s1.mass;
sig2 = s2.u_sim(4,:) .* mass_sim;

% Thrust
fig = figure('Position',[697,439,407,407]);

sgt = sgtitle('Thrust');
sgt.FontSize = 20;

subplot(1,2,1)
plot(prb.tau,prb.rho1*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on
plot(prb.tau,prb.rho2*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plt1 = plot(s1.tau,nrm_thrust1,'-k');
plt2 = plot(s1.tau,sig1,'--','Color',[0.7,0.7,0.7]);
plot(prb.tau,nrm_thrust_bar1,'.k');
ylabel('[N]');
xlim([0,s1.tvec(end)]);
ylim([0.9*prb.rho1,1.05*prb.rho2]);
yticks([3000,6000,9000,12000]);
if interpreter == "latex"
    leg = legend([plt1,plt2],{'$\|T(t)\|$','$\sigma(t)$'});
    xlabel('$t$ [s]');
elseif interpreter == "tex"
    leg = legend([plt1,plt2],{char(8214)+"{\itT{\rm(}t{\rm)}}"+char(8214),'\sigma({\itt})'});
    xlabel('{\itt} [s]');
end
set(leg,'LineWidth',1,'Position',[0.2357,0.5995,0.2334,0.1143],'FontSize',20);

ax = gca;
ax.FontSize = 20;
ax.Position = [0.2047,0.1317,0.3454,0.7915];

subplot(1,2,2)
plot(prb.tau,prb.rho1*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on
plot(prb.tau,prb.rho2*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plot(s2.tau_sim,nrm_thrust2,'-k');
plot(s2.tau_sim,sig2,'--','Color',[0.7,0.7,0.7]);
plot(prb.tau,nrm_thrust_bar2,'.k');
if interpreter == "latex"
    xlabel('$t$ [s]');
elseif interpreter == "tex"
    xlabel('{\it t} [s]')
end
xlim([0,s1.tvec(end)]);
ylim([0.9*prb.rho1,1.05*prb.rho2]);

ax = gca;
ax.Position = [0.5703,0.1312,0.3347,0.7938];
ax.FontSize = 20;
ax.YTickLabel = {};

%%% Magnified inset
annotation(fig,"rectangle",[0.7346,0.1506,0.1057,0.0361],'LineWidth',1);
annotation(fig,"arrow",'Position',[0.808,0.1966,0.0323,0.1081],...
           'LineWidth',1,'HeadStyle','plain','HeadLength',7,'HeadWidth',4);

axes('Position',[0.8203,0.317,0.1698,0.0982]);
plot(prb.tau,prb.rho1*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on
plot(prb.tau,prb.rho2*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plt1 = plot(s2.tau_sim,nrm_thrust2,'-k');
plt2 = plot(s2.tau_sim,sig2,'--','Color',[0.7,0.7,0.7]);
plot(prb.tau,nrm_thrust_bar2,'.k');
ylim([4850,5050]);
xlim([41,67]);
yticks([4850,5050]);

ax = gca;
ax.FontSize = 16;
set(ax,'LineWidth',1);
%%%

if save_figures
    exportgraphics(fig,'thrust.pdf','ContentType','vector');
    savefig(fig,'thrust.fig');
end