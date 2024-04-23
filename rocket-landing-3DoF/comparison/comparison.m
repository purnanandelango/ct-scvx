clc
clearvars
close all

interpreter = "tex";
% interpreter = "latex";

pxaxwidth = 1110; % [pixels]

save_figures = true;

% Load solutions
s1 = load('solution_eps_1e-5');

s2 = load('solution_node-only');

prb = s1.prb;

% Position
pxaxheight = pxaxwidth/2;
fig = figure('Position',[-1375,64,pxaxwidth*1.05,pxaxheight*1.05]);

plot(linspace(0,4000),tand(prb.gam_gs)*linspace(0,4000),'-','LineWidth',4.5,'Color',[1,0.5,0.5])
hold on
plot(s1.x(1,:),s1.x(3,:),'-k');
thrust = 30;
for k=1:prb.K    
    quiver(s1.xbar(1,k),s1.xbar(3,k),...
           -thrust*s1.ubar(1,k),-thrust*s1.ubar(3,k),...
           'Color',[213,114,14]/255,'LineWidth',5.5);
    quiver(s1.xbar(1,k),s1.xbar(3,k),...
           -0.6*thrust*s1.ubar(1,k),-0.6*thrust*s1.ubar(3,k),...
           'Color',[255,215,0]/255,'LineWidth',2.5);
end
plot(s1.xbar(1,:),s1.xbar(3,:),'.k');
plot(s2.x_sim(1,:),s2.x_sim(3,:),'-','Color',[0,0,1,0.5]);
plot(s2.x(1,:),s2.x(3,:),'.','Color',[0.5,0.5,1]);
ylim([0,2000]);
xlim([0,4000]);
xticks([500,1500,2500,3500]);
ylabel('Altitude [m]');
xlabel('Downrange [m]');
% title('Position');

ax = gca;
% ax.PlotBoxAspectRatio = [1,1,1];
% ax.DataAspectRatio = [1,1,1];
ax.Box = 'off';
ax.Units = "pixels";
ax.OuterPosition = [0, 0,  pxaxwidth, pxaxheight];

%%% Magnified inset
annotation(fig,"rectangle",[0.5713    0.2367    0.0913    0.0926],'LineWidth',1);
% annotation(fig,"arrow",'Position',[0.6661    0.3362    0.0506    0.2144],...
%            'LineWidth',2,'HeadStyle','plain','HeadLength',9,'HeadWidth',6,'Color',[0.2,0.2,0.2]);

axes('Position',[0.6011    0.5851    0.2500    0.2500])
plot(linspace(0,4000),tand(prb.gam_gs)*linspace(0,4000),'-','LineWidth',4.5,'Color',[1,0.5,0.5])
hold on
plot(s1.x(1,:),s1.x(3,:),'-k');
plot(s1.xbar(1,:),s1.xbar(3,:),'.k');
plot(s2.x_sim(1,:),s2.x_sim(3,:),'-','Color',[0,0,1,0.5]);
plot(s2.x(1,:),s2.x(3,:),'.','Color',[0.5,0.5,1]);
ylim([157,457]);
xlim([2420,2920]);
xticks([2500,2800]);
yticks([200,400]);

ax = gca;
% ax.PlotBoxAspectRatio = [1,1,1];
% ax.DataAspectRatio = [1,1,1];
ax.FontSize = 30;
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
pxaxheight = pxaxwidth/2;
fig = figure('Position',[-1375,64,pxaxwidth*1.05,pxaxheight*1.05]);

% sgt = sgtitle('Thrust');
% sgt.FontSize = 42;

subplot(1,2,1)
plot(prb.tau,prb.rho1*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
hold on
plot(prb.tau,prb.rho2*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
plot(s1.tau,nrm_thrust1,'-k');
plot(s1.tau,sig1,'--','Color',[0.5,0.5,1]);
plot(prb.tau,nrm_thrust_bar1,'.k');
ylabel('[N]');
xlim([0,s1.tvec(end)]);
ylim([0.9*prb.rho1,1.05*prb.rho2]);
xticks([0,25,50,75]);
yticks([3000,6000,9000,12000]);
if interpreter == "latex"
    xlabel('$t$ [s]');
elseif interpreter == "tex"
    xlabel('{\itt} [s]');
end

ax = gca;
ax.Box = "off";
ax.Units = "pixels";
ax.OuterPosition = [0, 0,  pxaxwidth/2+30, pxaxheight];

subplot(1,2,2)
plot(prb.tau,prb.rho1*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
hold on
plot(prb.tau,prb.rho2*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
plt1 = plot(s2.tau_sim,nrm_thrust2,'-k');
plt2 = plot(s2.tau_sim,sig2,'--','Color',[0.5,0.5,1]);
plot(prb.tau,nrm_thrust_bar2,'.k');
if interpreter == "latex"
    leg = legend([plt1,plt2],{'$\|T(t)\|$','$\sigma(t)$'});
    xlabel('$t$ [s]');
elseif interpreter == "tex"
    leg = legend([plt1,plt2],{char(8214)+"{\itT{\rm(}t{\rm)}}"+char(8214),'\sigma({\itt})'});    
    xlabel('{\it t} [s]')
end
set(leg,'LineWidth',1,'Position',[0.7377    0.7162    0.1162    0.1517],'FontSize',30);
xlim([0,s1.tvec(end)]);
xticks([0,25,50,75]);
ylim([0.9*prb.rho1,1.05*prb.rho2]);

ax = gca;
ax.Box = "off";
ax.Units = "pixels";
ax.OuterPosition = [pxaxwidth/2-30, 0, pxaxwidth/2+30, pxaxheight];
ax.YTickLabel = {};

%%% Magnified inset
annotation(fig,"rectangle",[0.6900    0.1935    0.1057    0.0361],'LineWidth',1);
% annotation(fig,"arrow",'Position',[0.7554    0.2316   -0.0189    0.1544],...
%            'LineWidth',2,'HeadStyle','plain','HeadLength',9,'HeadWidth',6,'Color',[0.2,0.2,0.2]);

axes('Position',[0.6117    0.4148    0.1698    0.0982]);
plot(prb.tau,prb.rho1*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
hold on
plot(prb.tau,prb.rho2*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
plt1 = plot(s2.tau_sim,nrm_thrust2,'-k');
plt2 = plot(s2.tau_sim,sig2,'--','Color',[0.5,0.5,1]);
plot(prb.tau,nrm_thrust_bar2,'.k');
ylim([4850,5050]);
xlim([41,67]);
xticks([45,55,65]);
yticks([4850,5050]);

ax = gca;
ax.FontSize = 30;
set(ax,'LineWidth',1);
%%%

if save_figures
    exportgraphics(fig,'thrust.pdf','ContentType','vector');
    savefig(fig,'thrust.fig');
end