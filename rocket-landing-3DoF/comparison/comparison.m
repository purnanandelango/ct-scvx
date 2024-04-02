clc
clearvars
close all

interpreter = "tex";
% interpreter = "latex";

s1 = load('../ctcs/recent_solution.mat');

s2 = load('../node-only-cnstr/recent_solution_cvx.mat');

prb = s1.prb;

fig = figure('Position',[214,19,612,340]);

plot3(linspace(0,4000),(430/2600)*linspace(0,4000),tand(prb.gam_gs)*linspace(0,4000),'-','LineWidth',3.5,'Color',[1,0.5,0.5])
hold on
plot3(s1.x(1,:),s1.x(2,:),s1.x(3,:),'-k');
plot3(s1.xbar(1,:),s1.xbar(2,:),s1.xbar(3,:),'.k');
plot3(s2.x_sim(1,:),s2.x_sim(2,:),s2.x_sim(3,:),'-','Color',[0.6,0.6,0.6]);
plot3(s2.x(1,:),s2.x(2,:),s2.x(3,:),'.','Color',[0.6,0.6,0.6]);

view(180,0);
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];
ax.ZLim = [0,2000];
ax.XLim = [0,4000];
zlabel('Altitude [m]');
xlabel('Downrange [m]');
ax.Box = 'off';
exportgraphics(fig,'traj.pdf','ContentType','vector');

mass_bar1 = exp(s1.xbar(7,:));
mass_bar2 = exp(s2.x(7,:));
mass_sim = exp(s2.x_sim(7,:));

nrm_thrust1 = misc.compute_vec_norm(s1.u(1:3,:)) .* s1.mass;
nrm_thrust2 = misc.compute_vec_norm(s2.u_sim(1:3,:)) .* mass_sim;

nrm_thrust_bar1 = misc.compute_vec_norm(s1.ubar(1:3,:)) .* mass_bar1;
nrm_thrust_bar2 = misc.compute_vec_norm(s2.u(1:3,:)) .* mass_bar2;

sig1 = s1.u(4,:) .* s1.mass;
sig2 = s2.u_sim(4,:) .* mass_sim;

fig = figure('Position',[697,439,612,577]);
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
if interpreter == "latex"
    legend([plt1,plt2],{'$\|T(t)\|$','$\sigma(t)$'});
    xlabel('$t$ [s]');
elseif interpreter == "tex"
    legend([plt1,plt2],{char(8214)+"{\itT{\rm(}t{\rm)}}"+char(8214),'\sigma({\itt})'});
    xlabel('{\itt} [s]');
end
ax = gca;
ax.Box = 'off';

subplot(1,2,2)
plot(prb.tau,prb.rho1*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on
plot(prb.tau,prb.rho2*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plt1 = plot(s2.tau_sim,nrm_thrust2,'-k');
plt2 = plot(s2.tau_sim,sig2,'--','Color',[0.7,0.7,0.7]);
plot(prb.tau,nrm_thrust_bar2,'.k');
% ylabel('[N]');
if interpreter == "latex"
    xlabel('$t$ [s]');
elseif interpreter == "tex"
    xlabel('{\it t} [s]')
end
xlim([0,s1.tvec(end)]);
ylim([0.9*prb.rho1,1.05*prb.rho2]);
ax = gca;
ax.Box = 'off';
plt.inset.MagInset(fig,ax,[42,67,4800,5100],[20,65,7000,8250],{'SW','SW';'NE','NE'});
ax = gca;
ax.XTickLabel = {};
ax.YTickLabel = {};
% exportgraphics(fig,'thrust.pdf','ContentType','vector');