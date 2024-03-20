clc
clearvars
close all

s1 = load('recent_solution.mat');

s2 = load('recent_solution_cvx.mat');

prb = s1.prb;

figure('Color',[1,1,1])

plot3(linspace(0,4000),(430/2600)*linspace(0,4000),tand(prb.gam_gs)*linspace(0,4000),'-','LineWidth',3.5,'Color',[1,0.5,0.5])
hold on
plot3(s1.x(1,:),s1.x(2,:),s1.x(3,:),'-k');
plot3(s2.x_sim(1,:),s2.x_sim(2,:),s2.x_sim(3,:),'-','Color',[0,0,0,0.5]);
plot3(s1.xbar(1,:),s1.xbar(2,:),s1.xbar(3,:),'.k');
plot3(s2.x(1,:),s2.x(2,:),s2.x(3,:),'.','Color',[0.5,0.5,0.5]);

view(180,0);
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];
ax.ZLim = [0,2000];
ax.XLim = [0,4000];
zlabel('Altitude [m]');
xlabel('Downrange [m]');
exportgraphics(ax,'traj.pdf','BackgroundColor','none');

mass_bar1 = exp(s1.xbar(7,:));
mass_bar2 = exp(s2.x(7,:));
mass_sim = exp(s2.x_sim(7,:));

nrm_thrust1 = misc.compute_vec_norm(s1.u(1:3,:)) .* s1.mass;
nrm_thrust2 = misc.compute_vec_norm(s2.u_sim(1:3,:)) .* mass_sim;

nrm_thrust_bar1 = misc.compute_vec_norm(s1.ubar(1:3,:)) .* mass_bar1;
nrm_thrust_bar2 = misc.compute_vec_norm(s2.u(1:3,:)) .* mass_bar2;

sig1 = s1.u(4,:) .* s1.mass;
sig2 = s2.u_sim(4,:) .* mass_sim;


figure
subplot(1,2,1)
plot(prb.tau,prb.rho1*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on
plot(prb.tau,prb.rho2*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plot(prb.tau,nrm_thrust_bar1,'.k');
plt1 = plot(s1.tau,nrm_thrust1,'-k');
plt2 = plot(s1.tau,sig1,'--','Color',[0.8,0.8,1]);
% title("Thrust magnitude [N]");
xlabel('$t$ [s]');
legend([plt1,plt2],{'$T(t)$','$\Gamma(t)$'});

subplot(1,2,2)
plot(prb.tau,prb.rho1*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on
plot(prb.tau,prb.rho2*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plot(prb.tau,nrm_thrust_bar2,'.k');
plt1 = plot(s2.tau_sim,nrm_thrust2,'-k');
plt2 = plot(s2.tau_sim,sig2,'--','Color',[0.8,0.8,1]);
% title("Thrust magnitude [N]");
xlabel('$t$ [s]');
% legend([plt1,plt2],{'$T(t)$','$\Gamma(t)$'});