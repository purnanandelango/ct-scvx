% clc
% close all
clearvars

load recent_solution_cvx

figure
plot3(x_sim(1,:),x_sim(2,:),x_sim(3,:),'-m');
hold on
plot3(x(1,:),x(2,:),x(3,:),'om','MarkerSize',13);
plot3(0*ones(1,100),linspace(0,4000),tand(prb.gam_gs)*linspace(0,4000),'-k')
view(-90,-1);
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];
ax.ZLim = [0,2000];
ax.YLim = [0,4000];

mass_sim = exp(x_sim(7,:));
mass = exp(x(7,:));
nrm_thrust_sim = misc.compute_vec_norm(u_sim(1:3,:)) .* mass_sim;
nrm_thrust = misc.compute_vec_norm(u(1:3,:)) .* mass;
sig = u_sim(4,:) .* mass_sim;

nrm_v = misc.compute_vec_norm(x(4:6,:));
nrm_v_sim = misc.compute_vec_norm(x_sim(4:6,:));

figure
subplot(1,2,1)
plot(prb.tau,nrm_thrust,'ob');
hold on
plot(tau_sim,nrm_thrust_sim,'-b');
plot(tau_sim,sig,'--r');
plot(prb.tau,prb.rho1*ones(1,prb.K),'-k');
plot(prb.tau,prb.rho2*ones(1,prb.K),'-k');

subplot(1,2,2)
plot(prb.tau,nrm_v,'or');
hold on
plot(tau_sim,nrm_v_sim,'-r');
plot(prb.tau,prb.vmax*ones(1,prb.K),'-k');
ylim([0,1.1*prb.vmax])


