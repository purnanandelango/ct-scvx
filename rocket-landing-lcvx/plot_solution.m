% clc
close all
clearvars

load recent_solution

figure
plot3(x(1,:),x(2,:),x(3,:),'-m');
hold on
plot3(xbar(1,:),xbar(2,:),xbar(3,:),'om','MarkerSize',13);
plot3(linspace(0,4000),(430/2600)*linspace(0,4000),tand(prb.gam_gs)*linspace(0,4000),'-k')
view(180,0);
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];
ax.ZLim = [0,2000];
ax.XLim = [0,4000];
zlabel('Altitude [m]');
xlabel('Downrange [m]');

mass_bar = exp(xbar(7,:));
nrm_thrust = misc.compute_vec_norm(u(1:3,:)) .* mass;
nrm_thrust_bar = misc.compute_vec_norm(ubar(1:3,:)) .* mass_bar;
sig = u(4,:) .* mass;

nrm_v = misc.compute_vec_norm(x(4:6,:));
nrm_v_bar = misc.compute_vec_norm(xbar(4:6,:));

figure
subplot(1,2,1)
plot(prb.tau,nrm_thrust_bar,'ob');
hold on
plt1 = plot(tau,nrm_thrust,'-b');
plt2 = plot(tau,sig,'--r');
plot(prb.tau,prb.rho1*ones(1,prb.K),'-k');
plot(prb.tau,prb.rho2*ones(1,prb.K),'-k');
title("Thrust magnitude [N]");
xlabel('$t$ [s]');
legend([plt1,plt2],{'$T(t)$','$\Gamma(t)$'});

subplot(1,2,2)
plot(prb.tau,nrm_v_bar,'or');
hold on
plot(tau,nrm_v,'-r');
plot(prb.tau,prb.vmax*ones(1,prb.K),'-k');
ylim([0,1.1*prb.vmax])
title("Speed [m/s]");
xlabel('$t$ [s]');