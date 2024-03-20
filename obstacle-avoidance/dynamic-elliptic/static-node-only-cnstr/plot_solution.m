clearvars
close all

load recent_solution

fig = figure('Position',[138,100,1034,437], ...
             'Color',[1,1,1]);

th = linspace(0,2*pi);
for j = 1:prb.n_obs
    pos_obs = prb.q_obs(:,j) + prb.H_obs{j}\[cos(th);sin(th)];
    plot(pos_obs(1,:),pos_obs(2,:),'-r');
    hold on
end

plot(r(1,:),r(2,:),'-k');
plot(xbar(1,:),xbar(2,:),'.k');    
title('Position [m]');

ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];
ax.XLim = [-65,65];
ax.YLim = [-35,35];

figure

Kfine = length(tau);

nrm_T(Kfine) = 0;
nrm_v(Kfine) = 0;
for j = 1:Kfine
    nrm_T(j) = norm(u(1:prb.n,j));
    nrm_v(j) = norm(v(1:prb.n,j));
end
nrm_Tbar(prb.K) = 0;
nrm_vbar(prb.K) = 0;
for j = 1:prb.K
    nrm_Tbar(j) = norm(ubar(1:prb.n,j));
    nrm_vbar(j) = norm(xbar(prb.n+1:2*prb.n,j));
end

subplot(2,2,2)
plot(tvec,nrm_v,'-k');
hold on 
plot(tvecbar,nrm_vbar,'.k');
plot(tvecbar,prb.vmax*ones(1,prb.K),'-r');
title('Speed')
xlabel('$t$');
xlim([0,tvec(end)])
ylim([0,1.1*prb.vmax])

subplot(2,2,3)
plot(tvecbar,prb.Tmin*ones(1,prb.K),'-r');
hold on 
plot(tvecbar,prb.Tmax*ones(1,prb.K),'-r');
plot(tvec,nrm_T,'-k');
plot(tvecbar,nrm_Tbar,'.k');
title('Thrust');
xlabel('$t$');
xlim([0,tvec(end)])
ylim([0,1.1*prb.Tmax])

subplot(2,2,4)
hold on
plot(prb.tau,prb.smin*ones(1,prb.K),'-r');
plot(prb.tau,prb.smax*ones(1,prb.K),'-r');
plot(prb.tau,ubar(prb.n+1,:),'.b');
plot(prb.tau,tvecbar,'.k');
p1 = plot(tau,tvec,'-k');
p2 = plot(tau,u(prb.n+1,:),'-b');
legend([p1,p2],{'$t(\tau)$','$s(\tau)$'})
xlabel('$\tau$');
title('Time \& Dilation')