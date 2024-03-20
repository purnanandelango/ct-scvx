clc 
clearvars
close all

s1 = load('ctcs/recent_solution_static.mat');

s2 = load('static-node-only-cnstr/recent_solution');

prb = s1.prb;
tau = s1.tau;
tvecbar = s1.tvecbar;
tvec = s1.tvec;

fig = figure('Position',[138,100,1034,437], ...
             'Color',[1,1,1]);

th = linspace(0,2*pi);
for j = 1:prb.n_obs
    pos_obs = s2.prb.q_obs(:,j) + prb.H_obs{j}\[cos(th);sin(th)];
    plot(pos_obs(1,:),pos_obs(2,:),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
    hold on
end

plot(s1.r(1,:),s1.r(2,:),'-k');
plot(s1.xbar(1,:),s1.xbar(2,:),'.k');    
plot(s2.r(1,:),s2.r(2,:),'-','Color',[0,0,0,0.5]);
plot(s2.xbar(1,:),s2.xbar(2,:),'.','Color',[0.5,0.5,0.5]);    
title('Position [m]');

ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];
ax.XLim = [-65,65];
ax.YLim = [-35,35];

ticklab = ax.XTickLabel;
for j = 1:length(ticklab)
    ticklab{j} = horzcat('$',ticklab{j},'$');
end
ax.XTickLabel = ticklab;
ticklab = ax.YTickLabel;
for j = 1:length(ticklab)
    ticklab{j} = horzcat('$',ticklab{j},'$');
end
ax.YTickLabel = ticklab;
exportgraphics(ax,'traj.pdf','BackgroundColor','none');

figure

Kfine = length(tau);

nrm_T1(Kfine) = 0;
nrm_v1(Kfine) = 0;
nrm_T2(Kfine) = 0;
nrm_v2(Kfine) = 0;
for j = 1:Kfine
    nrm_T1(j) = norm(s1.u(1:prb.n,j));
    nrm_v1(j) = norm(s1.v(1:prb.n,j));
    nrm_T2(j) = norm(s2.u(1:prb.n,j));
    nrm_v2(j) = norm(s2.v(1:prb.n,j));
end
nrm_Tbar1(prb.K) = 0;
nrm_vbar1(prb.K) = 0;
nrm_Tbar2(prb.K) = 0;
nrm_vbar2(prb.K) = 0;
for j = 1:prb.K
    nrm_Tbar1(j) = norm(s1.ubar(1:prb.n,j));
    nrm_vbar1(j) = norm(s1.xbar(prb.n+1:2*prb.n,j));
    nrm_Tbar2(j) = norm(s2.ubar(1:prb.n,j));
    nrm_vbar2(j) = norm(s2.xbar(prb.n+1:2*prb.n,j));    
end

subplot(2,2,2)
plot(tvecbar,prb.vmax*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on 
plot(tvec,nrm_v1,'-k');
plot(tvec,nrm_v2,'-','Color',[0,0,0,0.5]);
plot(tvecbar,nrm_vbar1,'.k');
plot(tvecbar,nrm_vbar2,'.','Color',[0.5,0.5,0.5]);
title('Speed [m/s]')
xlabel('$t$ [s]');
xlim([0,tvec(end)])
ylim([0,1.1*prb.vmax])
ax = gca;
exportgraphics(ax,'speed.pdf','BackgroundColor','none');

subplot(2,2,3)
plot(tvecbar,prb.Tmin*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on 
plot(tvecbar,prb.Tmax*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plot(tvec,nrm_T1,'-k');
plot(tvecbar,nrm_Tbar1,'.k');
plot(tvec,nrm_T2,'-','Color',[0,0,0,0.5]);
plot(tvecbar,nrm_Tbar2,'.','Color',[0.5,0.5,0.5]);
title('Acceleration [m/s$^2$]');
xlabel('$t$ [s]');
xlim([0,tvec(end)])
ylim([0,1.1*prb.Tmax])
ax = gca;
exportgraphics(ax,'acceleration.pdf','BackgroundColor','none');

subplot(2,2,4)
hold on
plot(prb.tau,prb.smin*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plot(prb.tau,prb.smax*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
% plot(prb.tau,tvecbar,'.k');
% plot(tau,tvec,'-k');
plot(prb.tau,s1.ubar(prb.n+1,:),'.k');
plot(tau,s1.u(prb.n+1,:),'-k');
plot(s2.prb.tau,s2.ubar(prb.n+1,:),'.','Color',[0.5,0.5,0.5]);
plot(tau,s2.u(prb.n+1,:),'-','Color',[0,0,0,0.5]);
% legend([p1,p2],{'$t(\tau)$','$s(\tau)$'})
xlabel('$\tau$');
title('Dilation Factor');
xlim([0,1])
ylim([0,1.1*prb.smax])
ax = gca;
exportgraphics(ax,'dilation.pdf','BackgroundColor','none');