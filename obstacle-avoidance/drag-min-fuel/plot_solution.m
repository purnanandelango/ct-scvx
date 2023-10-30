clearvars
close all

load recent_solution

figure

if prb.n == 2
    plot(r(1,:),r(2,:),'-k');
    hold on 
    plot(xbar(1,:),xbar(2,:),'.k','MarkerSize',15);    

    th = linspace(0,2*pi);
    hold on
    for j = 1:prb.nobs
        pobs = prb.qobs(:,j) + prb.Hobs{j}\[cos(th);sin(th)];
        plot(pobs(1,:),pobs(2,:),'-b','LineWidth',2);
    end

    ax = gca;
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

else
    plot3(r(1,:),r(2,:),r(3,:),'-k');
    hold on 
    plot3(xbar(1,:),xbar(2,:),xbar(3,:),'.k','MarkerSize',15);
    grid on
    
    [X1,Y1,Z1] = ellipsoid(0,0,0,1/0.3,1/0.1,1/0.3);
    X1 = X1 + prb.qobs(1,1);
    Y1 = Y1 + prb.qobs(2,1);
    Z1 = Z1 + prb.qobs(3,1);
    ellip1 = surf(X1,Y1,Z1,'EdgeColor',[1,0,0],'EdgeAlpha',0.5,'FaceColor',[1,0,0],'FaceAlpha',0.2);
    rotate(ellip1,[0,0,1],0);

    [X2,Y2,Z2] = ellipsoid(0,0,0,1/0.1,1/0.3,1/0.3);
    X2 = X2 + prb.qobs(1,2);
    Y2 = Y2 + prb.qobs(2,2);
    Z2 = Z2 + prb.qobs(3,2);
    ellip2 = surf(X2,Y2,Z2,'EdgeColor',[0,0,1],'EdgeAlpha',0.5,'FaceColor',[0,0,1],'FaceAlpha',0.2);
    rotate(ellip2,[0,0,1],0);    

    view(-6,1);
end
title('Position [m]');

ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];

ax = gca;
exportgraphics(ax,'results_3D/position.pdf','BackgroundColor','none');

figure
nrm_T(prb.Kfine) = 0;
nrm_v(prb.Kfine) = 0;
for j = 1:prb.Kfine
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
plot(tvec,nrm_v,'-m');
hold on 
plot(tvecbar,nrm_vbar,'.m');
plot(tvecbar,prb.vmax*ones(1,prb.K),'-r');
title('Speed [m s$^{-1}$]')
xlabel('$t$ [s]');
xlim([0,tvec(end)])
ylim([0,1.1*prb.vmax]);
ax = gca;
exportgraphics(ax,'results_3D/speed.pdf','BackgroundColor','none');

subplot(2,2,3)
plot(tvec,nrm_T,'-b');
hold on 
plot(tvecbar,nrm_Tbar,'.b');
plot(tvecbar,prb.umin*ones(1,prb.K),'-r');
plot(tvecbar,prb.umax*ones(1,prb.K),'-r');
title('Thrust [m s$^{-2}$]');
xlabel('$t$ [s]');
xlim([0,tvec(end)])
ylim([0,1.1*prb.umax])
ax = gca;
exportgraphics(ax,'results_3D/thrust.pdf','BackgroundColor','none');

subplot(2,2,4)
hold on
plot(prb.tau,tvecbar,'.k');
p1 = plot(tau,tvec,'-k');
plot(prb.tau,ubar(prb.n+1,:),'.g');
p2 = plot(tau,u(prb.n+1,:),'-g');
plot(prb.tau,prb.smin*ones(1,prb.K),'-r','LineWidth',1);
plot(prb.tau,prb.smax*ones(1,prb.K),'-r','LineWidth',1);
legend([p1,p2],{'$t(\tau)$','$s(\tau)$'})
xlabel('$\tau$');
ylim([0,1.1*prb.smax]);
title('Time \& Dilation')