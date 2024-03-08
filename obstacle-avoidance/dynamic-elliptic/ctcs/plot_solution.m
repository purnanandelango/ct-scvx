clearvars
close all

load recent_solution.mat
fig = figure('Position',[138,100,1034,437], ...
             'Color',[1,1,1]);

ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];        
ax.XLim = [-80,80];
ax.YLim = [-30,30];
ax.Units = 'pixels';
% pos = ax.Position;
% ti = ax.TightInset;
% rect = [0, 0, pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];

th = linspace(0,2*pi);

if prb.dyn_obs
    F(prb.Kfine) = struct('cdata',[],'colormap',[]);
    for k = 1:prb.Kfine
        plot(x_grid(1,:),x_grid(2,:),'-','Color',[0,0,0,0.3]);    
        plot(x_grid(1,k),x_grid(2,k),'.','Color',[1,0,0]);    
        hold on
        for j = 1:prb.n_obs
            pos_obs = prb.q_obs(t_grid(k),j) + prb.H_obs{j}\[cos(th);sin(th)];
            plot(pos_obs(1,:),pos_obs(2,:),'-r');
        end  
        ax = gca;
        ax.DataAspectRatio = [1,1,1];
        ax.PlotBoxAspectRatio = [1,1,1];        
        ax.XLim = [-80,80];
        ax.YLim = [-30,30];
        
        title('Position [m]');

        F(k) = getframe(fig);
        cla(ax);
    end

    % vid = VideoWriter('anim.avi','Motion JPEG AVI');
    % open(vid)
    % writeVideo(vid,F);
    % close(vid)

    for k = 1:prb.Kfine
        [A,map] = rgb2ind(frame2im(F(k)),256);
        if k==1
            imwrite(A,map,'anim.gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,'anim.gif','WriteMode','append','DelayTime',0.1);
        end
    end
else
    plot(r(1,:),r(2,:),'-b');
    plot(xbar(1,:),xbar(2,:),'.b','MarkerSize',25);    
    rinit = plot(r(1,1),r(2,1),'bo');
    rfinal = plot(r(1,end),r(2,end),'bx');
    legend([rinit,rfinal],{'$r_{\mathrm{i}}$','$r_{\mathrm{f}}$'},'Location','east');
    
    ax = gca;
    ax.DataAspectRatio = [1,1,1];
    ax.PlotBoxAspectRatio = [1,1,1];
    
    ax.XLim = [-21,2];
    ax.YLim = [-5,30];
    
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
    
    title('Position [m]');    
end

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
plot(tvecbar,prb.vmax*ones(1,prb.K),'-r');
hold on 
plot(tvec,nrm_v,'-k');
plot(tvecbar,nrm_vbar,'.k');
title('Speed [m s$^{-1}$]')
xlabel('$t$ [s]');
xlim([0,tvec(end)])
ylim([0,1.1*prb.vmax]);

subplot(2,2,3)
plot(tvecbar,prb.Tmin*ones(1,prb.K),'-r');
hold on 
plot(tvecbar,prb.Tmax*ones(1,prb.K),'-r');
plot(tvec,nrm_T,'-k');
plot(tvecbar,nrm_Tbar,'.k','MarkerSize',25);
title('Thrust [m s$^{-2}$]');
xlabel('$t$ [s]');
xlim([0,tvec(end)])
ylim([0,1.1*prb.Tmax])

subplot(2,2,4)
plot(prb.tau,prb.smin*ones(1,prb.K),'-r');
hold on
plot(prb.tau,prb.smax*ones(1,prb.K),'-r');
plot(prb.tau,tvecbar,'.k','MarkerSize',25);
p1 = plot(tau,tvec,'-k');
plot(prb.tau,ubar(prb.n+1,:),'.b','MarkerSize',25);
p2 = plot(tau,u(prb.n+1,:),'-b');
legend([p1,p2],{'$t(\tau)$','$s(\tau)$'})
xlabel('$\tau$');
ylim([0,1.1*prb.smax]);
title('Time \& Dilation')