clearvars
close all

file_name = "recent_solution";
% file_name = "recent_solution_static";

load(file_name);

fig = figure('Position',[138,100,1070,490], ...
             'Color',[1,1,1]);

ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];        
ax.XLim = [-80,80];
ax.YLim = [-30,30];

th = linspace(0,2*pi);

if prb.dyn_obs
    m = 15;
    F(1+(prb.Kfine-1)/m) = struct('cdata',[],'colormap',[]);
    cntr = 1;
    for k = 1:m:prb.Kfine
        plot(x_grid(1,:),x_grid(2,:),'-','Color',[0,0,0,0.3]);    
        hold on
        for j = 1:prb.n_obs
            pos_obs = prb.q_obs(t_grid(k),j) + prb.H_obs{j}\[cos(th);sin(th)];
            plot(pos_obs(1,:),pos_obs(2,:),'-r');
        end  
        plot(x_grid(1,k),x_grid(2,k),'.','Color',[0,0,0]);    
        ax = gca;
        ax.DataAspectRatio = [1,1,1];
        ax.PlotBoxAspectRatio = [1,1,1];        
        ax.XLim = [-80,80];
        ax.YLim = [-30,30];
        ax.YTickLabel = strrep(ax.YTickLabel,'-',char(8722));
        ax.XTickLabel = strrep(ax.XTickLabel,'-',char(8722));

        title('Position');
        ylabel('[m]');
        xlabel('[m]');

        F(cntr) = getframe(fig);
        cntr = cntr + 1; 
        cla(ax);
    end

    vid = VideoWriter('anim.avi','Motion JPEG AVI');
    open(vid)
    writeVideo(vid,F);
    close(vid)

    for k = 1:1+(prb.Kfine-1)/m
        [A,map] = rgb2ind(frame2im(F(k)),256);
        if k==1
            imwrite(A,map,'anim.gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,'anim.gif','WriteMode','append','DelayTime',0.1);
        end
    end
else

    th = linspace(0,2*pi);
    for j = 1:prb.n_obs
        pos_obs = prb.q_obs(t_grid(1),j) + prb.H_obs{j}\[cos(th);sin(th)];
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

n=2;
cnstr_fun_ndim       = @(rvpt,T) [-norm(prb.H_obs{ 1}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 1))) + 1; % Obstacle avoidance
                                 -norm(prb.H_obs{ 2}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 2))) + 1; % Obstacle avoidance
                                 -norm(prb.H_obs{ 3}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 3))) + 1; % Obstacle avoidance
                                 -norm(prb.H_obs{ 4}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 4))) + 1; % Obstacle avoidance
                                 -norm(prb.H_obs{ 5}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 5))) + 1; % Obstacle avoidance
                                 -norm(prb.H_obs{ 6}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 6))) + 1; % Obstacle avoidance
                                 -norm(prb.H_obs{ 7}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 7))) + 1; % Obstacle avoidance
                                 -norm(prb.H_obs{ 8}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 8))) + 1; % Obstacle avoidance
                                 -norm(prb.H_obs{ 9}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 9))) + 1; % Obstacle avoidance
                                 -norm(prb.H_obs{10}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2),10))) + 1; % Obstacle avoidance
                                  norm(rvpt(n+1:2*n))/prb.vmax - 1;                              % Speed upperbound
                                  norm(T)/prb.Tmax - 1;                                          % Thrust upperbound
                                 -norm(T)/prb.Tmin + 1;                                          % Thrust lowerbound
                                             ];

cnstr_viol(13,prb.Kfine) = 0;
for k = 1:prb.Kfine
    cnstr_viol(:,k) = arrayfun(@(y) max(0,y), cnstr_fun_ndim(x(:,k),u(1:2,k)));
end
figure
% semilogy(tvec,cnstr_viol(1,:),'DisplayName','Obstacle 1','Color',[0.5,0,0]);
semilogy(tvec,max(cnstr_viol(1:10,:),[],1),'DisplayName','Obstacle','Color',[0.5,0,0]);
hold on
% semilogy(tvec,cnstr_viol(2,:),'DisplayName','Obstacle 2','Color',[0,0,0.5]);
% semilogy(tvec,cnstr_viol(3,:),'DisplayName','Obstacle 3','Color',[0,0.5,0]);
% semilogy(tvec,cnstr_viol(4,:),'DisplayName','Obstacle 4','Color',[0.5,0.5,0]);
% semilogy(tvec,cnstr_viol(5,:),'DisplayName','Obstacle 5','Color',[0.5,0,0.5]);
% semilogy(tvec,cnstr_viol(6,:),'DisplayName','Obstacle 6','Color',[0,0.5,0.5]);
% semilogy(tvec,cnstr_viol(7,:),'DisplayName','Obstacle 7','Color',[0.5,0.7,0.5]);
% semilogy(tvec,cnstr_viol(8,:),'DisplayName','Obstacle 8','Color',[0.7,0.5,0.5]);
% semilogy(tvec,cnstr_viol(9,:),'DisplayName','Obstacle 9','Color',[0.7,0.5,0.7]);
% semilogy(tvec,cnstr_viol(10,:),'DisplayName','Obstacle 10','Color',[0.7,0.7,0.5]);
% semilogy(tvec,cnstr_viol(11,:),'DisplayName','Speed','Color',[0.5,0.7,0.7]);
% semilogy(tvec,cnstr_viol(12,:),'DisplayName','Accl. Upper','Color',[0.7,0.7,0.7]);
semilogy(tvec,cnstr_viol(13,:),'DisplayName','Accl. Lower','Color',[0.0,0.0,0.0]);
title("Constraint Violation")
legend('Location','best');
ylim([1e-4,1e0])

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

save(file_name,"cnstr_viol","-append");