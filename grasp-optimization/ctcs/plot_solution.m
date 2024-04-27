clc
close all
clearvars

load recent_solution

figure
plt_ctscvx = plot3(x(1,:),x(2,:),x(3,:),'--k');
hold on
plot3(xbar(1,:),xbar(2,:),xbar(3,:),'.k','MarkerSize',13);
axis('manual')

lw = 1.5;
scl = 1;
scl_box = 4;
for j=1:4:prb.K-1

    start_pos = xbar(1:3,j) + scl_box*prb.contact1; 
    end_pos = - scl*ubar(1:3,j);
    quiver3(start_pos(1),start_pos(2),start_pos(3),end_pos(1),end_pos(2),end_pos(3),'-','Color',[0,0,1],'LineWidth',2,'ShowArrowHead','off');

    start_pos = xbar(1:3,j) + scl_box*prb.contact2; 
    end_pos = - scl*ubar(4:6,j);
    quiver3(start_pos(1),start_pos(2),start_pos(3),end_pos(1),end_pos(2),end_pos(3),'-','Color',[1,0,0],'LineWidth',2,'ShowArrowHead','off');

    start_pos = xbar(1:3,j) + scl_box*prb.contact3; 
    end_pos = - scl*ubar(7:9,j);
    quiver3(start_pos(1),start_pos(2),start_pos(3),end_pos(1),end_pos(2),end_pos(3),'-','Color',[0,0.7,0.7],'LineWidth',2,'ShowArrowHead','off');
    
    cube_origin = xbar(1:3,j)' - scl_box*prb.box_width*ones(1,3);
    plotcube(scl_box*2*prb.box_width*ones(1,3),cube_origin,0.3,[1,0.8,0])
    
end
xlabel('$x$ [m]')
ylabel('$y$ [m]');
zlabel('$z$ [m]');
grid on
% view(0,90);

axis equal 

load('../node-only-cnstr/recent_solution_cvx');

plt_cvx = plot3(x_sim(1,:),x_sim(2,:),x_sim(3,:),'-m');
hold on
plot3(x(1,:),x(2,:),x(3,:),'om','MarkerSize',13);

legend([plt_ctscvx,plt_cvx],{'Prox-Linear','One-shot Convex Solve'},'Location','north')

ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];

% exportgraphics(ax,'results/position.pdf');

U = {zeros(3,prb.K),zeros(3,prb.K),zeros(3,prb.K)};  
nrmU = {zeros(1,prb.K),zeros(1,prb.K),zeros(1,prb.K)};
for k = 1:prb.K
    U{1}(:,k) = ubar(1:3,k);
    U{2}(:,k) = ubar(4:6,k);
    U{3}(:,k) = ubar(7:9,k);

    nrmU{1}(:,k) = norm(ubar(1:3,k));
    nrmU{2}(:,k) = norm(ubar(4:6,k));
    nrmU{3}(:,k) = norm(ubar(7:9,k));
end

figure
pltF1 = stairs(prb.tau,nrmU{1},'Color',[0,0,1]);
hold on
pltF2 = stairs(prb.tau,nrmU{2},'Color',[1,0,0]);
pltF3 = stairs(prb.tau,nrmU{3},'Color',[0.7,0.7,0]);
xlim([0,prb.tau(end)]);
legend([pltF1,pltF2,pltF3],{'Finger 1','Finger 2','Finger 3'});
title('Grasping Force [kg m s$^{-2}$]');
xlabel('$t$ [s]');

ax = gca;
% exportgraphics(ax,'results/grasp_force.pdf');