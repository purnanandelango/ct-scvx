clc
close all
clearvars

load recent_solution

figure
plt_ptr = plot3(xbar(1,:),xbar(2,:),xbar(3,:),'.-k','MarkerSize',13);
hold on
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
    quiver3(start_pos(1),start_pos(2),start_pos(3),end_pos(1),end_pos(2),end_pos(3),'-','Color',[0.7,0.7,0],'LineWidth',2,'ShowArrowHead','off');

    start_pos = xbar(1:3,j) + scl_box*prb.contact4; 
    end_pos = - scl*ubar(10:12,j);
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

load recent_solution_cvx

plt_cvx = plot3(x(1,:),x(2,:),x(3,:),'.-m','MarkerSize',13);

legend([plt_ptr,plt_cvx],{'Prox-Linear','One-shot Convex Solve'})

ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];
