% clc
% close all
clearvars

load recent_solution_cvx

figure
plot3(x(1,:),x(2,:),x(3,:),'om','MarkerSize',13);
hold on
plot3(x_sim(1,:),x_sim(2,:),x_sim(3,:),'-m');
axis('manual')

lw = 1.5;
scl = 1;
scl_box = 4;
for j=1:4:prb.K-1

    start_pos = x(1:3,j) + scl_box*prb.contact1; 
    end_pos = - scl*u(1:3,j);
    quiver3(start_pos(1),start_pos(2),start_pos(3),end_pos(1),end_pos(2),end_pos(3),'-','Color',[0,0,1],'LineWidth',2,'ShowArrowHead','off');

    start_pos = x(1:3,j) + scl_box*prb.contact2; 
    end_pos = - scl*u(4:6,j);
    quiver3(start_pos(1),start_pos(2),start_pos(3),end_pos(1),end_pos(2),end_pos(3),'-','Color',[1,0,0],'LineWidth',2,'ShowArrowHead','off');

    start_pos = x(1:3,j) + scl_box*prb.contact3; 
    end_pos = - scl*u(7:9,j);
    quiver3(start_pos(1),start_pos(2),start_pos(3),end_pos(1),end_pos(2),end_pos(3),'-','Color',[0,0.7,0.7],'LineWidth',2,'ShowArrowHead','off');
    
    cube_origin = x(1:3,j)' - scl_box*prb.box_width*ones(1,3);
    plotcube(scl_box*2*prb.box_width*ones(1,3),cube_origin,0.3,[1,0.8,0])
    
end
xlabel('$x$ [m]')
ylabel('$y$ [m]');
zlabel('$z$ [m]');
grid on

axis equal 

ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];


U = {zeros(3,prb.K),zeros(3,prb.K),zeros(3,prb.K)};  
nrmU = {zeros(1,prb.K),zeros(1,prb.K),zeros(1,prb.K)};
for k = 1:prb.K
    U{1}(:,k) = u(1:3,k);
    U{2}(:,k) = u(4:6,k);
    U{3}(:,k) = u(7:9,k);

    nrmU{1}(:,k) = norm(u(1:3,k));
    nrmU{2}(:,k) = norm(u(4:6,k));
    nrmU{3}(:,k) = norm(u(7:9,k));
end

nrm_v = misc.compute_vec_norm(x(4:6,:));
nrm_v_sim = misc.compute_vec_norm(x_sim(4:6,:));

figure
subplot(2,1,1)
stairs(prb.tau,nrmU{1},'Color',[0,0,1]);
hold on
stairs(prb.tau,nrmU{2},'Color',[1,0,0]);
stairs(prb.tau,nrmU{3},'Color',[0.8,0.6,0]);

subplot(2,1,2)
plot(prb.tau,nrm_v,'ro');
hold on
plot(prb.tau,prb.vmax*ones(1,prb.K),'--k');
plot(tau_sim',nrm_v_sim,'-r');
