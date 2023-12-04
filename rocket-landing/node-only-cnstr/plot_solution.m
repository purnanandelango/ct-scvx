clearvars
close all

load recent_solution

figure
% (mplot,nplot)
mplot = 3;
nplot = 4;

% Trajectory, vehicle axis, thrust vector and drag
subplot(mplot,nplot,[1,2,5,6])
n = 5;
plant.rocket6DoF.plot_vehicle_forces(u(1:3,1:n:end),rI(:,1:n:end),vI(:,1:n:end),qBI(:,1:n:end),0.4,0.2,struct('scl',0.4,'rho',prb.rho,'SA',prb.SA,'CA',prb.CA),{[2,3,1],{'y','z','x'},'x'});
plot3(xbar(3,:),xbar(4,:),xbar(2,:),'o','Color',[0.7,0.1,0.7],'LineWidth',0.5,'DisplayName','SCP solution');
legend('Position',[0.279,0.442,0.236,0.167])
grid on

% Thrust
subplot(mplot,nplot,[3,4])
plt.plot_vec_nrm(tvec,u(1:3,:),[prb.TBmin,prb.TBmax],2,'$\|T_{\mathcal{B}}(t)\|_2$','linear','blue');
nrm_Tbar = misc.compute_vec_norm(ubar(1:3,:));
hold on 
plot(tvecbar,nrm_Tbar,'ob');

nrm_vI = misc.compute_vec_norm(vI);
[~,idx] = min(abs(nrm_vI-prb.vmax_stc));
t_trig = tvec(min(idx));    % Trigger time

% Speed
subplot(mplot,nplot,[7,8])
plt.plot_vec_nrm(tvec,vI,prb.vmax,2,'Speed: $\|v_{\mathcal{I}}(t)\|_2$','linear','blue',"Single shot");
nrm_vIbar = misc.compute_vec_norm(xbar(5:7,:));
hold on 
legend('AutoUpdate','on');
plot(tvecbar,nrm_vIbar,'ob','DisplayName','SCP');
legend('AutoUpdate','off','Position',[0.773,0.507,0.122,0.075]);
plot(t_trig*ones(1,100),linspace(0,prb.vmax),'-k');
plot(tvec,prb.vmax_stc*ones(1,prb.Kfine),'-k');
ylabel("[L T$^{-1}]$",'FontSize',18);

% Dilation factors
subplot(mplot,nplot,[9,10])
plot(tau,u(4,:),'-b','DisplayName','Dilation Factors');
hold on
plot(tau,tvec,'-r','DisplayName','Time');
legend('Location','best','AutoUpdate','off');
plot(taubar,ubar(4,:),'ob');
plot(taubar,tvecbar,'or');
xlim([0,1]);

% Angle of attack
AoA = misc.compute_vec_norm(plant.rocket6DoF.compute_AoA(vI,qBI));
AoAbar = misc.compute_vec_norm(plant.rocket6DoF.compute_AoA(xbar(5:7,:),xbar(8:11,:)));

subplot(mplot,nplot,[11,12])
plot(tvec,AoA,'-b');
hold on
plot(tvecbar,AoAbar,'ob');
title("Angle of attack [${}^\circ$]");

% plot(t_trig*ones(1,100),linspace(0.9*min(AoA),1.1*max(AoA)),'-k');
plot(tvec,acosd(prb.cosaoamax)*ones(1,prb.Kfine),'-k');
xlim([0,tvec(end)]);
ylim([0,62]);
plot(t_trig*ones(1,100),linspace(0,62),'-k');
% ylim([0.9*min(AoA),1.1*max(AoA)]);

cnstr_viol(8,prb.Kfine) = 0;
stc_viol(2,prb.Kfine) = 0;
for k = 1:prb.Kfine
    cnstr_viol(:,k) = arrayfun(@(y) max(0,y), prb.cnstr_fun(x(:,k),u(1:3,k)));
    [~,~,stc_trig,~,stc_cnstr,~] = plant.rocket6DoF.q_aoa_cnstr(x(5:7,k),x(8:11,k),prb.vmax_stc,prb.cosaoamax,"v1");
    stc_viol(:,k) = [min(0,stc_trig)^2; max(0,stc_cnstr)^2];
end
figure
semilogy(tvec,cnstr_viol(1,:),'DisplayName','Mass');
hold on
semilogy(tvec,cnstr_viol(2,:),'DisplayName','Glide-slope');
semilogy(tvec,cnstr_viol(3,:),'DisplayName','Speed');
semilogy(tvec,cnstr_viol(4,:),'DisplayName','Tilt');
semilogy(tvec,cnstr_viol(5,:),'DisplayName','Angular Vel.');
semilogy(tvec,cnstr_viol(6,:),'DisplayName','Thrust gimbal');
semilogy(tvec,cnstr_viol(7,:),'DisplayName','Thrust upper');
semilogy(tvec,cnstr_viol(8,:),'DisplayName','Thrust lower');
title("Constraint Violation")
semilogy(tvec,stc_viol(1,:),"DisplayName","Trigger")
semilogy(tvec,stc_viol(2,:),"DisplayName","STC")
legend('Location','best');

figure
% Angular velocity
subplot(2,2,1)
plt.plot_vec_nrm(tvec,omgB*180/pi,prb.omgmax*180/pi,2,'Angular speed: $\|\omega_{\mathcal{B}}\|_2$','linear','blue');
nrm_omgbar = misc.compute_vec_norm(xbar(12:14,:))*180/pi;
hold on 
plot(tvecbar,nrm_omgbar,'ob');
ylim([0,35])
ylabel("[deg T$^{-1}$]",'FontSize',18);
xlabel("[T]",'FontSize',18);

% Tilt
subplot(2,2,2)
% legend('AutoUpdate','on');
plot(tvecbar,ones(1,prb.K)*prb.thetmax*180/pi+1,'--r');
% legend('AutoUpdate','off');
hold on
tilt_cnstr = 2*asind(misc.compute_vec_norm(prb.Hthet*x(8:11,:)));
plot(tvec,tilt_cnstr,'-b');
tilt_cnstr_bar = 2*asind(misc.compute_vec_norm(prb.Hthet*xbar(8:11,:)));
plot(tvecbar,tilt_cnstr_bar,'ob');
title('Tilt angle: $2\arcsin\|H_{\theta}q_{\mathcal{B}}\|_2$')
ylabel("[deg]",'FontSize',18);
ylim([0,50])
xlabel("[T]",'FontSize',18);

% Thrust
subplot(2,2,3)
plt.plot_vec_nrm(tvec,u(1:3,:),[prb.TBmin,prb.TBmax],2,'Thrust: $\|T_{\mathcal{B}}\|_2$','linear','blue');
nrm_Tbar = misc.compute_vec_norm(ubar(1:3,:));
hold on 
plot(tvecbar,nrm_Tbar,'ob');
ylabel("[M L T$^{-2}]$",'FontSize',18);
xlabel("[T]",'FontSize',18);

% Speed
subplot(2,2,4)
nrm_vI = misc.compute_vec_norm(vI);
plt.plot_vec_nrm(tvec,vI,prb.vmax,2,'Speed: $\|v_{\mathcal{I}}\|_2$','linear','blue');
nrm_vIbar = misc.compute_vec_norm(xbar(5:7,:));
hold on 
plot(tvecbar,nrm_vIbar,'ob');
ylabel("[L T$^{-1}]$",'FontSize',18);
xlabel("[T]",'FontSize',18);

%%% Diagnostics %%%
% figure
% subplot(2,2,1)
% plot3(rI(1,:),rI(2,:),rI(3,:),'-r');
% hold on
% plot3(x2(2,:),x2(3,:),x2(4,:),'--b');
% title("Position");
% axis equal
% 
% subplot(2,2,2)
% plot3(vI(1,:),vI(2,:),vI(3,:),'-r');
% hold on
% plot3(x2(5,:),x2(6,:),x2(7,:),'--b');
% title("Velocity");
% axis equal
% 
% subplot(2,2,3)
% plot3(qBI(1,:),qBI(2,:),qBI(3,:),'-r');
% hold on
% plot3(x2(8,:),x2(9,:),x2(10,:),'--b');
% title("Quaternion [1:3]");
% axis equal
%%%%%%