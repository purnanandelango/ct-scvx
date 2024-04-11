clc
clearvars
close all

interpreter = "tex";
% interpreter = "latex";

s1 = load('../ctcs/recent_solution');

s2 = load('../node-only-cnstr/recent_solution');

prb = s1.prb;

fig = figure('Position',[214,19,407,407]);
n = 10;
hold on
plot3(s1.x(3,:),s1.x(4,:),s1.x(2,:),'-k');
plant.rocket6DoF.plot_vehicle_forces(s1.u(1:3,1:n:end),s1.rI(:,1:n:end),s1.vI(:,1:n:end),s1.qBI(:,1:n:end),25,0.003,struct('scl',0.4,'rho',prb.rho,'SA',prb.SA,'CA',prb.CA),{[2,3,1],{'y','z','x'},'x'});
plot3(s1.xbar(3,:),s1.xbar(4,:),s1.xbar(2,:),'.k');
plot3(s2.x(3,:),s2.x(4,:),s2.x(2,:),'--','Color',[0.7,0.7,0.7]);
plot3(s2.xbar(3,:),s2.xbar(4,:),s2.xbar(2,:),'.','Color',[0.7,0.7,0.7]);
legend('off');
title('');
zlabel('Altitude [m]');
ylabel('Downrange [m]');
view(-90,0);
ax = gca;
if interpreter == "latex"
for k = 1:length(ax.YTickLabel)
    ax.YTickLabel{k} = "$"+ax.YTickLabel{k}+"$";
end
elseif interpreter == "tex"
    ax.YTickLabel = strrep(ax.YTickLabel,'-',char(8722));
end
ax.PlotBoxAspectRatio = [1,1,1];
ax.DataAspectRatio = [1,1,1];
ax.Box = 'off';
exportgraphics(fig,'traj.pdf','ContentType','vector');

% Thrust
fig = figure('Position',[215,669,407,228]);
nrm_T1 = misc.compute_vec_norm(s1.u(1:3,:));
nrm_Tbar1 = misc.compute_vec_norm(s1.ubar(1:3,:));
nrm_T2 = misc.compute_vec_norm(s2.u(1:3,:));
nrm_Tbar2 = misc.compute_vec_norm(s2.ubar(1:3,:));
hold on 
plot(s1.tvec,prb.TBmax*ones(1,length(s1.tvec)),'-','Color',[1,0.5,0.5],'LineWidth',3.5);
plot(s1.tvec,prb.TBmin*ones(1,length(s1.tvec)),'-','Color',[1,0.5,0.5],'LineWidth',3.5);
plot(s1.tvecbar,nrm_Tbar1,'.k');
plot(s1.tvec,nrm_T1,'-k');
plot(s2.tvecbar,nrm_Tbar2,'.','Color',[0.7,0.7,0.7]);
plot(s2.tvec,nrm_T2,'--','Color',[0.7,0.7,0.7]);
xlim([0,s1.tvec(end)]);
ylim([0.7*prb.TBmin,1.05*prb.TBmax]);
if interpreter == "tex"
    xlabel('{\it t} [s]');
elseif interpreter == "latex"
    xlabel('$t$ [s]');
end
ylabel('[N]');
ax = gca;
ax.Box = 'off';

annotation(fig,"rectangle",[0.1572,0.2807,0.0639,0.0570],'LineWidth',1);
annotation(fig,"arrow",'Position',[0.226,0.34,0.2948,0.1750],...
           'LineWidth',1,'HeadStyle','plain','HeadLength',7,'HeadWidth',4);

axes('Position',[0.5452,0.5553,0.2,0.2]);
plot(s1.tvec,prb.TBmin*ones(1,length(s1.tvec)),'-','Color',[1,0.5,0.5],'LineWidth',3.5);
hold on
plot(s1.tvecbar,nrm_Tbar1,'.k');
plot(s1.tvec,nrm_T1,'-k');
plot(s2.tvecbar,nrm_Tbar2,'.','Color',[0.7,0.7,0.7]);
plot(s2.tvec,nrm_T2,'--','Color',[0.7,0.7,0.7]);
ax = gca;
ax.XLim = [0,7];
ax.YLim = [4000,6000];
ax.XTick = [0,3,7];
ax.YTick = [4000,6000];
ax.FontSize = 16;
set(ax,'LineWidth',1);

exportgraphics(fig,'thrust.pdf','ContentType','vector');

% Angular speed
fig = figure('Position',[215,669,407,228]);
nrm_omg1 = misc.compute_vec_norm(s1.x(12:14,:))*180/pi;
nrm_omgbar1 = misc.compute_vec_norm(s1.xbar(12:14,:))*180/pi;
nrm_omg2 = misc.compute_vec_norm(s2.x(12:14,:))*180/pi;
nrm_omgbar2 = misc.compute_vec_norm(s2.xbar(12:14,:))*180/pi;
hold on 
plot(s1.tvec,prb.omgmax*ones(1,length(s1.tvec))*180/pi,'-','Color',[1,0.5,0.5],'LineWidth',3.5);
plot(s1.tvecbar,nrm_omgbar1,'.k');
plot(s1.tvec,nrm_omg1,'-k');
plot(s2.tvecbar,nrm_omgbar2,'.','Color',[0.7,0.7,0.7]);
plot(s2.tvec,nrm_omg2,'--','Color',[0.7,0.7,0.7]);
xlim([0,s1.tvec(end)]);
if interpreter == "tex"
    xlabel('{\it t} [s]');
    ylabel("[deg s^{"+char(8722)+"1}]");
elseif interpreter == "latex"
    xlabel('$t$ [s]');
    ylabel('[deg s$^{-1}$]')
end
ax = gca;
ax.Box = 'off';
exportgraphics(fig,'angspeed.pdf','ContentType','vector');

% Tilt angle
fig = figure('Position',[215,669,407,228]);
nrm_tilt1 = 2*asind(misc.compute_vec_norm(prb.Hthet*s1.x(8:11,:)));
nrm_tiltbar1 = 2*asind(misc.compute_vec_norm(prb.Hthet*s1.xbar(8:11,:)));
nrm_tilt2 = 2*asind(misc.compute_vec_norm(prb.Hthet*s2.x(8:11,:)));
nrm_tiltbar2 = 2*asind(misc.compute_vec_norm(prb.Hthet*s2.xbar(8:11,:)));
hold on 
plot(s1.tvec,prb.thetmax*ones(1,length(s1.tvec))*180/pi,'-','Color',[1,0.5,0.5],'LineWidth',3.5);
plot(s1.tvecbar,nrm_tiltbar1,'.k');
plot(s1.tvec,nrm_tilt1,'-k');
plot(s2.tvecbar,nrm_tiltbar2,'.','Color',[0.7,0.7,0.7]);
plot(s2.tvec,nrm_tilt2,'--','Color',[0.7,0.7,0.7]);
xlim([0,s1.tvec(end)]);
if interpreter == "tex"
    xlabel('{\it t} [s]');
elseif interpreter == "latex"
    xlabel('$t$ [s]');
end
ylabel('[deg]');
ax = gca;
ax.Box = 'off';
ax.YLim = [0,70];
exportgraphics(fig,'tilt.pdf','ContentType','vector');
