clc
clearvars
close all

interpreter = "tex";
% interpreter = "latex";

pxaxwidth = 1110; % [pixels]

% Load solutions
s1 = load('solution_eps_1e-4');

s2 = load('solution_node-only');

prb = s1.prb;

% Position
pxaxheight = pxaxwidth/1.75;
fig = figure('Position',[-1375,64,pxaxwidth*1.05,pxaxheight*1.05]);
hold on
plot(s1.x(4,:),s1.x(2,:),'-k');
xBI = plant.rocket6DoF.compute_bodyaxis(s1.qBI);
uBI = plant.rocket6DoF.compute_thrustdir(s1.u,s1.qBI);
bdscl = 25;
thrust = 0.003;
n = 40;
for k=1:n:prb.Kfine
    
    % bdscl = 6*(xBI(3,k)+xBI(1,k));
    bdscl = 25*(1.2*(1-k/prb.Kfine) + 3*(k/prb.Kfine))/1.5;

    quiver(s1.x(4,k),s1.x(2,k),...
           bdscl*xBI(3,k),bdscl*xBI(1,k),...
           'Color',[60,179,113]/255,'LineWidth',5.5);
    quiver(s1.x(4,k),s1.x(2,k),...
           -thrust*uBI(3,k),-thrust*uBI(1,k),...
           'Color',[213,114,14]/255,'LineWidth',5.5);
    quiver(s1.x(4,k),s1.x(2,k),...
           -0.5*thrust*uBI(3,k),-0.5*thrust*uBI(1,k),...
           'Color',[255,215,0]/255,'LineWidth',2.5);    
end
plot(s1.xbar(4,:),s1.xbar(2,:),'.k');
plot(s2.x(4,:),s2.x(2,:),'-','Color',[0,0,1,0.5]);
plot(s2.xbar(4,:),s2.xbar(2,:),'.','Color',[0.5,0.5,1]);
legend('off');
ylabel('Altitude [m]');
xlabel('Downrange [m]');
ylim([0,520]);
% title('Position');

ax = gca;
ax.XDir = "reverse";
if interpreter == "latex"
for k = 1:length(ax.YTickLabel)
    ax.YTickLabel{k} = "$"+ax.YTickLabel{k}+"$";
end
elseif interpreter == "tex"
    ax.YTickLabel = strrep(ax.YTickLabel,'-',char(8722));
end
ax.Box = "off";
ax.Units = "pixels";
ax.OuterPosition = [0, 0,  pxaxwidth, pxaxheight];

exportgraphics(fig,'traj.pdf','ContentType','vector');
savefig(fig,'traj.fig');

% Thrust
pxaxheight = pxaxwidth/2.5;
fig = figure('Position',[-1375,64,pxaxwidth*1.05,pxaxheight*1.05]);
nrm_T1 = misc.compute_vec_norm(s1.u(1:3,:));
nrm_Tbar1 = misc.compute_vec_norm(s1.ubar(1:3,:));
nrm_T2 = misc.compute_vec_norm(s2.u(1:3,:));
nrm_Tbar2 = misc.compute_vec_norm(s2.ubar(1:3,:));
hold on 
plot(s1.tvec,prb.TBmax*ones(1,length(s1.tvec)),'-','Color',[1,0.5,0.5],'LineWidth',4.5);
plot(s1.tvec,prb.TBmin*ones(1,length(s1.tvec)),'-','Color',[1,0.5,0.5],'LineWidth',4.5);
plot(s1.tvecbar,nrm_Tbar1,'.k');
plot(s1.tvec,nrm_T1,'-k');
plot(s2.tvecbar,nrm_Tbar2,'.','Color',[0.5,0.5,1]);
plot(s2.tvec,nrm_T2,'-','Color',[0,0,1,0.5]);
% title('Thrust');
if interpreter == "tex"
    xlabel('{\it t} [s]');
elseif interpreter == "latex"
    xlabel('$t$ [s]');
end
ylabel('[N]');
xlim([0,s1.tvec(end)]);
ylim([0.7*prb.TBmin,1.05*prb.TBmax]);
ax = gca;
ax.Box = "off";
ax.Units = "pixels";
ax.OuterPosition = [0, 0, pxaxwidth, pxaxheight];

%%% Magnified inset
annotation(fig,"rectangle",[0.1281    0.2414    0.0639    0.0570],'LineWidth',1);
% annotation(fig,"arrow",'Position',[ 0.1965    0.3003    0.3132    0.1330],...
%            'LineWidth',2,'HeadStyle','plain','HeadLength',9,'HeadWidth',6,'Color',[0.2,0.2,0.2]);

axes('Position',[0.5366    0.4794    0.2484    0.2605]);
plot(s1.tvec,prb.TBmin*ones(1,length(s1.tvec)),'-','Color',[1,0.5,0.5],'LineWidth',4.5);
hold on
plot(s1.tvecbar,nrm_Tbar1,'.k');
plot(s1.tvec,nrm_T1,'-k');
plot(s2.tvecbar,nrm_Tbar2,'.','Color',[0.5,0.5,1]);
plot(s2.tvec,nrm_T2,'-','Color',[0,0,1,0.5]);
xlim([0,7]);
ylim([4000,6000]);
xticks([0,3,6]);
yticks([4000,6000]);

ax = gca;
ax.FontSize = 30;
set(ax,'LineWidth',1);
%%%

exportgraphics(fig,'thrust.pdf','ContentType','vector');
savefig(fig,'thrust.fig');

% Angular speed
fig = figure('Position',[215,669,407,228]);
nrm_omg1 = misc.compute_vec_norm(s1.x(12:14,:))*180/pi;
nrm_omgbar1 = misc.compute_vec_norm(s1.xbar(12:14,:))*180/pi;
nrm_omg2 = misc.compute_vec_norm(s2.x(12:14,:))*180/pi;
nrm_omgbar2 = misc.compute_vec_norm(s2.xbar(12:14,:))*180/pi;
hold on 
plot(s1.tvec,prb.omgmax*ones(1,length(s1.tvec))*180/pi,'-','Color',[1,0.5,0.5],'LineWidth',4.5);
plot(s1.tvecbar,nrm_omgbar1,'.k');
plot(s1.tvec,nrm_omg1,'-k');
plot(s2.tvecbar,nrm_omgbar2,'.','Color',[0.5,0.5,1]);
plot(s2.tvec,nrm_omg2,'-','Color',[0,0,1,0.5]);
xlim([0,s1.tvec(end)]);
% title('Angular speed');
if interpreter == "tex"
    xlabel('{\it t} [s]');
    ylabel("[deg s^{"+char(8722)+"1}]");
elseif interpreter == "latex"
    xlabel('$t$ [s]');
    ylabel('[deg s$^{-1}$]')
end
ax = gca;
ax.FontSize = 20;
% exportgraphics(fig,'angspeed.pdf','ContentType','vector');

% Tilt angle
pxaxheight = pxaxwidth/2.5;
fig = figure('Position',[-1375,64,pxaxwidth*1.05,pxaxheight*1.05]);
nrm_tilt1 = 2*asind(misc.compute_vec_norm(prb.Hthet*s1.x(8:11,:)));
nrm_tiltbar1 = 2*asind(misc.compute_vec_norm(prb.Hthet*s1.xbar(8:11,:)));
nrm_tilt2 = 2*asind(misc.compute_vec_norm(prb.Hthet*s2.x(8:11,:)));
nrm_tiltbar2 = 2*asind(misc.compute_vec_norm(prb.Hthet*s2.xbar(8:11,:)));
hold on 
plot(s1.tvec,prb.thetmax*ones(1,length(s1.tvec))*180/pi,'-','Color',[1,0.5,0.5],'LineWidth',4.5);
plot(s1.tvecbar,nrm_tiltbar1,'.k');
plot(s1.tvec,nrm_tilt1,'-k');
plot(s2.tvecbar,nrm_tiltbar2,'.','Color',[0.5,0.5,1]);
plot(s2.tvec,nrm_tilt2,'-','Color',[0,0,1,0.5]);
xlim([0,s1.tvec(end)]);
if interpreter == "tex"
    xlabel('{\it t} [s]');
elseif interpreter == "latex"
    xlabel('$t$ [s]');
end
ylabel('[deg]');
ylim([0,70]);
% title('Tilt angle');
ax = gca;
ax.Box = "off";
ax.Units = "pixels";
ax.OuterPosition = [0, 0, pxaxwidth, pxaxheight];

%%% Magnified inset
% annotation(fig,"rectangle",[0.1348    0.6800    0.1672    0.1544],'LineWidth',1);
annotation(fig,"rectangle",[0.1348    0.7165    0.1672    0.1544],'LineWidth',1); % no title
% annotation(fig,"arrow",'Position',[0.3063    0.7722    0.1956   -0.2381],...
%            'LineWidth',2,'HeadStyle','plain','HeadLength',9,'HeadWidth',6,'Color',[0.2,0.2,0.2]);
axes('Position',[0.5159    0.4045    0.2975    0.3132]);

plot(s1.tvec,prb.thetmax*ones(1,length(s1.tvec))*180/pi,'-','Color',[1,0.5,0.5],'LineWidth',4.5);
hold on
plot(s1.tvecbar,nrm_tiltbar1,'.k');
plot(s1.tvec,nrm_tilt1,'-k');
plot(s2.tvecbar,nrm_tiltbar2,'.','Color',[0.5,0.5,1]);
plot(s2.tvec,nrm_tilt2,'-','Color',[0,0,1,0.5]);
ylim([54,70]);
xlim([1,17]);
yticks([55,65]);
xticks([5,15]);

ax = gca;
ax.FontSize = 30;
set(ax,'LineWidth',1);
%%%

exportgraphics(fig,'tilt.pdf','ContentType','vector');
savefig(fig,'tilt.fig');