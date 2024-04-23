clc 
clearvars
close all

interpreter = "tex";
% interpreter = "latex";

pxaxwidth = 1110; % [pixels]

% Load solutions
s1 = load('solution_eps_1e-5.mat');

s2 = load('solution_node-only.mat');

prb = s1.prb;
tau = s1.tau;
tvecbar = s1.tvecbar;
tvec = s1.tvec;

% Position
pxaxheight = pxaxwidth/2;
fig = figure('Position',[-1375,64,pxaxwidth*1.05,pxaxheight*1.05]);

th = linspace(0,2*pi);
for j = 1:prb.n_obs
    pos_obs = s2.prb.q_obs(:,j) + prb.H_obs{j}\[cos(th);sin(th)];
    plot(pos_obs(1,:),pos_obs(2,:),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
    hold on
end

plot(s1.r(1,:),s1.r(2,:),'-k');
plot(s1.xbar(1,:),s1.xbar(2,:),'.k');    
plot(s2.r(1,:),s2.r(2,:),'-','Color',[0,0,1,0.5]);
plot(s2.xbar(1,:),s2.xbar(2,:),'.','Color',[0.5,0.5,1]);    
ylabel('[m]');
xlabel('[m]');
xlim([-30,30]);
ylim([-30,30]);
xticks([-30,-10,10,30]);
yticks([-20,0,20]);
% title('Position');

ax = gca;
ax.Box = 'off';
ax.Units = "pixels";
ax.OuterPosition = [0, 0,  pxaxwidth, pxaxheight];


if interpreter == "latex"
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
elseif interpreter == "tex"
    ax.XTickLabel = strrep(ax.XTickLabel,'-',char(8722));
    ax.YTickLabel = strrep(ax.YTickLabel,'-',char(8722));
end

exportgraphics(fig,'traj.pdf','ContentType','vector');
savefig(fig,'traj.fig');

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

% Acceleration
pxaxheight = pxaxwidth/2.5;
fig = figure('Position',[-1375,64,pxaxwidth*1.05,pxaxheight*1.05]);
plot(tvecbar,prb.Tmin*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
hold on 
plot(tvecbar,prb.Tmax*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
plot(tvec,nrm_T1,'-k');
plot(tvecbar,nrm_Tbar1,'.k');
plot(tvec,nrm_T2,'-','Color',[0,0,1,0.5]);
plot(tvecbar,nrm_Tbar2,'.','Color',[0.5,0.5,1]);
if interpreter == "latex"
    ylabel('[m s$^{-2}$]');
    xlabel('$t$ [s]');
elseif interpreter == "tex"
    ylabel("[m s^{"+char(8722)+"2}]");
    xlabel('{\it t} [s]');
end
xlim([0,tvec(end)])
ylim([-0.25,1.1*prb.Tmax])
% title('Acceleration');
ax = gca;
ax.Box = "off";
ax.Units = "pixels";
ax.OuterPosition = [0, 0, pxaxwidth, pxaxheight];

%%% Magnified inset
annotation(fig,"rectangle",[0.3820    0.2479    0.1768    0.1255],'LineWidth',1);
% annotation(fig,"arrow",'Position',[0.5638    0.3359    0.0808    0.1684],...
%            'LineWidth',2,'HeadStyle','plain','HeadLength',9,'HeadWidth',6,'Color',[0.2,0.2,0.2]);

axes('Position',[0.5208,0.5269,0.3273,0.2822]);
plot(tvecbar,prb.Tmin*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
hold on 
plot(tvecbar,prb.Tmax*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
plot(tvec,nrm_T1,'-k');
plot(tvecbar,nrm_Tbar1,'.k');
plot(tvec,nrm_T2,'-','Color',[0,0,1,0.5]);
plot(tvecbar,nrm_Tbar2,'.','Color',[0.5,0.5,1]);
xlim([7.6,13.2]);
ylim([0,1.5]);
xticks([8,12]);
yticks([0.2,1]);

ax = gca;
ax.FontSize = 30;
set(ax,'LineWidth',1);
%%%

exportgraphics(fig,'acceleration.pdf','ContentType','vector');
savefig(fig,'acceleration.fig');



% Speed
fig = figure('Position',[215,669,407,228]);
plot(tvecbar,prb.vmax*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
hold on 
plot(tvec,nrm_v1,'-k');
plot(tvecbar,nrm_vbar1,'.k');
plot(tvec,nrm_v2,'-','Color',[0,0,1,0.5]);
plot(tvecbar,nrm_vbar2,'.','Color',[0.5,0.5,1]);
if interpreter == "latex"
    ylabel('[m s$^{-1}$]');
    xlabel('$t$ [s]');
elseif interpreter == "tex"
    ylabel("[m s^{"+char(8722)+"1}]");
    xlabel('{\it t} [s]');
end
% title('Speed');
xlim([0,tvec(end)])
ylim([0,1.1*prb.vmax])
ax = gca;
ax.FontSize = 20;
% exportgraphics(fig,'speed.pdf','ContentType','vector');

% Dilation and time
fig = figure('Position',[215,669,407,228]);
hold on
plot(prb.tau,prb.smin*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
plot(prb.tau,prb.smax*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
plot(prb.tau,s1.ubar(prb.n+1,:),'.k');
plot(tau,s1.u(prb.n+1,:),'-k');
plot(s2.prb.tau,s2.ubar(prb.n+1,:),'.','Color',[0.5,0.5,1]);
plot(tau,s2.u(prb.n+1,:),'-','Color',[0,0,1,0.5]);
if interpreter == "latex"
    xlabel('$\tau$');
elseif interpreter == "tex"
    xlabel('\tau');
end
% title('Dilation factor');
xlim([0,1])
ylim([-2,1.05*prb.smax])
ax = gca;
ax.FontSize = 20;
% exportgraphics(fig,'dilation.pdf','ContentType','vector');