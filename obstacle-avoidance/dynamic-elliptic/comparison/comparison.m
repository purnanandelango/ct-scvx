clc 
clearvars
close all

interpreter = "tex";
% interpreter = "latex";

s1 = load('../ctcs/recent_solution_static.mat');

s2 = load('../static-node-only-cnstr/recent_solution');

prb = s1.prb;
tau = s1.tau;
tvecbar = s1.tvecbar;
tvec = s1.tvec;

fig = figure('Position',[214,19,407,407]);

th = linspace(0,2*pi);
for j = 1:prb.n_obs
    pos_obs = s2.prb.q_obs(:,j) + prb.H_obs{j}\[cos(th);sin(th)];
    plot(pos_obs(1,:),pos_obs(2,:),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
    hold on
end

plot(s1.r(1,:),s1.r(2,:),'-k');
plot(s1.xbar(1,:),s1.xbar(2,:),'.k');    
plot(s2.r(1,:),s2.r(2,:),'--','Color',[0.7,0.7,0.7]);
plot(s2.xbar(1,:),s2.xbar(2,:),'.','Color',[0.7,0.7,0.7]);    
ylabel('[m]');
xlabel('[m]');

ax = gca;
ax.Box = 'off';
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];
ax.XLim = [-30,30];
ax.YLim = [-35,35];

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

fig = figure('Position',[215,669,407,228]);
plot(tvecbar,prb.vmax*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on 
plot(tvec,nrm_v1,'-k');
plot(tvecbar,nrm_vbar1,'.k');
plot(tvec,nrm_v2,'--','Color',[0.7,0.7,0.7]);
plot(tvecbar,nrm_vbar2,'.','Color',[0.7,0.7,0.7]);
if interpreter == "latex"
    ylabel('[m s$^{-1}$]');
    xlabel('$t$ [s]');
elseif interpreter == "tex"
    ylabel("[m s^{"+char(8722)+"1}]");
    xlabel('{\it t} [s]');
end
xlim([0,tvec(end)])
ylim([0,1.1*prb.vmax])
ax = gca;
ax.Box = 'off';
exportgraphics(fig,'speed.pdf','ContentType','vector');

fig = figure('Position',[215,669,407,228]);
plot(tvecbar,prb.Tmin*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on 
plot(tvecbar,prb.Tmax*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plot(tvec,nrm_T1,'-k');
plot(tvecbar,nrm_Tbar1,'.k');
plot(tvec,nrm_T2,'--','Color',[0.7,0.7,0.7]);
plot(tvecbar,nrm_Tbar2,'.','Color',[0.7,0.7,0.7]);
if interpreter == "latex"
    ylabel('[m s$^{-2}$]');
    xlabel('$t$ [s]');
elseif interpreter == "tex"
    ylabel("[m s^{"+char(8722)+"2}]");
    xlabel('{\it t} [s]');
end
xlim([0,tvec(end)])
ylim([-0.25,1.1*prb.Tmax])
ax = gca;
ax.Box = 'off';

annotation(fig,"rectangle",[0.3966,0.2807,0.2054,0.136],'LineWidth',1);
annotation(fig,"arrow",'Position',[0.609,0.364,0.0786,0.2149],...
           'LineWidth',1,'HeadStyle','plain','HeadLength',7,'HeadWidth',4);

axes('Position',[0.5377,0.6123,0.3737,0.3439]);
plot(tvecbar,prb.Tmin*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
hold on 
plot(tvecbar,prb.Tmax*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plot(tvec,nrm_T1,'-k');
plot(tvecbar,nrm_Tbar1,'.k');
plot(tvec,nrm_T2,'--','Color',[0.7,0.7,0.7]);
plot(tvecbar,nrm_Tbar2,'.','Color',[0.7,0.7,0.7]);
ax = gca;
ax.XLim = [7.6,13.2];
ax.YLim = [0,1.5];
ax.XTick = [8,13];
ax.YTick = [0,1.5];
ax.FontSize = 16;
set(ax,'LineWidth',1);

exportgraphics(fig,'acceleration.pdf','ContentType','vector');

fig = figure('Position',[215,669,407,228]);
hold on
plot(prb.tau,prb.smin*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plot(prb.tau,prb.smax*ones(1,prb.K),'-','LineWidth',3.5,'Color',[1,0.5,0.5]);
plot(prb.tau,s1.ubar(prb.n+1,:),'.k');
plot(tau,s1.u(prb.n+1,:),'-k');
plot(s2.prb.tau,s2.ubar(prb.n+1,:),'.','Color',[0.7,0.7,0.7]);
plot(tau,s2.u(prb.n+1,:),'--','Color',[0.7,0.7,0.7]);
if interpreter == "latex"
    xlabel('$\tau$');
elseif interpreter == "tex"
    xlabel('\tau');
end
xlim([0,1])
ylim([-2,1.05*prb.smax])
ax = gca;
ax.Box = 'off';
exportgraphics(fig,'dilation.pdf','ContentType','vector');