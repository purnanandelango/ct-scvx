clearvars
close all

load recent_solution

figure
subplot(2,2,1)
if prb.n == 3
    plot3(r(1,:),r(2,:),r(3,:),'-b');
    hold on 
    plot3(xbar(1,:),xbar(2,:),xbar(3,:),'ob');
else
    plot(r(1,:),r(2,:),'-b');
    hold on 
    plot(xbar(1,:),xbar(2,:),'ob');    
end 
title('Position');

if isfield(prb,'robs')
    if prb.n == 2
        th = linspace(0,2*pi);
        for j = 1:prb.nobs
            pobs = prb.robs(:,j) + prb.qobs(j)*[cos(th);sin(th)];
            plot(pobs(1,:),pobs(2,:),'-k');
        end
    end
end

ax = gca;
ax.DataAspectRatio = [1,1,1];
ax.PlotBoxAspectRatio = [1,1,1];

Kfine = length(tau);
if ~isequal(prb.disc,"Impulse")
    nrm_T(Kfine) = 0;
    for j = 1:Kfine
        nrm_T(j) = norm(u(1:prb.n,j));
    end
end
nrm_v(Kfine) = 0;
for j = 1:Kfine
    nrm_v(j) = norm(v(1:prb.n,j));
end
nrm_Tbar(prb.K) = 0;
nrm_vbar(prb.K) = 0;
for j = 1:prb.K
    nrm_Tbar(j) = norm(ubar(1:prb.n,j));
    nrm_vbar(j) = norm(xbar(prb.n+1:2*prb.n,j));
end

subplot(2,2,2)
plot(tvec,nrm_v,'-m');
hold on 
plot(tvecbar,nrm_vbar,'om');
plot(tvecbar,prb.vmax*ones(1,prb.K),'-r');
title('Speed')
xlabel('$t$');
xlim([0,tvec(end)])
ylim([0,1.1*prb.vmax])

subplot(2,2,3)
if ~isequal(prb.disc,"Impulse")
    plot(tvec,nrm_T,'-b');
end
hold on 
plot(tvecbar,nrm_Tbar,'ob');
plot(tvecbar,prb.Tmin*ones(1,prb.K),'-r');
plot(tvecbar,prb.Tmax*ones(1,prb.K),'-r');
title('Thrust');
xlabel('$t$');
xlim([0,tvec(end)])
ylim([0,1.1*prb.Tmax])