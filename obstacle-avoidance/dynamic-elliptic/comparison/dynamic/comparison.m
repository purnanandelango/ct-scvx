clc
clearvars
close all

interpreter = "latex";
% interpreter = "tex";

pxaxwidth = 1110; % [pixels]
pxaxheight = pxaxwidth/2.5; % [pixels]

save_figures = true;

% Load solutions
for k = 6:-1:1
    s(k) = load("solution_eps_1e-"+num2str(k)); 
end

fig = figure('Position',[-1249,425,pxaxwidth*1.01,pxaxheight*1.01]);
semilogy(s(6).tvec,max(s(5).cnstr_viol(1:10,:)*100,[],1),'-','DisplayName','10^{-6}','Color',[55, 129, 100, 175]/255);
hold on
% semilogy(s(5).tvec,max(s(5).cnstr_viol(1:10,:)*100,[],1),'-','DisplayName','10^{-5}','Color',[1,0.5,0.5]);
semilogy(s(4).tvec,max(s(4).cnstr_viol(1:10,:)*100,[],1),'-','DisplayName','10^{-4}','Color',[224, 75, 0, 175]/255);
% semilogy(s(3).tvec,max(s(3).cnstr_viol(1:10,:)*100,[],1),'-','DisplayName','10^{-3}','Color',[0.5,0.5,1]);
semilogy(s(2).tvec,max(s(2).cnstr_viol(1:10,:)*100,[],1),'-','DisplayName','10^{-2}','Color',[101, 51, 113, 175]/255);
% title("Constraint violation");
leg = legend();
leg.FontSize = 30; % normal scale
% leg.FontSize = 45; % small scale
if interpreter == "tex"
    leg.Title.String = "\epsilon";
    leg.String = strrep(leg.String,'-',char(8722));
elseif interpreter == "latex"
    leg.Title.String = "$\epsilon$";
    for k = 1:length(leg.String)
        leg.String{k} = horzcat('$',leg.String{k},'$');
    end
end
leg.Position = [0.25,0.468,0.089,0.288];
leg.LineWidth = 1;
xlim([7,25.5]);
ylim([1e-2,1e+1]);
yticks([1e-1,1e1]);
if interpreter == "tex"
    xlabel('{\itt} [s]');
    ylabel('%');
else
    xlabel('$t$ [s]');
    ylabel('$\%$');
end
ax = gca;
ax.Box = "off";
ax.Units = "pixels";
ax.OuterPosition = [0, 0, pxaxwidth, pxaxheight];
ax.YTickLabel = strrep(ax.YTickLabel,'-',char(8722));

if save_figures
    file_name = "cnstr-viol-eps-sweep-notitle-latex";
    exportgraphics(fig,file_name+".pdf",'ContentType','vector');
    savefig(fig,file_name+".fig");
end

prb = s(5).prb;
tau = s(5).tau;
tvec = s(5).tvec;
tvecbar = s(5).tvecbar;

Kfine = length(tau);

nrm_T1(Kfine) = 0;
for j = 1:Kfine
    nrm_T1(j) = norm(s(5).u(1:prb.n,j));
end
nrm_Tbar1(prb.K) = 0;
for j = 1:prb.K
    nrm_Tbar1(j) = norm(s(5).ubar(1:prb.n,j));
end

% Acceleration
pxaxheight = pxaxwidth/2.5;
fig = figure('Position',[-1375,64,pxaxwidth*1.05,pxaxheight*1.05]);
plot(tvecbar,prb.Tmin*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
hold on 
plot(tvecbar,prb.Tmax*ones(1,prb.K),'-','LineWidth',4.5,'Color',[1,0.5,0.5]);
plot(tvec,nrm_T1,'-k');
plot(tvecbar,nrm_Tbar1,'.k');
if interpreter == "tex"
    ylabel("[m s^{"+char(8722)+"2}]");
    xlabel('{\itt} [s]');
else
    ylabel("[m s$^{-2}$]");
    xlabel('$t$ [s]');    
end
xlim([0,tvec(end)])
ylim([-0.25,1.1*prb.Tmax])
% title('Acceleration');
ax = gca;
ax.Box = "off";
% ax.DataAspectRatioMode = "manual";
ax.Units = "pixels";
ax.OuterPosition = [0, 0, pxaxwidth, pxaxheight];

if save_figures
    file_name = "acceleration-dyn-notitle-latex";
    exportgraphics(fig,file_name+".pdf",'ContentType','vector');
    savefig(fig,file_name+".fig");
end