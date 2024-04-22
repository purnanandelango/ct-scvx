clc
clearvars
close all

pxaxwidth = 1110; % [pixels]
pxaxheight = pxaxwidth/2.5; % [pixels]

% Load solutions
for k = 6:-1:1
    s(k) = load("solution_eps_1e-"+num2str(k)); 
end

fig = figure('Position',[-1249,425,pxaxwidth*1.01,pxaxheight*1.01]);
semilogy(s(6).tvec,max(s(5).cnstr_viol(1:10,:)*100,[],1),'-','DisplayName','10^{-6}','Color',[0,0,0]);
hold on
% semilogy(s(5).tvec,max(s(5).cnstr_viol(1:10,:)*100,[],1),'-','DisplayName','10^{-5}','Color',[1,0.5,0.5]);
semilogy(s(4).tvec,max(s(4).cnstr_viol(1:10,:)*100,[],1),'-','DisplayName','10^{-4}','Color',[60,179,113,200]/255);
% semilogy(s(3).tvec,max(s(3).cnstr_viol(1:10,:)*100,[],1),'-','DisplayName','10^{-3}','Color',[0.5,0.5,1]);
semilogy(s(2).tvec,max(s(2).cnstr_viol(1:10,:)*100,[],1),'-','DisplayName','10^{-2}','Color',[255,99,71,175]/255);
title("Constraint violation");
leg = legend();
leg.FontSize = 30;
leg.Title.String = "\epsilon";
leg.String = strrep(leg.String,'-',char(8722));
leg.Position = [0.25,0.468,0.089,0.288];
leg.LineWidth = 1;
xlim([7,25.5]);
ylim([1e-2,1e+1]);
yticks([1e-1,1e1]);
xlabel('{\itt} [s]');
ylabel('%');
ax = gca;
ax.Box = "off";
ax.Units = "pixels";
ax.OuterPosition = [0, 0, pxaxwidth, pxaxheight];
ax.YTickLabel = strrep(ax.YTickLabel,'-',char(8722));

exportgraphics(fig,"cnstr-viol-eps-sweep.pdf","ContentType","vector");