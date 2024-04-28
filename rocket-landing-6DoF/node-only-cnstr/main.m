clearvars
% clc

% prb = problem_data(07,  ...         % K
%                    100,  ...        % scp_iters
%                    5e1, ...         % w_ep
%                    1.00, ...        % w_px
%                    0.01);           % cost_factor

prb = problem_data_lunar(05,  ...         % K
                         050,  ...        % scp_iters
                         2e1, ...         % w_ep
                         1.00, ...        % w_px
                         0.01);           % cost_factor

load('recent_solution','xbar','ubar','taubar');
[xbar,ubar] = misc.create_initialization(prb,1, ...
                                         xbar,ubar,taubar);

[xbar,ubar] = scp.ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost);
taubar = prb.tau;
tvecbar = prb.time_grid(prb.tau,xbar,ubar);

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc);%,prb.ode_solver);
tvec = prb.time_grid(tau,x,u);

m = x(1,:);
rI = x(2:4,:);
vI = x(5:7,:);
qBI = x(8:11,:);
omgB = x(12:14,:);

fprintf('\nFinal position error: %.3f m\nFinal velocity error: %.3f m/s\n',norm(rI(:,end)-prb.rIK),norm(vI(:,end)-prb.vIK));
fprintf("Fuel consumed: %.2f kg\n",x(1,1)-x(1,end));

save('recent_solution','m','rI','vI','qBI','omgB','tvec','tau','u','x','prb',...
                       'xbar','ubar','tvecbar','taubar');

xbar = [xbar;zeros(1,prb.K)];
save('recent_solution_guess','xbar','ubar','taubar');

% plot_solution_lunar;