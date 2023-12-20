clearvars
% clc

prb = problem_data(05,  ...         % K
                   100,  ...        % scp_iters
                   0.5, ...         % wvc
                   0.05, ...        % wtr
                   -0.0);           % cost_factor

% load("../node-only-cnstr/recent_solution.mat",'xbar','ubar','taubar');
% xbar = [xbar;zeros(1,prb.K)];

load('recent_solution','xbar','ubar','taubar');

[xbar,ubar] = misc.create_initialization(prb,1, ...
                                         xbar,ubar,taubar);

% YALMIP
% [xbar,ubar] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);

% Handparsed
% [xbar,ubar] = scp.run_ptr_dvar_handparse_noparam(xbar,ubar,prb);

% Prototype
[xbar,ubar] = run_ptr_dvar_noparam_mod(xbar,ubar,prb,@sys_cnstr_cost);

%%

taubar = prb.tau;
tvecbar = prb.time_grid(prb.tau,xbar,ubar);

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc);%,prb.ode_solver);
tvec = prb.time_grid(tau,x,u);

% Simulate on phyiscal time grid
[~,x2,~] = disc.simulate_dyn(xbar(:,1),{tvec,[u(1:3,:);ones(1,prb.Kfine)]},@(t,x,u) prb.dyn_func(t,x,u),[0,tvec(end)],prb.Kfine,prb.disc);%,prb.ode_solver);

m = x(1,:);
rI = x(2:4,:);
vI = x(5:7,:);
qBI = x(8:11,:);
omgB = x(12:14,:);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(rI(:,end)-prb.rIK),norm(vI(:,end)-prb.vIK));

save('recent_solution','m','rI','vI','qBI','omgB','tvec','tau','u','x','x2','prb',...
                       'xbar','ubar','tvecbar','taubar');

plot_solution;