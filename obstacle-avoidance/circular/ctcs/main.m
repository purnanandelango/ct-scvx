clearvars
clc

% cnstr_type = "exclusive-integrator-states";
cnstr_type = "single-integrator-state";

prb = problem_data(07, ...          % K
                   100, ...         % scp_iters
                   4e1, ...         % w_ep
                   1.00, ...        % w_px
                   0.70,...         % cost_factor
                   cnstr_type);


load('recent_solution','xbar','ubar');
[xbar,ubar] = misc.create_initialization(prb,1, ...
                                         xbar,ubar,[]);

[xbar,ubar,cost_val] = scp.ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost);
tvecbar = prb.time_grid(prb.tau,xbar,ubar);

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc);
tvec = prb.time_grid(tau,x,u);

% Simulate on phyiscal time grid
[~,x2,~] = disc.simulate_dyn(xbar(:,1),{tvec,[u(1:prb.n,:);ones(1,prb.Kfine)]},@(t,x,u) prb.dyn_func(t,x,u),[0,tvec(end)],prb.Kfine,prb.disc);

r = x(1:prb.n,:);
v = x(prb.n+1:2*prb.n,:);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(r(:,end)-prb.rK),norm(v(:,end)-prb.vK));

save('recent_solution','r','v','tvec','tau','u','x','x2','prb',...
                       'xbar','ubar','tvecbar','cost_val');

plot_solution;