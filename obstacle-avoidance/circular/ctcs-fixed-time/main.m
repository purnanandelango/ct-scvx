clearvars
clc

% cnstr_type = "exclusive-integrator-states";
cnstr_type = "single-integrator-state";

prb = problem_data(10, ...          % tau_f
                   12, ...          % K
                   025, ...         % scp_iters
                   1e2, ...         % wvc
                   10.00, ...       % wtr
                   0.01,...         % cost_factor
                   cnstr_type);


load('recent_solution','xbar','ubar');
[xbar,ubar] = misc.create_initialization(prb,2, ...
                                         xbar,ubar,[]);

[xbar,ubar,cost_val] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);
tvecbar = prb.time_grid(prb.tau,xbar,ubar);

% Simulate solution on fine grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,prb.tau(end)],prb.Kfine,prb.disc,prb.Eu2x);
tvec = prb.time_grid(tau,x,u);

r = x(1:prb.n,:);
v = x(prb.n+1:2*prb.n,:);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(r(:,end)-prb.rK),norm(v(:,end)-prb.vK));

save('recent_solution','r','v','tvec','tau','u','x','prb',...
                       'xbar','ubar','tvecbar','cost_val');

plot_solution;