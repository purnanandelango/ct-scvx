clearvars
clc

% cnstr_type = "exclusive-integrator-states";
cnstr_type = "single-integrator-state";

prb = problem_data(10, ...          % tau_f
                   06, ...          % K
                   100, ...         % scp_iters
                   5e1, ...         % wvc
                   2.00, ...        % wtr
                   0.001,...         % cost_factor
                   cnstr_type);


load('recent_solution','xbar','ubar');
[xbar,ubar] = misc.create_initialization(prb,1, ...
                                         xbar,ubar,[]);

scp.diagnose_ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost,{"","handparse"},{"",[]})

[xbar,ubar,cost_val] = scp.ctscvx_handparse_noparam(xbar,ubar,prb);

tvecbar = prb.time_grid(prb.tau,xbar,ubar);

% Simulate solution on fine grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,prb.tau(end)],prb.Kfine,prb.disc,prb.ode_solver);
% [tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,prb.tau(end)],prb.Kfine,prb.disc,prb.Eu2x,prb.ode_solver);
% [tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,prb.tau(end)],prb.Kfine,prb.disc,prb.t_burn,prb.ode_solver);
tvec = prb.time_grid(tau,x,u);

r = x(1:prb.n,:);
v = x(prb.n+1:2*prb.n,:);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(r(:,end)-prb.rK),norm(v(:,end)-prb.vK));

save('recent_solution','r','v','tvec','tau','u','x','prb',...
                       'xbar','ubar','tvecbar','cost_val');

% plot_solution;