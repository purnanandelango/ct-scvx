clearvars
clc

prb = problem_data_2D(10, ...          % K
                      30, ...          % scp_iters
                      5e1, ...         % wvc
                      1e0, ...         % wtr
                      1e-2);           % cost_factor

% prb = problem_data_3D(10, ...          % K
%                       30, ...          % scp_iters
%                       5e1, ...         % wvc
%                       1e0, ...         % wtr
%                       1e-2);           % cost_factor


load('recent_solution','xbar','ubar','taubar');
[xbar,ubar] = misc.create_initialization(prb,2, ...
                                         xbar,ubar,taubar);

[xbar,ubar] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);
tvecbar = prb.time_grid(prb.tau,xbar,ubar);
taubar = prb.tau;

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc,prb.ode_solver);
tvec = prb.time_grid(tau,x,u);

r = x(1:prb.n,:);
v = x(prb.n+1:2*prb.n,:);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(r(:,end)-prb.rK),norm(v(:,end)-prb.vK));

save('recent_solution','r','v','x','u','tvec','tau', ...
                       'prb', ...
                       'xbar','ubar','tvecbar','taubar');

% plot_solution;