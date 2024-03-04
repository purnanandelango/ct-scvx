clearvars
clc

prb = problem_data(07, ...          % K
                   12.1, ...        % T
                   010, ...         % scp_iters
                   5e1, ...         % wvc
                   1.00, ...        % wtr
                   0.001);           % cost_factor

% load('recent_solution','xbar','ubar','taubar');
% [xbar,ubar] = misc.create_initialization(prb,2, ...
%                                          xbar,ubar,taubar);

load('recent_solution_cvx','x_guess','u');
[xbar,ubar] = misc.create_initialization(prb,2, ...
                                         x_guess,u,[]);

[xbar,ubar] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);
tvecbar = prb.time_grid(prb.tau,xbar,ubar);
taubar = prb.tau;

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),taubar([1,end]),prb.Kfine,prb.disc,prb.ode_solver);
tvec = prb.time_grid(tau,x,u);

r = x(1:3,:);
v = x(4:6,:);

cost_val = x(7,end);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(r(:,end)-prb.rK),norm(v(:,end)-prb.vK));

save('recent_solution','r','v','x','u','tvec','tau', ...
                       'prb', ...
                       'xbar','ubar','tvecbar','taubar','cost_val');

% plot_solution;                                         