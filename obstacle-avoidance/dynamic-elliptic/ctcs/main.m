clearvars
% clc

prb = problem_data(10, ...          % K
                   100, ...         % scp_iters
                   1e2, ...         % wvc
                   1.00, ...        % wtr
                   0.0001);          % cost_factor

load('recent_solution','xbar','ubar','taubar');
[xbar,ubar] = misc.create_initialization(prb,2, ...
                                         xbar,ubar,taubar);

% scp.diagnose_ptr_handparse(xbar,ubar,prb,@sys_cnstr_cost,'affine-var')

% [xbar,ubar] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);
[xbar,ubar] = scp.run_ptr_handparse_noparam(xbar,ubar,prb);

tvecbar = prb.time_grid(prb.tau,xbar,ubar);
taubar = prb.tau;

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1), {prb.tau,ubar}, @(t,x,u) prb.dyn_func(t,x,u), [0,1], prb.Kfine, prb.disc, prb.ode_solver);
tvec = prb.time_grid(tau,x,u);

t_grid = linspace(0, tvec(end), prb.Kfine);
x_grid = interp1(tvec',x',t_grid')';
u_grid = interp1(tvec',u',t_grid')';

r = x(1:prb.n,:);
v = x(prb.n+1:2*prb.n,:);
p = x(2*prb.n+1,:);
t = x(2*prb.n+2,:);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(r(:,end)-prb.rK),norm(v(:,end)-prb.vK));

save('recent_solution','r','v','x','u','tvec','tau', ...
                       'prb', ...
                       'xbar','ubar','tvecbar','taubar',...
                       't_grid', 'x_grid', 'u_grid');

% plot_solution;