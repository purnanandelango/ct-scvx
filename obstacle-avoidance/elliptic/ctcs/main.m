clearvars
% clc

% prb = problem_data_2D(05, ...          % K
%                       100, ...         % scp_iters
%                       4e1, ...         % wvc
%                       1.00, ...        % wtr
%                       0.01);           % cost_factor

prb = problem_data_3D(05, ...          % K
                      010, ...         % scp_iters
                      4e1, ...         % wvc
                      2.00, ...        % wtr
                      0.70);           % cost_factor


load('recent_solution','xbar','ubar','taubar');
[xbar,ubar] = misc.create_initialization(prb,1, ...
                                         xbar,ubar,taubar);

% scp.diagnose_ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost,{"","handparse"},{"",[]})
% scp.diagnose_ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost,{"dvar_","handparse"},{"dvar_",[]})
scp.diagnose_ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost,{"",[]},{"dvar_",[]})

% [xbar,ubar] = scp.ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost);
[xbar,ubar] = scp.ctscvx_dvar_handparse_noparam(xbar,ubar,prb);

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

plot_solution;