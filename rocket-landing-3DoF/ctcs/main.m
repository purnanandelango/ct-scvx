clearvars
clc

% Initialize from straight-line

prb = problem_data(08, ...          % K
                   84, ...          % T
                   300, ...         % scp_iters
                   2e1, ...         % wvc
                   1.00, ...        % wtr
                   0.0200);         % cost_factor (0.0200,0.0002)

load('recent_solution','xbar','ubar','taubar');
[xbar,ubar] = misc.create_initialization(prb,1, ...
                                         xbar,ubar,taubar);

% Warm-start from node-only-cnstr solution

% prb = problem_data(08, ...          % K
%                    84, ...          % T
%                    300, ...         % scp_iters
%                    2e1, ...         % wvc
%                    1.00, ...        % wtr
%                    0.003);           % cost_factor
% 
% load('../node-only-cnstr/recent_solution_cvx','x_guess','u');
% [xbar,ubar] = misc.create_initialization(prb,2, ...
%                                          x_guess,u,[]);

% scp.diagnose_ptr_handparse(xbar,ubar,prb,@sys_cnstr_cost,"affine-var")

[xbar,ubar] = scp.run_ptr_handparse_noparam(xbar,ubar,prb);

tvecbar = prb.time_grid(prb.tau,xbar,ubar);
taubar = prb.tau;

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),taubar([1,end]),prb.Kfine,prb.disc,prb.ode_solver);
tvec = prb.time_grid(tau,x,u);

r = x(1:3,:);
v = x(4:6,:);
mass = exp(x(7,:));
thrust = u(1:3,:) .* mass;
sig = u(4,:) .* mass;

cost_val = mass(end);

fprintf('\nFinal position error: %.3f m\nFinal velocity error: %.3f m/s\n',norm(r(:,end)-prb.rK),norm(v(:,end)-prb.vK));
fprintf('Fuel consumed: %.2f kg\n',mass(1)-mass(end));

save('recent_solution','r','v','mass','thrust','sig','x','u','tvec','tau', ...
                       'prb', ...
                       'xbar','ubar','tvecbar','taubar','cost_val');

plot_solution;                                         