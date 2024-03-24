clearvars
clc

% prb = problem_data(06,  ...         % K
%                    120,  ...        % scp_iters
%                    2e1, ...         % wvc
%                    1.00, ...        % wtr
%                    0.01);           % cost_factor

prb = problem_data_lunar(06,  ...         % K
                         100,  ...        % scp_iters
                         5e1, ...         % wvc
                         1.00, ...        % wtr
                         0.01);           % cost_factor (0.01,0.10)

load('recent_solution','xbar','ubar','taubar');

% load('../node-only-cnstr/recent_solution_guess','xbar','ubar','taubar');

[xbar,ubar] = misc.create_initialization(prb,2, ...
                                         xbar,ubar,taubar);

% scp.diagnose_ptr_handparse(xbar,ubar,prb,@sys_cnstr_cost,'affine-var')

% [xbar,ubar] = scp.run_ptr_noparam(xbar,ubar,prb,@sys_cnstr_cost);

[xbar,ubar] = scp.run_ptr_handparse_noparam(xbar,ubar,prb);

taubar = prb.tau;
tvecbar = prb.time_grid(prb.tau,xbar,ubar);

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc,prb.ode_solver);
tvec = prb.time_grid(tau,x,u); % 

% Simulate on phyiscal time grid
[~,x2,~] = disc.simulate_dyn(xbar(:,1),{tvec,[u(1:3,:);ones(1,prb.Kfine)]},@(t,x,u) prb.dyn_func(t,x,u),[0,tvec(end)],prb.Kfine,prb.disc,prb.ode_solver);

m = x(1,:);
rI = x(2:4,:);
vI = x(5:7,:);
qBI = x(8:11,:);
omgB = x(12:14,:);

fprintf('\nFinal position error: %.3f m\nFinal velocity error: %.3f m/s\n',norm(rI(:,end)-prb.rIK),norm(vI(:,end)-prb.vIK));
fprintf("Fuel consumed: %.2f kg\n",x(1,1)-x(1,end));

save('recent_solution','m','rI','vI','qBI','omgB','tvec','tau','u','x','x2','prb',...
                       'xbar','ubar','tvecbar','taubar');

plot_solution_lunar;