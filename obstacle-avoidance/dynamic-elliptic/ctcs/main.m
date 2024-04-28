clearvars
clc

prb = problem_data(10, ...          % K
                   100, ...         % scp_iters
                   2e1, ...         % w_ep
                   1.00, ...        % w_px
                   0.003);           % cost_factor

load('recent_solution','xbar','ubar','taubar');
[xbar,ubar] = misc.create_initialization(prb,1, ...
                                         xbar,ubar,taubar);

% Initial guess for SCvxGEN
writematrix(xbar',"initialguess.csv")
writematrix((diag([1,1,1/(prb.K-1)])*ubar)',"initialguess.csv","WriteMode","append")

% scp.diagnose_ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost,{"","handparse"},{"",[]})
% scp.diagnose_ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost,{"","handparse"},{"dvar_","handparse"})

% [xbar,ubar] = scp.ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost);
[xbar,ubar] = scp.ctscvx_handparse_noparam(xbar,ubar,prb);

tvecbar = prb.time_grid(prb.tau,xbar,ubar);
taubar = prb.tau;

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1), {prb.tau,ubar}, @(t,x,u) prb.dyn_func(t,x,u), [0,1], prb.Kfine, prb.disc);%, prb.ode_solver);
tvec = prb.time_grid(tau,x,u);

t_grid = linspace(0, tvec(end), prb.Kfine);
x_grid = interp1(tvec',x',t_grid')';
u_grid = interp1(tvec',u',t_grid')';

r = x(1:prb.n,:);
v = x(prb.n+1:2*prb.n,:);
p = x(2*prb.n+1,:);
t = x(2*prb.n+2,:);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(r(:,end)-prb.rK),norm(v(:,end)-prb.vK));

if prb.dyn_obs
    file_name = "recent_solution";
else
    file_name = "recent_solution_static";
end

save(file_name,'r','v','x','u','tvec','tau', ...
                       'prb', ...
                       'xbar','ubar','tvecbar','taubar',...
                       't_grid', 'x_grid', 'u_grid');

% Compare against SCvxGEN solution
solution_file = "../SCvxGEN/data/solution.csv";
if exist(solution_file)
    xx = readmatrix(solution_file,"Range",[1 1 prb.K prb.nx])';
    uu = diag([1,1,prb.K-1])*readmatrix("solution.csv","Range",[prb.K+1 1 2*prb.K prb.nu])';
    norm(xx - xbar)/norm(xbar)
    norm(uu - ubar)/norm(ubar)
end

% plot_solution;