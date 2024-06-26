clearvars
clc

% prb = problem_data(06,  ...         % K
%                    120,  ...        % scp_iters
%                    2e1, ...         % w_ep
%                    1.00, ...        % w_px
%                    0.01);           % cost_factor

prb = problem_data_lunar(05,  ...        % K
                         300,  ...       % scp_iters
                         1e1, ...        % w_ep
                         1.00, ...       % w_px
                         0.005);         % cost_factor (0.005,0.050)

load('recent_solution','xbar','ubar','taubar');

% load('../node-only-cnstr/recent_solution_guess','xbar','ubar','taubar');

[xbar,ubar] = misc.create_initialization(prb,1, ...
                                         xbar,ubar,taubar);

% Initial guess for SCvxGEN
writematrix(xbar',"initialguess.csv")
writematrix((diag([1,1,1,1/(prb.K-1)])*ubar)',"initialguess.csv","WriteMode","append")

% scp.diagnose_ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost,{"",[]},{"","handparse"})
% scp.diagnose_ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost,{"","handparse"},{"dvar_","handparse"})
% scp.diagnose_ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost,{"",[]},{"dvar_",[]})

[xbar,ubar] = scp.ctscvx_handparse_noparam(xbar,ubar,prb);
% [xbar,ubar] = scp.ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost);

taubar = prb.tau;
tvecbar = prb.time_grid(prb.tau,xbar,ubar);

% Simulate solution on fine grid

% Simulate on [0,1] grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),[0,1],prb.Kfine,prb.disc);%,prb.ode_solver);
tvec = prb.time_grid(tau,x,u);

m = x(1,:);
rI = x(2:4,:);
vI = x(5:7,:);
qBI = x(8:11,:);
omgB = x(12:14,:);

fprintf('\nFinal position error: %.3f m\nFinal velocity error: %.3f m/s\n',norm(rI(:,end)-prb.rIK),norm(vI(:,end)-prb.vIK));
fprintf("Fuel consumed: %.2f kg\n",x(1,1)-x(1,end));

save('recent_solution','m','rI','vI','qBI','omgB','tvec','tau','u','x','prb',...
                       'xbar','ubar','tvecbar','taubar');

% Compare against SCvxGEN solution
solution_file = "../SCvxGEN/data/solution.csv";
if exist(solution_file)
    xx = readmatrix(solution_file,"Range",[1 1 prb.K prb.nx])';
    uu = diag([1,1,1,prb.K-1])*readmatrix("solution.csv","Range",[prb.K+1 1 2*prb.K prb.nu])';
    norm(xx - xbar)/norm(xbar)
    norm(uu - ubar)/norm(ubar)
end

% plot_solution_lunar;