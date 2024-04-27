clearvars
clc

prb = problem_data(10,...
                   07, ...
                   100, ...
                   5e1, ...
                   1.00, ...
                   1.00);

load('recent_solution','xbar','ubar');
[xbar,ubar] = misc.create_initialization(prb,1, ...
                                         xbar,ubar,[]);

[xbar,ubar] = scp.ctscvx_noparam(xbar,ubar,prb,@sys_cnstr_cost);
tvecbar = prb.time_grid(prb.tau,xbar,ubar);

% Simulate solution on fine grid
[tau,x,u] = disc.simulate_dyn(xbar(:,1),{prb.tau,ubar},@(t,x,u) prb.dyn_func(t,x,u),prb.tau([1,end]),prb.Kfine,prb.disc,prb.Eu2x);
tvec = prb.time_grid(tau,x,u);

r = x(1:prb.n,:);
v = x(prb.n+1:2*prb.n,:);

fprintf('\nFinal position error: %.3f\nFinal velocity error: %.3f\n',norm(r(:,end)-prb.rK),norm(v(:,end)-prb.vK));

save('recent_solution','r','v','tvec','tau','u','x','prb',...
                       'xbar','ubar','tvecbar');

plot_solution;