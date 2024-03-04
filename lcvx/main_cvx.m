clc
close all
clearvars

prb = problem_data_cvx(08,48);

yalmip clear

x_scl = sdpvar(prb.nx,prb.K);
u_scl = sdpvar(prb.nu,prb.K);

x = prb.Sx*x_scl + repmat(prb.cx,[1,prb.K]);
u = prb.Su*u_scl + repmat(prb.cu,[1,prb.K]);

cnstr = [];
cost_fun = prb.dtau*sum(u(4,:));
% cost_fun = -x(7,end);
for k = 1:prb.K
    
    if k < prb.K
        cnstr = [cnstr;
                 x(:,k+1) == prb.Ad*x(:,k) + prb.Bd*u(:,k) + prb.wd;
                 prb.tp_vec'*u(:,k) >= 0;
                 norm(u(1:3,k)) <= u(4,k);
                 % prb.thrst_bnd_mat{k}*[x(7,k);u(4,k)] <= prb.thrst_bnd_vec{k};
                 u(4,k) >= prb.mu_1(k)*(1 - (x(7,k) - prb.p_0(k)) + 0.5*(x(7,k) - prb.p_0(k))^2);
                 u(4,k) <= prb.mu_2(k)*(1 - (x(7,k) - prb.p_0(k)));
                 ];
    end
    cnstr = [cnstr;
             norm(x(4:6,k)) <= prb.vmax;
             % norm(x(1:2,k)) <= cotd(prb.gam_gs)*x(3,k);         
             x(3,k) - x(3,prb.K) >= tand(prb.gam_gs)*norm(x(1:2,k)-x(1:2,prb.K));
             x(7,k) <= prb.p_1(k);
             x(7,k) >= prb.p_0(k);
            ];

end
cnstr = [cnstr;
         u(:,prb.K) == u(:,prb.K-1)];
cnstr = [cnstr;
         x(1:3,1)     == prb.r1;
         x(1:3,prb.K) == prb.rK;
         x(4:6,1)     == prb.v1;
         x(4:6,prb.K) == prb.vK;
         x(7,1)       == prb.p1;
         x(7,prb.K)   >= prb.pK;
         ];

optimize(cnstr,cost_fun,prb.solver_settings);

x = value(x);
u = value(u);
cost_val = x(7,end);

x_guess = [x;
           zeros(1,prb.K)];

[tau_sim,x_sim] = ode45(@(tau,x) prb.dyn_func(tau,x,disc.u_zoh(tau,u,prb.tau)),[0,prb.tau(end)],x(:,1),odeset('RelTol',1e-5,'AbsTol',1e-7));
x_sim = x_sim';
u_sim = disc.u_zoh(tau_sim,u,prb.tau);

save('recent_solution_cvx','tau_sim','x_sim','u_sim','x','u','prb','cost_val','x_guess');

plot_solution_cvx;