function prb = problem_data_3D(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    % Dimension of system
    prb.n = 3; n = prb.n;

    prb.nx = 2*n+1+1;
    prb.nu = n+1;
    
    % Generate grid in [0,1]
    prb.tau = grid.generate_grid(0,1,K,'uniform'); 
    prb.dtau = diff(prb.tau); min_dtau = min(prb.dtau);
    
    prb.h = (1/10)*min_dtau;                    % Step size for integration that computes FOH matrices
    prb.Kfine = 1+50*round(1/min_dtau);         % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.accl        = [0;0;-1];                   % External acceleration    
    prb.c_d         = 0.01;                      % Drag coefficient
    
    % Bounds

    prb.rmax        = 40;
    prb.vmax        = 06;
    prb.pmin        = 0;
    prb.pmax        = 20;

    prb.ymin        = 0;
    prb.ymax        = 1;

    prb.umax        = 08;
    prb.umin        = 0.5;

    prb.ehat        = [0;0;1];
    prb.deltamax    = 85;                       % [deg] 

    prb.smin        = 0.1;
    prb.smax        = 15;
    prb.ToFguess    = 15;

    % Obstacle avoidance

    prb.nobs        = 2;

    % Centers
    prb.qobs        = [-5 -10;                  
                        6  20;
                       10  10];

    % Shape matrices
    prb.Hobs        = {blkdiag(diag([0.3,0.1])*geom.rot_mat_2D(0),0.3), ...
                       blkdiag(diag([0.1,0.3])*geom.rot_mat_2D(0),0.3)};  

    % Boundary conditions

    prb.r1 = [0;0;10];
    prb.v1 = [0;0;0];
    prb.p1 = 0;
    prb.y1 = 0;
    
    prb.rK = [-15;28;10];
    prb.vK = [0;0;0];

    % Initialization generator

    prb.x1 = [prb.r1; prb.v1+1; prb.p1; prb.y1];
    prb.xK = [prb.rK; prb.vK+1; prb.p1; prb.y1];    
    prb.u1 = [0.5*(prb.umax+prb.umin)*ones(n,1); prb.ToFguess];
    prb.uK = [0.5*(prb.umax+prb.umin)*ones(n,1); prb.ToFguess];

    % Scaling parameters
    % xmin = [-0.5*prb.rmax*ones(n,1); -0.5*prb.vmax*ones(n,1); prb.pmin; prb.ymin];
    xmin = zeros(prb.nx,1);
    xmax = [ 0.5*prb.rmax*ones(n,1);  0.5*prb.vmax*ones(n,1); prb.pmax; prb.ymax];
    
    % umin = [-0.5*prb.umax*ones(n,1); prb.smin+5];
    umin = zeros(prb.nu,1);
    umax = [ 0.5*prb.umax*ones(n,1); prb.smax-5];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    cnstr_scl = diag([...
                      1;
                      1;
                      0.1;
                      1; 
                      1;
                      0.1;
                      ]);
    cnstr_buffer = [0;
                    0;
                    0;
                    0;
                    0.01;
                    0;
                    ];

    % Path constraints

    prb.cnstr_fun       = @(x,u) cnstr_scl*[-norm(prb.Hobs{1}*(x(1:n)-prb.qobs(:,1)))^2 + 1;                          % Obstacle avoidance
                                            -norm(prb.Hobs{2}*(x(1:n)-prb.qobs(:,2)))^2 + 1;                          % Obstacle avoidance
                                             norm(x(n+1:2*n))^2 - prb.vmax^2;                                         % Speed upperbound
                                             norm(u(1:n)) - prb.umax;                                                 % Thrust upperbound
                                            -norm(u(1:n)) + prb.umin;                                                 % Thrust lowerbound
                                             norm(u(1:n)) - secd(prb.deltamax)*prb.ehat'*u(1:n);                      % Thrust pointing
                                             ] ... 
                                             + cnstr_buffer;

    prb.cnstr_fun_jac_x = @(x,u) cnstr_scl*[-2*(x(1:n)-prb.qobs(:,1))'*prb.Hobs{1}'*prb.Hobs{1}, zeros(1,n),    0;
                                            -2*(x(1:n)-prb.qobs(:,2))'*prb.Hobs{2}'*prb.Hobs{2}, zeros(1,n),    0;
                                             zeros(1,n),                                         2*x(n+1:2*n)', 0;
                                             zeros(1,n),                                         zeros(1,n),    zeros(1,1);
                                             zeros(1,n),                                         zeros(1,n),    zeros(1,1);
                                             zeros(1,n),                                         zeros(1,n),    zeros(1,1);
                                             ];

    prb.cnstr_fun_jac_u = @(x,u) cnstr_scl*[ zeros(3,n);
                                             u(1:n)'/norm(u(1:n));
                                            -u(1:n)'/norm(u(1:n));
                                             u(1:n)'/norm(u(1:n)) - secd(prb.deltamax)*prb.ehat';
                                             ];

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3";
    prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-7)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    % prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','mosek','verbose',false,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-9,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-9);
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.abstol',1e-8,'ecos.reltol',1e-8);
    prb.solver_settings = sdpsettings('solver','quadprog','verbose',false);
    % prb.solver_settings = sdpsettings('solver','osqp','verbose',false,'osqp.eps_abs',1e-7,'osqp.eps_rel',1e-7,'osqp.max_iter',5e4);        
   

    % prb.tr_norm = 2;
    % prb.tr_norm = inf;   
    prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-8;
    prb.epstr = 1e-3;

    % Takes in unscaled data
    prb.time_of_maneuver = @(z,u) disc.time_of_maneuver(prb.disc,prb.tau,u(n+1,:));    
    prb.time_grid = @(tau,z,u) disc.time_grid(prb.disc,tau,u(n+1,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(tau,xtil,util)             evaluate_dyn_func     (xtil,util,n,prb.c_d,prb.accl,prb.cnstr_fun);
    prb.dyn_func_linearize = @(tau,xtil,util)   evaluate_linearization(xtil,util,n,prb.c_d,prb.accl,prb.cnstr_fun,...
                                                                       prb.cnstr_fun_jac_x,prb.cnstr_fun_jac_u);

end

function F = evaluate_dyn_func(xtil,util,n,c_d,accl,cnstr_fun)

    x = xtil(1:2*n+1);
    v = x(n+1:2*n);
    u = util(1:n);
    s = util(n+1);

    cnstr_val = cnstr_fun(x,u);

    f = [v;
         u + accl - c_d*norm(v)*v;
         norm(u)];

    F = s*[f;
           sum( arrayfun(@(y) max(0,y)^2,cnstr_val) )];
end

function [A,B,w] = evaluate_linearization(xtil,util,n,c_d,accl,cnstr_fun, ...
                                          cnstr_fun_jac_x,cnstr_fun_jac_u)

    x = xtil(1:2*n+1);
    v = x(n+1:2*n);
    u = util(1:n);
    s = util(n+1);

    cnstr_val = cnstr_fun(x,u);
    abs_cnstr_val = arrayfun(@(y) max(0,y),cnstr_val);
    cnstr_val_jac_x = cnstr_fun_jac_x(x,u);
    cnstr_val_jac_u = cnstr_fun_jac_u(x,u);

    if norm(v) < 1e-7
        term_v = 0;
    else
        term_v = v*v'/norm(v); 
    end

    f = [v;
         u + accl - c_d*norm(v)*v;
         norm(u)];

    dfdx = [zeros(n), eye(n), zeros(n,1);
            zeros(n), -c_d*(norm(v)*eye(n)+term_v), zeros(n,1);
            zeros(1,2*n+1)];

    dfdu = [zeros(n);
            eye(n);
            u'/norm(u)];

    A = s*[dfdx,                             zeros(2*n+1,1);
           2*abs_cnstr_val'*cnstr_val_jac_x, 0];
    
    B = [s*dfdu,                             f;
         2*s*abs_cnstr_val'*cnstr_val_jac_u, sum( arrayfun(@(y) max(0,y)^2,cnstr_val) )];
    
    F = evaluate_dyn_func(xtil,util,n,c_d,accl,cnstr_fun);
    
    w = F - A*xtil - B*util;
end