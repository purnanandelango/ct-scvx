function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    % Dimension of system
    prb.n = 2; n = prb.n;

    prb.nx = 2*n+1;
    prb.nu = n+1;

    nx = prb.nx;
    nu = prb.nu;
    N = K;

    % Initialize PIPG solution
    prb.dx = zeros(nx, N);
    prb.du = zeros(nu, N);
    prb.phi = zeros(nx, N-1);
    prb.psi = zeros(nx, N-1);
    prb.w = zeros(nx, N-1);
    prb.v = zeros(1, N-1);
    
    % Generate grid in [0,1]
    prb.tau = grid.generate_grid(0,1,K,'uniform');
    prb.dtau = diff(prb.tau); min_dtau = min(prb.dtau);
    
    prb.h = (1/199)*min_dtau; % (1/10)*min_dtau;                    % Step size for integration that computes FOH matrices
    prb.Kfine = 1+50*round(1/min_dtau);         % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.mass        = 0.7;                      % Vehicle mass

    % Bounds

    prb.rmax        = 30;
    prb.vmax        = 7;

    prb.ymin        = 0;
    prb.ymax        = 0.1;

    prb.umax        = 7;
    prb.umin        = 1;

    % prb.ehat        = [0;1];
    % prb.deltamax    = 85;                       % [deg] 

    prb.ToFmax      = 20;
    prb.ToFguess    = 15;

    prb.dtmin       = 0.2*prb.ToFguess/(K-1);
    prb.dtmax       = 2.0*prb.ToFguess/(K-1);

    % Obstacle avoidance

    prb.nobs        = 2;

    % Centers
    prb.qobs        = [-5   5;                  
                        6  20];

    % Shape matrices
    prb.Hobs        = {1.5*diag([0.1,0.3])*geom.rot_mat_2D(0), ...
                       1.5*diag([0.3,0.1])*geom.rot_mat_2D(90)};  

    % Boundary conditions

    prb.r1 = [0;0];
    prb.v1 = [0.001;0.001];
    prb.p1 = 0;
    prb.y1 = 0;
    
    prb.rK = [0;30];
    prb.vK = [0.001;0.001];

    % Initialization generator

    prb.x1 = [prb.r1; prb.v1; prb.y1];
    prb.xK = [prb.rK; prb.vK; prb.y1];    
    prb.u1 = [0.5*(prb.umax+prb.umin)*ones(n,1); prb.ToFguess];
    prb.uK = [0.5*(prb.umax+prb.umin)*ones(n,1); prb.ToFguess];
    % prb.u1 = [1; 1; prb.ToFguess];
    % prb.uK = [4; 4; prb.ToFguess];

    % Scaling parameters
    xmin = zeros(prb.nx, 1);
    xmax = [prb.rmax*ones(n,1); prb.vmax*ones(n,1); prb.ymax];
    
    umin = zeros(prb.nu, 1);
    umax = [prb.umax*ones(n,1); (prb.K-1)*prb.dtmax];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % Scaling matrices
    Pr = prb.Sx(1:prb.n,1:prb.n);
    Pv = prb.Sx(prb.n+1:2*prb.n,prb.n+1:2*prb.n);
    Py = prb.Sx(2*prb.n+1,2*prb.n+1);
    Pu = prb.Su(1:prb.n,1:prb.n);
    Ps = prb.Su(prb.n+1,prb.n+1);

    % Include the scaling matrices in the struct for convenience
    prb.Pr = Pr;
    prb.Pv = Pv;
    prb.Py = Py;
    prb.Pu = Pu;
    prb.Ps = Ps;

    nc = 5; % number of constraints
    cnstr_indices = [4, 5]; % indices of the constraints to be imposed

    % 0: no constraints
    % 1: obstacle 1
    % 2: obstacle 2
    % 3: maximum speed
    % 4: acceleration upper bound
    % 5: acceleration lower bound

    I_nc = eye(nc+1);
    if cnstr_indices == 0
        cnstr_selector = I_nc(nc+1, :);
    else
        cnstr_selector = [];
        for i = cnstr_indices
            cnstr_selector = [cnstr_selector; I_nc(i, :)];
        end
    end

    cnstr_scl = diag([ ...
                      1; ...
                      1; ...
                      1/(sqrt(3) * 7); ...
                      1/(sqrt(3) * 7); ...
                      1/(sqrt(3) * 1); ...
                      0
                      ]);

    cnstr_buffer = cnstr_selector * [ ...
                                        0;
                                        0;
                                        0;
                                        0;
                                        0;
                                        0
                                        ];

    % Path constraints

    prb.cnstr_fun       = @(x,u) cnstr_selector * cnstr_scl * [ ...
                                                                -norm(prb.Hobs{1}*(x(1:n)-prb.qobs(:,1)))^2 + 1;                          % Obstacle avoidance
                                                                -norm(prb.Hobs{2}*(x(1:n)-prb.qobs(:,2)))^2 + 1;                          % Obstacle avoidance
                                                                 norm(x(n+1:2*n)) - prb.vmax;                                         % Speed upperbound
                                                                 norm(u(1:n)) - prb.umax;                                                 % Thrust upperbound
                                                                -norm(u(1:n)) + prb.umin;                                                 % Thrust lowerbound
                                                                 0
                                                                 ] ... 
                                                                 + cnstr_buffer;

    prb.cnstr_fun_jac_x = @(x,u) cnstr_selector * cnstr_scl * [ ...
                                                                -2*(x(1:n)-prb.qobs(:,1))'*prb.Hobs{1}'*prb.Hobs{1}, zeros(1,n),    0;
                                                                -2*(x(1:n)-prb.qobs(:,2))'*prb.Hobs{2}'*prb.Hobs{2}, zeros(1,n),    0;
                                                                 zeros(1,n),                                         x(n+1:2*n)'/norm(x(n+1:2*n)), 0;
                                                                 zeros(1,n),                                         zeros(1,n),    zeros(1,1);
                                                                 zeros(1,n),                                         zeros(1,n),    zeros(1,1);
                                                                 zeros(1,n),                                         zeros(1,n),    zeros(1,1);
                                                                 ];

    prb.cnstr_fun_jac_u = @(x,u) cnstr_selector * cnstr_scl * [ ...
                                                                 zeros(1,n);
                                                                 zeros(1,n);
                                                                 zeros(1,n);
                                                                 u(1:n)'/norm(u(1:n));
                                                                -u(1:n)'/norm(u(1:n));
                                                                 zeros(1,n);
                                                                 ];

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3"; % "v3";
    % prb.ode_solver = {'ode45',odeset('RelTol',1e-7,'AbsTol',1e-9)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    % prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','mosek','verbose',false,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-9,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-9);
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',false);
    prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.abstol',1e-10,'ecos.reltol',1e-10,'ecos.feastol',1e-10);
    % prb.solver_settings = sdpsettings('solver','quadprog','verbose',false,'quadprog.OptimalityTolerance',1e-9);
    % prb.solver_settings = sdpsettings('solver','osqp','verbose',false,'osqp.eps_abs',1e-10,'osqp.eps_rel',1e-10);

    % prb.tr_norm = 2;
    % prb.tr_norm = inf;   
    prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-8;
    prb.epstr = 5e-3;

    % Takes in unscaled data
    prb.time_of_maneuver = @(z,u) disc.time_of_maneuver(prb.disc,prb.tau,u(n+1,:));    
    prb.time_grid = @(tau,z,u) disc.time_grid(prb.disc,tau,u(n+1,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(tau,xtil,util)             evaluate_dyn_func(xtil,util,n,prb.cnstr_fun);
    prb.dyn_func_linearize = @(tau,xtil,util)   evaluate_linearization(xtil,util,n,prb.cnstr_fun,...
                                                                       prb.cnstr_fun_jac_x,prb.cnstr_fun_jac_u);

end

function F = evaluate_dyn_func(xtil,util,n,cnstr_fun)

    x = xtil(1:2*n+1);
    v = x(n+1:2*n);
    u = util(1:n);
    s = util(n+1);

    cnstr_val = cnstr_fun(x,u);

    f = [v;
         u];

    F = s*[f;
           sum( arrayfun(@(y) max(0,y),cnstr_val) .^ 2 )];
end

function [A,B,w] = evaluate_linearization(xtil,util,n,cnstr_fun, ...
                                          cnstr_fun_jac_x,cnstr_fun_jac_u)

    x = xtil(1:2*n+1);
    v = x(n+1:2*n);
    u = util(1:n);
    s = util(n+1);

    cnstr_val = cnstr_fun(x,u);
    abs_cnstr_val = arrayfun(@(y) max(0,y),cnstr_val);
    cnstr_val_jac_x = cnstr_fun_jac_x(x,u);
    cnstr_val_jac_u = cnstr_fun_jac_u(x,u);

    f = [v;
         u];

    dfdx = [zeros(n), eye(n);
            zeros(n), zeros(n)];

    dfdu = [zeros(n);
            eye(n)];

    A = s*[dfdx,                             zeros(2*n,1);
           2*abs_cnstr_val'*cnstr_val_jac_x];
    
    B = [s*dfdu,                             f;
         2*s*abs_cnstr_val'*cnstr_val_jac_u, sum( arrayfun(@(y) max(0,y),cnstr_val) .^ 2 )];
    
    F = evaluate_dyn_func(xtil,util,n,cnstr_fun);
    
    w = F - A*xtil - B*util;
end