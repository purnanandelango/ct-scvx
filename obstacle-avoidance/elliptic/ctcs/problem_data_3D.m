function prb = problem_data_3D(K,scp_iters,w_ep,w_px,cost_factor)
    
    prb.K = K;

    % Dimension of system
    prb.n = 3; n = prb.n;

    prb.nx = 2*n+1+1;
    prb.nu = n+1;
    
    % Generate grid in [0,1]
    prb.tau = grid.generate_grid(0,1,K,'uniform'); 
    prb.dtau = diff(prb.tau); min_dtau = min(prb.dtau);
    
    prb.h = (1/19)*min_dtau;                    % Step size for integration that computes FOH matrices
    prb.Kfine = 1+20*round(1/min_dtau);         % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.accl        = [0;0;-1];                 % External acceleration    
    prb.c_d         = 0.01;                     % Drag coefficient
    
    % Bounds

    prb.rmax        = 40;
    prb.vmax        = 05;
    prb.pmin        = 0;
    prb.pmax        = 10;

    prb.ymin        = 0;
    prb.ymax        = 0.1;

    prb.Tmin        = 1.5;
    prb.Tmax        = 08;    

    prb.ehat        = [0;0;1];
    prb.deltamax    = 85;                       % [deg] 

    prb.smin        = 01;
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

    prb.Ey = [zeros(1,2*n+1),1];

    prb.Ei      = eye(prb.nx);
    prb.zi      = [prb.r1;prb.v1;prb.p1;prb.y1];
    prb.i_idx   = 1:prb.nx;
    prb.Ef      = prb.Ei(1:2*n,:);
    prb.zf      = [prb.rK;prb.vK];
    prb.f_idx   = 1:2*prb.n;    
    prb.term_cost_vec = [zeros(2*n,1);1;0];

    % Initialization generator

    prb.x1 = [prb.r1; prb.v1; prb.p1; prb.y1];
    prb.xK = [prb.rK; prb.vK; prb.p1; prb.y1];    
    prb.u1 = [0.5*(prb.Tmax+prb.Tmin)*ones(n,1); prb.ToFguess];
    prb.uK = [0.5*(prb.Tmax+prb.Tmin)*ones(n,1); prb.ToFguess];

    % Scaling parameters
    xmin = [-0.5*prb.rmax*ones(n,1); -0.5*prb.vmax*ones(n,1); prb.pmin; prb.ymin];
    xmax = [ 0.5*prb.rmax*ones(n,1);  0.5*prb.vmax*ones(n,1); prb.pmax; prb.ymax];
    
    umin = [-prb.Tmax*ones(n,1); prb.smin]; prb.umin = umin;
    umax = [ prb.Tmax*ones(n,1); prb.smax]; prb.umax = umax;

    prb.scl_bnd = [0,1];
    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},prb.scl_bnd);
    % [Sz,cz] = misc.generate_scaling({[zeros(prb.nx,1),xmax],[zeros(prb.nu,1),umax]},prb.scl_bnd);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    prb.eps_cnstr = 1e-4;

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

    % Constraint parameters

    prb.cnstr_fun       = @(rvp,T) cnstr_scl*[-norm(prb.Hobs{1}*(rvp(1:n)-prb.qobs(:,1)))^2 + 1;            % Obstacle avoidance
                                              -norm(prb.Hobs{2}*(rvp(1:n)-prb.qobs(:,2)))^2 + 1;            % Obstacle avoidance
                                               norm(rvp(n+1:2*n))^2 - prb.vmax^2;                           % Speed upperbound
                                               norm(T) - prb.Tmax;                                          % Thrust upperbound
                                              -norm(T) + prb.Tmin;                                          % Thrust lowerbound
                                               norm(T) - secd(prb.deltamax)*prb.ehat'*T;                    % Thrust pointing
                                             ] ... 
                                             + cnstr_buffer;

    prb.cnstr_fun_jac_rvp = @(rvp,T) cnstr_scl*[-2*(rvp(1:n)-prb.qobs(:,1))'*prb.Hobs{1}'*prb.Hobs{1}, zeros(1,n),      0;
                                                -2*(rvp(1:n)-prb.qobs(:,2))'*prb.Hobs{2}'*prb.Hobs{2}, zeros(1,n),      0;
                                                 zeros(1,n),                                           2*rvp(n+1:2*n)', 0;
                                                 zeros(1,n),                                           zeros(1,n),      zeros(1,1);
                                                 zeros(1,n),                                           zeros(1,n),      zeros(1,1);
                                                 zeros(1,n),                                           zeros(1,n),      zeros(1,1);
                                             ];

    prb.cnstr_fun_jac_T = @(rvp,T) cnstr_scl*[ zeros(3,n);
                                               T'/norm(T);
                                              -T'/norm(T);
                                               T'/norm(T) - secd(prb.deltamax)*prb.ehat';
                                             ];

    % SCP parameters

    prb.disc = "ZOH";
    prb.foh_type = "v3";

    prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-7)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    % prb.solver_settings = sdpsettings('solver','quadprog','verbose',false);    
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.abstol',1e-8,'ecos.reltol',1e-8);    
    prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','mosek','verbose',false,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-9,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-9);
    % prb.solver_settings = sdpsettings('solver','osqp','verbose',false,'osqp.eps_abs',1e-9,'osqp.eps_rel',1e-9,'osqp.max_iter',1e5);        
   
    % prb.solver = struct('name',"quadprog",'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9,'Display','none');
    % prb.solver = struct('name',"piqp",'verbose',0,'eps_abs',1e-8,'eps_rel',1e-8,'eps_duality_gap_rel',1e-8,'eps_duality_gap_abs',1e-8);
    % prb.solver = struct('name',"ecos",'verbose',false,'abstol',1e-8,'reltol',1e-8);
    % prb.solver = struct('name',"gurobi",'verbose',0,'OptimalityTol',1e-9,'FeasibilityTol',1e-9);
    % prb.solver = struct('name',"scs",'eps_abs',1e-9,'eps_rel',1e-9,'verbose',false);
    % prb.solver = struct('name',"mosek",'MSK_DPAR_INTPNT_QO_TOL_PFEAS',1e-9,'MSK_DPAR_INTPNT_QO_TOL_DFEAS',1e-9,'MSK_DPAR_INTPNT_QO_TOL_REL_GAP',1e-9);
    % prb.solver = struct('name',"osqp",'eps_abs',1e-8,'eps_rel',1e-8,'verbose',0,'max_iter',5e4);
    prb.solver = struct('name',"pipg",'eps_abs',1e-12,'verbose',0,'max_iter',7e4,'rho',1.75,'lambda',1,'omega',80,'test_termination',200);
    
    % prb.px_norm = 2;
    % prb.px_norm = inf;   
    prb.px_norm = 'quad';
    
    prb.w_ep = w_ep;
    prb.w_px = w_px;
    prb.cost_factor = cost_factor;
    
    prb.eps_ep = 1e-7;
    prb.eps_px = 1e-3;

    % Time of maneuver and time grid
    prb.time_of_maneuver =     @(x,u) disc.time_of_maneuver(prb.disc,prb.tau,u(n+1,:));    
    prb.time_grid        = @(tau,x,u)        disc.time_grid(prb.disc,    tau,u(n+1,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func           = @(tau,x,u)   evaluate_dyn_func     (x,u,n,prb.c_d,prb.accl,prb.cnstr_fun);
    prb.dyn_func_linearize = @(tau,x,u)   evaluate_linearization(x,u,n,prb.c_d,prb.accl,prb.cnstr_fun,...
                                                                                        prb.cnstr_fun_jac_rvp,prb.cnstr_fun_jac_T);

end

function f = evaluate_dyn_func(x,u,n,c_d,accl,cnstr_fun)

    rvp = x(1:2*n+1);
    v = x(n+1:2*n);
    T = u(1:n);
    s = u(n+1);

    cnstr_val = cnstr_fun(rvp,T);

    F = [v;
         T + accl - c_d*norm(v)*v;
         norm(T)];

    f = s*[F;
           sum( arrayfun(@(y) max(0,y)^2,cnstr_val) )];
end

function [A,B,w] = evaluate_linearization(x,u,n,c_d,accl,cnstr_fun, ...
                                                         cnstr_fun_jac_rvp,cnstr_fun_jac_T)

    rvp = x(1:2*n+1);
    v = x(n+1:2*n);
    T = u(1:n);
    s = u(n+1);

    cnstr_val         = cnstr_fun(rvp,T);
    abs_cnstr_val     = arrayfun(@(y) max(0,y),cnstr_val);
    abs2_cnstr_val    = arrayfun(@(y) max(0,y)^2,cnstr_val);  
    cnstr_val_jac_rvp = cnstr_fun_jac_rvp(rvp,T);
    cnstr_val_jac_T   = cnstr_fun_jac_T(rvp,T);

    if norm(v) < 1e-7
        term_v = 0;
    else
        term_v = v*v'/norm(v); 
    end

    F = [v;
         T + accl - c_d*norm(v)*v;
         norm(T)];

    dFdrvp = [zeros(n),  eye(n),                      zeros(n,1);
              zeros(n), -c_d*(norm(v)*eye(n)+term_v), zeros(n,1);
              zeros(1,2*n+1)];

    dFdT = [zeros(n);
            eye(n);
            T'/norm(T)];

    A = s*[dFdrvp,                             zeros(2*n+1,1);
           2*abs_cnstr_val'*cnstr_val_jac_rvp, 0];
    
    B = [s*dFdT,                             F;
         2*s*abs_cnstr_val'*cnstr_val_jac_T, sum(abs2_cnstr_val)];
    
    f = evaluate_dyn_func(x,u,n,c_d,accl,cnstr_fun);
    
    w = f - A*x - B*u;
end