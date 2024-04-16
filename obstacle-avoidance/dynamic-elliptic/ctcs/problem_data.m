function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    % Dimension of system
    prb.n = 2; n = prb.n;

    prb.nx = 2*n+1+1+1;
    prb.nu = n+1;
    
    % Generate grid in [0,1]
    prb.tau = grid.generate_grid(0,1,K,'uniform'); 
    prb.dtau = diff(prb.tau); min_dtau = min(prb.dtau);
    
    prb.h = (1/49)*min_dtau;                    % Step size for integration that computes discretization matrices
    prb.Kfine = 1+100*round(1/min_dtau);         % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.accl        = [0;0];                   % External acceleration    
    prb.c_d         = 0.01;                    % Drag coefficient
    
    % Bounds

    prb.rmax        = 40;
    prb.vmax        = 06;
    % prb.pmin        = 0;
    prb.pmax        = 100;
    % prb.tmin        = 1;
    % prb.tmax        = 50;  

    prb.ymin        = 0;
    prb.ymax        = 1;

    prb.Tmin        = 0.5;
    prb.Tmax        = 6.0;    

    prb.smin        = 01;
    prb.smax        = 60;
    prb.ToFguess    = 30;

    % Obstacle avoidance

    % Activate dynamic obstacles
    prb.dyn_obs = 1;
    
    % Centers
    prb.psi_obs     = [
                       34 -32 42 -24 34 -32  42 -24  34 -32 ;
                       20  20 10  10  0   0 -10 -10 -20 -20 ;
                      ];

    [~, prb.n_obs]  = size(prb.psi_obs);

    % Shape matrices
    prb.H_obs       = {diag([0.9/2, 1/15/2])*geom.rot_mat_2D(90), ...
                       diag([0.9/2, 1/15/2])*geom.rot_mat_2D(90), ...
                       diag([0.9/2, 1/15/2])*geom.rot_mat_2D(90), ...
                       diag([0.9/2, 1/15/2])*geom.rot_mat_2D(90), ...
                       diag([0.9/2, 1/15/2])*geom.rot_mat_2D(90), ...
                       diag([0.9/2, 1/15/2])*geom.rot_mat_2D(90), ...
                       diag([0.9/2, 1/15/2])*geom.rot_mat_2D(90), ...
                       diag([0.9/2, 1/15/2])*geom.rot_mat_2D(90), ...
                       diag([0.9/2, 1/15/2])*geom.rot_mat_2D(90), ...
                       diag([0.9/2, 1/15/2])*geom.rot_mat_2D(90), ...
                       };  

    prb.delpsi_obs  = prb.dyn_obs*[10,10,10,10,10,10,10,10,10,10];
    prb.phi_obs     = 0.5*pi*[1, 1, 0, 0, 1, 1, 0, 0, 1, 1];
    prb.freq_obs    = pi*[1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20, 1/20];

    prb.q_obs       = @(t, idx) prb.psi_obs(:,idx) + [prb.delpsi_obs(idx) * sin( prb.phi_obs(idx) + prb.freq_obs(idx) * t );
                                                      0];
    prb.q_obs_jac_t = @(t, idx) [prb.delpsi_obs(idx) * cos( prb.phi_obs(idx) + prb.freq_obs(idx) * t ) * prb.freq_obs(idx);
                                 0];

    prb.obs_inflate = 1.00;

    % Boundary conditions

    prb.r1 = [0;-28];
    prb.v1 = [0.1;0];
    prb.p1 = 0;
    prb.t1 = 0;
    prb.y1 = 0;

    
    prb.rK = [0; 28];
    prb.vK = [0.1;0];
    % prb.pK = prb.pmax;
    prb.pK = 0;
    % prb.tK = prb.ToFguess;
    prb.tK = 0;

    prb.Ey = [zeros(1,2*n+1+1),1];

    prb.Ei = eye(prb.nx);
    prb.zi = [prb.r1;prb.v1;prb.p1;prb.t1;prb.y1];
    prb.Ef = prb.Ei(1:2*n,:);
    prb.zf = [prb.rK;prb.vK];
    prb.term_cost_vec = [zeros(2*n,1);1;0;0]; 

    % Initialization generator

    prb.x1 = [prb.r1; prb.v1; prb.p1; prb.t1; prb.y1];
    prb.xK = [prb.rK; prb.vK; prb.pK; prb.tK; prb.y1];

    prb.u1 = [prb.Tmin*ones(n,1);     prb.ToFguess];
    prb.uK = [prb.Tmin*ones(n,1);     prb.ToFguess];

    % Scaling parameters
    xmin =   [-0.5*prb.rmax*ones(n,1); -0.5*prb.vmax*ones(n,1);        0;            0; prb.ymin]; 
    xmax =   [ 0.5*prb.rmax*ones(n,1);  0.5*prb.vmax*ones(n,1); prb.pmax; prb.ToFguess; prb.ymax];
    
    umin =   [-prb.Tmax*ones(n,1); 0];          prb.umin = umin;
    umax =   [ prb.Tmax*ones(n,1); prb.smax];   prb.umax = umax;

    prb.scl_bnd = [0,1];
    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},prb.scl_bnd);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    prb.eps_cnstr = 1e-5;    

    cnstr_scl = 1.0*diag([...
                         1.0*[1.0; ...
                              1.0; ...
                              1.0; ...
                              1.0; ...
                              1.0; ...
                              1.0; ...
                              1.0; ...
                              1.0; ...
                              1.0; ...
                              1.0]; ...
                          0.01; ...
                          0.1; ...
                          1.0; ...
                          ]);

    cnstr_buffer = 0.000*[1.00;
                          1.00;
                          1.00;
                          1.00;
                          1.00;
                          1.00;
                          1.00;
                          1.00;
                          1.00;
                          1.00;
                          1.00;
                          1.00;
                          1.00];

    % Constraint parameters

    for j = 1:prb.n_obs
        prb.H_obs{j} = prb.H_obs{j} / prb.obs_inflate;
    end    

    prb.cnstr_fun       = @(rvpt,T) cnstr_scl*[-norm(prb.H_obs{ 1}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 1)))^2 + 1; % Obstacle avoidance
                                               -norm(prb.H_obs{ 2}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 2)))^2 + 1; % Obstacle avoidance
                                               -norm(prb.H_obs{ 3}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 3)))^2 + 1; % Obstacle avoidance
                                               -norm(prb.H_obs{ 4}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 4)))^2 + 1; % Obstacle avoidance
                                               -norm(prb.H_obs{ 5}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 5)))^2 + 1; % Obstacle avoidance
                                               -norm(prb.H_obs{ 6}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 6)))^2 + 1; % Obstacle avoidance
                                               -norm(prb.H_obs{ 7}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 7)))^2 + 1; % Obstacle avoidance
                                               -norm(prb.H_obs{ 8}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 8)))^2 + 1; % Obstacle avoidance
                                               -norm(prb.H_obs{ 9}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 9)))^2 + 1; % Obstacle avoidance
                                               -norm(prb.H_obs{10}*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2),10)))^2 + 1; % Obstacle avoidance
                                                norm(rvpt(n+1:2*n))^2 - prb.vmax^2;                              % Speed upperbound
                                                norm(T)^2 - prb.Tmax^2;                                          % Thrust upperbound
                                               -norm(T)^2 + prb.Tmin^2;                                          % Thrust lowerbound
                                             ] ... 
                                             + cnstr_buffer;

    prb.cnstr_fun_jac_rvpt = @(rvpt,T) cnstr_scl*[-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 1))'*prb.H_obs{ 1}'*prb.H_obs{ 1}, zeros(1,n),       zeros(1,1), -(-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 1))'*prb.H_obs{ 1}'*prb.H_obs{ 1}) * prb.q_obs_jac_t(rvpt(2*n+2),  1); 
                                                  -2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 2))'*prb.H_obs{ 2}'*prb.H_obs{ 2}, zeros(1,n),       zeros(1,1), -(-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 2))'*prb.H_obs{ 2}'*prb.H_obs{ 2}) * prb.q_obs_jac_t(rvpt(2*n+2),  2); 
                                                  -2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 3))'*prb.H_obs{ 3}'*prb.H_obs{ 3}, zeros(1,n),       zeros(1,1), -(-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 3))'*prb.H_obs{ 3}'*prb.H_obs{ 3}) * prb.q_obs_jac_t(rvpt(2*n+2),  3); 
                                                  -2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 4))'*prb.H_obs{ 4}'*prb.H_obs{ 4}, zeros(1,n),       zeros(1,1), -(-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 4))'*prb.H_obs{ 4}'*prb.H_obs{ 4}) * prb.q_obs_jac_t(rvpt(2*n+2),  4); 
                                                  -2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 5))'*prb.H_obs{ 5}'*prb.H_obs{ 5}, zeros(1,n),       zeros(1,1), -(-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 5))'*prb.H_obs{ 5}'*prb.H_obs{ 5}) * prb.q_obs_jac_t(rvpt(2*n+2),  5); 
                                                  -2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 6))'*prb.H_obs{ 6}'*prb.H_obs{ 6}, zeros(1,n),       zeros(1,1), -(-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 6))'*prb.H_obs{ 6}'*prb.H_obs{ 6}) * prb.q_obs_jac_t(rvpt(2*n+2),  6); 
                                                  -2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 7))'*prb.H_obs{ 7}'*prb.H_obs{ 7}, zeros(1,n),       zeros(1,1), -(-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 7))'*prb.H_obs{ 7}'*prb.H_obs{ 7}) * prb.q_obs_jac_t(rvpt(2*n+2),  7); 
                                                  -2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 8))'*prb.H_obs{ 8}'*prb.H_obs{ 8}, zeros(1,n),       zeros(1,1), -(-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 8))'*prb.H_obs{ 8}'*prb.H_obs{ 8}) * prb.q_obs_jac_t(rvpt(2*n+2),  8); 
                                                  -2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 9))'*prb.H_obs{ 9}'*prb.H_obs{ 9}, zeros(1,n),       zeros(1,1), -(-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2), 9))'*prb.H_obs{ 9}'*prb.H_obs{ 9}) * prb.q_obs_jac_t(rvpt(2*n+2),  9); 
                                                  -2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2),10))'*prb.H_obs{10}'*prb.H_obs{10}, zeros(1,n),       zeros(1,1), -(-2*(rvpt(1:n)-prb.q_obs(rvpt(2*n+2),10))'*prb.H_obs{10}'*prb.H_obs{10}) * prb.q_obs_jac_t(rvpt(2*n+2), 10); 
                                                   zeros(1,n),                                                            2*rvpt(n+1:2*n)', zeros(1,1),  zeros(1,1);
                                                   zeros(1,n),                                                            zeros(1,n),       zeros(1,1),  zeros(1,1);
                                                   zeros(1,n),                                                            zeros(1,n),       zeros(1,1),  zeros(1,1);
                                                  ];

    prb.cnstr_fun_jac_T = @(rvp,T) cnstr_scl*[ zeros(10+1,n);
                                               2*T';
                                              -2*T';
                                              ];

    for j = 1:prb.n_obs
        prb.H_obs{j} = prb.H_obs{j} * prb.obs_inflate;
    end        

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3_parallel";

    % prb.ode_solver = {'ode45',odeset('RelTol',1e-3,'AbsTol',1e-4)};
    % prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-6)};    
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    % prb.solver_settings = sdpsettings('solver','quadprog','verbose',false);    
    prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.abstol',1e-8,'ecos.reltol',1e-8);    
    % prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','mosek','verbose',false,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-9,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-9);
    % prb.solver_settings = sdpsettings('solver','osqp','verbose',false,'osqp.eps_abs',1e-8,'osqp.eps_rel',1e-8,'osqp.max_iter',5e4);        
   
    % prb.solver = struct('name',"quadprog",'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9,'Display','none');
    prb.solver = struct('name',"piqp",'verbose',0,'eps_abs',1e-8,'eps_rel',1e-8,'eps_duality_gap_rel',1e-8,'eps_duality_gap_abs',1e-8);
    % prb.solver = struct('name',"ecos",'verbose',false,'abstol',1e-8,'reltol',1e-8);
    % prb.solver = struct('name',"gurobi",'verbose',0,'OptimalityTol',1e-9,'FeasibilityTol',1e-9);
    % prb.solver = struct('name',"scs",'eps_abs',1e-9,'eps_rel',1e-9,'verbose',false);
    % prb.solver = struct('name',"mosek",'MSK_DPAR_INTPNT_QO_TOL_PFEAS',1e-9,'MSK_DPAR_INTPNT_QO_TOL_DFEAS',1e-9,'MSK_DPAR_INTPNT_QO_TOL_REL_GAP',1e-9);
    % prb.solver = struct('name',"osqp",'eps_abs',1e-8,'eps_rel',1e-8,'verbose',0,'max_iter',5e4);
    % prb.solver = struct('name',"pipg",'eps_abs',1e-9,'verbose',0,'max_iter',5e4,'rho',1.75,'lambda',1,'omega',20,'test_termination',500);

    % prb.tr_norm = 2;
    % prb.tr_norm = inf;   
    prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-7;
    prb.epstr = 7.5e-4;

    % Time of maneuver and time grid
    prb.time_of_maneuver =     @(x,u) x(end-1,end);    
    prb.time_grid        = @(tau,x,u) disc.time_grid(prb.disc,    tau,u(n+1,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(tau,x,u)             evaluate_dyn_func     (x,u,n,prb.c_d,prb.accl,prb.cnstr_fun);
    prb.dyn_func_linearize = @(tau,x,u)   evaluate_linearization(x,u,n,prb.c_d,prb.accl,prb.cnstr_fun,...
                                                                                        prb.cnstr_fun_jac_rvpt,prb.cnstr_fun_jac_T);

end

function f = evaluate_dyn_func(x,u,n,c_d,accl,cnstr_fun)

    rvpt = x(1:2*n+2);
    v = x(n+1:2*n);
    % t = x(2*n+2);
    T = u(1:n);
    s = u(n+1);

    cnstr_val = cnstr_fun(rvpt,T);

    F = [v;
         T + accl - c_d*norm(v)*v;
         norm(T)^2;
         1];

    f = s*[F;
           sum( arrayfun(@(y) max(0,y)^2,cnstr_val) )];
end

function [A,B,w] = evaluate_linearization(x,u,n,c_d,accl,cnstr_fun, ...
                                                         cnstr_fun_jac_rvpt,cnstr_fun_jac_T)

    rvpt = x(1:2*n+2);
    % t = x(2*n+2);
    v = x(n+1:2*n);
    T = u(1:n);
    s = u(n+1);

    cnstr_val         = cnstr_fun(rvpt,T);
    abs_cnstr_val     = arrayfun(@(y) max(0,y),cnstr_val);
    abs2_cnstr_val    = arrayfun(@(y) max(0,y)^2,cnstr_val);
    cnstr_val_jac_rvp = cnstr_fun_jac_rvpt(rvpt,T);
    cnstr_val_jac_T   = cnstr_fun_jac_T(rvpt,T);

    if norm(v) < 1e-7
        term_v = 0;
    else
        term_v = v*v'/norm(v); 
    end

    F = [v;
         T + accl - c_d*norm(v)*v;
         norm(T)^2;
         1];

    dFdrvp = [zeros(n),  eye(n),                      zeros(n,1), zeros(n,1);
              zeros(n), -c_d*(norm(v)*eye(n)+term_v), zeros(n,1), zeros(n,1);
              zeros(1,2*n+1), zeros(1,1);
              zeros(1,2*n+1), zeros(1,1)];

    dFdT = [zeros(n);
            eye(n);
            2*T';
            zeros(1,n)];

    A = s*[dFdrvp,                             zeros(2*n+2,1); 
           2*abs_cnstr_val'*cnstr_val_jac_rvp, 0];
    
    B = [s*dFdT,                             F;
         2*s*abs_cnstr_val'*cnstr_val_jac_T, sum(abs2_cnstr_val)];
    
    f = evaluate_dyn_func(x,u,n,c_d,accl,cnstr_fun);
    
    w = f - A*x - B*u;
end