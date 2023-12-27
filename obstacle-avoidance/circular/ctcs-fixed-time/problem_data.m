function prb = problem_data(tau_f,K,scp_iters,wvc,wtr,cost_factor,cnstr_type)
    
    prb.K = K;

    % Number of state and input constraints embedded in dynamics
    prb.m = 3;

    % Dimension of double integrator
    prb.n = 2;
    % nx is defined later
    prb.nu = prb.n;
    
    prb.tau = grid.generate_grid(0,tau_f,K,'uniform');  % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    prb.h = (1/999)*min(prb.dtau);                   % Step size for integration that computes FOH matrices
    prb.Kfine = 1+1e3*round(1/min(prb.dtau));        % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.accl = [0;0];
    
    prb.c_d = 0.01; 
    
    % Bounds

    prb.rmax = 40;
    prb.vmax = 10;
    prb.pmin = 0;    

    % FOH/ZOH
    % prb.pmax = 20;
    % prb.Tmax = 07;       

    % FBP
    prb.pmax = 200;
    prb.Tmax = 90;

    % Impulse
    % prb.pmax = 20;
    % prb.Tmax = 10;    

    % Obstacle avoidance
    prb.nobs = 2;

    prb.robs = [-5 -10;
                 6  20];
    prb.qobs = [6 7];

    % Boundary conditions

    prb.r1 = [0;0];
    prb.v1 = 0.1*[1;0];
    prb.p1 = 0;
    
    prb.rK = [-15;28];
    prb.vK = 0.1*[-1;0];    

    prb.cnstr_type = cnstr_type;
    if cnstr_type == "exclusive-integrator-states"
        prb.nx              = 2*prb.n+1+prb.m; 

        prb.ymin            = 0*ones(prb.m,1);
        prb.ymax            = 1*ones(prb.m,1);        
        
        prb.y1              = zeros(prb.m,1);            
        prb.yK              = zeros(prb.m,1);   
        
        prb.eps_cnstr       = 1e-4;

        prb.Ey              = [zeros(prb.m,2*prb.n+1),eye(prb.m)];        
        prb.term_cost_vec   = [zeros(2*prb.n,1);1;zeros(prb.m,1)]; 

        prb.Eu2x            = [zeros(prb.n);
                               eye(prb.n);
                               zeros(1+prb.m,prb.n)];        
    elseif cnstr_type == "single-integrator-state"
        prb.nx              = 2*prb.n+1+1;         
        
        prb.ymin            = 0;
        prb.ymax            = 1;        
        
        prb.y1              = 0;    
        prb.yK              = 0;
        
        prb.eps_cnstr       = 1e-5;

        prb.Ey              = [zeros(1,2*prb.n+1),1];        
        prb.term_cost_vec   = [zeros(2*prb.n,1);1;0];        

        prb.Eu2x            = [zeros(prb.n);
                               eye(prb.n);
                               zeros(2,prb.n)];        
    end

    prb.Ei = eye(prb.nx);
    prb.zi = [prb.r1;prb.v1;prb.p1;prb.y1];
    prb.Ef = prb.Ei(1:2*prb.n,:);
    prb.zf = [prb.rK;prb.vK];    

    % Straight-line initialization

    prb.x1 = [prb.r1;prb.v1;prb.p1;prb.y1];
    prb.xK = [prb.rK;prb.vK;prb.p1;prb.yK];    
    prb.u1 = 0.5*prb.Tmax*ones(prb.n,1);
    prb.uK = 0.5*prb.Tmax*ones(prb.n,1);

    % Scaling parameters

    xmin = [-0.5*prb.rmax*ones(prb.n,1); -0.5*prb.vmax*ones(prb.n,1); prb.pmin; prb.ymin];
    xmax = [ 0.5*prb.rmax*ones(prb.n,1);  0.5*prb.vmax*ones(prb.n,1); prb.pmax; prb.ymax];
    
    umin = -prb.Tmax*ones(prb.n,1); prb.umin = umin;
    umax =  prb.Tmax*ones(prb.n,1); prb.umax = umax;
    
    prb.scl_bnd = [0,1];
    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},prb.scl_bnd);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % Constraint parameters

    cnstr_scl = diag([ ...
                      1.0;
                      5.0;
                      1.0;
                      ]);
    % cnstr_scl = eye(prb.m);

    cnstr_buffer = [...
                    0.0;
                    0.0;
                    0.0;
                    ];

    prb.cnstr_fun       = @(rvp,T) cnstr_scl\[-norm(rvp(1:prb.n)-prb.robs(:,1)) + prb.qobs(1);
                                              -norm(rvp(1:prb.n)-prb.robs(:,2)) + prb.qobs(2);
                                               norm(rvp(prb.n+1:2*prb.n))^2 - prb.vmax^2;
                                             ] + cnstr_buffer;

    prb.cnstr_fun_jac_rvp = @(rvp,T) cnstr_scl\[-(rvp(1:prb.n)-prb.robs(:,1))'/norm(rvp(1:prb.n)-prb.robs(:,1)) zeros(1,prb.n+1);
                                                -(rvp(1:prb.n)-prb.robs(:,2))'/norm(rvp(1:prb.n)-prb.robs(:,2)) zeros(1,prb.n+1);
                                                 zeros(1,prb.n) 2*rvp(prb.n+1:2*prb.n)' 0];

    prb.cnstr_fun_jac_T = @(rvp,T) cnstr_scl\[ zeros(1,prb.n);
                                               zeros(1,prb.n);
                                               zeros(1,prb.n)];

    % SCP parameters

    % prb.disc = "FOH";
    % prb.foh_type = "v3_parallel";

    % prb.disc = "ZOH";
    % prb.foh_type = "v3";    

    % prb.disc = "Impulse";
    % prb.impulse_type = "v3";    

    prb.disc = "FBP";
    prb.fbp_type = "v3_parallel";
    prb.t_burn = 0.05*prb.dtau(1);

    prb.scp_iters = scp_iters; % Maximum SCP iterations

    prb.ode_solver = {'ode45',odeset('AbsTol',1e-7,'RelTol',1e-5)};

    % prb.solver_settings = sdpsettings('solver','quadprog','verbose',false);    
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.abstol',1e-8,'ecos.reltol',1e-8);    
    prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
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
    prb.epstr = 1e-3;

    % Time grid and time of manuever
    prb.time_of_maneuver = @(x,u)     tau_f;    
    prb.time_grid        = @(tau,x,u) linspace(0,tau_f,length(tau));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func           = @(tau,x,u)      evaluate_dyn_func(x,u,prb.n,prb.c_d,prb.accl,prb.cnstr_fun,cnstr_type);
    prb.dyn_func_linearize = @(tau,x,u) evaluate_linearization(x,u,prb.n,prb.c_d,prb.accl,prb.cnstr_fun,...
                                                                                          prb.cnstr_fun_jac_rvp,prb.cnstr_fun_jac_T,cnstr_type);

end

function dx = evaluate_dyn_func(x,u,n,c_d,accl,cnstr_fun,cnstr_type)

    rvp         = x(1:2*n+1);
    T           = u(1:n);
    cnstr_val   = cnstr_fun(rvp,T);
    drvp        = [plant.doubleint.dyn_func(rvp(1:end-1),T,1.0,n,c_d,accl);
                   T'*T];

    if cnstr_type == "exclusive-integrator-states"
        dy      = (arrayfun(@(z) max(0,z),cnstr_val) .^ 2);
    elseif cnstr_type == "single-integrator-state" 
        dy      = sum( arrayfun(@(z) max(0,z),cnstr_val) .^ 2 );
    else
        error("Invalid constraint formulation");
    end

    dx          = [drvp;
                   dy];
end

function [A,B,w] = evaluate_linearization(x,u,n,c_d,accl,cnstr_fun,cnstr_fun_jac_rv,cnstr_fun_jac_T,cnstr_type)

    rvp                 = x(1:2*n+1);
    T                   = u(1:n);
    [dFrv_drv,dFrv_dT]  = plant.doubleint.compute_linearization(rvp(1:end-1),T,1.0,n,c_d,accl);
    dF_drvp             = [dFrv_drv zeros(2*n,1);
                           zeros(1,2*n+1)];
    dF_dT               = [dFrv_dT;
                           2*T'];

    cnstr_val           = cnstr_fun(rvp,T);
    cnstr_val_jac_rvp   = cnstr_fun_jac_rv(rvp,T);
    cnstr_val_jac_T     = cnstr_fun_jac_T(rvp,T);

    m                   = length(cnstr_val);
    absg2               = zeros(m,1);
    absg2_jac_rvp       = zeros(m,2*n+1);
    absg2_jac_T         = zeros(m,n);
    for j = 1:m
        absg2(j)            = max(0,cnstr_val(j))^2;
        absg2_jac_rvp(j,:)  = 2*max(0,cnstr_val(j))*cnstr_val_jac_rvp(j,:);
        absg2_jac_T(j,:)    = 2*max(0,cnstr_val(j))*cnstr_val_jac_T(j,:);
    end

    if cnstr_type == "exclusive-integrator-states"
        A = [dF_drvp       zeros(2*n+1,m);
             absg2_jac_rvp zeros(m,m)];
        B = [dF_dT;
             absg2_jac_T];    
    elseif cnstr_type == "single-integrator-state"
        A = [dF_drvp              zeros(2*n+1,1);
             sum(absg2_jac_rvp,1) zeros(1,1)];
        B = [dF_dT;
             sum(absg2_jac_T,1)];    
    else
        error("Invalid constraint formulation");
    end

    dx = evaluate_dyn_func(x,u,n,c_d,accl,cnstr_fun,cnstr_type);
    w = dx - A*x - B*u;
end