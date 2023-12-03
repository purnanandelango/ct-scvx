function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor,cnstr_type)
    
    prb.K = K;

    % Number of state and input constraints embedded in dynamics
    prb.m = 4;

    % Dimension of double integrator
    prb.n = 2;
    % nx is defined later
    prb.nu = prb.n+1;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform');  % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    prb.h = (1/19)*min(prb.dtau);                   % Step size for integration that computes FOH matrices
    prb.Kfine = 1+20*round(1/min(prb.dtau));        % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.accl = [0;-1];
    
    prb.c_d = 0.1; 
    
    % Bounds

    prb.rmax = 40;
    prb.vmax = 6;
    prb.Tmax = 6;
    prb.Tmin = 2;

    prb.smin = 1;
    prb.smax = 15;
    prb.dtmin = 1;
    prb.dtmax = 3;
    prb.ToFmax = 20;
    prb.snom = [1,10];
    prb.ToFguess = 15;

    % Obstacle avoidance
    prb.nobs = 2;

    prb.robs = [-5 -10;
                 6  20];
    prb.qobs = [6 7];

    % Boundary conditions

    prb.r1 = [0;0];
    prb.v1 = 0.1*[1;0];
    
    prb.rK = [-15;28];
    prb.vK = 0.1*[-1;0];

    prb.cnstr_type = cnstr_type;
    if cnstr_type == "exclusive-integrator-states"
        prb.nx          = 2*prb.n+prb.m; 

        prb.ymin        = 0*ones(prb.m,1);
        prb.ymax        = 1*ones(prb.m,1);        
        
        prb.y1          = zeros(prb.m,1);            
        prb.yK          = zeros(prb.m,1);   
        
        prb.eps_cnstr   = 1e-4;
    elseif cnstr_type == "single-integrator-state"
        prb.nx          = 2*prb.n+1;         
        
        prb.ymin        = 0;
        prb.ymax        = 1;        
        
        prb.y1          = 0;    
        prb.yK          = 0;
        
        prb.eps_cnstr   = 5e-4;
    end

    % Straight-line initialization

    prb.x1 = [prb.r1;prb.v1;prb.y1];
    prb.xK = [prb.rK;prb.vK;prb.yK];    
    prb.u1 = [0.5*prb.Tmax*ones(prb.n,1);prb.ToFguess];
    prb.uK = [0.5*prb.Tmax*ones(prb.n,1);prb.ToFguess];

    % Scaling parameters

    xmin = [-0.5*prb.rmax*ones(prb.n,1); -0.5*prb.vmax*ones(prb.n,1); prb.ymin];
    xmax = [ 0.5*prb.rmax*ones(prb.n,1);  0.5*prb.vmax*ones(prb.n,1); prb.ymax];
    
    umin = [zeros(prb.n,1);         prb.snom(1)];
    umax = [prb.Tmax*ones(prb.n,1); prb.snom(2)];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % Constraint parameters

    % cnstr_scl = diag([ 0.25*prb.rmax, ...
    %                    0.25*prb.rmax, ...
    %                    1.00*prb.vmax, ...
    %                    1.00*prb.Tmax]);
    cnstr_scl = eye(prb.m);

    cnstr_buffer = [...
                    0.01;
                    0.01;
                    0.0;
                    0.05;
                    ];

    prb.cnstr_fun       = @(rv,T) cnstr_scl\[-norm(rv(1:prb.n)-prb.robs(:,1)) + prb.qobs(1);
                                             -norm(rv(1:prb.n)-prb.robs(:,2)) + prb.qobs(2);
                                              norm(rv(prb.n+1:2*prb.n))^2 - prb.vmax^2;
                                             -norm(T) + prb.Tmin;
                                             ] + cnstr_buffer;

    prb.cnstr_fun_jac_rv = @(rv,T) cnstr_scl\[-(rv(1:prb.n)-prb.robs(:,1))'/norm(rv(1:prb.n)-prb.robs(:,1)) zeros(1,prb.n);
                                              -(rv(1:prb.n)-prb.robs(:,2))'/norm(rv(1:prb.n)-prb.robs(:,2)) zeros(1,prb.n);
                                               zeros(1,prb.n) 2*rv(prb.n+1:2*prb.n)'
                                               zeros(1,2*prb.n)];

    prb.cnstr_fun_jac_T = @(rv,T) cnstr_scl\[ zeros(1,prb.n);
                                              zeros(1,prb.n);
                                              zeros(1,prb.n);
                                             -T'/norm(T)];

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3";
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.AbsTol',1e-8,'ecos.RelTol',1e-8);
    % prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);

    % prb.tr_norm = 2;
    % prb.tr_norm = inf;    
    prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-7;
    prb.epstr = 1e-3;

    % Time grid and time of manuever
    prb.time_of_maneuver = @(x,u)     disc.time_of_maneuver(prb.disc,prb.tau,u(prb.n+1,:));    
    prb.time_grid        = @(tau,x,u)        disc.time_grid(prb.disc,    tau,u(prb.n+1,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func           = @(tau,x,u)      evaluate_dyn_func(x,u,prb.n,prb.c_d,prb.accl,prb.cnstr_fun,cnstr_type);
    prb.dyn_func_linearize = @(tau,x,u) evaluate_linearization(x,u,prb.n,prb.c_d,prb.accl,prb.cnstr_fun,...
                                                                                          prb.cnstr_fun_jac_rv,prb.cnstr_fun_jac_T,cnstr_type);

end

function dx = evaluate_dyn_func(x,u,n,c_d,accl,cnstr_fun,cnstr_type)

    rv          = x(1:2*n);
    T           = u(1:n);
    s           = u(n+1);
    cnstr_val   = cnstr_fun(rv,T);
    drv         = plant.doubleint.dyn_func(rv,T,s,n,c_d,accl);

    if cnstr_type == "exclusive-integrator-states"
        dy      = s*(arrayfun(@(z) max(0,z),cnstr_val) .^ 2);
    elseif cnstr_type == "single-integrator-state" 
        dy      = s*sum( arrayfun(@(z) max(0,z),cnstr_val) .^ 2 );
    else
        error("Invalid constraint formulation");
    end

    dx          = [drv;
                   dy];
end

function [A,B,w] = evaluate_linearization(x,u,n,c_d,accl,cnstr_fun,cnstr_fun_jac_rv,cnstr_fun_jac_T,cnstr_type)

    rv                  = x(1:2*n);
    T                   = u(1:n);
    s                   = u(n+1);
    [dFdrv,dFdT,dFds]   = plant.doubleint.compute_linearization(rv,T,s,n,c_d,accl);
    dFdu                = [dFdT,dFds];

    cnstr_val           = cnstr_fun(rv,T);
    cnstr_val_jac_rv    = cnstr_fun_jac_rv(rv,T);
    cnstr_val_jac_T     = cnstr_fun_jac_T(rv,T);

    m                   = length(cnstr_val);
    absg2               = zeros(m,1);
    absg2_jac_rv        = zeros(m,2*n);
    absg2_jac_T         = zeros(m,n);
    for j = 1:m
        absg2(j)            = max(0,cnstr_val(j))^2;
        absg2_jac_rv(j,:)   = 2*max(0,cnstr_val(j))*cnstr_val_jac_rv(j,:);
        absg2_jac_T(j,:)    = 2*max(0,cnstr_val(j))*cnstr_val_jac_T(j,:);
    end

    if cnstr_type == "exclusive-integrator-states"
        A = [dFdrv          zeros(2*n,m);
             s*absg2_jac_rv zeros(m,m)];
        B = [dFdu;
             s*absg2_jac_T absg2];    
    elseif cnstr_type == "single-integrator-state"
        A = [dFdrv                 zeros(2*n,1);
             s*sum(absg2_jac_rv,1) zeros(1,1)];
        B = [dFdu;
             s*sum(absg2_jac_T,1)  sum(absg2)];    
    else
        error("Invalid constraint formulation");
    end

    dx = evaluate_dyn_func(x,u,n,c_d,accl,cnstr_fun,cnstr_type);
    w = dx - A*x - B*u;
end