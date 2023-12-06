function prb = problem_data(n,K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    % Dimension of double integrator
    prb.n = n;

    % Number of integrated constraints
    prb.m = 4;

    prb.nx = 2*prb.n + prb.m;
    prb.nu = prb.n+1;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/19)*min(prb.dtau);            % Step size for integration that computes FOH matrices
    prb.Kfine = 1+20*round(1/min(prb.dtau));    % Size of grid on which SCP solution is simulated
    
    % System parameters

    % accl defined later
    
    prb.c_d = 0.1; 
    
    % Bounds

    prb.rmax = 40;
    prb.vmax = 7;
    prb.Tmax = 7;
    prb.Tmin = 1.5;

    prb.smin = 1;
    prb.smax = 20;
    prb.snom = [5,10];
    prb.dtmin = 0.5;
    prb.dtmax = 3;
    prb.ToFmax = 20;
    prb.ToFguess = 10;

    prb.ymin = 0*ones(prb.m,1);
    prb.ymax = 1*ones(prb.m,1);

    prb.eps_cnstr = 1e-4;

    % Obstacle avoidance
    prb.nobs = 2;
    prb.Hobs = cell(1,prb.nobs);
    prb.hobs = cell(1,prb.nobs);

    if prb.n == 2

        [prb.Hobs{1},prb.hobs{1}] = geom.construct_warped_box(3,[-9;27],2,2);       
        [prb.Hobs{2},prb.hobs{2}] = geom.construct_warped_box(4,[-9;7],5,5);
    
        % Boundary conditions
    
        prb.r1 = [0;0];
        prb.v1 = [0;0];
        
        prb.rK = [-15;28];
        prb.vK = -0.1*[1;0];
    
        prb.accl = [0;-1];

    elseif prb.n == 3

        [prb.Hobs{1},prb.hobs{1}] = geom.construct_warped_box(5,[-5;10;0],2,2);       
        [prb.Hobs{2},prb.hobs{2}] = geom.construct_warped_box(3,[-10;22;0],5,5);
    
        % Boundary conditions
    
        prb.r1 = [0;0;0];
        prb.v1 = [0;0;0];
        
        prb.rK = [-15;28;0];
        prb.vK = -0.1*[1;0;1];  

        prb.accl = [0;0;-1];        

    end

    % Straight-line initialization
    prb.y1 = zeros(prb.m,1);    
    prb.yK = zeros(prb.m,1); % Not imposed as constraint; used to generate initialization    

    prb.x1 = [prb.r1;prb.v1;prb.y1];
    prb.xK = [prb.rK;prb.vK;prb.yK];    
    prb.u1 = [0.5*prb.Tmax*ones(prb.n,1);prb.ToFguess];
    prb.uK = [0.5*prb.Tmax*ones(prb.n,1);prb.ToFguess];

    % Scaling parameters
    xmin = [-0.5*prb.rmax*ones(prb.n,1); -0.5*prb.vmax*ones(prb.n,1); prb.ymin];
    xmax = [ 0.5*prb.rmax*ones(prb.n,1);  0.5*prb.vmax*ones(prb.n,1); prb.ymax];
    
    umin = [prb.Tmin*ones(prb.n,1); prb.snom(1)];
    umax = [prb.Tmax*ones(prb.n,1); prb.snom(2)];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % Constraint parameters
    prb.sgndst_min = 0.01;

    prb.cnstr_fun       = @(rv,T,sgndst) [-sgndst' + prb.sgndst_min;
                                           norm(rv(prb.n+1:2*prb.n))^2 - prb.vmax^2;
                                          -norm(T) + prb.Tmin];

    prb.cnstr_fun_jac_rv = @(rv,T,sgndst,sgndst_jac) [-sgndst_jac'        zeros(prb.nobs,prb.n);
                                                       zeros(1,prb.n)     2*rv(prb.n+1:2*prb.n)'
                                                       zeros(1,prb.n)     zeros(1,prb.n)];

    prb.cnstr_fun_jac_T = @(rv,T) [ zeros(prb.nobs,prb.n);
                                    zeros(1,prb.n);
                                   -T'/norm(T)];    

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3_parallel";
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

    % Takes in unscaled data
    prb.time_of_maneuver =     @(x,u) disc.time_of_maneuver(prb.disc,prb.tau,u(prb.n+1,:));    
    prb.time_grid        = @(tau,x,u)        disc.time_grid(prb.disc,    tau,u(prb.n+1,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func           = @(tau,x,u)      evaluate_dyn_func(x,u,prb.n,prb.c_d,prb.accl,prb.Hobs,prb.hobs,prb.cnstr_fun);
    prb.dyn_func_linearize = @(tau,x,u) evaluate_linearization(x,u,prb.n,prb.c_d,prb.accl,prb.Hobs,prb.hobs,prb.cnstr_fun,...
                                                                                                            prb.cnstr_fun_jac_rv,prb.cnstr_fun_jac_T);

end

function dx = evaluate_dyn_func(x,u,n,c_d,accl,Hobs,hobs,cnstr_fun)
    sgndst = zeros(1,length(Hobs));
    for j = 1:length(Hobs)
        [~,sgndst(j)] = geom.sign_dist_polyhed(x(1:n),Hobs{j},hobs{j});
    end
    rv = x(1:2*n);
    T = u(1:n);
    s = u(n+1);
    cnstr_val = cnstr_fun(rv,T,sgndst);
    drv = plant.doubleint.dyn_func(rv,T,s,n,c_d,accl);
    dx = [drv;
          s*(arrayfun(@(y) max(0,y),cnstr_val) .^ 2)];
end

function [A,B,w] = evaluate_linearization(x,u,n,c_d,accl,Hobs,hobs,cnstr_fun,cnstr_fun_jac_rv,cnstr_fun_jac_T)
    sgndst = zeros(1,length(Hobs));
    sgndst_jac = zeros(n,length(Hobs));
    for j = 1:length(Hobs)
        [~,sgndst(j),sgndst_jac(:,j)] = geom.sign_dist_polyhed(x(1:n),Hobs{j},hobs{j});
    end
    rv = x(1:2*n);
    T = u(1:n);
    s = u(n+1);
    [dFdrv,dFdT,dFds] = plant.doubleint.compute_linearization(rv,T,s,n,c_d,accl);
    dFdu = [dFdT,dFds];

    cnstr_val           = cnstr_fun(rv,T,sgndst);
    cnstr_val_jac_rv    = cnstr_fun_jac_rv(rv,T,sgndst,sgndst_jac);
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

    A = [dFdrv          zeros(2*n,m);
         s*absg2_jac_rv zeros(m,m)];
    B = [dFdu;
         s*absg2_jac_T absg2];

    dx = evaluate_dyn_func(x,u,n,c_d,accl,Hobs,hobs,cnstr_fun);
    w = dx - A*x - B*u;
end