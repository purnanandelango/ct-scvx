function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    % Dimension of double integrator
    prb.n = 2;

    prb.nx = 2*prb.n;
    prb.nu = prb.n+1;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/19)*min(prb.dtau);               % Step size for integration that computes discretization matrices
    prb.Kfine = 1+20*round(1/min(prb.dtau));    % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.accl = [0;0];    
    prb.c_d = 0.01; 
    
    % Bounds

    prb.rmax = 40;  
    prb.vmax = 06;
    prb.Tmax = 6.0;
    prb.Tmin = 0.5;

    prb.smin = 01;
    prb.smax = 60;
    prb.ToFguess = 30;

    % Obstacle avoidance
    
    % Centers
    prb.q_obs   = [
                   34 -32 42 -24 34 -32  42 -24  34 -32 ;
                   20  20 10  10  0   0 -10 -10 -20 -20 ;
                   ];
    [~,prb.n_obs] = size(prb.q_obs);

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

    % Boundary conditions

    prb.r1 = [0;-28];
    prb.v1 = [0.1; 0];
    
    prb.rK = [0; 28];
    prb.vK = [0.1; 0];

    % Straight-line initialization    

    prb.x1 = [prb.r1;prb.v1];
    prb.xK = [prb.rK;prb.vK];    
    prb.u1 = [prb.Tmin*ones(prb.n,1);prb.ToFguess];
    prb.uK = [prb.Tmin*ones(prb.n,1);prb.ToFguess];

    % Scaling parameters

    xmin = [-0.5*prb.rmax*ones(prb.n,1); -0.5*prb.vmax*ones(prb.n,1)];
    xmax = [ 0.5*prb.rmax*ones(prb.n,1);  0.5*prb.vmax*ones(prb.n,1)];
    
    umin = [zeros(prb.n,1); 1];
    umax = [prb.Tmax*ones(prb.n,1); prb.ToFguess];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3_parallel";

    prb.ode_solver = {'ode45',odeset('AbsTol',1e-5,'RelTol',1e-6)};
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
    prb.epstr = 7.5e-4;

    % Time grid and time of manuever
    prb.time_of_maneuver = @(x,u)     disc.time_of_maneuver(prb.disc,prb.tau,u(prb.n+1,:));    
    prb.time_grid        = @(tau,x,u)        disc.time_grid(prb.disc,    tau,u(prb.n+1,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func           = @(tau,x,u) plant.doubleint.dyn_func(x,u(1:prb.n),u(prb.n+1),prb.n,prb.c_d,prb.accl);
    prb.dyn_func_linearize = @(tau,x,u)   evaluate_linearization(x,u,prb.n,prb.c_d,prb.accl);

end
function [A,B,w] = evaluate_linearization(x,u,n,c_d,accl)
    [A,B,S,w] = plant.doubleint.compute_linearization(x,u(1:n),u(n+1),n,c_d,accl);
    B = [B,S];
end