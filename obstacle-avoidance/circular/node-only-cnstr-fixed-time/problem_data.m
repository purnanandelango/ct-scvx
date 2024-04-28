function prb = problem_data(tau_f,K,scp_iters,w_ep,w_px,cost_factor)
    
    prb.K = K;

    % Dimension of double integrator
    prb.n = 2;

    prb.nx = 2*prb.n;
    prb.nu = prb.n;
    
    prb.tau = grid.generate_grid(0,tau_f,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/19)*min(prb.dtau);               % Step size for integration that computes discretization matrices
    prb.Kfine = 1+100*round(1/min(prb.dtau));    % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.accl = [0;-1];
    
    prb.c_d = 0.1; 
    
    % Bounds

    prb.rmax = 40;  
    prb.vmax = 6;
    prb.Tmax = 7;
    prb.Tmin = 0.5;

    % Obstacle avoidance
    prb.nobs = 2;

    prb.robs = [-5 -10;
                 6  20];
    prb.aobs = [6 7];

    % Boundary conditions

    prb.r1 = [0;0];
    prb.v1 = 0.1*[1;0];
    
    prb.rK = [-15;28];
    prb.vK = 0.1*[-1;0];

    % Straight-line initialization    

    prb.x1 = [prb.r1;prb.v1];
    prb.xK = [prb.rK;prb.vK];    
    prb.u1 = ones(prb.n,1);
    prb.uK = ones(prb.n,1);

    % Scaling parameters

    xmin = [-0.5*prb.rmax*ones(prb.n,1); -0.5*prb.vmax*ones(prb.n,1)];
    xmax = [ 0.5*prb.rmax*ones(prb.n,1);  0.5*prb.vmax*ones(prb.n,1)];
    
    umin = prb.Tmin*ones(prb.n,1);
    umax = prb.Tmax*ones(prb.n,1);

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % SCP parameters

    % prb.disc = "FOH";
    % prb.foh_type = "v3";

    prb.disc = "Impulse";
    prb.Eu2x = [zeros(prb.n);
                eye(prb.n)];

    prb.scp_iters = scp_iters; % Maximum SCP iterations

    prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.AbsTol',1e-8,'ecos.RelTol',1e-8);
    % prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);    
    
    % prb.px_norm = 2;
    % prb.px_norm = inf;
    prb.px_norm = 'quad';
    
    prb.w_ep = w_ep;
    prb.w_px = w_px;
    prb.cost_factor = cost_factor;
    
    prb.eps_ep = 1e-7;
    prb.eps_px = 1e-3;

    % Time grid and time of manuever
    prb.time_of_maneuver = @(x,u)      tau_f;    
    prb.time_grid        = @(tau,x,u)  linspace(0,tau_f,length(tau));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func           = @(tau,x,u) plant.doubleint.dyn_func(x,u,1.0,prb.n,prb.c_d,prb.accl);
    prb.dyn_func_linearize = @(tau,x,u)   evaluate_linearization(x,u,prb.n,prb.c_d,prb.accl);

end
function [A,B,w] = evaluate_linearization(x,u,n,c_d,accl)
    [A,B,~,~] = plant.doubleint.compute_linearization(x,u,1.0,n,c_d,accl);
    f = plant.doubleint.dyn_func(x,u,1.0,n,c_d,accl);
    w = f - A*x - B*u;
end