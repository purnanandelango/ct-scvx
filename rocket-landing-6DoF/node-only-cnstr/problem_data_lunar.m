function prb = problem_data_lunar(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    prb.nx = 14;
    prb.nu = 4;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/19)*prb.dtau;                    % Step size for integration that computes discretization matrices
    prb.Kfine = 1+round(50/min(prb.dtau));      % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.g0 = 1.61;
    
    % Not used
    prb.c_ax = 0; prb.c_ayz = 0;
    prb.rho = 1;
    prb.SA = 0.5;    
    prb.CA = diag([prb.c_ax,prb.c_ayz,prb.c_ayz]);
    prb.rcpB = 0.05*[1;0;0];
    
    prb.JB = diag([19150,13600,13600]);
    prb.gI = [-prb.g0;0;0];
    
    prb.rTB = -0.25*[1;0;0];
    
    prb.Isp = 225;
    ge = 9.806;
    prb.alphmdt = 1/(prb.Isp*ge);
    prb.betmdt = 0.00; % Not used
    
    % Bounds

    prb.thetmax = 60*pi/180;    prb.sinthetmaxby2 = sin(prb.thetmax/2);     % Vehicle tilt
    prb.gamgs   = 85*pi/180;    prb.cotgamgs = cot(prb.gamgs);              % Glide-slope
    
    prb.omgmax  = 10*pi/180;                                                % Angular velocity (inf-norm 28.6)
    prb.delmax  = 45*pi/180;    prb.cosdelmax = cos(prb.delmax);            % Gimbal
    
    prb.Hgam    = [0,1,0;0,0,1];
    prb.Hthet   = [0,1,0,0;0,0,1,0];
    
    prb.TBmin    = 5000;
    prb.TBmax    = 22500;
    
    prb.vmax     = 50;
    
    prb.vmax_stc    = 25;
    prb.cosaoamax   = cosd( 10 ); 
    prb.stc_flag    = "v1";

    prb.mdry    = 2100;
    prb.mwet    = 3250;

    prb.smin    = 1;
    prb.smax    = 150;
   
    prb.ToFguess= 55;   
    
    % Boundary conditions

    prb.rI1     = [433;0;250];
    prb.vI1     = [10;0;-30];
    
    prb.rIK     = [30;0;0];
    prb.vIK     = [-1;0;0];
    prb.omgB1   = zeros(3,1);
    prb.omgBK   = zeros(3,1);
    prb.q1      = [0;0;0;1];

    % Straight-line initialization
    prb.x1      = [prb.mwet;prb.rI1;prb.vI1;prb.q1;prb.omgB1];
    prb.xK      = [prb.mdry;prb.rIK;prb.vIK;prb.q1;prb.omgBK];    
    prb.u1      = [-2.0*prb.mwet*prb.gI;prb.ToFguess];
    prb.uK      = [-2.0*prb.mdry*prb.gI;prb.ToFguess];

    % Scaling parameters
    xmin = [prb.mdry; -400*ones(3,1); -100*ones(3,1); -ones(4,1); -prb.omgmax*ones(3,1)];
    xmax = [prb.mwet;  400*ones(3,1);  100*ones(3,1);  ones(4,1);  prb.omgmax*ones(3,1)];
    
    umin = [-prb.TBmax*ones(3,1); prb.smin];
    umax = [ prb.TBmax*ones(3,1); prb.smax];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3_parallel";
    prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-7)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations
    
    prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.AbsTol',1e-8,'ecos.RelTol',1e-8,'ecos.FeasTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','mosek','verbose',false);

    % prb.tr_norm = 2;
    % prb.tr_norm = inf;
    prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-7;
    prb.epstr = 1e-3;

    % Time of maneuver and time grid
    prb.time_of_maneuver = @(x,u) disc.time_of_maneuver(prb.disc,prb.tau,u(4,:));    
    prb.time_grid = @(tau,x,u) disc.time_grid(prb.disc,tau,u(4,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func           = @(tau,x,u) plant.rocket6DoF.dyn_func(x,u(1:3),u(4),prb.c_ax,prb.c_ayz,diag(prb.JB),prb.gI,prb.rho,prb.SA,prb.rTB,prb.rcpB,prb.alphmdt,prb.betmdt);
    prb.dyn_func_linearize = @(tau,x,u)    evaluate_linearization(x,u,          prb.c_ax,prb.c_ayz,diag(prb.JB),prb.gI,prb.rho,prb.SA,prb.rTB,prb.rcpB,prb.alphmdt,prb.betmdt);

end

function [A,B,w] = evaluate_linearization(x,u,c_ax,c_ayz,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt)
    [A,B,S,w] = plant.rocket6DoF.compute_linearization(x,u(1:3),u(4),c_ax,c_ayz,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
    B = [B,S];
end