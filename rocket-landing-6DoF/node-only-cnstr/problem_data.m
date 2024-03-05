function prb = problem_data(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    prb.nx = 14;
    prb.nu = 4;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/19)*prb.dtau;                    % Step size for integration that computes discretization matrices
    prb.Kfine = 1+round(50/min(prb.dtau));      % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.g0 = 1;
    
    prb.c_ax = 0.5; prb.c_ayz = 1;
    % prb.c_ax = 0; prb.c_ayz = 0;
    
    prb.CA = diag([prb.c_ax,prb.c_ayz,prb.c_ayz]);
    
    prb.JB = 0.168*diag([2e-2,1,1]);
    prb.gI = [-prb.g0;0;0];
    prb.rho = 1;
    prb.SA = 0.5;
    
    prb.rTB = -0.25*[1;0;0];
    prb.rcpB = 0.05*[1;0;0];
    
    prb.Isp = 30;
    prb.alphmdt = 1/(prb.Isp*prb.g0);
    prb.betmdt = 0.01;
    
    % Bounds

    prb.thetmax = 75*pi/180;    prb.sinthetmaxby2 = sin(prb.thetmax/2);     % Vehicle tilt
    prb.gamgs   = 75*pi/180;    prb.cotgamgs = cot(prb.gamgs);              % Glide-slope
    
    prb.omgmax  = 21.5*pi/180;                                              % Angular velocity (inf-norm 28.6)
    prb.delmax  = 20*pi/180;    prb.cosdelmax = cos(prb.delmax);            % Gimbal
    
    prb.Hgam    = [0,1,0;0,0,1];
    prb.Hthet   = [0,1,0,0;0,0,1,0];
    
    prb.TBmin    = 1.5;
    prb.TBmax    = 6.5;
    
    prb.vmax     = 3;
    
    prb.vmax_stc    = 1.5;
    prb.cosaoamax   = cosd( 10 ); 
    prb.stc_flag    = "v1";

    prb.mdry    = 1;
    prb.mwet    = 2;

    prb.smin    = 1.0;
    prb.smax    = 50;
    prb.dtmin   = 0.1;
    prb.dtmax   = 10;
    prb.ToFmax  = 20;
   
    prb.snom    = [1,15];
    prb.ToFguess= 10;    
    
    % Boundary conditions

    prb.rI1     = [7.5;4.5;2];
    prb.vI1     = [-0.5;-2.8;0];
    % prb.vI1    = [-2.5;-0.5;0];
    
    prb.rIK     = zeros(3,1);
    prb.vIK     = -0.1*[1;0;0];
    prb.omgB1   = zeros(3,1);
    prb.omgBK   = prb.omgB1;
    prb.q1      = [0;0;0;1];

    % Straight-line initialization
    prb.x1      = [prb.mwet;prb.rI1;prb.vI1;prb.q1;prb.omgB1];
    prb.xK      = [prb.mdry;prb.rIK;prb.vIK;prb.q1;prb.omgBK];    
    prb.u1      = [-prb.mwet*prb.gI;prb.ToFguess];
    prb.uK      = [-prb.mdry*prb.gI;prb.ToFguess];

    % Scaling parameters
    xmin = [prb.mdry; -10*ones(3,1); -1*ones(3,1); -ones(4,1); -prb.omgmax*ones(3,1)];
    xmax = [prb.mwet;  10*ones(3,1);  1*ones(3,1);  ones(4,1);  prb.omgmax*ones(3,1)];
    
    umin = [prb.TBmin*ones(3,1); prb.snom(1)];
    umax = [prb.TBmax*ones(3,1); prb.snom(2)];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % Constraint parameters


    prb.cnstr_fun       = @(xi,TB)               [prb.mdry - xi(1);
                                                  norm(prb.Hgam*xi(2:4)) - xi(2)/prb.cotgamgs;
                                                  norm(xi(5:7))^2 - prb.vmax^2;
                                                  norm(prb.Hthet*xi(8:11))^2 - prb.sinthetmaxby2^2;
                                                  norm(xi(12:14))^2 - prb.omgmax^2;
                                                  norm(TB) - TB(1)/prb.cosdelmax;
                                                  norm(TB) - prb.TBmax;
                                                  prb.TBmin - norm(TB)];
    
    prb.cnstr_fun_jac_xi = @(xi,TB)               [-1, zeros(1,13);
                                                    0, (prb.Hgam'*prb.Hgam*xi(2:4))'/norm(prb.Hgam*xi(2:4))+[-1/prb.cotgamgs 0 0], zeros(1,10);
                                                    0, zeros(1,3), 2*xi(5:7)', zeros(1,7);
                                                    0, zeros(1,6), 2*(prb.Hthet'*prb.Hthet*xi(8:11))', zeros(1,3);
                                                    0, zeros(1,10), 2*xi(12:14)';
                                                    zeros(1,14);
                                                    zeros(1,14);
                                                    zeros(1,14)];
    
    prb.cnstr_fun_jac_TB = @(xi,TB)               [ zeros(1,3);
                                                    zeros(1,3);
                                                    zeros(1,3);
                                                    zeros(1,3);
                                                    zeros(1,3);
                                                    TB'/norm(TB)+[-1/prb.cosdelmax 0 0];
                                                    TB'/norm(TB);
                                                   -TB'/norm(TB)];

    prb.stc_fun = @(x,u) plant.rocket6DoF.q_aoa_cnstr(x(5:7),x(8:11),prb.vmax_stc,prb.cosaoamax,prb.stc_flag);    

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3";
    prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-7)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations
    
    % prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.AbsTol',1e-8,'ecos.RelTol',1e-8,'ecos.FeasTol',1e-9);

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