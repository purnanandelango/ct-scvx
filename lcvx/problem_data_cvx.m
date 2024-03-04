function prb = problem_data_cvx(K,T)

    prb.K = K;

    prb.nx = 7;
    prb.nu = 4;

    prb.tau = linspace(0,T,K);
    prb.dtau = T/(K-1);

    % System parameters
    prb.mwet = 1905;             % Wet mass [kg]
    prb.mdry = 1505;             % Dry mass [kg]

    prb.gvec = [0;0;-3.71];     % Acceleration due to gravity [m/s^2]

    % Bounds
    prb.rmax     = 4000;                            % Position [m]
    prb.Tmax     = 24000;                           % Maximum thrust [N]
    prb.vmax     = 130;                             % Speed [m/s]
    prb.rho1     = 0.2*prb.Tmax*cosd(0);         % Thrust lower-bound [N]        
    prb.rho2     = 0.8*prb.Tmax*cosd(0);         % Thrust upper-bound [N]  
    prb.gam_gs   = 18;                               % Minimum glide-slope angle [deg]
    prb.thet_tp  = 50;                              % Thrust pointing [deg]
    prb.tp_vec   = [0; 
                    0; 
                   -1/cosd(prb.thet_tp); 
                    1];
    prb.alpha    = 0.0005;

    prb.mu_1     = prb.rho1 ./ (prb.mwet - prb.alpha*prb.rho2*prb.tau);
    prb.mu_2     = prb.rho2 ./ (prb.mwet - prb.alpha*prb.rho2*prb.tau);
    prb.p_0      = log(prb.mwet - prb.alpha*prb.rho2*prb.tau);
    prb.p_1      = log(prb.mwet - prb.alpha*prb.rho1*prb.tau);
    prb.p_min    = max(log(prb.mdry),prb.p_0);

    prb.thrst_bnd_mat           = cell(prb.K-1,1);
    prb.thrst_bnd_vec           = cell(prb.K-1,1);
    for k = 1:prb.K-1    
        prb.thrst_bnd_mat{k}    = [-prb.mu_1(k) -1;
                                    prb.mu_2(k)  1];
        prb.thrst_bnd_vec{k}    = [-prb.mu_1(k)*(1+prb.p_0(k));
                                    prb.mu_2(k)*(1+prb.p_0(k))];
    end    
    
    % Boundary conditions

    prb.r1 = [0;2000;1400];             % [m]
    prb.v1 = [0;100;-50];               % [m/s]
    prb.p1 = log(prb.mwet); 

    prb.rK = [0;0;0];                   % [m]
    prb.vK = [0;0;0];                   % [m/s]
    prb.pK = log(prb.mdry); 

    prb.xK = [prb.rK;prb.vK];

    % Scaling parameters

    xmin =   [-0.5*prb.rmax*ones(3,1);-0.5*prb.vmax*ones(3,1); 0];
    xmax =   [ 0.5*prb.rmax*ones(3,1); 0.5*prb.vmax*ones(3,1); prb.p1];

    umin =   [ prb.rho1*ones(3,1);
               prb.rho1;
              ];
    umax =   [ prb.rho2*ones(3,1);
               prb.rho2;
              ];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};    

    % System matrices
    prb.Ac = [zeros(3), eye(3), zeros(3,1);
              zeros(4,7)];
    prb.Bc = [zeros(3,4);
              eye(3), zeros(3,1);
              zeros(1,3), -prb.alpha];
    prb.wc = [zeros(3,1);prb.gvec;0];


    M = expm([prb.Ac, prb.Bc, prb.wc;
              zeros(5,12)]*prb.dtau);
    prb.Ad = M(1:7,1:7);
    prb.Bd = M(1:7,7+1:7+4);
    prb.wd = M(1:7,7+5);

    prb.solver_settings = sdpsettings('solver','gurobi','verbose',true,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',true,'ecos.abstol',1e-8,'ecos.reltol',1e-8);
    % prb.solver_settings = sdpsettings('solver','mosek');

    prb.dyn_func = @(tau,x,u) prb.Ac*x + prb.Bc*u + prb.wc;    
end