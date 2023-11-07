function prb = problem_data_cvx(K,T)

    prb.K = K;

    prb.nx = 7;
    prb.nu = 4;

    prb.tau = linspace(0,T,K);
    prb.dtau = T/(K-1);

    % System parameters
    prb.mwet = 2000;            % Wet mass [kg]
    prb.mdry = 150;             % Dry mass [kg]

    prb.gvec = [0;0;-3.71];     % Acceleration due to gravity [m/s^2]

    % Bounds
    prb.rmax     = 4000;                % Position [m]
    prb.Tmax     = 10;                  % Maximum thrust [N]
    prb.vmax     = 300;                 % Speed [m/s]
    prb.rho1     = 0.2*Tmax;            % Thrust lower-bound [N]        
    prb.rho2     = 0.8*Tmax;            % Thrust upper-bound [N]  
    prb.tangamgs = tand( 85 );          % Tangent of maximum glide-slope angle
    
    % Boundary conditions

    prb.r1 = [0;4000;1500];             % [m]
    prb.v1 = [0;100;-75];               % [m/s]
    prb.p1 = log(prb.mwet); 

    prb.rK = [0;0;0];                   % [m]
    prb.vK = [0;0;0];                   % [m/s]
    prb.pK = log(prb.mdry); 

    prb.xK = [prb.rK;prb.vK];

    % Scaling parameters

    xmin =   0*[-0.5*prb.rmax*ones(3,1);-0.5*prb.vmax*ones(3,1)];
    xmax =   [ 0.5*prb.rmax*ones(3,1); 0.5*prb.vmax*ones(3,1)];

    umin =   [
              ];
    umax =   [ 0.5*prb.F1max*ones(3,1);
               0.5*prb.F2max*ones(3,1);
               0.5*prb.F3max*ones(3,1);
               0.5*prb.F4max*ones(3,1);
              ];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};    

    % System matrices
    prb.Ac = [zeros(3), eye(3);
              zeros(3,6)];
    prb.Bc = repmat(eye(3),[1,4])/prb.mass;
    prb.wc = [zeros(3,1);0;0;-prb.accl];


    prb.Ad = [eye(3),   prb.dtau*eye(3);
              zeros(3), eye(3)];
    
    prb.Bd = [repmat(prb.dtau*prb.dtau*0.5*eye(3)/prb.mass,[1,4]);
              repmat(prb.dtau*eye(3)/prb.mass,[1,4])];
      
    prb.wd = [0;
              0;
             -0.5*prb.dtau*prb.dtau*prb.accl;
              0;
              0;
             -prb.dtau*prb.accl];

    prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',true,'ecos.abstol',1e-8,'ecos.reltol',1e-8);    
end