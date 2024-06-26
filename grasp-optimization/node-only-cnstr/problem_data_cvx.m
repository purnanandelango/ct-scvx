function prb = problem_data_cvx(K,T)

    prb.K = K;

    prb.nx = 6;
    prb.nu = 9;

    prb.tau = linspace(0,T,K);
    prb.dtau = T/(K-1);

    % System parameters
    prb.box_width   = 0.1;                          % Width of box [m]
    prb.mass        = 0.4;                          % Mass [kg]
    prb.accl        = 9.806;                        % Acceleration due to gravity [kg m/s^2]
    prb.mu          = [1.5;                         % Friction coefficient
                       1.7;
                       3.1];

    % Input indices    
    prb.idx = {1:3,4:6,7:9};


    % Position vector of contact points in block frame
    prb.contact1 = [-prb.box_width;  0;               0];
    prb.contact2 = [ prb.box_width;  0;               0];
    prb.contact3 = [ 0;              prb.box_width;   0];

    % Matrix in linear equality constraints due to rotational equlibrium
    prb.contact_mat = [qlib.skew(prb.contact1), ...
                       qlib.skew(prb.contact2), ...
                       qlib.skew(prb.contact3)];

    % Bounds
    prb.rmax  = 10;                 % Position
    prb.vmax  = 1.3;                % Speed
    prb.F1max = 2.7;                % First finger        
    prb.F2max = 2.7;                % Second finger         
    prb.F3max = 3;                  % Third finger
    prb.F1min = 0.5;                % Minimum normal force
    prb.F2min = 0.5;                % Minimum normal force
    prb.F3min = 2.0;                % Minimum normal force

    % Boundary conditions

    prb.r1 = [20;-2;10];
    prb.v1 = [0;0;0];

    prb.rK = [15;10;10];
    prb.vK = [0;0;0];

    prb.xK = [prb.rK;prb.vK];

    % Scaling parameters

    xmin =   [-0.5*prb.rmax*ones(3,1);-0.5*prb.vmax*ones(3,1)];
    xmax =   [ 0.5*prb.rmax*ones(3,1); 0.5*prb.vmax*ones(3,1)];

    umin =   [-prb.F1max*ones(3,1);
              -prb.F2max*ones(3,1);
              -prb.F3max*ones(3,1);
              ];
    umax =   [ prb.F1max*ones(3,1);
               prb.F2max*ones(3,1);
               prb.F3max*ones(3,1);
              ];

    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},[0,1]);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};    

    % System matrices
    prb.Ac = [zeros(3), eye(3);
              zeros(3,6)];
    prb.Bc = [zeros(3,9);repmat(eye(3),[1,3])/prb.mass];
    prb.wc = [zeros(3,1);0;0;-prb.accl];


    prb.Ad = [eye(3),   prb.dtau*eye(3);
              zeros(3), eye(3)];
    
    prb.Bd = [repmat(prb.dtau*prb.dtau*0.5*eye(3)/prb.mass,[1,3]);
              repmat(prb.dtau*eye(3)/prb.mass,[1,3])];
      
    prb.wd = [0;
              0;
             -0.5*prb.dtau*prb.dtau*prb.accl;
              0;
              0;
             -prb.dtau*prb.accl];

    prb.solver_settings = sdpsettings('solver','gurobi','verbose',true,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',true,'ecos.abstol',1e-8,'ecos.reltol',1e-8); 

    prb.dyn_func = @(tau,x,u) prb.Ac*x + prb.Bc*u + prb.wc;
end