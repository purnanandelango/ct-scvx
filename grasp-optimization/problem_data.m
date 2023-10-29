function prb = problem_data(K,T, ...
                            scp_iters,wvc,wtr,cost_factor)

    prb.K = K;

    prb.nx = 6+1+1;
    prb.nu = 12;

    prb.tau = linspace(0,T,K);
    prb.dtau = T/(K-1);
        
    prb.h = (1/10)*prb.dtau;
    prb.Kfine = 20*K;

    % System parameters
    prb.box_width   = 0.1;                          % Width of box [m]
    prb.mass        = 0.4;                          % Mass [kg]
    prb.accl        = 9.806;                        % Acceleration due to gravity [kg m/s^2]
    prb.mu          = [1.1;                         % Friction coefficient
                       1.7;
                       1.8;
                       1.9];

    % Input indices
    prb.idx = {1:3,4:6,7:9,10:12};

    % Position vector of contact points in block frame
    prb.contact1 = [-prb.box_width;  0;               0];
    prb.contact2 = [ prb.box_width;  0;               0];
    prb.contact3 = [ 0;             -prb.box_width;   0];
    prb.contact4 = [ 0;              prb.box_width;   0];

    % Matrix in linear equality constraints due to rotational equlibrium
    prb.contact_mat = [qlib.skew(prb.contact1), ...
                       qlib.skew(prb.contact2), ...
                       qlib.skew(prb.contact3), ...
                       qlib.skew(prb.contact4)];

    % Bounds
    prb.rmax  = 10;                  % Position
    prb.vmax  = 2;                   % Speed
    prb.pmax  = 50;                  % Cost 
    prb.F1max = 2;                   % First finger        
    prb.F2max = 2;                   % Second finger         
    prb.F3max = 3;                   % Third finger        
    prb.F4max = 3;                   % Fourth finger 
    prb.Fmin = 0.2;                  % Minimum normal force

    prb.ymin = 0;
    prb.ymax = 0.1;

    % Boundary conditions

    prb.r1 = [20;-2;10];
    prb.v1 = [0;0;0];
    prb.p1 = 0;
    prb.y1 = 0;

    prb.rK = [20;10;20];
    prb.vK = [0;0;0];

    % Initialization generator

    prb.x1 = [prb.r1; 0.5*prb.vmax*ones(3,1); 0.5*prb.pmax; prb.y1];
    prb.xK = [prb.rK; 0.5*prb.vmax*ones(3,1); 0.5*prb.pmax; prb.y1];    

    prb.u1 = min([prb.F1max,prb.F2max,prb.F3max,prb.F4max])*ones(prb.nu,1);
    prb.uK = min([prb.F1max,prb.F2max,prb.F3max,prb.F4max])*ones(prb.nu,1);

    % Scaling parameters

    xmin = 0*[-0.5*prb.rmax*ones(3,1);-0.5*prb.vmax*ones(3,1); 0;            prb.ymin];
    xmax =   [ 0.5*prb.rmax*ones(3,1); 0.5*prb.vmax*ones(3,1); 0.5*prb.pmax; prb.ymax];

    umin =   [-0.5*prb.F1max*ones(3,1);
              -0.5*prb.F2max*ones(3,1);
              -0.5*prb.F3max*ones(3,1);
              -0.5*prb.F4max*ones(3,1);
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
    
    cnstr_scl = diag([0.05; ...
                      1*ones(3,1); ...
                      1; ...
                      1; ...
                      1; ...
                      1; ...
                      1; ...
                      1; ...
                      1; ...
                      1;
                      ]);
    cnstr_buffer = zeros(12,1);

    % Path constraints

    prb.cnstr_fun = @(x,u) cnstr_scl*[norm(x(4:6))^2 - prb.vmax^2;
                                      prb.contact_mat*u;
                                      norm(u(prb.idx{1}(2:3))) - prb.mu(1)*u(prb.idx{1}(1));
                                      norm(u(prb.idx{2}(2:3))) + prb.mu(2)*u(prb.idx{2}(1));
                                      norm(u(prb.idx{3}([1,3]))) - prb.mu(3)*u(prb.idx{3}(2));
                                      norm(u(prb.idx{4}([1,3]))) + prb.mu(4)*u(prb.idx{4}(2));
                                      norm(u(prb.idx{1})) - prb.F1max;
                                      norm(u(prb.idx{2})) - prb.F2max;
                                      norm(u(prb.idx{3})) - prb.F3max;
                                      norm(u(prb.idx{4})) - prb.F4max;
                                      ] + cnstr_buffer;

    prb.cnstr_fun_jac_x = @(x,u) cnstr_scl*[zeros(1,3), 2*x(4:6)', 0;
                                            zeros(11,7);
                                            ];

    prb.cnstr_fun_jac_u = @(x,u) cnstr_scl*[zeros(1,prb.nu);
                                            prb.contact_mat;
                                            -prb.mu(1), u(2:3)'/norm(u(2:3)), zeros(1,9);
                                            zeros(1,3), prb.mu(2), u(5:6)'/norm(u(5:6)), zeros(1,6);
                                            zeros(1,6), u(7)/norm(u([7,9])), -prb.mu(3), u(9)/norm(u([7,9])), zeros(1,3);
                                            zeros(1,9), u(10)/norm(u([10,12])),  prb.mu(4), u(12)/norm(u([10,12]));
                                            u(1:3)'/norm(u(1:3)), zeros(1,9);
                                            zeros(1,3), u(4:6)'/norm(u(4:6)), zeros(1,6);
                                            zeros(1,6), u(7:9)'/norm(u(7:9)), zeros(1,3);
                                            zeros(1,9), u(10:12)'/norm(u(10:12));
                                            ];

% SCP parameters

    prb.disc = "ZOH";
    prb.foh_type = "v3";
    prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-7)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','mosek','verbose',false,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-9,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-9);
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.abstol',1e-8,'ecos.reltol',1e-8);   
    % prb.solver_settings = sdpsettings('solver','quadprog','verbose',false);
    % prb.solver_settings = sdpsettings('solver','osqp','verbose',false,'osqp.eps_abs',1e-7,'osqp.eps_rel',1e-7,'osqp.max_iter',5e4);        

    prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-8;
    prb.epstr = 1e-3;

    % Takes in unscaled data
    prb.time_of_maneuver = @(z,u) T;    
    prb.time_grid = @(tau,z,u) prb.tau;    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(tau,xtil,u)             evaluate_dyn_func     (xtil,u,prb.mass,prb.accl,prb.cnstr_fun);
    prb.dyn_func_linearize = @(tau,xtil,u)   evaluate_linearization(xtil,u,prb.mass,prb.accl,prb.cnstr_fun,...
                                                                    prb.cnstr_fun_jac_x,prb.cnstr_fun_jac_u);    

end
function F = evaluate_dyn_func(xtil,u,mass,accl,cnstr_fun)
   x = xtil(1:7);
   v = x(4:6);
   
   cnstr_val = cnstr_fun(x,u);

   Bsys = repmat(eye(3),[1,4]);
   
   f = [v;
        Bsys*u/mass + [0;0;-accl];
        % norm(u(1:3)) + norm(u(4:6)) + norm(u(7:9)) + norm(u(10:12));
        norm(u(1:3))^2 + norm(u(4:6))^2 + norm(u(7:9))^2 + norm(u(10:12))^2;
        ];
    
   F = [f;
        sum([arrayfun(@(y) max(0,y)^2,cnstr_val(1)); ...
             arrayfun(@(y) y^2,       cnstr_val(2:4)); ...
             arrayfun(@(y) max(0,y)^2,cnstr_val(5:12))])];
end
function [A,B,w] = evaluate_linearization(xtil,u,mass,accl,cnstr_fun,cnstr_fun_jac_x,cnstr_fun_jac_u)
   x = xtil(1:7);
   % v = x(4:6);

   cnstr_val = cnstr_fun(x,u);
   abs_cnstr_val = [arrayfun(@(y) max(0,y),cnstr_val(1));
                    cnstr_val(2:4);
                    arrayfun(@(y) max(0,y),cnstr_val(5:12))];
   cnstr_val_jac_x = cnstr_fun_jac_x(x,u);
   cnstr_val_jac_u = cnstr_fun_jac_u(x,u);

   Bsys = repmat(eye(3),[1,4]);

   % f = [v;
   %      Bsys*u/mass + [0;0;-accl];
   %      % norm(u(1:3)) + norm(u(4:6)) + norm(u(7:9)) + norm(u(10:12));
   %      norm(u(1:3))^2 + norm(u(4:6))^2 + norm(u(7:9))^2 + norm(u(10:12))^2;
   % ];   

   dfdx = [zeros(3,3), eye(3), zeros(3,1);
           zeros(3,7);
           zeros(1,7)];

   dfdu = [zeros(3,12);
           Bsys/mass;
           % u(1:3)'/norm(u(1:3)), u(4:6)'/norm(u(4:6)), u(7:9)'/norm(u(7:9)), u(10:12)'/norm(u(10:12));
           2*u(1:3)', 2*u(4:6)', 2*u(7:9)', 2*u(10:12)';
           ];

   A = [dfdx, zeros(7,1);
        2*abs_cnstr_val'*cnstr_val_jac_x, 0];

   B = [dfdu;
        2*abs_cnstr_val'*cnstr_val_jac_u];
    
   F = evaluate_dyn_func(xtil,u,mass,accl,cnstr_fun);
    
   w = F - A*xtil - B*u;   
end