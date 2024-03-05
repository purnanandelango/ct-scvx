function prb = problem_data(K,T, ...
                            scp_iters,wvc,wtr,cost_factor)

    prb.K = K;

    prb.nx = 6+1+1;
    prb.nu = 4;

    prb.tau = linspace(0,T,K);
    prb.dtau = T/(K-1);
        
    prb.h = (1/10)*prb.dtau;
    prb.Kfine = 20*K;

    % System parameters
    prb.mwet = 1905;             % Wet mass [kg]
    prb.mdry = 1505;             % Dry mass [kg]

    prb.gvec = [0;0;-3.71];     % Acceleration due to gravity on Mars [m/s^2]

    % Bounds
    prb.rmax     = 4000;                            % Position [m]
    prb.Tmax     = 16572;                           % Maximum thrust [N]
    prb.vmax     = 500*5/18;                        % Maximum speed [m/s]
    prb.rho1     = 0.3*prb.Tmax*cosd(0);            % Thrust lower-bound [N]        
    prb.rho2     = 0.8*prb.Tmax*cosd(0);            % Thrust upper-bound [N]  
    prb.gam_gs   = 6;                               % Minimum glide-slope angle [deg]
    prb.thet_tp  = 40;                              % Maximum thrust-pointing angle [deg]
    prb.tp_vec   = [0; 
                    0; 
                    1/cosd(prb.thet_tp); 
                   -1];
    Isp          = 225;                             % Rocket engine specific impulse [s]
    ge           = 9.806;                           % Earth's gravity at sea-level [m/s^2] 
    prb.alpha    = 1/(Isp*ge);

    mu_1  = @(tau) prb.rho1/(prb.mwet - prb.alpha*prb.rho2*tau);
    mu_2  = @(tau) prb.rho2/(prb.mwet - prb.alpha*prb.rho2*tau);    
    p_0   = @(tau) log(prb.mwet - prb.alpha*prb.rho2*tau);
    p_1   = @(tau) log(prb.mwet - prb.alpha*prb.rho1*tau);
    p_min = @(tau) max(log(prb.mdry),p_0(tau));
    a_1   = @(tau) mu_1(tau)*(1+p_0(tau)+0.5*p_0(tau)^2);
    b_1   = @(tau) -mu_1(tau)*(1+p_0(tau));
    c_1   = @(tau) 0.5*mu_1(tau);
    a_2   = @(tau) -mu_2(tau)*(1+p_0(tau));
    b_2   = @(tau) mu_2(tau);

    prb.ymin = 0;
    prb.ymax = 1;

    prb.eps_cnstr = 1e-6;    

    % Boundary conditions

    prb.r1 = [2000;0;1500];             % [m]
    prb.v1 = [288;108;-270]*5/18;       % [m/s]
    prb.p1 = log(prb.mwet); 
    prb.y1 = 0;

    prb.rK = [0;0;0];                   % [m]
    prb.vK = [0;0;0];                   % [m/s]
    prb.pK = log(prb.mdry);

    prb.Ey = [zeros(1,prb.nx-1),1];

    Inx = eye(prb.nx);
    prb.Ei      = Inx;
    prb.zi      = [prb.r1;prb.v1;prb.p1;prb.y1];
    prb.i_idx   = 1:prb.nx;
    prb.Ef      = Inx(1:6,:);
    prb.zf      = [prb.rK;prb.vK];
    prb.f_idx   = 1:6;    
    prb.term_cost_vec = [zeros(6,1);-1;0];        

    % Initialization generator

    prb.x1 = [prb.r1; 0.5*prb.vmax*ones(3,1); prb.p1; prb.y1];
    prb.xK = [prb.rK; 0.5*prb.vmax*ones(3,1); prb.pK; prb.y1];    

    prb.u1 = [-0.5*prb.gvec;
              (prb.rho1+prb.rho2)/(prb.mwet+prb.mdry)];
    prb.uK = [-0.5*prb.gvec;
              (prb.rho1+prb.rho2)/(prb.mwet+prb.mdry)];

    % Scaling parameters

    xmin =   [-0.5*prb.rmax*ones(3,1);-0.5*prb.vmax*ones(3,1); prb.pK; prb.ymin];
    xmax =   [ 0.5*prb.rmax*ones(3,1); 0.5*prb.vmax*ones(3,1); prb.p1; prb.ymax];

    prb.umin = [-prb.rho2*ones(3,1)/prb.mdry;0];
    prb.umax =  prb.rho2*ones(4,1)/prb.mdry;

    prb.scl_bnd = [0,1];
    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[prb.umin,prb.umax]},prb.scl_bnd);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};    
    
    cnstr_scl = diag([1; ...
                      1; ...
                      0.01; ...
                      1; ...
                      1; ...
                      1; ...
                      100; ...
                      1; ...
                      1000; ...
                      1000;
                      ]);
    cnstr_buffer = zeros(10,1);

    % Path constraints

    prb.cnstr_fun = @(tau,x,u) cnstr_scl*[norm(tand(prb.gam_gs)*x(1:2))^2 - x(3)^2;
                                          -x(3);
                                          norm(x(4:6))^2 - prb.vmax^2;
                                          x(7) - p_1(tau);
                                          -x(7) + p_min(tau);
                                          -prb.tp_vec'*u;
                                          norm(u(1:3))^2 - u(4)^2;
                                          -u(4);
                                          a_1(tau) + b_1(tau)*x(7) + c_1(tau)*x(7)^2 - u(4);
                                          b_2(tau)*x(7) + u(4) + a_2(tau);
                                          ] + cnstr_buffer;

    prb.cnstr_fun_jac_x = @(tau,x,u) cnstr_scl*[2*(tand(prb.gam_gs)^2)*x(1:2)', -2*x(3), zeros(1,4);
                                                0, 0, -1, zeros(1,4);
                                                zeros(1,3), 2*x(4:6)', 0;
                                                zeros(1,6), 1;
                                                zeros(1,6), -1;
                                                zeros(3,7);
                                                zeros(1,6), b_1(tau) + 2*c_1(tau)*x(7);
                                                zeros(1,6), b_2(tau);
                                                ];

    prb.cnstr_fun_jac_u = @(tau,x,u) cnstr_scl*[zeros(5,4);
                                                -prb.tp_vec';
                                                2*u(1:3)',-2*u(4);
                                                zeros(1,3),-1;
                                                zeros(1,3),-1;
                                                zeros(1,3), 1;
                                                ];

% SCP parameters

    prb.disc = "ZOH";
    prb.zoh_type = "v3_parallel";
    prb.ode_solver = {'ode45',odeset('RelTol',1e-4,'AbsTol',1e-5)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations

    prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','mosek','verbose',false,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-9,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-9);
    % prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.abstol',1e-8,'ecos.reltol',1e-8);   
    % prb.solver_settings = sdpsettings('solver','quadprog','verbose',false);
    % prb.solver_settings = sdpsettings('solver','osqp','verbose',false,'osqp.eps_abs',1e-7,'osqp.eps_rel',1e-7,'osqp.max_iter',5e4);        

    % prb.solver = struct('name',"quadprog",'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9,'Display','none');
    prb.solver = struct('name',"piqp",'verbose',0,'eps_abs',1e-8,'eps_rel',1e-8,'eps_duality_gap_rel',1e-8,'eps_duality_gap_abs',1e-8);
    % prb.solver = struct('name',"ecos",'verbose',false,'abstol',1e-8,'reltol',1e-8);
    % prb.solver = struct('name',"gurobi",'verbose',0,'OptimalityTol',1e-9,'FeasibilityTol',1e-9);
    % prb.solver = struct('name',"scs",'eps_abs',1e-9,'eps_rel',1e-9,'verbose',false);
    % prb.solver = struct('name',"mosek",'MSK_DPAR_INTPNT_QO_TOL_PFEAS',1e-9,'MSK_DPAR_INTPNT_QO_TOL_DFEAS',1e-9,'MSK_DPAR_INTPNT_QO_TOL_REL_GAP',1e-9);
    % prb.solver = struct('name',"osqp",'eps_abs',1e-8,'eps_rel',1e-8,'verbose',0,'max_iter',5e4);
    % prb.solver = struct('name',"pipg",'eps_abs',1e-9,'verbose',0,'max_iter',5e4,'rho',1.5,'lambda',0.05,'omega',100,'test_termination',500);    

    prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-7;
    prb.epstr = 5e-4;

    % Takes in unscaled data
    prb.time_of_maneuver = @(z,u) T;    
    prb.time_grid = @(tau,z,u) prb.tau;    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func = @(tau,xtil,u)             evaluate_dyn_func     (tau,xtil,u,prb.alpha,prb.gvec,prb.cnstr_fun);
    prb.dyn_func_linearize = @(tau,xtil,u)   evaluate_linearization(tau,xtil,u,prb.alpha,prb.gvec,prb.cnstr_fun,...
                                                                    prb.cnstr_fun_jac_x,prb.cnstr_fun_jac_u);    

end
function F = evaluate_dyn_func(tau,xtil,u,alph,gvec,cnstr_fun)
   x = xtil(1:7);
   v = x(4:6);
   
   cnstr_val = cnstr_fun(tau,x,u);
   
   f = [v;
        u(1:3) + gvec;
        -alph*u(4);
        ];
    
   F = [f;
        sum(arrayfun(@(y) max(0,y)^2,cnstr_val))];
end
function [A,B,w] = evaluate_linearization(tau,xtil,u,alph,gvec,cnstr_fun,cnstr_fun_jac_x,cnstr_fun_jac_u)
   x = xtil(1:7);
   % v = x(4:6);

   cnstr_val = cnstr_fun(tau,x,u);
   abs_cnstr_val = arrayfun(@(y) max(0,y),cnstr_val);
   cnstr_val_jac_x = cnstr_fun_jac_x(tau,x,u);
   cnstr_val_jac_u = cnstr_fun_jac_u(tau,x,u);

   dfdx = [zeros(3,3), eye(3), zeros(3,1);
           zeros(3,7);
           zeros(1,7)];

   dfdu = [zeros(3,4);
           eye(3) zeros(3,1);
           zeros(1,3) -alph;
           ];

   A = [dfdx, zeros(7,1);
        2*abs_cnstr_val'*cnstr_val_jac_x, 0];

   B = [dfdu;
        2*abs_cnstr_val'*cnstr_val_jac_u];
    
   F = evaluate_dyn_func(tau,xtil,u,alph,gvec,cnstr_fun);
    
   w = F - A*xtil - B*u;   
end