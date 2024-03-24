function prb = problem_data_lunar(K,scp_iters,wvc,wtr,cost_factor)
    
    prb.K = K;

    prb.ny = 1;
    prb.nx = 14+prb.ny;
    prb.nu = 4;
    
    prb.tau = grid.generate_grid(0,1,K,'uniform'); % Generate grid in [0,1]

    prb.dtau = diff(prb.tau);
    
    prb.h = (1/19)*prb.dtau;                 % Step size for integration that computes discretization matrices
    prb.Kfine = 1+round(20/min(prb.dtau));   % Size of grid on which SCP solution is simulated
    
    % System parameters

    prb.g0 = 1.61;
    
    % Not used
    prb.c_ax = 0; prb.c_ayz = 0;    
    prb.CA = diag([prb.c_ax,prb.c_ayz,prb.c_ayz]);
    prb.rho = 1;
    prb.SA = 0.5;    
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
    
    prb.vmax_stc      = 25;
    prb.vmax_stc_aug  = 25-3;     % Trigger conservatively for better numerics
    prb.cosaoamax     = cosd( 10 ); 
    prb.cosaoamax_aug = cosd( 7 );
    prb.stc_flag      = "v1";
    prb.trig_scl      = 0.1;

    prb.mdry    = 2100;
    prb.mwet    = 3250;

    prb.smin    = 1;
    prb.smax    = 60;
    % prb.dtmin   = 0.1;
    % prb.dtmax   = 10;
    % prb.ToFmax  = 30;
   
    prb.snom    = [1,15];
    prb.ToFguess= 30;   

    prb.ymin = 0*ones(prb.ny,1);
    prb.ymax = 1*ones(prb.ny,1);

    prb.eps_cnstr = 1e-4;
    
    % Boundary conditions

    prb.rI1     = [433;0;250];
    prb.vI1     = [10;0;-30];
    
    prb.rIK     = [30;0;0];
    prb.vIK     = [-1;0;0];
    prb.omgB1   = zeros(3,1);
    prb.omgBK   = zeros(3,1);
    prb.q1      = [0;0;0;1];

    prb.y1      = zeros(prb.ny,1);
    prb.yK      = zeros(prb.ny,1);

    prb.Ey      = [zeros(prb.ny,prb.nx-prb.ny),eye(prb.ny)];

    Inx = eye(prb.nx);
    prb.Ei      = Inx([1:7,12:14+prb.ny],:);
    prb.zi      = [prb.mwet;prb.rI1;prb.vI1;prb.omgB1;prb.y1];
    prb.i_idx   = [1:7,12:14+prb.ny];
    prb.Ef      = Inx(2:14,:);
    prb.zf      = [prb.rIK;prb.vIK;prb.q1;prb.omgBK];
    prb.f_idx   = 2:14;    
    prb.term_cost_vec = -[1;zeros(13+prb.ny,1)];    

    % Straight-line initialization
    prb.x1      = [prb.mwet;prb.rI1;prb.vI1;prb.q1;prb.omgB1;prb.y1];
    prb.xK      = [prb.mdry;prb.rIK;prb.vIK;prb.q1;prb.omgBK;prb.yK];    
    prb.u1      = [-prb.mwet*prb.gI;prb.ToFguess];
    prb.uK      = [-prb.mdry*prb.gI;prb.ToFguess];

    % Scaling parameters
    xmin = [prb.mdry; -400*ones(3,1); -100*ones(3,1);  -ones(4,1); -prb.omgmax*ones(3,1); prb.ymin];
    xmax = [prb.mwet;  400*ones(3,1);  100*ones(3,1);   ones(4,1);  prb.omgmax*ones(3,1); prb.ymax];
    
    umin = [-prb.TBmax*ones(3,1); prb.smin]; prb.umin = umin;
    umax = [ prb.TBmax*ones(3,1); prb.smax]; prb.umax = umax;  

    prb.scl_bnd = [0,1];
    [Sz,cz] = misc.generate_scaling({[xmin,xmax],[umin,umax]},prb.scl_bnd);
    % [Sz,cz] = misc.generate_scaling({[zeros(prb.nx,1),xmax],[zeros(prb.nu,1),umax]},prb.scl_bnd);

    prb.Sx = Sz{1}; prb.invSx = inv(Sz{1});
    prb.Su = Sz{2}; prb.invSu = inv(Sz{2});
    prb.cx = cz{1};
    prb.cu = cz{2};

    % Constraint parameters

    prb.cnstr_scl = diag([ ...
                          1e-3;
                          1e-5;
                          1e-2;
                          1e-5;
                          1000;
                          1;
                          1e-5;
                          1e-3;
                          1e-5;
                          1e-5;
                          ]);
    prb.cnstr_buffer = [...
                        0;
                        0;
                        0;
                        0;
                        0;
                        0;
                        0;
                        0;
                        0;
                        0;
                        ];

    prb.cnstr_fun       = @(xi,TB) prb.cnstr_scl*[ prb.mdry - xi(1);
                                                   norm(prb.cotgamgs*prb.Hgam*xi(2:4))^2 - xi(2)^2;
                                                  -xi(2);
                                                   norm(xi(5:7))^2 - prb.vmax^2;
                                                   norm(prb.Hthet*xi(8:11))^2 - prb.sinthetmaxby2^2;
                                                   norm(xi(12:14))^2 - prb.omgmax^2;
                                                   norm(prb.cosdelmax*TB)^2 - TB(1)^2;
                                                  -TB(1);
                                                   norm(TB)^2 - prb.TBmax^2;
                                                   prb.TBmin^2 - norm(TB)^2] ...
                                                 + prb.cnstr_buffer;

    prb.cnstr_fun_jac_xi = @(xi,TB) prb.cnstr_scl*[-1, zeros(1,13);
                                                    0, 2*(prb.cotgamgs^2)*(prb.Hgam'*prb.Hgam*xi(2:4))'-[2*xi(2) 0 0], zeros(1,10);
                                                    0, -1, zeros(1,12);
                                                    0, zeros(1,3),  2*xi(5:7)', zeros(1,7);
                                                    0, zeros(1,6),  2*(prb.Hthet'*prb.Hthet*xi(8:11))', zeros(1,3);
                                                    0, zeros(1,10), 2*xi(12:14)';
                                                    zeros(1,14);
                                                    zeros(1,14);
                                                    zeros(1,14);
                                                    zeros(1,14)];

    prb.cnstr_fun_jac_TB = @(xi,TB) prb.cnstr_scl*[ zeros(1,3);
                                                    zeros(1,3);
                                                    zeros(1,3);
                                                    zeros(1,3);
                                                    zeros(1,3);
                                                    zeros(1,3);
                                                    2*(prb.cosdelmax^2)*TB'-[2*TB(1) 0 0];
                                                    -1, 0, 0;
                                                    2*TB';
                                                   -2*TB'];

    prb.stc_fun = @(xi,TB) plant.rocket6DoF.q_aoa_cnstr(xi(5:7),xi(8:11),prb.vmax_stc_aug,prb.cosaoamax_aug,prb.stc_flag);

    % SCP parameters

    prb.disc = "FOH";
    prb.foh_type = "v3_parallel";
    prb.ode_solver = {'ode45',odeset('RelTol',1e-5,'AbsTol',1e-7)};
    prb.scp_iters = scp_iters; % Maximum SCP iterations
    
    % prb.solver_settings = sdpsettings('solver','quadprog','verbose',false);    
    prb.solver_settings = sdpsettings('solver','ecos','verbose',false,'ecos.abstol',1e-8,'ecos.reltol',1e-8);    
    % prb.solver_settings = sdpsettings('solver','gurobi','verbose',false,'gurobi.OptimalityTol',1e-9,'gurobi.FeasibilityTol',1e-9);
    % prb.solver_settings = sdpsettings('solver','mosek','verbose',false,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-9,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-9);
    % prb.solver_settings = sdpsettings('solver','osqp','verbose',false,'osqp.eps_abs',1e-8,'osqp.eps_rel',1e-8,'osqp.max_iter',5e4);        
   
    % prb.solver = struct('name',"quadprog",'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9,'Display','none');
    prb.solver = struct('name',"piqp",'verbose',0,'eps_abs',1e-8,'eps_rel',1e-8,'eps_duality_gap_rel',1e-8,'eps_duality_gap_abs',1e-8);
    % prb.solver = struct('name',"ecos",'verbose',false,'abstol',1e-8,'reltol',1e-8);
    % prb.solver = struct('name',"gurobi",'verbose',0,'OptimalityTol',1e-9,'FeasibilityTol',1e-9);
    % prb.solver = struct('name',"scs",'eps_abs',1e-9,'eps_rel',1e-9,'verbose',false);
    % prb.solver = struct('name',"mosek",'MSK_DPAR_INTPNT_QO_TOL_PFEAS',1e-9,'MSK_DPAR_INTPNT_QO_TOL_DFEAS',1e-9,'MSK_DPAR_INTPNT_QO_TOL_REL_GAP',1e-9);
    % prb.solver = struct('name',"osqp",'eps_abs',1e-8,'eps_rel',1e-8,'verbose',0,'max_iter',5e4);
    % prb.solver = struct('name',"pipg",'eps_abs',1e-9,'verbose',0,'max_iter',5e4,'rho',1.5,'lambda',0.05,'omega',100,'test_termination',500);

    % prb.tr_norm = 2;
    % prb.tr_norm = inf;
    prb.tr_norm = 'quad';
    
    prb.wvc = wvc;
    prb.wtr = wtr;
    prb.cost_factor = cost_factor;
    
    prb.epsvc = 1e-7;
    prb.epstr = 1e-3;

    % Takes in unscaled data
    prb.time_of_maneuver =     @(x,u) disc.time_of_maneuver(prb.disc,prb.tau,u(4,:));    
    prb.time_grid        = @(tau,x,u)        disc.time_grid(prb.disc,    tau,u(4,:));    
    
    % Convenient functions for accessing RHS of nonlinear and linearized ODE
    prb.dyn_func           = @(tau,x,u)      evaluate_dyn_func(x,u,prb.cnstr_fun,                                          prb.stc_fun,prb.trig_scl,prb.c_ax,prb.c_ayz,diag(prb.JB),prb.gI,prb.rho,prb.SA,prb.rTB,prb.rcpB,prb.alphmdt,prb.betmdt);
    prb.dyn_func_linearize = @(tau,x,u) evaluate_linearization(x,u,prb.cnstr_fun,prb.cnstr_fun_jac_xi,prb.cnstr_fun_jac_TB,prb.stc_fun,prb.trig_scl,prb.c_ax,prb.c_ayz,diag(prb.JB),prb.gI,prb.rho,prb.SA,prb.rTB,prb.rcpB,prb.alphmdt,prb.betmdt);
end

function f = evaluate_dyn_func(x,u,cnstr_fun,stc_fun,trig_scl,c_ax,c_ayz,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt)

    xi = x(1:14);
    TB = u(1:3);
    s  = u(4);

    % \derv{\xi} = sF(\xi,T_B)
    sF = plant.rocket6DoF.dyn_func(xi,TB,s,c_ax,c_ayz,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
    
    cnstr_val = cnstr_fun(xi,TB);

    f = [sF;
         s*sum( arrayfun(@(y) max(0,y),cnstr_val) .^ 2 )];

    if length(x) == 16
        [~,~,stc_trig,~,stc_cnstr] = stc_fun(xi,TB);
        min_stc_trig = min(0,stc_trig);
        max_stc_cnstr = max(0,stc_cnstr);

        stc_val = s*((trig_scl*min_stc_trig*max_stc_cnstr)^2);
        % stc_val = s*( (min_stc_trig^2)*(max_stc_cnstr^2)/(trig_scl+min_stc_trig^2) );

        f(end+1) = stc_val;
    end
end

function [A,B,w] = evaluate_linearization(x,u,cnstr_fun,cnstr_fun_jac_xi,cnstr_fun_jac_TB,stc_fun,trig_scl,c_ax,c_ayz,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt)

    xi = x(1:14);
    TB = u(1:3);
    s  = u(4);
    
    [dsF_dxi,dsF_dTB,dsF_ds] = plant.rocket6DoF.compute_linearization(xi,TB,s,c_ax,c_ayz,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
    dsF_dTB = [dsF_dTB,dsF_ds];

    cnstr_val = cnstr_fun(xi,TB);
    cnstr_val_jac_xi = cnstr_fun_jac_xi(xi,TB);
    cnstr_val_jac_TB = cnstr_fun_jac_TB(xi,TB);
    [~,~,stc_trig,dstc_trig_dxi,stc_cnstr,dstc_cnstr_dxi] = stc_fun(xi,TB);
    min_stc_trig = min(0,stc_trig);
    max_stc_cnstr = max(0,stc_cnstr);

    ny = length(cnstr_val);
    max_cnstr_2_jac_xi = zeros(ny,14);
    max_cnstr_2_jac_TB = zeros(ny,3);
    for j = 1:ny        
        max_cnstr_2_jac_xi(j,:) = 2*max(0,cnstr_val(j))*cnstr_val_jac_xi(j,:);
        max_cnstr_2_jac_TB(j,:) = 2*max(0,cnstr_val(j))*cnstr_val_jac_TB(j,:);
    end

    if length(x) == 15

        A = [dsF_dxi                        zeros(14,1);
             s*sum(max_cnstr_2_jac_xi,1)    0];
        B = [dsF_dTB;
             s*sum(max_cnstr_2_jac_TB,1), sum( arrayfun(@(y) max(0,y)^2, cnstr_val) )];

    else

        dstc_dxi = (trig_scl^2)*s*( 2*min_stc_trig*(max_stc_cnstr^2)*dstc_trig_dxi + 2*(min_stc_trig^2)*max_stc_cnstr*dstc_cnstr_dxi );
        dstc_dTB = [zeros(1,3), (trig_scl*min_stc_trig*max_stc_cnstr)^2];
        
        % dstc_dxi = s*( 2*trig_scl*min_stc_trig*(max_stc_cnstr^2)*dstc_trig_dxi/(trig_scl+min_stc_trig^2)^2 + 2*(min_stc_trig^2)*max_stc_cnstr*dstc_cnstr_dxi/(trig_scl+min_stc_trig^2) );
        % dstc_dTB = [zeros(1,3), (min_stc_trig^2)*(max_stc_cnstr^2)/(trig_scl+min_stc_trig^2) ];
    
        A = [dsF_dxi                        zeros(14,2);
             s*sum(max_cnstr_2_jac_xi,1)    zeros(1,2);
             dstc_dxi                       zeros(1,2)];
        B = [dsF_dTB;
             s*sum(max_cnstr_2_jac_TB,1), sum( arrayfun(@(y) max(0,y)^2, cnstr_val) );
             dstc_dTB];

    end

    f = evaluate_dyn_func(x,u,cnstr_fun,stc_fun,trig_scl,c_ax,c_ayz,JB,gI,rho,SA,rTB,rcpB,alphmdt,betmdt);
    w = f - A*x - B*u;    
end