clc
close all
clearvars

prb = problem_data_cvx(20,09);

yalmip clear

x_scl = sdpvar(prb.nx,prb.K);
u_scl = sdpvar(prb.nu,prb.K);

x = prb.Sx*x_scl + repmat(prb.cx,[1,prb.K]);
u = prb.Su*u_scl + repmat(prb.cu,[1,prb.K]);

cnstr = [];
cost_fun = 0;
for k = 1:prb.K-1

    cost_fun = cost_fun + norm(u(prb.idx{1},k))^2 ...
                            + norm(u(prb.idx{2},k))^2 ...
                            + norm(u(prb.idx{3},k))^2 ...
                            + norm(u(prb.idx{4},k))^2;    
    
    cnstr = [cnstr;
             x(:,k+1) == prb.Ad*x(:,k) + prb.Bd*u(:,k) + prb.wd;
             norm(x(4:6,k)) <= prb.vmax;
             prb.contact_mat*u(:,k) == 0;
             norm(u(prb.idx{1}(2:3),k))   <=  prb.mu(1)*u(prb.idx{1}(1),k);
             norm(u(prb.idx{2}(2:3),k))   <= -prb.mu(2)*u(prb.idx{2}(1),k);
             norm(u(prb.idx{3}([1,3]),k)) <=  prb.mu(3)*u(prb.idx{3}(2),k);
             norm(u(prb.idx{4}([1,3]),k)) <= -prb.mu(4)*u(prb.idx{4}(2),k);
             norm(u(prb.idx{1},k)) <= prb.F1max;
             norm(u(prb.idx{2},k)) <= prb.F2max;
             norm(u(prb.idx{3},k)) <= prb.F3max;
             norm(u(prb.idx{4},k)) <= prb.F4max;
             u(prb.idx{1}(1),k) >= prb.Fmin;
            -u(prb.idx{2}(1),k) >= prb.Fmin;
             u(prb.idx{3}(2),k) >= prb.Fmin;
            -u(prb.idx{4}(2),k) >= prb.Fmin;
             ];

end
cnstr = [cnstr;
         u(:,prb.K) == u(:,prb.K-1)];
cnstr = [cnstr;
         x(1:3,1)     == prb.r1;
         x(1:3,prb.K) == prb.rK;
         x(4:6,1)     == prb.v1;
         x(4:6,prb.K) == prb.vK;
         ];

optimize(cnstr,cost_fun,prb.solver_settings);

x = value(x);
u = value(u);
cost_val = value(cost_fun)*prb.dtau;

cum_cost(1,prb.K) = 0;
for k = 1:prb.K-1
    cum_cost(k+1) = cum_cost(k) + prb.dtau*(norm(u(prb.idx{1},k))^2 ...
                                            + norm(u(prb.idx{2},k))^2 ...
                                            + norm(u(prb.idx{3},k))^2 ...
                                            + norm(u(prb.idx{4},k))^2);
end

x_guess = [x;
           cum_cost;
           zeros(1,prb.K)];

RotEq(1,prb.K-1) = 0;
Fcone(4,prb.K-1) = 0;
Fmag(4,prb.K-1) = 0;
Vmag(1,prb.K-1) = 0;
for j=1:prb.K-1
    RotEq(j)    = norm(prb.contact_mat*u(:,k),Inf);                                 % Rotational equilibrium
    Fcone(1,j)  = norm(u(prb.idx{1}(2:3),k))   - prb.mu(1)*u(prb.idx{1}(1),k);              % Friction cone
    Fcone(2,j)  = norm(u(prb.idx{2}(2:3),k))   - -prb.mu(2)*u(prb.idx{2}(1),k);             % Friction cone            
    Fcone(3,j)  = norm(u(prb.idx{3}([1,3]),k)) - prb.mu(3)*u(prb.idx{3}(2),k);              % Friction cone            
    Fcone(4,j)  = norm(u(prb.idx{4}([1,3]),k)) - -prb.mu(4)*u(prb.idx{4}(2),k);
    Fmag(1,j)   = norm(u(prb.idx{1},k)) - prb.F1max;                                    % Grasp force magnitude
    Fmag(2,j)   = norm(u(prb.idx{2},k)) - prb.F2max;                                    % Grasp force magnitude
    Fmag(3,j)   = norm(u(prb.idx{3},k)) - prb.F3max;                                    % Grasp force magnitude
    Fmag(4,j)   = norm(u(prb.idx{4},k)) - prb.F4max;                                    % Grasp force magnitude
    Vmag(j)     = norm(x(4:6,k));                                        % Velocity magnitude                  
end
fprintf("\nConstraint Activity:\n");
fprintf("Rotational eqb.: %.3f\n",max(abs(RotEq)));
fprintf("Friction cone 1: %.3f\n",max(Fcone(1,:)));
fprintf("Friction cone 2: %.3f\n",max(Fcone(2,:)));
fprintf("Friction cone 3: %.3f\n",max(Fcone(3,:)));
fprintf("Friction cone 4: %.3f\n",max(Fcone(4,:)));
fprintf("Friction mag. 1: %.3f\n",max(Fmag(1,:)));
fprintf("Friction mag. 2: %.3f\n",max(Fmag(2,:)));
fprintf("Friction mag. 3: %.3f\n",max(Fmag(3,:)));
fprintf("Friction mag. 4: %.3f\n",max(Fmag(4,:)));
fprintf("Velocity mag.  : %.3f\n",max(Vmag));

save('recent_solution_cvx','x','u','prb','cost_val','x_guess');

% plot_solution_cvx;

