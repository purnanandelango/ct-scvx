function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    xbar,ubar)

    K = prb.K;

    % Unscaled variables
    r    = x(1:prb.n,:);
    v    = x(prb.n+1:2*prb.n,:);
    T    = u(1:prb.n,:);
    s    = u(prb.n+1,:);

    % Obstacle avoidance buffer
    nu_ncvx = sdpvar(prb.n_obs+1,K);
    
    % Boundary conditions
    cnstr = [r(:,1) == prb.r1;
             v(:,1) == prb.v1;
             r(:,K) == prb.rK;
             v(:,K) == prb.vK];

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 norm(T(:,k)) <= prb.Tmax;              % Thrust magnitude upper bound
                 norm(v(:,k)) <= prb.vmax;              % Velocity magnitude upper bound
                 prb.smin <= s(k) <= prb.smax];         % Lower and upper bounds on dilation factor

        for j = 1:prb.n_obs
            cnstr = [cnstr;
                     norm(prb.H_obs{j}*(xbar(1:prb.n,k)-prb.q_obs(:,j)))^2 - 1 + ...
                     2*dot( (xbar(1:prb.n,k) - prb.q_obs(:, j))'*prb.H_obs{j}'*prb.H_obs{j}, ...
                     r(:,k) - xbar(1:prb.n,k) ) ...
                     + nu_ncvx(j,k) >= 0;
                     nu_ncvx(j,k) >= 0];
        end
        cnstr = [cnstr;
                 norm(ubar(1:prb.n,k)) - prb.Tmin + dot(ubar(1:prb.n,k),T(:,k)-ubar(1:prb.n,k))/norm(ubar(1:prb.n,k)) + nu_ncvx(end,k) >= 0;
                 nu_ncvx(end,k) >= 0];

        cost_fun = cost_fun + prb.cost_factor*norm(prb.invSu(1:prb.n,1:prb.n)*(T(:,k)-prb.cu(1:prb.n)))^2;
    
    end  

    vc_cnstr = sum(nu_ncvx(:));    

    cost_fun = cost_fun + prb.w_ep*vc_cnstr;

end