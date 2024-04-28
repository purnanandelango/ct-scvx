function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    xbar,ubar)

    K = prb.K;

    % Unscaled variables
    r    = x(1:prb.n,:);
    v    = x(prb.n+1:2*prb.n,:);
    T    = u(1:prb.n,:);

    % Obstacle avoidance buffer
    nu_ncvx = sdpvar(prb.nobs+1,K);
    
    % Boundary conditions
    cnstr = [r(:,1)   == prb.r1;
             v(:,1)   == prb.v1;
             r(:,K)   == prb.rK;
             v(:,K)   == prb.vK];

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 norm(T(:,k)) <= prb.Tmax;              % Thrust magnitude upper bound
                 norm(v(:,k)) <= prb.vmax];             % Velocity magnitude upper bound

        for j = 1:prb.nobs
            cnstr = [cnstr;
                     norm(xbar(1:prb.n,k)-prb.robs(:,j)) - prb.aobs(j) + dot(xbar(1:prb.n,k)-prb.robs(:,j),r(:,k)-xbar(1:prb.n,k))/norm(xbar(1:prb.n,k)-prb.robs(:,j)) + nu_ncvx(j,k) >= 0;
                     nu_ncvx(j,k) >= 0];
        end
        cnstr = [cnstr;
                norm(ubar(1:prb.n,k)) - prb.Tmin + dot(ubar(1:prb.n,k),T(:,k)-ubar(1:prb.n,k))/norm(ubar(1:prb.n,k)) + nu_ncvx(end,k) >= 0;
                nu_ncvx(end,k) >= 0];

    end  

    vc_cnstr = sum(nu_ncvx(:));    

    cost_fun = cost_fun ...
               + prb.cost_factor*prb.invSu(1,1)*norm(T(:)) ...
               + prb.w_ep*vc_cnstr;

end