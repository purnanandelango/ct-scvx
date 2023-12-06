function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    ~,~)

    K = prb.K;
    n = prb.n;

    % Unscaled variables
    r = x(1:n,:);
    v = x(n+1:2*n,:);
    y = x(2*n+1:2*n+prb.m,:);
    T = u(1:n,:);
    s = u(n+1,:);
    
    % Boundary conditions
    cnstr = [r(:,1) == prb.r1;
             v(:,1) == prb.v1;
             r(:,K) == prb.rK;
             v(:,K) == prb.vK;
             y(:,1) == zeros(prb.m,1)];

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 norm(T(:,k)) <= prb.Tmax;          % Thrust magnitude upper bound
                 prb.smin <= s(k) <= prb.smax];     % Lower and upper bounds on dilation factor

        if k < K
            cnstr = [cnstr;
                     y(:,k+1) <= y(:,k) + prb.eps_cnstr];
        end

    end  

    cost_fun = cost_fun ...
               + prb.cost_factor*prb.invSu(1,1)*norm(T(:));

    % Compute time of maneuver and constrain time step
    cnstr = [cnstr;
             misc.time_cnstr(s,prb.dtau,{prb.dtmin,prb.dtmax,prb.ToFmax},prb.disc)];

    vc_cnstr = 0;

end