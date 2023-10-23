function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(xtil,util,prb,...
                                                    ~,~)
% r    = xtil(1:n)
% v    = xtil(n+1:2*n)
% p    = xtil(2*n+1)
% bet  = xtil(2*n+2)
% s    = util(1:n)
% s    = util(n+1)

    K = prb.K;

    % Unscaled variables
    bet = sdpvar(1,K);
    u   = sdpvar(prb.n,K);
    s   = sdpvar(1,K);

    for k = 1:K
        bet(k)   = prb.Sx(2*prb.n+2,2*prb.n+2)              *xtil(2*prb.n+2,k)        + prb.cx(2*prb.n+2);
        u(:,k)   = prb.Su(1:prb.n,1:prb.n)                  *util(1:prb.n,k)          + prb.cu(1:prb.n);        
        s(k)     = prb.Su(prb.n+1,prb.n+1)                  *util(prb.n+1,k)          + prb.cu(prb.n+1); 
    end
    r1  = prb.Sx(1:prb.n,1:prb.n)                  *xtil(1:prb.n,1)          + prb.cx(1:prb.n);
    v1  = prb.Sx(prb.n+1:2*prb.n,prb.n+1:2*prb.n)  *xtil(prb.n+1:2*prb.n,1)  + prb.cx(prb.n+1:2*prb.n);
    rK  = prb.Sx(1:prb.n,1:prb.n)                  *xtil(1:prb.n,K)          + prb.cx(1:prb.n);
    vK  = prb.Sx(prb.n+1:2*prb.n,prb.n+1:2*prb.n)  *xtil(prb.n+1:2*prb.n,K)  + prb.cx(prb.n+1:2*prb.n);    
    p1  = prb.Sx(2*prb.n+1,2*prb.n+1)              *xtil(2*prb.n+1,1)        + prb.cx(2*prb.n+1);
    pK  = prb.Sx(2*prb.n+1,2*prb.n+1)              *xtil(2*prb.n+1,K)        + prb.cx(2*prb.n+1);
    
    % Boundary conditions
    cnstr = [r1     == prb.r1;
             v1     == prb.v1;
             rK     == prb.rK;
             vK     == prb.vK;
             p1     == prb.p1;   
             bet(1) == prb.bet1];

    cost_fun = prb.cost_factor*pK;

    for k = 1:K    
        
        cnstr = [cnstr;
                 -prb.umax <= u(:,k) <= prb.umax;   % Pointwise thrust bound
                 prb.smin <= s(k) <= prb.smax];     % Lower and upper bounds on dilation factor

        if k < K
            cnstr = [cnstr;
                     bet(k+1) <= bet(k) + 1e-5];
        end

    end

    % Compute time of maneuver and constrain time step
    % ToF = 0;
    % switch prb.disc
    %     case "ZOH"
    %         for k = 1:prb.K-1
    %             ToF = ToF + prb.dtau(k)*s(k);
    %             cnstr = [cnstr; prb.dtmin <= prb.dtau(k)*s(k) <= prb.dtmax]; 
    %         end
    %     case "FOH"
    %         for k = 1:prb.K-1
    %             ToF = ToF + 0.5*prb.dtau(k)*(s(k+1)+s(k));
    %             cnstr = [cnstr; prb.dtmin <= 0.5*prb.dtau(k)*(s(k+1)+s(k)) <= prb.dtmax];
    %         end
    % end    
    % 
    % % Time of maneuver upper bound
    % cnstr = [cnstr;ToF <= prb.ToFmax];

    vc_cnstr = 0;

end