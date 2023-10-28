function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(xtil,util,prb,...
                                                    ~,~)
% r    = xtil(1:n)
% v    = xtil(n+1:2*n)
% y    = xtil(2*n+1)
% s    = util(1:n)
% s    = util(n+1)

    K = prb.K;

    % Scaling matrices
    Pr = prb.Sx(1:prb.n,1:prb.n);
    Pv = prb.Sx(prb.n+1:2*prb.n,prb.n+1:2*prb.n);
    Py = prb.Sx(2*prb.n+1,2*prb.n+1);
    Pu = prb.Su(1:prb.n,1:prb.n);
    Ps = prb.Su(prb.n+1,prb.n+1);

    for k = 1:K
        y(k)   = Py * xtil(2*prb.n+1,k)        + prb.cx(2*prb.n+1);
        u(:,k) = Pu * util(1:prb.n,k)          + prb.cu(1:prb.n);        
        s(k)   = Ps * util(prb.n+1,k)          + prb.cu(prb.n+1);
    end
    r1         = Pr * xtil(1:prb.n,1)          + prb.cx(1:prb.n);
    v1         = Pv * xtil(prb.n+1:2*prb.n,1)  + prb.cx(prb.n+1:2*prb.n);
    rK         = Pr * xtil(1:prb.n,K)          + prb.cx(1:prb.n);
    vK         = Pv * xtil(prb.n+1:2*prb.n,K)  + prb.cx(prb.n+1:2*prb.n);    
    
    % Boundary conditions
    cnstr = [Pr \ r1     == Pr \ prb.r1;
             Pv \ v1     == Pv \ prb.v1;
             Pr \ rK     == Pr \ prb.rK;
             Pv \ vK     == Pv \ prb.vK;
             (1/Py)*y(1) == (1/Py)*prb.y1];

    ToFfinal = 0.0;
    for k = 1:K-1
        ToFfinal = ToFfinal + 0.5*prb.dtau(k)*(s(k) + s(k+1));
    end

    cost_fun = prb.cost_factor*(1/Ps)*(ToFfinal);

    for k = 1:K    
        
        cnstr = [cnstr;
                 -Pu \ (prb.umax*ones(prb.nu-1,1)) <= Pu \ u(:,k) <= Pu \ (prb.umax*ones(prb.nu-1,1)); % thrust bound
                 (1/Ps)*s(k) >= (1/Ps)*(K-1)*prb.dtmin % lower bound on the dilation factor
                ];

        if k < K
            cnstr = [cnstr;
                     (1/Ps)*(1/2)*prb.dtau(k)*(s(k) + s(k+1)) <= (1/Ps)*prb.dtmax % upper bound on the dilation factor
                    ];

            cnstr = [cnstr;
                     (1/Py)*y(k+1) <= (1/Py)*y(k) + (1/Py)*1e-5]; % relaxation
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