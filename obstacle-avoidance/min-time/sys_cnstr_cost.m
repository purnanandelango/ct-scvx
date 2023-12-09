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
        y(k)   = xtil(2*prb.n+1,k);
        u(:,k) = util(1:prb.n,k);        
        s(k)   = util(prb.n+1,k);
    end
    r1         = xtil(1:prb.n,1);
    v1         = xtil(prb.n+1:2*prb.n,1);
    rK         = xtil(1:prb.n,K);
    vK         = xtil(prb.n+1:2*prb.n,K);
    
    % Boundary conditions
    cnstr = [Pr \ r1     == Pr \ prb.r1;
             Pv \ v1     == Pv \ prb.v1;
             Pr \ rK     == Pr \ prb.rK;
             Pv \ vK     == Pv \ prb.vK;
             (1/Py)*y(1) == (1/Py)*prb.y1
             ];

    ToFfinal = 0.0;
    for k = 1:K-1
        ToFfinal = ToFfinal + 0.5*prb.dtau(k)*(s(k) + s(k+1));
    end

    cost_fun = prb.cost_factor*(1/Ps)*(ToFfinal);

    for k = 1:K   

        % Conservatively approximated constraints imposed in practice 
        % (the resulting solution will never violate the actual constraint)
        
        % Approximation here is to ensure that the dilation factor is
        % always nonnegative
        cnstr = [cnstr;
                 % -Pu \ (prb.Tmin*ones(prb.nu-1,1)) <= Pu \ u(:,k) <= Pu \ (prb.Tmax*ones(prb.nu-1,1)); % thrust bound
                 (1/Ps)*s(k) >= (1/Ps)*prb.smin % lower bound on the dilation factor
                ];

        % Approximation here is to avoid temporal entanglement of the
        % constraint sets
        cnstr = [cnstr;
                 (1/Ps)*s(k) <= (1/Ps)*prb.smax % upper bound on the dilation factor
                ];

    end

    for k = 1:K-1

        % Actual dilation constraints

        % cnstr = [cnstr;
        %          (1/Ps)*(1/2)*(s(k) + s(k+1)) <= (1/Ps)*(1/prb.dtau(k))*prb.dtmax % upper bound on the dilation factor
        %         ];

        % cnstr = [cnstr;
        %          (1/Ps)*(1/2)*(s(k) + s(k+1)) >= (1/Ps)*(1/prb.dtau(k))*prb.dtmin % lower bound on the dilation factor
        %         ];

        cnstr = [cnstr;
                 (1/Py)*y(k+1) <= (1/Py)*y(k) + (1/Py)*prb.eps_cnstr]; % relaxation

        % cnstr = [cnstr;
        %          (1/Py)*y(k+1) == (1/Py)*y(k)]; % no relaxation

    end

    vc_cnstr = 0;

end