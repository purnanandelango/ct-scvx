function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(xtil,util,prb,...
                                                    ~,~)
% r    = xtil(1:n)
% v    = xtil(n+1:2*n)
% p    = xtil(2*n+1)
% y    = xtil(2*n+2)
% s    = util(1:n)
% s    = util(n+1)

    K = prb.K;

    % Unscaled variables
    y   = sdpvar(1,K);
    u   = sdpvar(prb.n,K);
    s   = sdpvar(1,K);

    Sscl_r = prb.Sx(1:prb.n,1:prb.n);                  cscl_r = prb.cx(1:prb.n);
    Sscl_v = prb.Sx(prb.n+1:2*prb.n,prb.n+1:2*prb.n);  cscl_v = prb.cx(prb.n+1:2*prb.n);
    Sscl_p = prb.Sx(2*prb.n+1,2*prb.n+1);              cscl_p = prb.cx(2*prb.n+1);              
    Sscl_y = prb.Sx(2*prb.n+2,2*prb.n+2);              cscl_y = prb.cx(2*prb.n+2);   
    Sscl_u = prb.Su(1:prb.n,1:prb.n);                  cscl_u = prb.cu(1:prb.n);
    Sscl_s = prb.Su(prb.n+1,prb.n+1);                  cscl_s = prb.cu(prb.n+1);  

    for k = 1:K
        y(k)     = Sscl_y*xtil(2*prb.n+2,k) + cscl_y;
        u(:,k)   = Sscl_u*util(1:prb.n,k)   + cscl_u;        
        s(k)     = Sscl_s*util(prb.n+1,k)   + cscl_s; 
    end
    r1  = Sscl_r*xtil(1:prb.n,1)            + cscl_r;
    v1  = Sscl_v*xtil(prb.n+1:2*prb.n,1)    + cscl_v;
    rK  = Sscl_r*xtil(1:prb.n,K)            + cscl_r;
    vK  = Sscl_v*xtil(prb.n+1:2*prb.n,K)    + cscl_v;    
    p1  = Sscl_p*xtil(2*prb.n+1,1)          + cscl_p;
    pK_scl  = xtil(2*prb.n+1,K);
    
    % Boundary conditions
    cnstr = [Sscl_r \ (r1 - cscl_r)     == Sscl_r \ (prb.r1 - cscl_r);
             Sscl_v \ (v1 - cscl_v)     == Sscl_v \ (prb.v1 - cscl_v);
             Sscl_r \ (rK - cscl_r)     == Sscl_r \ (prb.rK - cscl_r);
             Sscl_v \ (vK - cscl_v)     == Sscl_v \ (prb.vK - cscl_v);
             Sscl_p \ (p1 - cscl_p)     == Sscl_p \ (prb.p1 - cscl_p);   
             Sscl_y \ (y(1) - cscl_y)   == Sscl_y \ (prb.y1 - cscl_y);];

    cost_fun = prb.cost_factor*pK_scl;

    for k = 1:K    
        
        cnstr = [ cnstr;
                 -Sscl_u \ (prb.umax - cscl_u) <= util(1:prb.n,k) <= Sscl_u \ (prb.umax - cscl_u);      % Pointwise thrust bound
                  Sscl_s \ (prb.smin - cscl_s) <= util(prb.n+1,k) <= Sscl_s \ (prb.smax - cscl_s)];     % Lower and upper bounds on dilation factor

        if k < K
            cnstr = [cnstr;
                     Sscl_y \ (y(k+1) - y(k)) <= 1e-5];
        end

    end

    vc_cnstr = 0;

end