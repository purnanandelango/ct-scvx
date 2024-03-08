function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    ~,~)

    K = prb.K;
    n = prb.n;

    r  = x(1:n,:);
    v  = x(n+1:2*n,:);
    p  = x(2*n+1,:);
    t  = x(2*n+2,:);
    y  = x(2*n+3,:);
    
    % Boundary conditions
    cnstr = [r(:,1) == prb.r1;
             v(:,1) == prb.v1;
             r(:,K) == prb.rK;
             v(:,K) == prb.vK;
             p(1)   == prb.p1;   
             t(1)   == prb.t1; 
             y(1)   == prb.y1;
             ];

    pK_scl = prb.invSx(2*n+1,2*n+1)*(p(K) - prb.cx(2*n+1));

    cost_fun = prb.cost_factor*pK_scl;

    for k = 1:K    
        
        cnstr = [ cnstr;
                  prb.umin <= u(:,k) <= prb.umax];

        if k < K
            cnstr = [cnstr;
                     y(k+1) - y(k) <= prb.eps_cnstr];
        end

    end

    vc_cnstr = 0;

end