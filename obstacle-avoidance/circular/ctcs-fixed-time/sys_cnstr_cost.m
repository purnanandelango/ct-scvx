function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                    ~,~)

    K = prb.K;
    n = prb.n;

    % Unscaled variables
    r = x(1:n,:);
    v = x(n+1:2*n,:);
    if prb.cnstr_type == "exclusive-integrator-states"
        y = x(2*n+1+1:2*n+1+prb.m,:);
    elseif prb.cnstr_type == "single-integrator-state" 
        y = x(2*n+1+1,:);
    end
    p1 = x(2*n+1,1);
    % T = u(1:n,:);

    pK_scl = prb.invSx(2*n+1,2*n+1)*(x(2*n+1,K)-prb.cx(2*n+1));
    
    % Boundary conditions
    cnstr = [r(:,1) == prb.r1;
             v(:,1) == prb.v1;
             p1 == prb.p1;
             r(:,K) == prb.rK;
             v(:,K) == prb.vK;
             y(:,1) == 0];

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 -prb.Tmax <= u(:,k) <= prb.Tmax];
       
        if k < K
            cnstr = [cnstr;
                     y(:,k+1) <= y(:,k) + prb.eps_cnstr];                
        end        

    end  

    cost_fun = cost_fun + prb.cost_factor*pK_scl;

    vc_cnstr = 0;    

end