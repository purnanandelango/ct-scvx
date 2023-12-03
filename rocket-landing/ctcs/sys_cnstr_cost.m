function [cnstr,cost_fun,vcvb_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                      ~,~)

    K = prb.K;
    ny = prb.ny;

    m     = x(1,:);
    rI    = x(2:4,:);
    vI    = x(5:7,:);
    qBI   = x(8:11,:); 
    omgB  = x(12:14,:);
    y     = x(14+1:14+ny,:);
    % TB    = u(1:3,:);
    % s     = u(4,:);

    % Boundary conditions
    cnstr = [m(1)      == prb.mwet;
             rI(:,1)   == prb.rI1;
             vI(:,1)   == prb.vI1;
             rI(:,K)   == prb.rIK;
             vI(:,K)   == prb.vIK;
             qBI(:,K)  == prb.q1;
             omgB(:,1) == prb.omgB1;
             omgB(:,K) == prb.omgBK;
             y(:,1)    == zeros(ny,1)];

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 prb.umin <= u(:,k) <= prb.umax;
                ];

        if k < K
            cnstr = [cnstr;
                     y(:,k+1) <= y(:,k) + prb.eps_cnstr];
        end        
    
    end  

    cost_fun = cost_fun + prb.cost_factor*prb.invSx(1,1)*(x(1,K)-prb.cx(1));    

    % Compute time of maneuver and constrain time step
    % cnstr = [cnstr;
    %          misc.time_cnstr(s,prb.dtau,{prb.dtmin,prb.dtmax,prb.ToFmax},prb.disc)];

    vcvb_cnstr = 0;

end