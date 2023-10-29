function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(xtil_scl,u_scl,prb,...
                                                    ~,~)

% r = xtil(1:3)
% v = xtil(4:6)
% p = xtil(7)
% y = xtil(8)

    K = prb.K;

    % Unscaled state control input
    r1      = prb.Sx(1:3,1:3)   *xtil_scl(1:3,1) + prb.cx(1:3);
    rK      = prb.Sx(1:3,1:3)   *xtil_scl(1:3,K) + prb.cx(1:3);
    v1      = prb.Sx(4:6,4:6)   *xtil_scl(4:6,1) + prb.cx(4:6);
    vK      = prb.Sx(4:6,4:6)   *xtil_scl(4:6,K) + prb.cx(4:6);
    p1      = prb.Sx(7,7)       *xtil_scl(7,1)   + prb.cx(7);
    pK_scl  = xtil_scl(7,K);
    y       = prb.Sx(8,8)       *xtil_scl(8,:)   + repmat(prb.cx(8),[1,K]);
    u       = prb.Su            *u_scl           + repmat(prb.cu,[1,K]);
    
    cnstr = [];
    cost_fun = prb.cost_factor*pK_scl;

    % Boundary conditions
    cnstr = [cnstr;
             r1     == prb.r1;
             rK     == prb.rK;
             v1     == prb.v1;
             vK     == prb.vK;
             p1     == prb.p1;
             y(1)   == prb.y1];

    for k=1:K
    
        % Face normal force lower-bound
        cnstr = [cnstr;
                 -u(prb.idx{1}(1),k) + prb.Fmin <= 0;
                  u(prb.idx{2}(1),k) + prb.Fmin <= 0;
                 -u(prb.idx{3}(2),k) + prb.Fmin <= 0;
                  u(prb.idx{4}(2),k) + prb.Fmin <= 0;
                ];

        if k < K
        cnstr = [cnstr;
                 y(k+1) - y(k) <= 1e-5];
        end
    
    end

    vc_cnstr = 0;

end