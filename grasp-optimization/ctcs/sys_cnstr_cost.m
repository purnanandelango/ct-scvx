function [cnstr,cost_fun,vc_cnstr] = sys_cnstr_cost(xtil,u,prb,...
                                                    ~,~)

% r = xtil(1:3)
% v = xtil(4:6)
% p = xtil(7)
% y = xtil(8)

    K = prb.K;

    % Unscaled state control input
    r1      = xtil(1:3,1);
    rK      = xtil(1:3,K);
    v1      = xtil(4:6,1);
    vK      = xtil(4:6,K);
    p1      = xtil(7,1);
    pK_scl  = prb.invSx(7,7)*(xtil(7,K)-prb.cx(7));
    y       = xtil(8,:);
    
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
        % cnstr = [cnstr;
        %          -u(prb.idx{1}(1),k) + prb.F1min <= 0;
        %           u(prb.idx{2}(1),k) + prb.F2min <= 0;
        %           u(prb.idx{3}(2),k) + prb.F3min <= 0;
        %         ];

        cnstr = [cnstr;
                 prb.umin <= u(:,k) <= prb.umax;
                ];

        if k < K
        cnstr = [cnstr;
                 y(k+1) - y(k) <= prb.eps_cnstr];
        end        
    
    end

    vc_cnstr = 0;

end