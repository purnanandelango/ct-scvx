function [cnstr,cost_fun,vcvb_cnstr] = sys_cnstr_cost(x,u,prb,...
                                                      xbar,ubar)

    K = prb.K;

    m     = x(1,:);
    rI    = x(2:4,:);
    vI    = x(5:7,:);
    qBI   = x(8:11,:); 
    omgB  = x(12:14,:);    
    TB    = u(1:3,:);
    s     = u(4,:);    
    
    % Boundary conditions
    cnstr = [m(1)      == prb.mwet;
             rI(:,1)   == prb.rI1;
             vI(:,1)   == prb.vI1;
             rI(:,K)   == prb.rIK;
             vI(:,K)   == prb.vIK;
             qBI(:,K)  == prb.q1;
             omgB(:,1) == prb.omgB1;
             omgB(:,K) == prb.omgBK];

    vb_TBmin = sdpvar(1,K);       
    % vb_stc  = sdpvar(1,K);

    cost_fun = 0;

    for k = 1:K    
        
        cnstr = [cnstr;
                 m(k) >= prb.mdry                                                                           % Vehicle mass lower bound
                 norm(prb.Hgam*rI(:,k)) <= rI(1,k)/prb.cotgamgs;                                            % Glide slope constraint
                 norm(prb.Hthet*qBI(:,k)) <= prb.sinthetmaxby2;                                            % Vehicle tilt angle constraint
                 norm(omgB(:,k)) <= prb.omgmax;                                                             % Angular velocity magnitude upper bound
                 norm(TB(:,k)) <= prb.TBmax;                                                                % Thrust magnitude upper bound
                 prb.cosdelmax*norm(TB(:,k)) <= TB(1,k);                                                    % Thrust pointing constraint
                 norm(vI(:,k)) <= prb.vmax;                                                                 % Velocity magnitude upper bound
                 prb.smin <= s(k) <= prb.smax;                                                              % Lower and upper bounds on dilation factor
                 (TB(:,k))'*(-ubar(1:3,k)/norm(ubar(1:3,k))) + prb.TBmin <= vb_TBmin(k);                    % Linearized thrust magnitude lower bound
                 vb_TBmin(k) >= 0];

        % Airspeed-triggered angle-of-attack STC
        % [stc,dstc_dx] = plant.rocket6DoF.q_aoa_cnstr(xbar(5:7,k),xbar(8:11,k),prb.vmax_stc,prb.cosaoamax,prb.stc_flag);
        % cnstr = [cnstr;
        %          stc + dstc_dx*(x(:,k)-xbar(:,k)) <= vb_stc(k);
        %          vb_stc(k) >= 0];
    
    end  

    cost_fun = cost_fun + prb.cost_factor*prb.invSx(1,1)*x(1,K);

    cnstr = [cnstr;
             misc.time_cnstr(s,prb.dtau,{prb.dtmin,prb.dtmax,prb.ToFmax},prb.disc)];

    vcvb_cnstr = 0;

    vcvb_cnstr = vcvb_cnstr + prb.wvc*sum(vb_TBmin(:));

    % vcvb_cnstr = vcvb_cnstr + prb.wvc*sum(vb_stc(:));

    cost_fun = cost_fun + vcvb_cnstr;

end