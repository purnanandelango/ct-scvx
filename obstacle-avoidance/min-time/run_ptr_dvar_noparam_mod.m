function [xbar,ubar,converged] = run_ptr_dvar_noparam_mod(xbar,ubar,prb,sys_constr_cost_fun,varargin)
    % PTR SCP without parameters as decision variables and ZOH/FOH discretization
    % Provision for updating problem parameters after each SCP iteration
    % Scaling terms cx and cu are not used
    % Exact penalty weight can be matrix-valued
    
        converged = false;
        K = prb.K;
    
        % Check if type of FOH computation is specified
        if isfield(prb,'foh_type')
            foh_type = string(prb.foh_type);
            assert(ismember(foh_type,["v1","v2","v3"]),"Incorrect type of FOH discretization.");        
        else
            foh_type = "v3"; % Default
        end
    
        % Exact penalty weight
        if isfield(prb,'wvc')
            expnwt =  prb.wvc;
        elseif isfield(prb,'Wvc')
            expnwt = prb.Wvc;
        end    
        
        fprintf("+------------------------------------------------------------------------------------------------------+\n");
        fprintf("|                                   ..:: Penalized Trust Region ::..                                   |\n");
        fprintf("+-------+------------+-----------+-----------+---------+---------+----------+---------+----------------+\n");
        fprintf("| Iter. | Prop. [ms] | Prs. [ms] | Slv. [ms] | log(TR) | log(VC) |   Cost   |   ToF   | log(VC cnstr.) |\n");
        fprintf("+-------+------------+-----------+-----------+---------+---------+----------+---------+----------------+\n");
        
        for j = 1:prb.scp_iters
            
            yalmip clear
    
            % Variables
            dx = sdpvar(prb.nx,K);
            du = sdpvar(prb.nu,K);
            vc_minus = sdpvar(prb.nx,K-1);
            vc_plus = sdpvar(prb.nx,K-1);
    
            % Unscaled state and control input
            x_unscl = sdpvar(prb.nx,K);
            u_unscl = sdpvar(prb.nu,K);
            for k = 1:K
                x_unscl(:,k) = prb.Sx*dx(:,k) + xbar(:,k);
                u_unscl(:,k) = prb.Su*du(:,k) + ubar(:,k);
            end
    
            Jvc = 0;
            cnstr = [];
            for k = 1:K-1
                % Virtual control penalty 
                Jvc = Jvc + sum(expnwt*(vc_plus(:,k) + vc_minus(:,k)));
                cnstr = [cnstr;
                         vc_plus(:,k) >= 0;
                         vc_minus(:,k) >= 0];
            end
    
            % Trust region penalty
            switch prb.tr_norm
                case {2,inf}
                    Jtr = sdpvar(1,K);        
                    for k = 1:K
                        cnstr = [cnstr; norm([dx(:,k);du(:,k)],prb.tr_norm) <= Jtr(k)]; 
                    end                
                case 'quad'
                    Jtr = 0;        
                    for k = 1:K
                        Jtr = Jtr + 0.5*([dx(:,k);du(:,k)])'*([dx(:,k);du(:,k)]);
                    end                
            end
    
            % Linearized dynamics constraint
            if prb.disc == "FOH"
                % Propagation
                tic
                if isfield(prb,'ode_solver')
                    [Ak,Bmk,Bpk,~,~,xbarprop] = feval("disc.compute_foh_noparam_"+foh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
                else
                    [Ak,Bmk,Bpk,~,~,xbarprop] = feval("disc.compute_foh_noparam_"+foh_type,prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
                end
                propagate_time = toc*1000;
                
                % Prescaling, preconditioning %
    
                nx = prb.nx;
                nu = prb.nu;
    
                Ek= zeros(nx,nx,K-1);
                dx_prop = zeros(nx,K);
    
                for k = 1:K-1
                    Ek(:,:,k) = prb.Sx;
                    dx_prop(:,k+1) = xbarprop(:,k+1) - xbar(:,k+1);
                end
    
                Ek_hat = zeros(nx,nx,K-1);
                Ak_hat = zeros(nx,nx,K-1);
                Bmk_hat = zeros(nx,nu,K-1);
                Bpk_hat = zeros(nx,nu,K-1);
                dx_prop_hat = zeros(nx,K);
    
                for k = 1:K-1
                    Ek_hat(:,:,k) = prb.invSx*Ek(:,:,k);
                    Ak_hat(:,:,k) = prb.invSx*Ak(:,:,k)*prb.Sx;
                    Bmk_hat(:,:,k) = prb.invSx*Bmk(:,:,k)*prb.Su;
                    Bpk_hat(:,:,k) = prb.invSx*Bpk(:,:,k)*prb.Su;
                    dx_prop_hat(:,k+1) = prb.invSx*dx_prop(:,k+1);
                end

                % Row-normalization %
    
                for k = 1:(K-1)
                    for l = 1:nx
                        row = [Ak_hat(l,:,k), -Ek_hat(l,:,k), Bmk_hat(l,:,k), Bpk_hat(l,:,k), dx_prop_hat(l,k+1)];
                        row_norm = norm(row, 'inf');
                        % fprintf("row_norm: %f\n", row_norm);
                        if row_norm > 1e-4
                            row_norm_inv = 1 / row_norm;
                            Ek_hat(l,:,k) = row_norm_inv*Ek_hat(l,:,k);
                            Ak_hat(l,:,k) = row_norm_inv*Ak_hat(l,:,k);
                            Bmk_hat(l,:,k) = row_norm_inv*Bmk_hat(l,:,k);
                            Bpk_hat(l,:,k) = row_norm_inv*Bpk_hat(l,:,k);
                            dx_prop_hat(l,k+1) = row_norm_inv*dx_prop_hat(l,k+1);
                        end
                    end
                end

                % %%%%%%%%%%%%%%%%%%%

                % Compute maximum singular value %

                % SVD %

                Amhat = [];
                Aphat = [];
                Bmhat = [];
                Bphat = [];

                for k = 1:prb.K-1
                    Amhat = blkdiag(Amhat, Ak_hat(:,:,k));
                    Aphat = blkdiag(Aphat, -Ek_hat(:,:,k));
                    Bmhat = blkdiag(Bmhat, Bmk_hat(:,:,k));
                    Bphat = blkdiag(Bphat, Bpk_hat(:,:,k));
                end  

                Hx = [Amhat, sparse(nx*(K-1),nx)] + [sparse(nx*(K-1),nx), Aphat];
                Hu = [Bmhat, sparse(nx*(K-1),nu)] + [sparse(nx*(K-1),nu), Bphat];
                Ey = [zeros(1, nx-1), 1];
                Hy = -[kron(speye(K-1),Ey), sparse((K-1),nx)] + [sparse((K-1),nx), kron(speye(K-1),Ey)];
                Hy_ = [Hy, sparse((K-1), nu*K + 2*nx*(K-1))];
                H = [Hx, Hu, speye(nx*(K-1)), -speye(nx*(K-1)); Hy_];

                sigs_svd = svd(full(H.'*H));
                sig_svd = max(sigs_svd);

                % Customized power iteration %

                phibar = ones(nx, K-1);
                psibar = ones(nx, K-1);

                sig = power_iteration_custom(Ak_hat, -Ek_hat, Bmk_hat, Bpk_hat, prb.invSx*xbar, prb.invSu*ubar, phibar, psibar);

                sig_rel_err = norm(sig - sig_svd) / sig_svd;

                fprintf("\nMaximum singular value: %8.8f\n", sig);
                fprintf("Error in the maximum singular value: %f\n", sig_rel_err);

                [dx_pipg, du_pipg, phi_pipg, psi_pipg, w_pipg, v_pipg] = pipg_custom(sig, Ak_hat, -Ek_hat, Bmk_hat, Bpk_hat, dx_prop_hat, prb.invSx*xbar, prb.invSu*ubar, prb);

                % Warm start %

                % prb.dx = dx_pipg;
                % prb.du = du_pipg;
                % prb.phi = phi_pipg;
                % prb.psi = psi_pipg;
                % prb.w = w_pipg;
                % prb.v = v_pipg;

                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                for k = 1:K-1
                    cnstr = [cnstr;
                             vc_minus(:,k) - vc_plus(:,k) == - Ek_hat(:,:,k)*dx(:,k+1) +...
                                                               Ak_hat(:,:,k)*dx(:,k) +...
                                                               Bmk_hat(:,:,k)*du(:,k) +...
                                                               Bpk_hat(:,:,k)*du(:,k+1) +...
                                                               dx_prop_hat(:,k+1)];
                end                        
            elseif prb.disc == "ZOH"
                % Propagation
                tic
                if isfield(prb,'ode_solver')
                    [Ak,Bk,~,~,xbarprop] = disc.compute_zoh_noparam(prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize,prb.ode_solver);
                else
                    [Ak,Bk,~,~,xbarprop] = disc.compute_zoh_noparam(prb.tau,xbar,ubar,prb.h,prb.dyn_func,prb.dyn_func_linearize);
                end
                propagate_time = toc*1000;
    
                for k = 1:K-1
                    cnstr = [cnstr;
                             vc_minus(:,k) - vc_plus(:,k) == - dx(:,k+1) +...
                                                               prb.invSx*Ak(:,:,k)*prb.Sx*dx(:,k) +...
                                                               prb.invSx*Bk(:,:,k)*prb.Su*du(:,k) +...
                                                               prb.invSx*(xbarprop(:,k+1) - xbar(:,k+1))];
                end
                cnstr = [cnstr; u(:,K) == u(:,K-1)];             
            end
            
            % Constraints
            [cnstr_sys,cost_fun,vc_constr_term] = sys_constr_cost_fun(x_unscl,u_unscl,prb,...
                                                                      xbar,ubar);
    
    
            cnstr = [cnstr;cnstr_sys];
            
            % Objective
            obj_fun = Jvc + prb.wtr*sum(Jtr) + cost_fun;            
            
            % Solve
            % Model = export(cnstr,obj_fun,prb.solver_settings); % Export input to solver from YALMIP
            yalmip_out = optimize(cnstr,obj_fun,prb.solver_settings);
            % assert(ismember(yalmip_out.problem,[0,3]),"Subproblem is unsolved.\nSolver message: %s",yalmiperror(yalmip_out.problem));
            if ~ismember(yalmip_out.problem,[0,4])
                fprintf("+------------------------------------------------------------------------------------------------------+\n");
                fprintf('Subproblem is unsolved. Returning the previous iterate.\n'); 
                break
            end
            if yalmip_out.problem == 4
                warning("Solver numerical issues.");
            end
            
            % Post process
            solve_time = yalmip_out.solvertime*1000;
            parse_time = yalmip_out.yalmiptime*1000;            
            dx = value(dx);
            du = value(du);
            x_unscl = value(x_unscl);
            u_unscl = value(u_unscl);        
            cost_val = value(cost_fun)/prb.cost_factor;
            vc_term = value(Jvc);
            vc_constr_term = value(vc_constr_term)/max(expnwt(:));

            fprintf("Error in `dx`: %f\n", norm(dx_pipg - dx) / norm(dx));
            fprintf("Error in `du`: %f\n\n", norm(du_pipg - du) / norm(dx));
    
            % Ensure that the TR value is always displayed consistently with infinity norm
            % Note that this is for display and for evaluation of termination criteria 
            Jtr_post_solve = zeros(1,K);        
            for k = 1:K
                Jtr_post_solve(k) = norm([dx(:,k);du(:,k)],'inf');
            end
            tr_term = max(Jtr_post_solve);
    
            % Update reference trajectory
            xbar = x_unscl;
            ubar = u_unscl;     

            % xbar = prb.Sx * dx_pipg + xbar;
            % ubar = prb.Su * du_pipg + ubar;
    
            ToF = prb.time_of_maneuver(xbar,ubar);        
            
            % Console output
            fprintf('|  %02d   |   %7.1e  |  %7.1e  |  %7.1e  | %5.1f   | %5.1f   | %8.1e | %8.5f |    %5.1f       |\n',j,propagate_time,parse_time,solve_time,log10(tr_term),log10(vc_term),cost_val,ToF,log10(vc_constr_term));
            
            if vc_term < prb.epsvc && tr_term < prb.epstr
                converged = true;
                fprintf("+------------------------------------------------------------------------------------------------------+\n")
                fprintf('Converged!\n')
                break
            end
    
            if nargin == 5 && j < prb.scp_iters % Update problem parameters
                prb = varargin{1}(prb,xbar,ubar);
            end
            
        end
    
end