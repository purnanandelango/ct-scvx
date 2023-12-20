function [dx, du, phi, psi, w, v, prb] = pipg_custom(sig, A1, A2, B1, B2, d, x, u, prb, verbose)

j_max = 1000000;
j_min = 100000;
j_check = 100000;
eps_abs = 1e-12;
eps_rel = 1e-12;
omg = 10.0;
rho = 1.85;

invPy = prb.invSx(end, end);

epsilon = invPy*prb.eps_cnstr;
lam = prb.wtr;

nx = size(x, 1);
nu = size(u, 1);
N = size(x, 2);

qx = zeros(nx, N);
qu = zeros(nu, N);
qphi = prb.wvc*ones(nx, N-1);
qpsi = prb.wvc*ones(nx, N-1);

dx = prb.dx;
du = prb.du;
phi = prb.phi;
psi = prb.psi;
w = prb.w;
v = prb.v;

dx_ = prb.dx;
du_ = prb.du;
phi_ = prb.phi;
psi_ = prb.psi;
w_ = prb.w;
v_ = prb.v;

zdx = prb.dx;
zdu = prb.du;
zphi = prb.phi;
zpsi = prb.psi;
eta = prb.w;
gam = prb.v;

zeros_nx_1 = zeros(nx-1, 1);

alf = 2 / (lam + sqrt(lam^2 + 4*omg*sig));
bet = omg*alf;

for j = 1:j_max

    if j >= j_min && mod(j, j_check) == 0
        dx_ = dx;
        du_ = du;
        phi_ = phi;
        psi_ = psi;
        w_ = w;
        v_ = v;
    end

    dx(:, 1) = zdx(:, 1) - alf*(lam*zdx(:, 1) + qx(:, 1) + A1(:, :, 1).'*eta(:, 1) - [zeros_nx_1; gam(:, 1)]);
    dx(:, 1) = proj_Dx1(dx(:, 1) + x(:, 1), prb) - x(:, 1); % shifted

    du(:, 1) = zdu(:, 1) - alf*(lam*zdu(:, 1) + qu(:, 1) + B1(:, :, 1).'*eta(:, 1));
    du(:, 1) = proj_Du(du(:, 1) + u(:, 1), prb) - u(:, 1); % shifted

    phi(:, 1) = zphi(:, 1) - alf*(qphi(:, 1) + eta(:, 1));
    phi(:, 1) = proj_Dphi(phi(:, 1));

    psi(:, 1) = zpsi(:, 1) - alf*(qpsi(:, 1) - eta(:, 1));
    psi(:, 1) = proj_Dpsi(psi(:, 1));

    for k = 2:N-1

        dx(:, k) = zdx(:, k) - alf*(lam*zdx(:, k) + qx(:, k) + A1(:, :, k).'*eta(:, k) + A2(:, :, k-1).'*eta(:, k-1) - [zeros_nx_1; gam(:, k)] + [zeros_nx_1; gam(:, k-1)]);
        % dx(:, k) = proj_Dxk(dx(:, k) + x(:, k), prb) - x(:, k); % shifted
    
        du(:, k) = zdu(:, k) - alf*(lam*zdu(:, k) + qu(:, k) + B1(:, :, k).'*eta(:, k) + B2(:, :, k-1).'*eta(:, k-1));
        du(:, k) = proj_Du(du(:, k) + u(:, k), prb) - u(:, k); % shifted
    
        phi(:, k) = zphi(:, k) - alf*(qphi(:, k) + eta(:, k));
        phi(:, k) = proj_Dphi(phi(:, k));
    
        psi(:, k) = zpsi(:, k) - alf*(qpsi(:, 1) - eta(:, k));
        psi(:, k) = proj_Dpsi(psi(:, k));

    end

    dx(:, N) = zdx(:, N) - alf*(lam*zdx(:, N) + qx(:, N) + A2(:, :, N-1).'*eta(:, N-1) + [zeros_nx_1; gam(:, N-1)]);
    dx(:, N) = proj_DxN(dx(:, N) + x(:, N), prb) - x(:, N); % shifted

    du(:, N) = zdu(:, N) - alf*(lam*zdu(:, N) + qu(:, N) + B2(:, :, N-1).'*eta(:, N-1));
    du(:, N) = proj_Du(du(:, N) + u(:, N), prb) - u(:, N); % shifted

    for k = 1:N-1

        w(:, k) = eta(:, k) + ... 
                  bet*( ...
                          A1(:, :, k)*(2*dx(:, k) - zdx(:, k)) + A2(:, :, k)*(2*dx(:, k+1) - zdx(:, k+1)) + ...
                          B1(:, :, k)*(2*du(:, k) - zdu(:, k)) + B2(:, :, k)*(2*du(:, k+1) - zdu(:, k+1)) + ...
                          (2*phi(:, k) - zphi(:, k)) - (2*psi(:, k) - zpsi(:, k)) + d(:, k+1) ...
                      );

        v(:, k) = max(0, gam(:, k) + ...
                  bet*( ...
                          (2*dx(end, k+1) - zdx(end, k+1)) - (2*dx(end, k) - zdx(end, k)) + (x(end, k+1) - x(end, k)) - epsilon ...
                      ));

    end

    zdx = (1 - rho)*zdx + rho*dx;
    zdu = (1 - rho)*zdu + rho*du;
    zphi = (1 - rho)*zphi + rho*phi;
    zpsi = (1 - rho)*zpsi + rho*psi;
    eta = (1 - rho)*eta + rho*w;
    gam = (1 - rho)*gam + rho*v;

    if (j >= j_min && mod(j, j_check) == 0)

        z_inf_jp1 = max([norm(dx(:), 'inf'), norm(du(:), 'inf'), norm(phi(:), 'inf'), norm(psi(:), 'inf')]);
        z_inf_j = max([norm(dx_(:), 'inf'), norm(du_(:), 'inf'), norm(phi_(:), 'inf'), norm(psi_(:), 'inf')]);
        z_inf_del_j = max([norm(dx(:) - dx_(:), 'inf'), norm(du(:) - du_(:), 'inf'), norm(phi(:) - phi_(:), 'inf'), norm(psi(:) - psi_(:), 'inf')]);

        r_inf_jp1 = max([norm(w(:), 'inf'), norm(v(:), 'inf')]);
        r_inf_j = max([norm(w_(:), 'inf'), norm(v_(:), 'inf')]);
        r_inf_del_j = max([norm(w(:) - w_(:), 'inf'), norm(v(:) - v_(:), 'inf')]);

        if z_inf_del_j <= eps_abs + eps_rel*max(z_inf_jp1, z_inf_j) && r_inf_del_j <= eps_abs + eps_rel*max(r_inf_jp1, r_inf_j)

            if verbose == 1

                if j == j_min
    
                    fprintf("\n");
                    fprintf("+-----------------------------------------------------------------------+\n");
                    fprintf("|                             ..:: PIPG ::..                            |\n");
                    fprintf("+-----------+-----------+------------------+--------------+-------------+\n");
                    fprintf("| Iteration | Objective | Constraint Viol. | Primal Diff. |  Dual Diff. |\n");
                    fprintf("+-----------+-----------+------------------+--------------+-------------+\n");
        
                end
    
                fprintf("| %9.2e | %9.2e |     %9.2e    |   %9.2e  |  %9.2e  |\n", ...
                j,NaN,NaN,z_inf_del_j,r_inf_del_j);

            end

            break

        end

        if verbose == 1

            if j == j_min
    
                fprintf("\n");
                fprintf("+-----------------------------------------------------------------------+\n");
                fprintf("|                             ..:: PIPG ::..                            |\n");
                fprintf("+-----------+-----------+------------------+--------------+-------------+\n");
                fprintf("| Iteration | Objective | Constraint Viol. | Primal Diff. |  Dual Diff. |\n");
                fprintf("+-----------+-----------+------------------+--------------+-------------+\n");
    
            end
    
            fprintf("| %9.2e | %9.2e |     %9.2e    |   %9.2e  |  %9.2e  |\n", ...
                j,NaN,NaN,z_inf_del_j,r_inf_del_j);

        end

    end

end

if verbose == 1

    fprintf("+-----------+-----------+------------------+--------------+-------------+\n\n");

end

end

%%%%%%%%%%%%%%%%%%%
%%% Projections %%%
%%%%%%%%%%%%%%%%%%%

function xp = proj_Dx1(x, prb)

    xp = x;

    Ei = prb.Ei;
    zi = prb.zi;
    Px_inv = prb.invSx;

    x1 = (Ei * Px_inv * Ei.') * zi;

    i_idx = prb.i_idx;

    xp(i_idx) = x1;

end

function xp = proj_Dxk(x, ~)
    
    xp = x;

end

function xp = proj_DxN(x, prb)
   
    xp = x;

    Ef = prb.Ef;
    zf = prb.zf;
    Px_inv = prb.invSx;

    xN = (Ef * Px_inv * Ef.') * zf;
    
    f_idx = prb.f_idx;

    xp(f_idx) = xN;

end

function up = proj_Du(u, prb)
    
    Pu_inv = prb.invSu;
    u_min = Pu_inv * prb.umin;
    u_max = Pu_inv * prb.umax;

    up = u;
    up(end) = min(max(u(end), u_min(end)), u_max(end));

end

function phip = proj_Dphi(phi)

    phip = max(0, phi);

end

function psip = proj_Dpsi(psi)

    psip = max(0, psi);

end