function sig = power_iteration_custom(A1, A2, B1, B2, x, u, phi, psi)
%POWER_ITERATION_CUSTOM
%   SIG = POWER_ITERATION_CUSTOM(A1, A2, B1, B2, X, U, PHI, PSI)

j_max = 10000;
eps_abs = 1e-16;
eps_rel = 1e-16;
eps_buff = 0;

nx = size(x, 1);
nu = size(u, 1);
N = size(x, 2);

x = ones(nx, N);
u = ones(nu, N);
phi = ones(nx, N-1);
psi = ones(nx, N-1);

w = ones(nx, N-1);
v = ones(1, N-1);
zeros_nx_1 = zeros(nx-1, 1);

sig = 0.0;

for k = 1:N-1
    sig = sig + x(:, k).'*x(:, k) + u(:, k).'*u(:, k) + phi(:, k).'*phi(:, k) + psi(:, k).'*psi(:, k);
end

sig = sig + x(:, N).'*x(:, N) + u(:, N).'*u(:, N);

sig = sqrt(sig);

for j = 1:j_max

    for k = 1:N-1
        w(:, k) = (1 / sig) * (A1(:, :, k)*x(:, k) + A2(:, :, k)*x(:, k+1) + B1(:, :, k)*u(:, k) + B2(:, :, k)*u(:, k+1) + phi(:, k) - psi(:, k));
        v(:, k) = (1 / sig) * (x(end, k+1) - x(end, k));
    end
    
    x(:, 1) = A1(:, :, 1).'*w(:, 1) - [zeros_nx_1; v(:, k)];
    u(:, 1) = B1(:, :, 1).'*w(:, 1);
    phi(:, 1) = w(:, 1);
    psi(:, 1) = -w(:, 1);
    
    for k = 2:N-1
        x(:, k) = A1(:, :, k).'*w(:, k) + A2(:, :, k-1).'*w(:, k-1) - [zeros_nx_1; v(:, k)] + [zeros_nx_1; v(:, k-1)];
        u(:, k) = B1(:, :, k).'*w(:, k) + B2(:, :, k-1).'*w(:, k-1);
        phi(:, k) = w(:, k);
        psi(:, k) = -w(:, k);
    end

    x(:, N) = A2(:, :, N-1).'*w(:, N-1) + [zeros_nx_1; v(:, N-1)];
    u(:, N) = B2(:, :, N-1).'*w(:, N-1);

    sig_star = 0.0;

    for k = 1:N-1
        sig_star = sig_star + x(:, k).'*x(:, k) + u(:, k).'*u(:, k) + phi(:, k).'*phi(:, k) + psi(:, k).'*psi(:, k);
    end

    sig_star = sig_star + x(:, N).'*x(:, N) + u(:, N).'*u(:, N);

    sig_star = sqrt(sig_star);

    if abs(sig_star - sig) <= eps_abs + eps_rel * max(sig_star, sig)
        break
    elseif j < j_max
        sig = sig_star;
    end

end

sig = (1 + eps_buff) * sig_star;

end