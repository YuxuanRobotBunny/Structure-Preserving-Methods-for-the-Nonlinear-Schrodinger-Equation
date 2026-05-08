
clear; close all;

L = 8; N = 64; T = 4*pi; tau = 0.01;
omega = 1; x0 = 1.5;          % HO frequency, initial displacement

h = 2*L/N;
x = (-L + (0:N-1)'*h);
k = (pi/L) * [(0:N/2-1)'; (-N/2:-1)'];
V = 0.5 * (omega*x).^2;

% Initial state: displaced ground state
psi0 = (omega/pi)^(1/4) * exp(-omega*(x - x0).^2 / 2);

% ---- build dense H once for Cayley LU and Magnus eigendecomposition ----
fprintf('Building dense H (N=%d)...\n', N);
D2 = zeros(N);
for j = 1:N
    e = zeros(N,1); e(j) = 1;
    D2(:,j) = ifft(-k.^2 .* fft(e));
end
D2 = (D2 + D2')/2;                       % symmetrize roundoff
H_full = -D2 + diag(V);
H_full = (H_full + H_full')/2;

fprintf('Eigendecomposing H...\n');
[U_eig, Lam] = eig(H_full, 'vector');
Lam = real(Lam);

% ---- adjust tau so Nsteps*tau == T exactly ----
Nsteps = round(T/tau); tau = T/Nsteps;

% ---- factor Cayley operator once (constant tau) ----
A_cay = eye(N) + 1i*tau/2 * H_full;
B_cay = eye(N) - 1i*tau/2 * H_full;
[Lc, Uc, Pc] = lu(A_cay);

% ---- stepper handles ----
step_rk4_lin     = @(p) rk4_step(p, H_full, tau);
step_strang_lin = @(p) strang_op_step(p, k, V, tau);
step_cayley     = @(p) Uc \ (Lc \ (Pc * (B_cay * p)));
step_magnus     = @(p) U_eig * (exp(-1i*tau*Lam) .* (U_eig' * p));

methods = { step_rk4_lin,    'RK4';
            step_strang_lin, 'Strang split-op';
            step_cayley,     'Cayley / CN';
            step_magnus,     'Magnus exact' };

% ---- run + record ----
t_grid = (0:Nsteps)' * tau;
norms  = zeros(Nsteps+1, size(methods,1));
energs = zeros(Nsteps+1, size(methods,1));

for j = 1:size(methods,1)
    fprintf('Running %s ...\n', methods{j,2});
    psi = psi0;
    norms(1,j)  = sqrt(h * sum(abs(psi).^2));
    energs(1,j) = h * real(psi' * (H_full * psi));
    for n = 1:Nsteps
        psi = methods{j,1}(psi);
        norms(n+1,j)  = sqrt(h * sum(abs(psi).^2));
        energs(n+1,j) = h * real(psi' * (H_full * psi));
    end
end

% ---- summary ----
fprintf('\n%-18s %14s %14s\n', 'Method', 'max |dN|/N0', 'max |dE|/|E0|');
fprintf('%s\n', repmat('-', 1, 50));
for j = 1:size(methods,1)
    dN = max(abs(norms(:,j) - norms(1,j))) / norms(1,j);
    dE = max(abs(energs(:,j) - energs(1,j))) / abs(energs(1,j));
    fprintf('%-18s %14.3e %14.3e\n', methods{j,2}, dN, dE);
end

% ---- plots ----
colors = {[0.85 0.33 0.10], [0.00 0.45 0.74], [0.49 0.18 0.56], [0.47 0.67 0.19]};

figure('Name','Linear Schrödinger - L^2 norm error','Position',[100 100 850 500]);
for j = 1:size(methods,1)
    relN = abs(norms(:,j) - norms(1,j)) / norms(1,j);
    relN(relN == 0) = eps;
    semilogy(t_grid, relN, 'LineWidth', 1.5, 'Color', colors{j}, ...
             'DisplayName', methods{j,2});
    hold on;
end
xlabel('t'); ylabel('| ||\psi||_2 - ||\psi_0||_2 | / ||\psi_0||_2');
title(sprintf('Linear Schrödinger (HO, \\omega=%g): unitarity (\\tau=%g)', omega, tau));
legend('Location','best','Interpreter','none'); grid on;
ylim([1e-16 1e2]);

figure('Name','Linear Schrödinger - energy error','Position',[120 120 850 500]);
for j = 1:size(methods,1)
    relE = abs(energs(:,j) - energs(1,j)) / abs(energs(1,j));
    relE(relE == 0) = eps;
    semilogy(t_grid, relE, 'LineWidth', 1.5, 'Color', colors{j}, ...
             'DisplayName', methods{j,2});
    hold on;
end
xlabel('t'); ylabel('|<H>(t) - <H>_0| / |<H>_0|');
title(sprintf('Linear Schrödinger (HO): energy expectation (\\tau=%g)', tau));
legend('Location','best','Interpreter','none'); grid on;
ylim([1e-16 1e2]);


% ============== local stepper functions ==============

function psi_next = rk4_step(psi, H, tau)
    F  = @(p) -1i * (H * p);
    k1 = F(psi);
    k2 = F(psi + tau/2 * k1);
    k3 = F(psi + tau/2 * k2);
    k4 = F(psi + tau   * k3);
    psi_next = psi + tau/6 * (k1 + 2*k2 + 2*k3 + k4);
end

function psi_next = strang_op_step(psi, k, V, tau)
    % H_T = -d_xx, Fourier symbol +k^2; H_V = V(x), pointwise
    psi = exp(-1i*tau/2 * V) .* psi;
    psi = ifft(exp(-1i*tau * k.^2) .* fft(psi));
    psi_next = exp(-1i*tau/2 * V) .* psi;
end
