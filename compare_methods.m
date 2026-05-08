clear; close all; clc;

% ---- problem ----
P = nlse_setup('lambda', 1, 'L', 20, 'N', 256, 'T', 2, 'tau', 1e-3, ...
               'eta', 1, 'xi', 1);

fprintf('NLSE structure-preservation comparison\n');
fprintf('  N = %d, L = %g, tau = %g, T = %g, Nsteps = %d, lambda = %g\n', ...
        P.N, P.L, P.tau, P.T, P.Nsteps, P.lambda);

% ---- initial states (all share psi0; Besse stores phi^{n-1/2}) ----
init_plain = struct('psi', P.psi0);
init_besse = struct('psi', P.psi0, 'phi', abs(P.psi0).^2);   % startup: phi^{-1/2} = |psi^0|^2


fprintf('\n[1/5] Explicit Euler...\n');
res_eu = run_method(@step_euler,  init_plain, P, 'Explicit Euler');
fprintf('[2/5] RK4...\n');
res_r4 = run_method(@step_rk4,    init_plain, P, 'RK4');
fprintf('[3/5] Strang SSFM...\n');
res_st = run_method(@step_strang, init_plain, P, 'Strang SSFM');
fprintf('[4/5] Crank-Nicolson (DFP)...\n');
res_cn = run_method(@step_cn,     init_plain, P, 'Crank-Nicolson');
fprintf('[5/5] Besse Relaxation...\n');
res_be = run_method(@step_besse,  init_besse, P, 'Besse Relaxation');

results = {res_eu, res_r4, res_st, res_cn, res_be};
colors  = {[0.85 0.33 0.10], [0.93 0.69 0.13], [0.00 0.45 0.74], ...
           [0.49 0.18 0.56], [0.47 0.67 0.19]};

fprintf('\n%-22s %14s %14s %12s\n', ...
        'Method', 'max |dM|/M0', 'max |dH|/|H0|', 'diverged');
fprintf('%s\n', repmat('-', 1, 65));
for i = 1:numel(results)
    r = results{i};
    if r.diverged
        fprintf('%-22s %14s %14s %12s\n', r.name, '   --', '   --', 'YES');
    else
        dM = max(abs(r.M - r.M(1))) / r.M(1);
        dH = max(abs(r.H - r.H(1))) / abs(r.H(1));
        fprintf('%-22s %14.3e %14.3e %12s\n', r.name, dM, dH, '');
    end
end

% ---- plot 1: relative mass error ----
figure('Name', 'Mass Error', 'Position', [80 80 800 500]);
for i = 1:numel(results)
    r = results{i};
    relM = abs(r.M - r.M(1)) / r.M(1);
    relM(relM == 0) = eps;            % avoid log(0) gaps
    semilogy(r.t, relM, 'LineWidth', 1.5, 'Color', colors{i}, ...
             'DisplayName', r.name);
    hold on;
end
xlabel('t'); ylabel('|M(t) - M(0)| / M(0)');
title(sprintf('Relative mass error  (\\tau = %g, N = %d)', P.tau, P.N));
legend('Location', 'best', 'Interpreter', 'none'); grid on;
ylim([1e-16 1e2]);

% ---- plot 2: relative energy error ----
figure('Name', 'Energy Error', 'Position', [120 120 800 500]);
for i = 1:numel(results)
    r = results{i};
    relH = abs(r.H - r.H(1)) / abs(r.H(1));
    relH(relH == 0) = eps;
    semilogy(r.t, relH, 'LineWidth', 1.5, 'Color', colors{i}, ...
             'DisplayName', r.name);
    hold on;
end
xlabel('t'); ylabel('|H(t) - H(0)| / |H(0)|');
title(sprintf('Relative energy error  (\\tau = %g, N = %d)', P.tau, P.N));
legend('Location', 'best', 'Interpreter', 'none'); grid on;
ylim([1e-16 1e2]);

% ---- plot 3: final density profiles ----

styles = {'-', '-', '--', '-.', ':'};  
widths = [1.6, 1.6, 1.8, 2.0, 2.2];     

figure('Name', 'Final Density |psi(x,T)|^2', 'Position', [160 160 900 500]);
plot(P.x, abs(P.psi0).^2, 'k:', 'LineWidth', 1.2, 'DisplayName', 'initial');
hold on;
for i = 1:numel(results)
    r = results{i};
    if ~r.diverged
        plot(P.x, abs(r.psi_T).^2, 'LineStyle', styles{i}, ...
             'LineWidth', widths(i), 'Color', colors{i}, ...
             'DisplayName', r.name);
    end
end
xlabel('x'); ylabel('|\psi|^2');
title(sprintf('Density at t = T = %g  (soliton velocity = 2\\xi = %g)', ...
              P.T, 2*P.xi));
legend('Location', 'best', 'Interpreter', 'none'); grid on;
xlim([-P.L P.L]);

% ---- plot 4: density error vs reference (CN) ----

ref_idx = 4;   % CN
if ~results{ref_idx}.diverged
    rho_ref = abs(results{ref_idx}.psi_T).^2;
    figure('Name', 'Density Error vs CN', 'Position', [200 200 900 500]);
    for i = 1:numel(results)
        if i == ref_idx, continue; end
        r = results{i};
        if r.diverged, continue; end
        err = abs(abs(r.psi_T).^2 - rho_ref);
        err(err < eps) = eps;
        semilogy(P.x, err, 'LineStyle', styles{i}, 'LineWidth', 1.5, ...
                 'Color', colors{i}, 'DisplayName', r.name);
        hold on;
    end
    xlabel('x'); ylabel('| |\psi_{method}|^2 - |\psi_{CN}|^2 |');
    title(sprintf('Pointwise density error vs CN at t = T = %g', P.T));
    legend('Location', 'best', 'Interpreter', 'none'); grid on;
    xlim([-P.L P.L]);
end
