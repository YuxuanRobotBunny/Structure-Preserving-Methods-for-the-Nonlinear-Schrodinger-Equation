
clear; close all;

P = nlse_setup('lambda',1, 'L',20, 'N',128, 'T',100, 'tau',0.01, ...
               'eta',1, 'xi',1);
fprintf('Long-time drift test: T=%g, tau=%g, Nsteps=%d, N=%d\n', ...
        P.T, P.tau, P.Nsteps, P.N);
fprintf('(CN with Picard inner iteration is the slow one; expect ~1-3 min total.)\n\n');

init_plain = struct('psi', P.psi0);
init_besse = struct('psi', P.psi0, 'phi', abs(P.psi0).^2);   % phi^{-1/2} startup

t_start = tic;
fprintf('[1/4] RK4...\n');
res_r4 = run_method(@step_rk4,    init_plain, P, 'RK4');
fprintf('  elapsed: %.1f s\n', toc(t_start));

t_start = tic;
fprintf('[2/4] Strang SSFM...\n');
res_st = run_method(@step_strang, init_plain, P, 'Strang SSFM');
fprintf('  elapsed: %.1f s\n', toc(t_start));

t_start = tic;
fprintf('[3/4] Crank-Nicolson (DFP)...\n');
res_cn = run_method(@step_cn,     init_plain, P, 'Crank-Nicolson');
fprintf('  elapsed: %.1f s\n', toc(t_start));

t_start = tic;
fprintf('[4/4] Besse Relaxation...\n');
res_be = run_method(@step_besse,  init_besse, P, 'Besse Relaxation');
fprintf('  elapsed: %.1f s\n', toc(t_start));

results = {res_r4, res_st, res_cn, res_be};
colors  = {[0.93 0.69 0.13], [0.00 0.45 0.74], [0.49 0.18 0.56], [0.47 0.67 0.19]};

% summary
fprintf('\n%-22s %14s %14s\n', 'Method', 'max |dM|/M0', 'max |dH|/|H0|');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:numel(results)
    r = results{i};
    if r.diverged
        fprintf('%-22s %14s %14s\n', r.name, '--', '--');
    else
        dM = max(abs(r.M - r.M(1))) / r.M(1);
        dH = max(abs(r.H - r.H(1))) / abs(r.H(1));
        fprintf('%-22s %14.3e %14.3e\n', r.name, dM, dH);
    end
end

% ---- plot 1: long-time mass drift ----
figure('Name','Long-time mass drift','Position',[100 100 900 500]);
for i = 1:numel(results)
    r = results{i};
    if r.diverged, continue; end
    relM = abs(r.M - r.M(1)) / r.M(1);
    relM(relM == 0) = eps;
    semilogy(r.t, relM, 'LineWidth', 1.4, 'Color', colors{i}, ...
             'DisplayName', r.name);
    hold on;
end
xlabel('t'); ylabel('|M(t) - M_0| / M_0');
title(sprintf('Long-time relative mass error  (T = %g, \\tau = %g)', P.T, P.tau));
legend('Location','best','Interpreter','none'); grid on;
ylim([1e-16 1e2]);

% ---- plot 2: long-time energy drift ----
figure('Name','Long-time energy drift','Position',[120 120 900 500]);
for i = 1:numel(results)
    r = results{i};
    if r.diverged, continue; end
    relH = abs(r.H - r.H(1)) / abs(r.H(1));
    relH(relH == 0) = eps;
    semilogy(r.t, relH, 'LineWidth', 1.4, 'Color', colors{i}, ...
             'DisplayName', r.name);
    hold on;
end
xlabel('t'); ylabel('|H(t) - H_0| / |H_0|');
title(sprintf('Long-time relative energy error  (T = %g, \\tau = %g)', P.T, P.tau));
legend('Location','best','Interpreter','none'); grid on;
ylim([1e-16 1e2]);
