
clear; close all;

eta = 1; xi = 1; lambda = 1;
L = 20; N = 256; T = 1.0;

tau_list = [2e-3, 1e-3, 5e-4, 2.5e-4];

tau_ref = tau_list(end) / 8;
fprintf('Computing reference: RK4 at tau_ref = %.3e (Nsteps=%d)\n', ...
        tau_ref, round(T/tau_ref));
P_ref = nlse_setup('lambda',lambda,'L',L,'N',N,'T',T,'tau',tau_ref, ...
                   'eta',eta,'xi',xi);
init_ref = struct('psi', P_ref.psi0);
r_ref = run_method(@step_rk4, init_ref, P_ref, 'reference');
psi_ref = r_ref.psi_T;

methods = { @step_rk4,    'plain', 'RK4';
            @step_strang, 'plain', 'Strang SSFM';
            @step_cn,     'plain', 'Crank-Nicolson';
            @step_besse,  'besse', 'Besse Relaxation' };

errors = NaN(numel(tau_list), size(methods,1));

for j = 1:size(methods,1)
    fprintf('\nMethod: %s\n', methods{j,3});
    for k = 1:numel(tau_list)
        P = nlse_setup('lambda',lambda,'L',L,'N',N,'T',T,'tau',tau_list(k), ...
                       'eta',eta,'xi',xi);
        if strcmp(methods{j,2}, 'besse')
            init = struct('psi', P.psi0, 'phi', abs(P.psi0).^2);   % phi^{-1/2} startup
        else
            init = struct('psi', P.psi0);
        end
        r = run_method(methods{j,1}, init, P, methods{j,3});
        if ~r.diverged
            errors(k,j) = sqrt(P.h * sum(abs(r.psi_T - psi_ref).^2));
        end
        fprintf('  tau = %.3e   ||err||_2 = %.3e\n', tau_list(k), errors(k,j));
    end
end

% empirical orders (regression in log-log)
fprintf('\nEmpirical orders (least-squares slope of log(err) vs log(tau)):\n');
log_tau = log(tau_list(:));
for j = 1:size(methods,1)
    e = errors(:,j);
    mask = isfinite(e) & e > 0;
    if sum(mask) >= 2
        p = polyfit(log_tau(mask), log(e(mask)), 1);
        fprintf('  %-22s  order ~ %.2f\n', methods{j,3}, p(1));
    end
end

% log-log plot
markers = {'o','s','^','d'};
colors  = {[0.93 0.69 0.13], [0.00 0.45 0.74], [0.49 0.18 0.56], [0.47 0.67 0.19]};

figure('Name','Convergence Order','Position',[100 100 850 600]);
for j = 1:size(methods,1)
    loglog(tau_list, errors(:,j), [markers{j} '-'], 'LineWidth', 1.6, ...
           'Color', colors{j}, 'MarkerSize', 8, ...
           'MarkerFaceColor', colors{j}, 'DisplayName', methods{j,3});
    hold on;
end

% reference slopes anchored at finest tau of order-2 and order-4 methods
tau_x = [tau_list(end), tau_list(1)];
ref2 = errors(end,2) * (tau_x/tau_list(end)).^2;
ref4 = errors(end,1) * (tau_x/tau_list(end)).^4;
loglog(tau_x, ref2, 'k--', 'LineWidth', 1, 'DisplayName','slope 2 (ref)');
loglog(tau_x, ref4, 'k:',  'LineWidth', 1, 'DisplayName','slope 4 (ref)');

xlabel('\tau'); ylabel('||\psi^N - \psi_{ref}||_2');
title(sprintf('Time-stepping convergence (T = %g, N = %d, L = %g)', T, N, L));
legend('Location','northwest','Interpreter','none');
grid on;
