clear; close all; clc;

lambda = 1; L = 20; N = 192; T = 1.0; tau = 1e-3;
eta = 1; xi = 1;

P0 = nlse_setup('lambda',lambda, 'L',L, 'N',N, 'T',T, 'tau',tau, ...
                'eta',eta, 'xi',xi, 'dealias',false);
P1 = nlse_setup('lambda',lambda, 'L',L, 'N',N, 'T',T, 'tau',tau, ...
                'eta',eta, 'xi',xi, 'dealias',true);

mode_hi = floor(0.44 * N);
packet = sech(0.35 * P0.x) .* exp(1i * (pi/L) * mode_hi * P0.x);
psi0 = P0.psi0 + 0.12 * packet;
P0.psi0 = psi0;
P1.psi0 = nlse_filter(psi0, P1);

fprintf('De-aliasing study: N=%d, tau=%g, T=%g, high mode=%d\n', ...
        N, P0.tau, P0.T, mode_hi);
fprintf('2/3 filter retains %d / %d Fourier modes.\n\n', ...
        nnz(P1.dealias_mask), N);

methods = { @step_rk4,    'plain', 'RK4';
            @step_strang, 'plain', 'Strang SSFM';
            @step_besse,  'besse', 'Besse Relaxation' };

fprintf('%-18s %-10s %14s %14s %14s\n', ...
        'Method', 'filter', 'max |dM|/M0', 'max |dH|/|H0|', 'tail energy');
fprintf('%s\n', repmat('-', 1, 78));

results = cell(size(methods,1), 2);
for j = 1:size(methods,1)
    for f = 1:2
        if f == 1
            P = P0;
            label = 'off';
        else
            P = P1;
            label = 'on';
        end
        if strcmp(methods{j,2}, 'besse')
            init = struct('psi', P.psi0, 'phi', abs(P.psi0).^2);
        else
            init = struct('psi', P.psi0);
        end
        r = run_method(methods{j,1}, init, P, methods{j,3});
        results{j,f} = r;
        if r.diverged
            fprintf('%-18s %-10s %14s %14s %14s\n', methods{j,3}, label, '--', '--', '--');
        else
            dM = max(abs(r.M - r.M(1))) / r.M(1);
            dH = max(abs(r.H - r.H(1))) / abs(r.H(1));
            tail = spectral_tail_energy(r.psi_T, P1.dealias_mask);
            fprintf('%-18s %-10s %14.3e %14.3e %14.3e\n', ...
                    methods{j,3}, label, dM, dH, tail);
        end
    end
end

figure('Name','De-aliasing Fourier spectrum','Position',[100 100 900 520]);
for j = 1:size(methods,1)
    for f = 1:2
        r = results{j,f};
        if r.diverged, continue; end
        spec = abs(fftshift(fft(r.psi_T))).^2 / N^2;
        modes = (-N/2:N/2-1)';
        if f == 1
            ls = '-';
            suffix = 'off';
        else
            ls = '--';
            suffix = 'on';
        end
        semilogy(modes, max(spec, eps), ls, 'LineWidth', 1.3, ...
                 'DisplayName', sprintf('%s (%s)', methods{j,3}, suffix));
        hold on;
    end
end
xline(-floor(N/3), 'k:', 'HandleVisibility','off');
xline( floor(N/3), 'k:', 'HandleVisibility','off');
xlabel('Fourier mode'); ylabel('|hat{\psi}_k(T)|^2');
title('Final Fourier spectrum with and without 2/3 de-aliasing');
legend('Location','best','Interpreter','none'); grid on;

figure('Name','De-aliasing energy error','Position',[140 140 900 520]);
colors = {[0.93 0.69 0.13], [0.00 0.45 0.74], [0.47 0.67 0.19]};
for j = 1:size(methods,1)
    for f = 1:2
        r = results{j,f};
        if r.diverged, continue; end
        relH = abs(r.H - r.H(1)) / abs(r.H(1));
        relH(relH == 0) = eps;
        if f == 1
            ls = '-';
            suffix = 'off';
        else
            ls = '--';
            suffix = 'on';
        end
        semilogy(r.t, relH, ls, 'Color', colors{j}, 'LineWidth', 1.3, ...
                 'DisplayName', sprintf('%s (%s)', methods{j,3}, suffix));
        hold on;
    end
end
xlabel('t'); ylabel('|H(t) - H(0)| / |H(0)|');
title('Energy error with and without 2/3 de-aliasing');
legend('Location','best','Interpreter','none'); grid on;

function tail = spectral_tail_energy(psi, mask)
    psi_hat = fft(psi);
    tail = sum(abs(psi_hat(~mask)).^2) / max(sum(abs(psi_hat).^2), eps);
end
