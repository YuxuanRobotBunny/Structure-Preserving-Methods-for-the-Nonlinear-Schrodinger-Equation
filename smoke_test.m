
P = nlse_setup('N', 128, 'L', 20, 'T', 0.5, 'tau', 1e-3, ...
               'eta', 1, 'xi', 1);
fprintf('Smoke test: N=%d  tau=%g  T=%g  Nsteps=%d\n', ...
        P.N, P.tau, P.T, P.Nsteps);

init_plain = struct('psi', P.psi0);
init_besse = struct('psi', P.psi0, 'phi', abs(P.psi0).^2);   % startup phi^{-1/2}

[M0, H0, ~] = nlse_diagnostics(P.psi0, P);
E0 = nlse_besse_energy(init_besse, P);
fprintf('Initial:  M = %.10f   H = %.10f   E_besse = %.10f\n\n', M0, H0, E0);

methods = { @step_euler,  init_plain, 'Explicit Euler';
            @step_rk4,    init_plain, 'RK4';
            @step_strang, init_plain, 'Strang SSFM';
            @step_cn,     init_plain, 'Crank-Nicolson';
            @step_besse,  init_besse, 'Besse Relaxation' };

fprintf('%-22s %14s %14s %14s %10s\n', ...
        'Method', 'max |dM|/M0', 'max |dH|/|H0|', 'max |dE|/|E0|', 'diverged');
fprintf('%s\n', repmat('-', 1, 82));
for i = 1:size(methods, 1)
    r = run_method(methods{i,1}, methods{i,2}, P, methods{i,3});
    if r.diverged
        fprintf('%-22s %14s %14s %14s %10s\n', r.name, '--', '--', '--', 'YES');
    else
        dM = max(abs(r.M - r.M(1))) / r.M(1);
        dH = max(abs(r.H - r.H(1))) / abs(r.H(1));
        if isfield(r, 'Emod')
            dE = max(abs(r.Emod - r.Emod(1))) / abs(r.Emod(1));
            fprintf('%-22s %14.3e %14.3e %14.3e %10s\n', r.name, dM, dH, dE, '');
        else
            fprintf('%-22s %14.3e %14.3e %14s %10s\n', r.name, dM, dH, '--', '');
        end
    end
end
fprintf('\nSmoke test complete.\n');
