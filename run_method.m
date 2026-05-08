function out = run_method(stepper, init_state, P, name)

    state = init_state;
    out.name     = name;
    out.t        = (0:P.Nsteps)' * P.tau;
    out.M        = zeros(P.Nsteps+1, 1);
    out.H        = zeros(P.Nsteps+1, 1);
    out.Pm       = zeros(P.Nsteps+1, 1);
    out.diverged = false;
    track_besse_energy = isfield(state, 'phi');
    if track_besse_energy
        out.Emod = zeros(P.Nsteps+1, 1);
    end

    [out.M(1), out.H(1), out.Pm(1)] = nlse_diagnostics(state.psi, P);
    if track_besse_energy
        out.Emod(1) = nlse_besse_energy(state, P);
    end

    for n = 1:P.Nsteps
        state = stepper(state, P);

        % divergence guard (forward Euler will trip this)
        if any(~isfinite(state.psi)) || max(abs(state.psi)) > 1e8
            warning('%s diverged at step %d / %d (t = %g)', ...
                    name, n, P.Nsteps, n*P.tau);
            out.diverged = true;
            out.M(n+1:end)  = NaN;
            out.H(n+1:end)  = NaN;
            out.Pm(n+1:end) = NaN;
            if track_besse_energy
                out.Emod(n+1:end) = NaN;
            end
            break;
        end

        [out.M(n+1), out.H(n+1), out.Pm(n+1)] = nlse_diagnostics(state.psi, P);
        if track_besse_energy
            out.Emod(n+1) = nlse_besse_energy(state, P);
        end
    end

    out.psi_T = state.psi;
end
