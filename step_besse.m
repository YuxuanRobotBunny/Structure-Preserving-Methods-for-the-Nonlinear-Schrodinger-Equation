function state = step_besse(state, P)

    psi = state.psi;
    phi = state.phi;                          % phi^{n-1/2}

    % auxiliary update (explicit, real-valued)
    phi_half = 2 * abs(psi).^2 - phi;         % phi^{n+1/2}
    phi_half = real(nlse_filter(phi_half, P));

    % build RHS:  [I + i tau/2 (d_xx + lambda phi_half)] psi
    rhs = apply_operator(psi, +P.tau/2, phi_half, P);

    % solve LHS u = rhs  with LHS = [I - i tau/2 (d_xx + lambda phi_half)]
    A = @(u) apply_operator(u, -P.tau/2, phi_half, P);
    [psi_next, flag] = gmres(A, rhs, 30, 1e-12, 50);
    if flag ~= 0
        warning('Besse: GMRES did not converge cleanly (flag = %d)', flag);
    end

    state.psi = nlse_filter(psi_next, P);
    state.phi = phi_half;                     % advance phi for next step
end

function y = apply_operator(u, sgn_tau_half, phi_half, P)
%   y = u + i*sgn_tau_half * (d_xx u + lambda * phi_half .* u)
%   sgn_tau_half = +tau/2 for RHS, -tau/2 for LHS.
    Du = ifft(-P.k2 .* fft(u));
    y  = u + 1i * sgn_tau_half * (Du + P.lambda * phi_half .* u);
end
