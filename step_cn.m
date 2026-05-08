function state = step_cn(state, P)

    psi_n   = state.psi;
    psi_old = psi_n;

    max_picard = 50;
    tol_picard = 1e-12;

    for k = 1:max_picard
        rho_half = 0.5 * (abs(psi_old).^2 + abs(psi_n).^2);
        rho_half = real(nlse_filter(rho_half, P));

        % RHS:  [I + i*tau/2 * (d_xx + lambda*rho)] psi^n
        rhs = apply_op(psi_n, +P.tau/2, rho_half, P);

        % LHS:  [I - i*tau/2 * (d_xx + lambda*rho)] psi^{n+1}  =  rhs
        A = @(u) apply_op(u, -P.tau/2, rho_half, P);
        [psi_new, flag] = gmres(A, rhs, 30, 1e-13, 50);
        if flag ~= 0
            warning('CN: inner GMRES flag = %d at Picard k = %d', flag, k);
        end

        if norm(psi_new - psi_old) / max(norm(psi_n), eps) < tol_picard
            break;
        end
        psi_old = psi_new;
    end

    if k == max_picard
        warning('CN: Picard did not converge to %g in %d iterations', ...
                tol_picard, max_picard);
    end

    state.psi = nlse_filter(psi_new, P);
end

function y = apply_op(u, sgn_tau_half, rho_half, P)
%APPLY_OP  y = u + i * sgn_tau_half * (d_xx u + lambda * rho_half .* u).
    Du = ifft(-P.k2 .* fft(u));
    y  = u + 1i * sgn_tau_half * (Du + P.lambda * rho_half .* u);
end
