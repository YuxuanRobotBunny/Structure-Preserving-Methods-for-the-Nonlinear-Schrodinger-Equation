function state = step_rk4(state, P)

    psi = state.psi;
    k1 = rhs(psi, P);
    k2 = rhs(psi + P.tau/2 * k1, P);
    k3 = rhs(psi + P.tau/2 * k2, P);
    k4 = rhs(psi + P.tau   * k3, P);
    state.psi = nlse_filter(psi + P.tau/6 * (k1 + 2*k2 + 2*k3 + k4), P);
end

function r = rhs(psi, P)
    Lpsi = ifft(-P.k2 .* fft(psi));
    r = 1i * Lpsi + 1i * P.lambda * nlse_nonlinear(psi, P);
end
