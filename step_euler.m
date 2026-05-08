function state = step_euler(state, P)

    psi  = state.psi;
    Lpsi = ifft(-P.k2 .* fft(psi));                       % psi_xx
    F    = 1i * Lpsi + 1i * P.lambda * nlse_nonlinear(psi, P);
    state.psi = nlse_filter(psi + P.tau * F, P);
end
