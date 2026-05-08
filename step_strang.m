function state = step_strang(state, P)

    psi = state.psi;

    % half step of linear part (Fourier diagonal)
    psi = ifft(exp(-1i * P.k2 * P.tau/2) .* fft(psi));

    % full step of nonlinear part (pointwise phase rotation)
    psi = exp(1i * P.lambda * abs(psi).^2 * P.tau) .* psi;
    psi = nlse_filter(psi, P);

    % half step of linear part again
    psi = ifft(exp(-1i * P.k2 * P.tau/2) .* fft(psi));

    state.psi = nlse_filter(psi, P);
end
