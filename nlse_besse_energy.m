function E = nlse_besse_energy(state, P)

    psi_x = ifft(1i * P.k .* fft(state.psi));
    phi_prev = state.phi;
    phi_next = 2 * abs(state.psi).^2 - phi_prev;

    E = P.h * sum(abs(psi_x).^2) ...
        - 0.5 * P.lambda * P.h * sum(phi_next .* phi_prev);
end
