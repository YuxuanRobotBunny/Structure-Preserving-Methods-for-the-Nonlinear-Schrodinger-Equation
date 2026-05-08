function [M, H, Pm] = nlse_diagnostics(psi, P)

    psi_x = ifft(1i * P.k .* fft(psi));

    M  = P.h * sum(abs(psi).^2);
    H  = P.h * sum(abs(psi_x).^2) - 0.5 * P.lambda * P.h * sum(abs(psi).^4);
    Pm = P.h * sum(imag(conj(psi) .* psi_x));
end
