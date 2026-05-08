function nonlin = nlse_nonlinear(psi, P)

    nonlin = abs(psi).^2 .* psi;
    nonlin = nlse_filter(nonlin, P);
end
