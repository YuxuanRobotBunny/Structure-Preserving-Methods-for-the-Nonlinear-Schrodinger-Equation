function u = nlse_filter(u, P)

    if P.dealias
        u = ifft(P.dealias_mask .* fft(u));
    end
end
