function X = fft_series_approx_2(x, Nu, Nv, threshold_dB)
    N = Nv * Nu;

    X = 1 / N * fftshift(fft2(x));
    X_dB = 20 * log10(abs(X));

    max_X_dB = max(X_dB, [], "all");

    X(X_dB < max_X_dB - threshold_dB) = 0;
    X = sparse(X);
end
