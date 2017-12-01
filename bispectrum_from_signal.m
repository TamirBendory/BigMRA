function B = bispectrum_from_signal(x)
% Given signals in x (a matrix whose columns are different signals),
% returns their bispectra in B (a 3D matrix whose slices are bispectra).

    [L, K] = size(x);
    fftx = fft(x);
    B  = zeros(L, L, K);
    for k = 1 : K
        B(:, :, k) = (fftx(:, k) * fftx(:, k)') .* circulant_AD(fftx(:, k));
    end

end
