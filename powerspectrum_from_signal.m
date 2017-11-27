function P = powerspectrum_from_signal(x)
% Given signals in x (a matrix whose columns are the vectors),
% returns their power spectra in P (a matrix of same size as x).

    fftx = fft(x);
    P = fftx .* conj(fftx);

end
