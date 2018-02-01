% NB: see notes Jan. 31, 2018 (notebook 33)

% This script shows precise relationaships between the moments we compute
% in BigMRA and monomials in Fourier domain for zero-padded signals.

% Pick a signal length and a true signal (the code works both for real and
% complex signals.)
L = 5;
% x = randn(L, 1);
x = randn(L, 1) + 1i*randn(L, 1);

% Pick a zero-padding length; the only requirement is W >= 2L-1, and we
% might as well take is as small as possible, so:
W = 2*L-1;

% Zero-pad the signal up to length W, then take its DFT
Fxz = fft(x, W);

% Power spectrum of the zero padded signal
Pxz = real(Fxz).^2 + imag(Fxz).^2; % Fxz .* conj(Fxz);

% Inverse DFT of Pxz, and extract the L leading entries
M2_Fourier = ifft(Pxz);
M2_Fourier = M2_Fourier(1:L);

% Do a direct computation of the second moment (s is the shift)
M2_direct = zeros(L, 1);
for s = 0 : (L-1)
    range1 = s : (L-1);
    range2 = range1 - s;
    M2_direct(1+s) = sum(x(1+range1) .* conj(x(1+range2)));
end

fprintf('Error on 2nd moment: %g\n', norm(M2_direct - M2_Fourier, 2)/norm(M2_direct, 2));

% Now for third order moments via Fourier: take the usual bispectrum of the
% zero padded signal
C = circulant(Fxz);
Bxz = (Fxz*Fxz') .* C;

% Inverse 2D DFT of Bxz, and extract leading LxL submatrix
M3_Fourier = ifft2(Bxz);
M3_Fourier = M3_Fourier(1:L, 1:L);

% Do a direct computation of the third moment (s1, s2 are shifts)
% Do a direct computation of the second moment (s is the shift)
M3_direct = zeros(L, L);
for s1 = 0 : (L-1)
    for s2 = 0 : (L-1)
        if s1 + s2 <= L-1 % other moments are 0
            range1 = s2 : (L-1-s1);
            range2 = range1 - s2;
            range3 = range1 + s1;
            M3_direct(1+s1, 1+s2) = sum(x(1+range1) .* conj(x(1+range2)) .* x(1+range3));
        end
    end
end

fprintf('Error on 3rd moment: %g\n', norm(M3_direct - M3_Fourier, 'fro')/norm(M3_direct, 'fro'));
