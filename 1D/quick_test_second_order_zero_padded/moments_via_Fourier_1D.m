% NB: see notes Jan. 31, 2018 (notebook 33)
clear;
close all;
clc;

% This script shows precise relationaships between the moments we compute
% in BigMRA and monomials in Fourier domain for zero-padded signals.

% Pick a micrograph length N and a signal length L, and a micrograph
% (the code works both for real and complex data.)
L = 100;
N = 2*L-1; % N = 2L-1 is what we would see for the optimization part
x = randn(N, 1); % + 1i*randn(N, 1);

tic();

% Since we are only interested in autocorrelations with maximum shift
% separation up to L-1, we only need to zero pad with L-1 elements. If we
% actually did this, the micrograph would end up with length N+L-1.
% The following call to fft implicitly zero pads x up to length N+L-1, then
% takes the DFT. The returned vector has length N+L-1.
Fxz = fft(x, N+L-1);

% Power spectrum of the zero padded signal
Pxz = real(Fxz).^2 + imag(Fxz).^2; % Fxz .* conj(Fxz);

% Inverse DFT of Pxz, and extract the L leading entries
M2_Fourier = ifft(Pxz);
M2_Fourier = M2_Fourier(1:L);

fprintf('M2 Fourier: %.2g [ms]\n', 1000*toc());

tic();

% Do a direct computation of the second moment (s is the shift)
M2_direct = zeros(L, 1);
for s = 0 : (L-1)
    range1 = s : (N-1);
    range2 = range1 - s;
    M2_direct(1+s) = sum(x(1+range1) .* conj(x(1+range2)));
end

fprintf('M2 direct:  %.2g [ms]\n', 1000*toc());

fprintf('Error on 2nd moment: %g\n', norm(M2_direct - M2_Fourier, 2)/norm(M2_direct, 2));

tic();

% Now for third order moments via Fourier: take the usual bispectrum of the
% zero padded signal: this might not be practical for large N.
C = circulant(Fxz);
Bxz = (Fxz*Fxz') .* C;

% Inverse 2D DFT of Bxz
M3_Fourier = ifft2(Bxz);

% Extract leading LxL submatrix and zero-out coefficients with separation
% more than L-1.
M3_Fourier = M3_Fourier(1:L, 1:L);
M3_Fourier = transpose(flipud((tril(rot90(M3_Fourier))))); % this has got to have a nicer expression

fprintf('M3 Fourier: %.2g [ms]\n', 1000*toc());

tic();

% Do a direct computation of the third moment (s1, s2 are shifts)
% Do a direct computation of the second moment (s is the shift)
M3_direct = zeros(L, L);
for s1 = 0 : (L-1)
    for s2 = 0 : (L-1)
        if s1 + s2 <= L-1 % other moments are 0
            range1 = s2 : (N-1-s1);
            range2 = range1 - s2;
            range3 = range1 + s1;
            M3_direct(1+s1, 1+s2) = sum(x(1+range1) .* conj(x(1+range2)) .* x(1+range3));
        end
    end
end

fprintf('M3 direct:  %.2g [ms]\n', 1000*toc());

fprintf('Error on 3rd moment: %g\n', norm(M3_direct - M3_Fourier, 'fro')/norm(M3_direct, 'fro'));
