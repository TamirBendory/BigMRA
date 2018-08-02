% Test file for bigmra_em_1D using a 1D micrograph
% NB, Aug. 2, 2018
clear all; close all; clc;

%%
L = 21;             % Length of the true signal
W = 1*L;            % Length of the observed windows (W >= L)
x = 1+randn(L, 1);  % True signal
sigma = 2.0;        % Standard deviation of the additive white Gaussian noise

N = 3e6;   % Length of the micrograph
m = 1e4;   % Desired number of occurrences of x in the micrograph
[y, placed] = generate_clean_micrograph_1D(x, W, N, m);
y = y + sigma*randn(size(y));
% Extract windows from the micrograph for EM
y = [y ; zeros((W+L)-mod(numel(y), W+L), 1)]; % ensure length of y can be divided by W+L
% y = y(1:(numel(y)-mod(numel(y), W+L)));
Y = reshape(y, [W+L, numel(y)/(W+L)]);
Y = Y(1:W, :);

%% Pick initial guesses for the signal and for alpha
x0 = randn(L, 1);
alpha0 = rand(1);

% Run EM
[x_est, alpha_est] = bigmra_em_1D(Y, sigma, x0, alpha0);

% Display the results
[~, shift] = max(real(ifft(conj(fft(x)).*fft(x_est))));
shift = rem(shift, L);
if shift > L/2
    shift = shift - L;
end
I = 0:(L-1);
shifted = I - shift + 1;
plot(I, x, '.-', shifted, x_est, 'o-');
legend('true x', 'estimated x');
title(sprintf('Estimated alpha: %g', alpha_est));
