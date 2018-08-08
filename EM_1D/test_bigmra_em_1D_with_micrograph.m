% Test file for bigmra_em_1D using a 1D micrograph
% NB, Aug. 2, 2018
clear all; close all; clc;

%rng(6875769);
rng(687);

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
skip = -(W-1); %typical values: 0, L, -(W-1)
if skip >= 0
    y = [y ; zeros((W+skip)-mod(numel(y), W+skip), 1)]; % ensure length of y can be divided by W+skip
    % y = y(1:(numel(y)-mod(numel(y), W+skip)));
    Y = reshape(y, [W+skip, numel(y)/(W+skip)]);
    Y = Y(1:W, :);
else
    num_windows = floor((length(y)-W)/(W+skip));
    Y = zeros(W, num_windows);
    for k = 1 : num_windows
        Y(:, k) = y((1:W)+(k-1)*(W+skip));
    end
end
fprintf('Done, with %d windows.\n', size(Y, 2));

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
