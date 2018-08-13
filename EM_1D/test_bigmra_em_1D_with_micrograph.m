% Test file for bigmra_em_1D using a 1D micrograph
% NB, Aug. 2, 2018
clear all; clf; clc;

% rng(6875769);

%%
L = 21;             % Length of the true signal
W = 1*L;            % Length of the observed windows (W >= L)
x = 0+randn(L, 1);  % True signal
sigma = 3.0;        % Standard deviation of the additive white Gaussian noise

N = 3e6;   % Length of the micrograph
m = 3e4;   % Desired number of occurrences of x in the micrograph
[y, placed] = generate_clean_micrograph_1D(x, W, N, m);
y = y + sigma*randn(size(y));
% Extract windows from the micrograph for EM
skip = L; %typical values: 0, L, -(W-1)
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

%%
% Run EM a first time
[x_est, alpha_est] = bigmra_em_1D(Y, sigma, x0, alpha0);
fprintf('Initial EM done.\n');

%% Compute some kind of log-likelihood for each cyclic shift of x_est, and
%  use the most promising one to warm-start a second EM run.
shiftQs = zeros(L, 1);
for k = 0 : L-1
    x_shifted = circshift(x_est, k);
    shiftQs(k+1) = bigmra_1D_likelihood_proxy(Y, sigma, x_shifted, alpha_est);
end

% this forces us to pick a non-trivial shift
shiftQs(1+0) = -inf;

[a, b] = max(shiftQs);

% Assuming odd L
subplot(2, 1, 1);
plot((0:L-1)-(L-1)/2, fftshift(shiftQs), '.-');
hold all;
plot([-10, 10], (.95*(max(shiftQs)-min(shiftQs))+min(shiftQs))*[1, 1], 'r--');
hold off;

fprintf('Attempting a shift of %d.\n', b-1);

[x_est2, alpha_est2] = bigmra_em_1D(Y, sigma, circshift(x_est, b-1), alpha_est);
fprintf('Secondary EM done.\n');

Q1 = bigmra_1D_likelihood_proxy(Y, sigma, x_est, alpha_est);
Q2 = bigmra_1D_likelihood_proxy(Y, sigma, x_est2, alpha_est2);

% This Q2 > Q1 criterion does not work ... I saw it fail pretty obviously.
% if Q2 > Q1
if norm(x_est2) > norm(x_est)
    x_est = x_est2;
    alpha_est = alpha_est2;
    fprintf('Picked the second one.\n');
else
    fprintf('Picked the first one.\n');
end

%%

%{
% Now, we're going to slide the first estimate L-1 times to the left and
% L-1 times to the right, and re-run EM everytime, with warm starts.
xs = zeros(L, 2*L-1);
xs(:, L) = x_est;
alphas = zeros(2*L-1, 1);
alphas(L) = alpha_est;

fprintf('Sliding right: ');
for right = L-1 : -1 : 1
    x0 = [randn(1); xs(1:end-1, right+1)];
    alpha0 = alphas(right+1);
    [xs(:, right), alphas(right)] = bigmra_em_1D(Y, sigma, x0, alpha0);
    fprintf('.');
end
fprintf('\nSliding left:  ');
for left = L+1 : +1 : 2*L-1
    x0 = [xs(2:end, left-1) ; randn(1)];
    alpha0 = alphas(left-1);
    [xs(:, left), alphas(left)] = bigmra_em_1D(Y, sigma, x0, alpha0);
    fprintf('.');
end
fprintf('\n');
%}

%%
%{
clf;

for plt = 1 : 2*L-1
    subplot(12, 7, plt);
    plot(0:L-1, xs(:, plt));
end

% Selection of the best one is a bit cheating: we compare to the true x
[err, best] = min(sqrt(sum((xs - repmat(x, [1, 2*L-1])).^2, 1)));
fprintf('Relative error of the best estimator of x: %g\n', err/norm(x));

subplot(12, 7, (43):(12*7));
plot(0:L-1, x, 'o-', 0:L-1, xs(:, best), '.-');
legend('true x', 'estimated x');
%}

%% Display the results (older code)

[~, shift] = max(real(ifft(conj(fft(x)).*fft(x_est))));
shift = rem(shift, L);
if shift > L/2
    shift = shift - L;
end
I = 0:(L-1);
shifted = I - shift + 1;
subplot(2, 1, 2);
plot(I, x, '.-', shifted, x_est, 'o-');
legend('true x', 'estimated x');
title(sprintf('Estimated alpha: %g', alpha_est));
