% Test file for bigmra_em_1D
% NB, Aug. 1, 2018

L = 21;             % Length of the true signal
W = 1*L;            % Length of the observed windows (W >= L)
x = 1+randn(L, 1);  % True signal
alpha = .50;        % Probability that a window doesn't contain any signal
sigma = 2.0;        % Standard deviation of the additive white Gaussian noise
N = 1e4;            % Number of observed windows

% Probability distribution of the shifts (0 to W+L-1) in the observed windows.
% No shift leads to no signal at all in the window, whereas all other
% shifts are equally likely and lead to at least some signal in the window.
rho = [alpha ; ones(W+L-1, 1)*(1-alpha)/(W+L-1)];

% The linear operators of the forward model (works on matrices too):
Z = @(u) [zeros(W, size(u, 2)) ; u];    % Zero pad with W zeros at beginning
R = @(k, u) circshift(u, k, 1);         % Circular shift by k entries
P = @(u) u(1:W, :);                     % Keep only W first entries

% Store the N observed windows as columns of W.
Y = zeros(W, N);
% Pick the N true shifts according to distribution rho.
shifts = randsample(0:(W+L-1), N, true, rho);
for k = 1 : N
    Y(:, k) = P(R(shifts(k), Z(x)));
end
% Add noise
Y = Y + sigma*randn(W, N);

%% Pick initial guesses for the signal and for alpha
x0 = randn(L, 1);
alpha0 = rand(1);

% Run EM
[x_est, alpha_est] = bigmra_em_1D(Y, sigma, x0, alpha0);

% Display the results
disp([x, x_est]);
fprintf('\n');
disp([alpha, alpha_est]);
[~, shift] = max(real(ifft(conj(fft(x)).*fft(x_est))));
shift = rem(shift, L);
I = 0:(L-1);
shifted = I - shift + 1;
plot(I, x, '.-', shifted, x_est, 'o-');
legend('true x', 'estimated x');
title(sprintf('True alpha: %g, estimated alpha: %g, empirical alpha: %g', ...
      alpha, alpha_est, nnz(shifts == 0)/N));
