clear;
close all;
clc;

%% Problem parameters
L = 11;             % length of signal
sigma = 0;          % noise level
W = 4*L-1;

%% Generate signals
x_true = randn(L, 1);

%% Generate data
N = 10*L;
y = zeros(N, 1);
y(L+1:2*L) = x_true;
% y(50+L+1:50+2*L) = x_true;
m = 1;

%%  Process the data: compute moments

% computing the empirical invariants of the data
[M1, M2, M3] = moments_from_data_no_debias(y, W);

%% Run the least squares

x_est = least_squares(M1, M2, M3, W, sigma, N, L, m);


%% Visualization
if length(x_est) == W
    x_true_zp = [x_true ; zeros(W-L, 1)];
    x_est = align_to_reference(x_est, x_true_zp);
    plot(1:W, x_true_zp, 'o-', 1:W, x_est, '.-');
    title(sprintf('Relative error: %g', norm(x_true_zp - x_est, 'fro')/norm(x_true_zp, 'fro')));
elseif length(x_est) == L
%     x_est = align_to_reference(x_est, x_true);
    plot(1:L, x_true, 'o-', 1:L, x_est, '.-');
    title(sprintf('Relative error: %g', norm(x_true - x_est, 'fro')/norm(x_true, 'fro')));
end