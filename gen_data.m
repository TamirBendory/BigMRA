function [y, y_clean, ind, class] = gen_data(x, N, m, sigma, W)
% Generate data following the Big-MRA model
%
% Input:
%
% x - the K signals to be hidden in noise (LxK matrix)
% N - total length of the observation
% m - vector of length K such that x(:, k) is repeated (nearly) m(k) times
% sigma - noise level
% W - length of intended sliding window: only used to ensure proper separation
%
% Output:
%
% y - noisy measurement of length N
% y_clean - clean (noiseless) measurement
% ind - locations of x in y_clean
% class - for each location in ind, specifies which signal (1..K) is there.

    y_clean = zeros(N, 1);
    [L, K] = size(x);
    ind = randi(N-L+1, sum(m), 1);
    ind = sort(ind);
    indn = ind(1);

    % removing adjacent signals (need to be rewritten)

    for i = 1:length(ind)-1
        if ind(i+1) - ind(i) > 2*W
            indn = [indn; ind(i+1)];
        end
    end
    ind = indn;
    % drawing classes
    class = randi([1,K], length(ind), 1);

    % generating noiseless measurement(depends on the window length)
    for i = 1:length(ind)
        y_clean( ind(i) : ind(i)+L-1 ) = y_clean( ind(i) : ind(i)+L-1 ) +  x(:, class(i));
    end

    % add noise
    y = y_clean + sigma*randn(N, 1);

end
