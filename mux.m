function mu = mux(x, m, N)
% Produce the sliding-window average mean one would see without noise
% 
% x - matrix of size L x K, where each column is a signal
% m - vector of length K, such that x(:, k) appears m(k) times
% N - length of the long observation
%
% mu is a scalar

    [L, K] = size(x); %#ok<ASGLU>

    m = m(:);
    assert(length(m) == K);

    mu = (sum(x, 1)*m)/N;

end
