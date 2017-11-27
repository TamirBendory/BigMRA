function B = bsx(x, m, N, W)
% Produce the sliding-window average bispectrum one would see without noise
% 
% x - matrix of size L x K, where each column is a signal
% m - vector of length K, such that x(:, k) appears m(k) times
% N - length of the long observation
% W - length of the sliding window
%
% B is a matrix of size W x W.

    [L, K] = size(x);

    m = m(:);
    assert(length(m) == K);

    B = zeros(W);

    for k = 1 : K

        B = B + m(k)*(W-L+1)*bispectrum_from_signal([x(:, k) ; zeros(W-L, 1)]);

        for ii = 0:(L-2)
            B = B + m(k)*bispectrum_from_signal([zeros(W-ii-1, 1) ; x(1:(ii+1), k)]);
        end

        for ii = 1:(L-1)
            B = B + m(k)*bispectrum_from_signal([x(ii+1:end, k) ; zeros(W-L+ii, 1)]);
        end

    end

    B = B/N;

end
