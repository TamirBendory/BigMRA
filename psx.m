function P = psx(x, m, N, W)
% Produce the sliding-window average power spectrum one would see w/o noise
% 
% x - matrix of size L x K, where each column is a signal
% m - vector of length K, such that x(:, k) appears m(k) times
% N - length of the long observation
% W - length of the sliding window
%
% P is a vector of length W

    [L, K] = size(x);

    m = m(:);
    assert(length(m) == K);

    P = zeros(W, 1);

    for k = 1 : K

        P = P + m(k)*(W-L+1)*powerspectrum_from_signal([x(:, k) ; zeros(W-L, 1)]);

        for ii = 0:(L-2)
            P = P + m(k)*powerspectrum_from_signal([zeros(W-ii-1, 1) ; x(1:(ii+1), k)]);
        end

        for ii = 1:(L-1)
            P = P + m(k)*powerspectrum_from_signal([x(ii+1:end, k) ; zeros(W-L+ii, 1)]);
        end

    end

    P = P/N;

end
