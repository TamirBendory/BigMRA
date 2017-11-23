function [mean_est, P_est, B_est] = invariants_from_data_no_debias(X, w)
% See invariants_from_data -- but no debiasing is attempted.
% The bispectrum is computed on the data directly (no centering).
% w is a weight vector: nonnegative entries, length = size(X, 2)

    [N, M] = size(X);
    
    if ~exist('w', 'var') || isempty(w)
        w = ones(M, 1);
    end
    w = w(:);
    assert(length(w) == M, 'w must have length M = size(X, 2)');
    assert(all(w >= 0), 'w must be nonnegative');
    w = w / sum(w);
    
    %% Estimate the mean and subtract it from the observations.
    %  This also gives an estimate of the first component of the fft of x.
    %  Centering the observations helps for the bispectrum estimation part.
    mean_est = mean(X, 1) * w;
    
    %% Prepare fft's
    X_fft = fft(X);

    %% Estimate the power spectrum (gives estimate of modulus of fft of x).
    P_est = abs(X_fft).^2 * w;

    %% Estimate the bispectrum
    if nargout >= 3
        
        B_est = zeros(N);
        parfor m = 1 : M
            xm_fft = X_fft(:, m);
            Bm = (xm_fft*xm_fft') .* circulant(xm_fft);
            B_est = B_est + w(m) * Bm;
        end
        
        % % Below is a way to compute the average bispectrum without loop,
        % % but it does not take advantage of parallelization and creates a
        % % potentially huge matrix as intermediate step (N^2 * M in size).
        % vec = @(A) A(:);
        % I1 = vec((1:N)' * ones(1, N));
        % I2 = vec(ones(N, 1) * (1:N));
        % I3 = vec(circulant(1:N));
        % B_est_2 = reshape((X_fft(I1, :) .* conj(X_fft(I2, :)) .* X_fft(I3, :)) * w, [N, N]);
        
    end

end
