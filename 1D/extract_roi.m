function X_short = extract_roi(X_long, L_short)
% Extract region of interest of K signals
%
% function X_short = extract_roi(X_long, L_short)
%
% Given K signals of length L_long as the columns of a matrix X_long,
% extracts from each signal a contiguous subsignal (with cyclic indexing)
% of length L_short with largest 2-norm (the region of interest), and
% returns these as the columns of X_short.

    [L_long, K] = size(X_long);
    
    assert(L_Long >= L_short, 'Can only extract a shorter region.');

    X_short = zeros(L_short, K);
    for k = 1 : K
        for s = 0 : (L_long - 1)
            x = circshift(X_long(:, k), s);
            x = x(1:L_short);
            if norm(x) > norm(X_short(:, k))
                X_short(:, k) = x;
            end
        end
    end
    
end
