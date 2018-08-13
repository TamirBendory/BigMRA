function Q = bigmra_1D_likelihood_proxy(Y, sigma, x, alpha)
% NB, Aug. 9, 2018 -- some type of proxy for the log likelihood of
% x, alpha given Y, sigma.
%
% It's unclear whether this function is reasonable.

    assert(isscalar(alpha));
    
    [W, N] = size(Y); %#ok<ASGLU>
    L = size(x, 1);
    rhofun = @(alpha) [alpha ; ones(W+L-1, 1)*(1-alpha)/(W+L-1)];
    rho = rhofun(alpha);
    
    % The following operators apply to matrices as well, column by column.
    %
    % Project to the first W entries of a vector of length W+L, and adjoint
    P = @(u) u(1:W, :);
    Pt = @(u) [u ; zeros(L, size(u, 2))];
    %
    % Prepend W zeros to a vector of length L, and adjoint
    Z = @(u) [zeros(W, size(u, 2)) ; u];
    Zt = @(u) u(W+1:end, :); %#ok<NASGU>
    %
    R = @(k, u) circshift(u, k, 1);
    
    % Precomputations with the observed windows in Y
    FPtY = fft(Pt(Y));

    sqnorms = zeros(W+L, 1);
    for ell = 0 : W+L-1
        sqnorms(ell+1) = norm(P(R(ell, Z(x))), 2)^2;
    end

    inner_products = real(ifft(bsxfun(@times, FPtY, conj(fft(Z(x))))));

    T = bsxfun(@minus, inner_products, sqnorms/2)/sigma^2;

    T = bsxfun(@minus, T, max(T, [], 1));

    WW = bsxfun(@times, exp(T), rho);

    WW = bsxfun(@times, WW, 1./sum(WW, 1));

    v = sum(WW, 2);
    
    inner = @(A, B) A(:)'*B(:);
    sqnormsY = sum(Y.^2, 1);
    allsqnorms = sqnorms + sqnormsY - 2*inner_products;
    
    Q = v'*log(rho) - inner(WW, allsqnorms)/(2*sigma^2);
    
end
