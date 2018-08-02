function [x, alpha] = bigmra_em_1D(Y, sigma, x, alpha)
% NB, Aug. 1st, 2018
    
    [W, N] = size(Y);
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
    Zt = @(u) u(W+1:end, :);
    %
    R = @(k, u) circshift(u, k, 1);
    
    % Precomputations with the observed windows in Y
    FPtY = fft(Pt(Y));

    stop = false;
    while ~stop
    
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

        % Update x by solving a linear system
        b = real(Zt(ifft(sum(FPtY .* conj(fft(WW)), 2))));
        A = zeros(L, L);
        for ell = 0 : W+L-1
            Aell = P(R(ell, Z(eye(L))));
            A = A + v(ell+1)*(Aell'*Aell);
        end
        x_new = A\b;

        % Update alpha
        alpha_new = v(1)/N;       % this is naturally in the interval [0, 1]
        
        % Stopping criterion (somewhat arbitrary..)
        if norm(x_new - x) < 1e-6*norm(x)
            stop = true;
        end
        
        x = x_new;
        alpha = alpha_new;
        rho = rhofun(alpha);
        
    end
    
end
