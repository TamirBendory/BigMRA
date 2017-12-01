function f = bigmra_lsq_cost(x, params)

    M1 = params.M1;
    M2 = params.M2;
    M3 = params.M3;
    sigma = params.sigma;
    N = params.N;
    m = params.m;
    W = params.W;
    
    M1_est = mux(x, m, N);
    M2_est = psx(x, m, N, W);
    M3_est = bsx(x, m, N, W);
    
    % Add bias on power spectrum
    M2_est = M2_est + sigma^2*W;

    % Add bias on bispectrum
    A = eye(W);
    A(:, 1) = A(:, 1) + 1;
    A(1, :) = A(1, :) + 1;
    M3_est = M3_est + sigma^2*W^2*M1*A;
    
    w1 = 1;
    w2 = 1/W;
    w3 = 1/(W^2);
    
    f = w1 * ( M1 - M1_est )^2;
    
    f = f + w2 * sum( ( M2 - M2_est ).^2 );
    
    M3_diff = M3(:) - M3_est(:);
    f = f + w3 * sum( M3_diff .* conj(M3_diff) );

end
