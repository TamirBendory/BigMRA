function f = least_squares_cost(x, params)

    M1 = params.M1;
    M2 = params.M2;
    M3 = params.M3;
    sigma = params.sigma;
    N = params.N;
    m = params.m;
    W = params.W;
    
    [M1x, M2x, M3x] = moments_from_signal_with_bias(x, W, N, m, sigma);
    
    w1 = 1;
    w2 = 1/W;
    w3 = 1/(W^2);
    
    f = w1 * ( M1 - M1x )^2;
    
    f = f + w2 * sum( ( M2 - M2x ).^2 );
    
    M3_diff = M3(:) - M3x(:);
    f = f + w3 * sum( M3_diff.^2 );

    f = f / m;
    
end
