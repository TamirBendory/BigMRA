function f = least_squares_cost_2D(X, params)

    N = params.N;
    m = params.m;
    W = params.W;
    M1 = params.M1;
    M2 = params.M2;
    M3 = params.M3;
    sigma = params.sigma;
    list2 = params.list2;
    list3 = params.list3;
    
    [M1x, M2x, M3x] = moments_from_signal_with_bias_2D(X, W, N, m, sigma, list2, list3);
    
    w1 = 1;
    w2 = 1/(W^2);
    w3 = 1/(W^4);
    
    f = w1 * ( M1 - M1x )^2;
    
    f = f + w2 * sum( ( M2 - M2x ).^2 );
    
    f = f + w3 * sum( ( M3 - M3x ).^2 );
    
    f = f / m;
    
end
