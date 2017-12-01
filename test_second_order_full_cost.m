function f = test_second_order_full_cost(x, params)

%     L = params.L;
    W = params.W;
    M = params.M;
    
    x_padded = [zeros(W-1, 1) ; x ; zeros(W-1, 1)];

    M_x = zeros(W, W);
    count = length(x_padded) - (W-1);
    y = zeros(W, 1); %#ok<PREALL> % for AD
    ffty = zeros(W, 1); %#ok<PREALL> % for AD
    for k = 1 : count
        y = x_padded(k + (0:(W-1)));
        ffty = fft(y);
        M_x = M_x + ffty * ffty';
    end
    M_x = M_x / count;
    
    M_minus_Mx = M - M_x;
    
    f = .5*sum( M_minus_Mx(:) .* conj(M_minus_Mx(:)) );

end
