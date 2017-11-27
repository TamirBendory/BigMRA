function f = test_second_order_cost(x, params)

%     L = params.L;
    W = params.W;
    M = params.M;
    
    x_padded = [zeros(W-1, 1) ; x ; zeros(W-1, 1)];

    M_x = zeros(W, 1);
    count = length(x_padded) - (W-1);
    y = zeros(W, 1); %#ok<PREALL> % for AD
    ffty = zeros(W, 1); %#ok<PREALL> % for AD
    for k = 1 : count
        y = x_padded(k + (0:(W-1)));
        ffty = fft(y);
        M_x = M_x + ffty .* conj(ffty);
    end
    M_x = M_x / count;
    
    f = .5*sum((M - M_x).^2);

end
