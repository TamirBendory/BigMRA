function [f, g] = least_squares_2D_cost_grad_new(X, params)

    N = params.N; %#ok<NASGU>
    m = params.m;
    W = params.W;
    M1data = params.M1;
    M2data = params.M2;
    M3data = params.M3;
    sigma = params.sigma; %#ok<NASGU>
    list2 = params.list2;
    list3 = params.list3;
    bias2 = params.bias2;
    bias3 = params.bias3;
    
    % This L is as defined by the size of the optimization variable: it
    % need not be equal to the L as defined by the size of the true signal.
    L = size(X, 1);
    assert(size(X, 2) == L, 'X must be square');
    
    f = 0;
    
    if nargout >= 2
        g = zeros(L, L);
    end
    
    w1 = 1;
    w2 = 1/(W^2);
    w3 = 1/(W^4);
    
    % First-order moment, forward model
    M1 = m*sum(X(:));
    R1 = M1 - M1data;
    f = f + .5*w1*R1^2/m;
    
    if nargout >= 2
        g = g + w1*R1*ones(L, L);
    end
    
    
    
    % Second-order moments, forward model
    n2 = size(list2, 1);
    for k = 1 : n2
        
        shift1 = list2(k, :);
        
        vals1 = [0, shift1(1)];
        range1 = (1+max(vals1)) : (L+min(vals1));
        vals2 = [0, shift1(2)];
        range2 = (1+max(vals2)) : (L+min(vals2));
        X1 = X(range1, range2); %#ok<PFBNS>
        X2 = X(range1-shift1(1), range2-shift1(2));
        
        M2k = m*sum(X1(:) .* X2(:)) + bias2(k);
        
        R2k = M2k - M2data(k);
        f = f + .5*w2*R2k^2/m;
        
        % This part of the code for gradient only
        if nargout >= 2
            T1 = zeros(L);
            T1(range1, range2) = X2;
            T2 = zeros(L);
            T2(range1-shift1(1), range2-shift1(2)) = X1;
            G = T1 + T2;
            g = g + w2*R2k*G;
        end
        
    end
    
    
    % Third-order moments, forward model
    n3 = size(list3, 1);
    for k = 1 : n3
        
        shifts = list3(k, :);
        shift1 = shifts([1, 2]);
        shift2 = shifts([3, 4]);

        vals1 = [0, shift1(1), shift2(1)];
        range1 = (1+max(vals1)) : (L+min(vals1));
        vals2 = [0, shift1(2), shift2(2)];
        range2 = (1+max(vals2)) : (L+min(vals2));
        X1 = X(range1, range2); %#ok<PFBNS>
        X2 = X(range1-shift1(1), range2-shift1(2));
        X3 = X(range1-shift2(1), range2-shift2(2));
        X1X2 = X1 .* X2;
        
        M3k = m*sum(X1X2(:) .* X3(:)) + bias3(k);
        
        R3k = M3k - M3data(k);
        f = f + .5*w3*R3k^2/m;
        
        % This part of the code for gradient only
        if nargout >= 2
            T1 = zeros(L);
            T1(range1, range2) = X2.*X3;
            T2 = zeros(L);
            T2(range1-shift1(1), range2-shift1(2)) = X1.*X3;
            T3 = zeros(L);
            T3(range1-shift2(1), range2-shift2(2)) = X1X2;
            G = T1 + T2 + T3;
            g = g + w3*R3k*G;
        end
        
    end
    
end
