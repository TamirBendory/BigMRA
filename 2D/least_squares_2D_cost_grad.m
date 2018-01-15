function [f, g] = least_squares_2D_cost_grad(X, params)

    N = params.N;
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
    g = zeros(L, L);
    
    % inner product over matrices
    inner = @(X, Y) X(:)'*Y(:);
    % Zero-padding operator
    ZP = @(X) [X zeros(L, W); zeros(W, L+W)];
    % Adjoint of the zero-padding operator
    ZPadj = @(Y) Y(1:L, 1:L);
    % 2D circular shift operator
    CS = @(X, k) circshift(X, k);
    % Adjoint of the 2D circular shift operator
    CSadj = @(X, k) circshift(X, -k);
    
    w1 = 1/m;
    w2 = 1/(W^2)/m;
    w3 = 1/(W^4)/m;
    
    % First-order moment, forward model
    M1 = m*sum(X(:));
    R1 = M1 - M1data;
    f = f + w1*R1^2;
    g = g + w1*2*R1*ones(L, L)*m;
    
    ZPX = ZP(X);
    
    % Second-order moments, forward model
    n2 = size(list2, 1);
    parfor k = 1 : n2
        
        shift = list2(k, :);
        CX = CS(ZPX, shift);
        M2k = m*inner(CX, ZPX);
        
        M2k = M2k + bias2(k);
        
        R2 = M2k - M2data(k);
        f = f + w2*R2^2;
        
        % remainder for gradient
        CaX = CSadj(ZPX, shift);
        G = ZPadj(CX + CaX);
        g = g + w2*2*R2*G*m;
        
    end
    
    % Third-order moments, forward model
    n3 = size(list3, 1);
    parfor k = 1 : n3
        
        line = list3(k, :);     % intermediate step to help parfor
        shift1 = line([1, 2]);
        shift2 = line([3, 4]);
        CX1 = CS(ZPX, shift1);
        CX2 = CS(ZPX, shift2);
        T1 = CX1.*CX2;
        M3k = m*inner(T1, ZPX);

        M3k = M3k + bias3(k);
        
        R3 = M3k - M3data(k);
        f = f + w3*R3^2;
        
        % remainder for gradient
        T2 = CSadj(CX1 .* ZPX, shift2);
        T3 = CSadj(CX2 .* ZPX, shift1);
        G = ZPadj(T1 + T2 + T3);
        g = g + w3*2*R3*G*m;
        
    end
    
end
