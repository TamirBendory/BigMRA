function x = NCN_cryo(x0, alpha, beta, tolGrad, tolEig, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda)

x = x0;
[f,g] = costgrad_fminunc(x, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda);
f0 = f

while norm(g) > tolGrad
    [f,g,h] = costgrad_fminunc(x, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda);
    disp(['f = ' num2str(f)])
    [V,D] = eig(h);
    d = diag(D);
    
    
    d(abs(d)>=tolEig) = abs(d(abs(d)>=tolEig));
%     d(abs(d)<tolEig) = tolEig;
    d(abs(d)<tolEig) = 0;
    d = d.^(-1);
    d(isinf(d)) = 0;
    d = real(V*diag(d)*V')*g;
    
%     d = h\g; % ~inv(h)*g

    eta = 1;
    xn = x - eta*d;
    fn = costgrad_fminunc(xn, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda);
    while fn > f - alpha*eta*(g.'*d)
        eta = eta*beta;
        xn = x - eta*d;
        fn = costgrad_fminunc(xn, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, gamma, lambda);
    end
    disp(['eta = ' num2str(eta)])
    x = xn;
end
    