function x_est = least_squares(M1, M2, M3, W, sigma, N, L, m)
% Still-Big-But-Not-MRA via least-squares solve. Homogeneous for now.

    manifold = euclideanfactory(L, 1);
    
    params.M1 = M1;
    params.M2 = M2;
    params.M3 = M3;
    params.W = W;
    params.sigma = sigma;
    params.N = N;
    params.m = m;
    
    problem = manoptAD(manifold, @least_squares_cost, params);
    
    warning('off', 'manopt:getHessian:approx');
    
    x_est = trustregions(problem);
    
    warning('on', 'manopt:getHessian:approx');
    
end
