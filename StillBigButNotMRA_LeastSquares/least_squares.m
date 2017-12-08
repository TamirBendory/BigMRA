function x_est = least_squares(M1, M2, M3, W, sigma, N, L, m)
% Still-Big-But-Not-MRA via least-squares solve. Homogeneous for now.

    manifold = euclideanfactory(W, 1); %% !!! W or L: chose wisely -- W seems to help
    
    params.M1 = M1;
    params.M2 = M2;
    params.M3 = M3;
    params.W = W;
    params.sigma = sigma;
    params.N = N;
    params.m = m;
    
    problem = manoptAD(manifold, @least_squares_cost, params);
    
    warning('off', 'manopt:getHessian:approx');
    
    [x_est, loss] = trustregions(problem); %#ok<ASGLU>
%     
%     % Trying to escape local optima...
%     x_base = x_est;
%     for shift = 1 : L-1
%         x_est_shifted = circshift_ad(x_base, shift);
%         [x_est_shifted, loss_shifted] = trustregions(problem, x_est_shifted);
%         if loss_shifted < loss
%             loss = loss_shifted;
%             x_est = x_est_shifted;
%         end
%     end
%     
    warning('on', 'manopt:getHessian:approx');
    
end
