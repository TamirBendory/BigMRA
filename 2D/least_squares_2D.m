function [X_est, problem] = least_squares_2D(M1, M2, M3, W, sigma, N, L, m, list2, list3, X0)

    manifold = euclideanfactory(W, W); % choose wisely if do LxL or WxW -- WxW seems necessary to avoid local opts.
    
    params.W = W;
    params.N = N;
    params.m = m;
    params.M1 = M1;
    params.M2 = M2;
    params.M3 = M3;
    params.list2 = list2;
    params.list3 = list3;
    params.sigma = sigma;
    
% % % % %     problem = manoptAD(manifold, @least_squares_cost_2D, params);

    problem.M = manifold;
	problem.costgrad = @(X) least_squares_2D_cost_grad(X, params);
% % %     checkgradient(problem); pause;
    
    warning('off', 'manopt:getHessian:approx');
    
%     keyboard;

    opts = struct();
    opts.maxtime = 240;
    [X_est, loss] = trustregions(problem, X0, opts); %#ok<ASGLU>

%     [X_est, loss] = conjugategradient(problem, X0, opts); %#ok<ASGLU>

    warning('on', 'manopt:getHessian:approx');
    
end
