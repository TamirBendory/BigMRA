function [X_est, problem, stats] = least_squares_2D(M1, M2, M3, W, sigma, N, L, m, list2, list3, X0,method,opt_window,max_iter) %#ok<INUSL>

    params.W = W;
    params.N = N;
    params.m = m;
    params.M1 = M1;
    params.M2 = M2;
    params.M3 = M3;
    params.list2 = list2;
    params.list3 = list3;
    params.sigma = sigma;

    %% Precompute biases once and for all
    n2 = size(list2, 1);
    bias2 = zeros(n2, 1);
    for k = 1 : n2
        shift = list2(k, :);
        if all(shift == [0 0])
            bias2(k) = N^2*sigma^2;
        end
    end
    
    n3 = size(list3, 1);
    bias3 = zeros(n3, 1);
    for k = 1 : n3
        shift1 = list3(k, [1, 2]);
        shift2 = list3(k, [3, 4]);
        if all(shift1 == [0, 0])
            bias3(k) = bias3(k) + M1*sigma^2; % this is the 'data' M1 (fixed)
        end
        if all(shift2 == [0, 0])
            bias3(k) = bias3(k) + M1*sigma^2;
        end
        if all(shift1 == shift2)
            bias3(k) = bias3(k) + M1*sigma^2;
        end
    end
    
    params.bias2 = sparse(bias2);
    params.bias3 = sparse(bias3);

    %% Setup Manopt problem

    % Choose wisely if optimize in LxL or WxW
    % It seems that WxW might be necessary to avoid local optima.
   % manifold = euclideanfactory(W, W);
     manifold = euclideanfactory(opt_window, opt_window);
    
%    manifold = euclideanfactory((W+1)/2, (W+1)/2);
    
    problem.M = manifold;
    problem.costgrad = @(X) least_squares_2D_cost_grad_parallel(X, params);

    if ~exist('X0', 'var')
        X0 = [];
    end

    if ~exist('method', 'var')
        method = 'trust_regions';
    end

    
    opts = struct();
    opts.tolgradnorm = 1e-5;
    opts.maxtime = 60*60*10; %10hours
    opts.maxiter = max_iter;

    %opts.tolcost = 1e-2;
    opts.verbosity = 2;
    if strcmp(method,'trust_regions')
    warning('off', 'manopt:getHessian:approx');
    %opts.maxiter = 100;
    [X_est, loss, stats] = trustregions(problem, X0, opts); %#ok<ASGLU>
    warning('on', 'manopt:getHessian:approx');
    elseif strcmp(method,'conjugate_gradient')
    %opts.maxiter = 20;
    [X_est, loss, stats] = conjugategradient(problem, X0, opts); %#ok<ASGLU>
    end        
end