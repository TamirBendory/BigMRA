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
    
    % Precompute biases once and for all
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
    params.bias2 = bias2;
    params.bias3 = bias3;
    
% % % % %     problem = manoptAD(manifold, @least_squares_cost_2D, params);

    problem.M = manifold;
	problem.costgrad = @(X) least_squares_2D_cost_grad(X, params);
% % %     checkgradient(problem); pause;
    
    warning('off', 'manopt:getHessian:approx');
    
%     keyboard;

    opts = struct();
%     opts.maxtime = 240;
    [X_est, loss] = trustregions(problem, X0, opts); %#ok<ASGLU>

%     [X_est, loss] = conjugategradient(problem, X0, opts); %#ok<ASGLU>

    warning('on', 'manopt:getHessian:approx');
    
end
