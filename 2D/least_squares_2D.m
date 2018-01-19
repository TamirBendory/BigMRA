function [X_est, problem, stats] = least_squares_2D(M1, M2, M3, W, sigma, N, L, m, list2, list3, X0)

    manifold = euclideanfactory(W, W); % choose wisely if do LxL or WxW or something else -- WxW seems useful to avoid local opts.
    
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
% % % 	problem.costgrad = @(X) least_squares_2D_cost_grad(X, params);
	problem.costgrad = @(X) least_squares_2D_cost_grad_new(X, params);
% % %     checkgradient(problem); pause;
    
    
%     keyboard;

    opts = struct();
%     opts.maxtime = 240;

%     [X_est, loss, stats] = conjugategradient(problem, X0, opts); %#ok<ASGLU>

%     [X_est, loss, stats] = barzilaiborwein(problem, X0, opts); %#ok<ASGLU>

    opts.maxiter = 1000;
%     problem.linesearch = @(in1, in2) 2; % optimism in BFGS linesearch -- not sure this is a good idea
    [X_est, loss, stats] = rlbfgs(problem, X0, opts); %#ok<ASGLU>
    warning('off', 'manopt:getHessian:approx');
    opts.maxiter = 200;
    [X_est, loss, stats] = trustregions(problem, X_est, opts); %#ok<ASGLU>
    warning('on', 'manopt:getHessian:approx');

% % %     problem.partialgrad = @(X, sample3) partialgrad(X, params, sample3);
% % %     
% % %     problem.costgrad = @(X) least_squares_2D_cost_grad_new(X, params);
% % %     
% % %     problem.ncostterms = size(list3, 1);
% % %     
% % %     opts.batchsize = size(list2, 1);
% % %     
% % %     opts.maxiter = 100000;
% % %     
% % %     metrics.cost = @(problem, x) getCost(problem, x);
% % %     metrics.gradnorm = @(problem, x) problem.M.norm(x, getGradient(problem, x));
% % %     opts.statsfun = statsfunhelper(metrics);
% % %     
% % %     opts.checkperiod = 1000; % for debugging
% % %     
% % %     if isempty(X0)
% % %         X0 = problem.M.rand();
% % %     end
% % %     opts.stepsize_init = 1/problem.M.norm(X0, getGradient(problem, X0)); % veeery careful here....
% % %     
% % %     [X_est, stats] = stochasticgradient(problem, X0, opts);
% % %     disp(stats)
    
end

function G = partialgrad(X, params, sample3)

    [~, G] = least_squares_2D_cost_grad_new(X, params, sample3);

end

