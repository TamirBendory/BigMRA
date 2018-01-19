function [X_est, problem, stats] = least_squares_2D(M1, M2, M3, W, sigma, N, L, m, list2, list3, X0)

<<<<<<< HEAD
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
=======
    manifold = euclideanfactory(W, W); % choose wisely if do LxL or WxW -- WxW seems useful to avoid local opts.
    
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
>>>>>>> 6bb6424567985d4d069416a830401112a2988a87
    end
    if all(shift1 == shift2)
        bias3(k) = bias3(k) + M1*sigma^2;
    end
end
params.bias2 = bias2;
params.bias3 = bias3;

% % % % %     problem = manoptAD(manifold, @least_squares_cost_2D, params);

problem.M = manifold;
% 	problem.costgrad = @(X) least_squares_2D_cost_grad(X, params);
<<<<<<< HEAD
problem.costgrad = @(X) least_squares_2D_cost_grad_new(X, params);
=======
% 	problem.costgrad = @(X) least_squares_2D_cost_grad_new(X, params);
	problem.costgrad = @(X) least_squares_2D_cost_grad_new_parallel(X, params);
>>>>>>> 6bb6424567985d4d069416a830401112a2988a87
% % %     checkgradient(problem); pause;


%     keyboard;

<<<<<<< HEAD
opts = struct();
=======
    if ~exists('X0', 'var')
        X0 = [];
    end

    opts = struct();
>>>>>>> 6bb6424567985d4d069416a830401112a2988a87
%     opts.maxtime = 240;

%     [X_est, loss] = conjugategradient(problem, X0, opts); %#ok<ASGLU>

%     [X_est, loss] = barzilaiborwein(problem, X0, opts); %#ok<ASGLU>

<<<<<<< HEAD
opts.maxiter = 1000;
%     problem.linesearch = @(in1, in2) 2; % optimism in BFGS linesearch -- not sure this is a good idea
 %opts.tolgradnorm = 1e-7;
[X_est, loss] = rlbfgs(problem, X0, opts); %#ok<ASGLU>
%warning('off', 'manopt:getHessian:approx');
%<<<<<<< HEAD
 opts.tolgradnorm = 1e-5;
%=======
%   opts.maxiter = 1000;
%>>>>>>> 465debbc500803ffe1d3ce4e9add830d9f1f8595
   [X_est, loss] = trustregions(problem, X_est, opts); %#ok<ASGLU>
%  warning('on', 'manopt:getHessian:approx');

=======
%     opts.maxiter = 1000;
%     problem.linesearch = @(in1, in2) 2; % optimism in BFGS linesearch -- not sure this is a good idea
%     [X0, loss] = rlbfgs(problem, X0, opts); %#ok<ASGLU>
    warning('off', 'manopt:getHessian:approx');
    opts.tolgradnorm = 1e-5;
    opts.maxiter = 1000;
    [X_est, loss, stats] = trustregions(problem, X0, opts); %#ok<ASGLU>
    warning('on', 'manopt:getHessian:approx');
    
>>>>>>> 6bb6424567985d4d069416a830401112a2988a87
end
