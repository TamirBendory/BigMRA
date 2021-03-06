function x_est = bigmra_lsq(M1, M2, M3, W, sigma, N, L, m)
% Big-MRA via least-squares solve.
    
    m = m(:);
    K = length(m);

    manifold = euclideanfactory(L, K);
    
    % Needed to add this in reverse mode (it solved the problem I had then)
    manifold.egrad2rgrad = @(x, g) real(g);
    
    params.M1 = M1;
    params.M2 = M2;
    params.M3 = M3;
    params.W = W;
    params.sigma = sigma;
    params.N = N;
    params.m = m;
    
    problem = manoptAD(manifold, @bigmra_lsq_cost, params);
    
    warning('off', 'manopt:getHessian:approx');
    
    [x_est, loss] = trustregions(problem);
    
    % Cyclic shifts of the true signal can create local minima, hard to
    % escape. Scan through all cyclic shifts of the computed estimator,
    % identify the most promising one, and re-optimize. This could be
    % iterated. Need to think how to do that for K > 1 at reasonable cost.
    %
    % Actually, this heuristic doesn't work too well: need to consider more
    % than just the most promising one (which tends to be the just-computed
    % estimated.) In the revised version, we pick the most promising shift,
    % discarding the non-shift, and pick the best of both obtained signals.
    %
    % So, that also didn't work. Seems like we need to reoptimize from all
    % possible shifts...
    if K == 1
        x_base = x_est;
        for shift = 1 : L-1
            x_est_shifted = circshift(x_base, shift);
            [x_est_shifted, loss_shifted] = trustregions(problem, x_est_shifted);
            if loss_shifted < loss
                loss = loss_shifted;
                x_est = x_est_shifted;
            end
        end
    end
    
    warning('on', 'manopt:getHessian:approx');
    
end
