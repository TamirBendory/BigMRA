function problem = manoptAD(M, costfun, param) %, parameter_positions, parameter_values)
% Assuming f = costfun(param, x).

    % Follow tips to make it faster here:
    % https://adimat.sc.informatik.tu-darmstadt.de/doc/adimat-16.html

    problem = struct();
    
    problem.M = M;
    
    problem.cost = @(x) costfun(param, x);
    
    problem.costgrad = @costgrad;

    % Could also define directional derivative for cheaper. think this
    % needs forward mode.
    
    % Might also want to generate the gradient function once, then place a
    % handle to that function in the problem structure directly to avoid
    % ADiMat's call overhead.
        
    ad_opts = admOptions(...
        'independents', 2, ...   % positions of inputs of f that are variables
        'functionResults', {1}); % any matrix of the same size as the output of f, in a cell (allows multiple outputs)

    
    function [f, g] = costgrad(x)

        [g, f] = admDiffRev(costfun, 1, param, x, ad_opts);
        g = reshape(g, size(x));
        g = M.egrad2rgrad(x, g);
        
    end

end
