function [M1, M2, M3] = moments_from_signal_with_bias_2D(X, W, N, m, sigma, list2, list3)
% Inputs:
%   X is an image of size LxL
%   W is the maximum shift used in computing correlations
%   N is the size of an observed image (NxN)
%   m is the number of times X is supposed to appear in an observation
%   sigma is the noise standard deviation (iid standard Gaussian, white noise)
    
    % Zero-pad X
    L = size(X, 1);
    assert(size(X, 2) == L, 'X must be square.');
    
%     X_zp = [X zeros(L, W-L); zeros(W-L, W)]; % triple check how much we need to shift here
    X_zp = [X zeros(L, W); zeros(W, L+W)]; % triple check how much we need to shift here
    
    % Get the moments as if the zero-padded X were an observation
    [M1, M2, M3] = moments_from_data_no_debias_2D(X_zp, list2, list3);
    
    % Correct scaling and bias to match what observed M{123} should look like.
    M1 = M1 * m;
    
    M2 = M2 * m;
    
    % This could be precomputed
    bias2 = zeros(size(M2));
    for k = 1 : size(M2, 1)
        if list2(k, 1) == 0 && list2(k, 2) == 0
            bias2(k) = N^2*sigma^2;
        end
    end
%     bias2 = sparse(bias2); % for the precomputed version later
    M2 = M2 + bias2;
    
    
    M3 = M3 * m;
    bias3 = zeros(size(M3));
    for k = 1 : size(M3, 1)
        if list3(k, 1) == 0 && list3(k, 2) == 0
            bias3(k) = bias3(k) + M1*sigma^2; % careful that M1 was already scaled
        end
        if list3(k, 3) == 0 && list3(k, 4) == 0
            bias3(k) = bias3(k) + M1*sigma^2; % careful that M1 was already scaled
        end
        if list3(k, 1) == list3(k, 3) && list3(k, 2) == list3(k, 4)
            bias3(k) = bias3(k) + M1*sigma^2; % careful that M1 was already scaled
        end
    end
    M3 = M3 + bias3;
    
end
