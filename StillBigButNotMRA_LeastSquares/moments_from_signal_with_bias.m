function [M1, M2, M3] = moments_from_signal_with_bias(x, W, N, m, sigma)
% Inputs:
%   x is a signal of length L
%   W is the window length (maximum shift used in computing correlations)
%   N is the length of a long observation y
%   m is the number of times x is supposed to appear in a long observation y
%   sigma is the noise standard deviation (iid standard Gaussian, white noise)
    
    % Zero-pad x
    x_zp = [x(:) ; zeros(W, 1)];
    
    % Get the moments as if the zero-padded x were a long observation
    [M1, M2, M3] = moments_from_data_no_debias(x_zp, W);
    
    % Correct scaling and bias to match what observed M{123} should look like.
    M1 = M1 * m;
    
    M2 = M2 * m;
    M2(1) = M2(1) + N*sigma^2;
    
    M3 = M3 * m;
    A = eye(W);
    A(:, 1) = A(:, 1) + 1;
    A(1, :) = A(1, :) + 1;
    M3 = M3 + (M1*sigma^2)*A; % careful that M1 was already scaled
    
end
