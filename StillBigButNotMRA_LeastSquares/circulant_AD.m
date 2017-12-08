function C = circulant_AD(v)
% Given a vector v, returns a circulant matrix C whose first row is v.
%
% Produces same output as gallery('circul', v) but faster because it avoids
% a succession of checks. This code is adapted from Matlab's circul.m.
%
% May 2017
% https://arxiv.org/abs/1705.00641
% https://github.com/NicolasBoumal/MRA

% % % % %     % Make v into a row vector
% % % % %     v = v(:).';
% % % % %     n = length(v);
% % % % % 
% % % % %     if n == 1
% % % % %         C = v;
% % % % %     else
% % % % %         C = toeplitz([v(1), v(n:-1:2)], v);
% % % % %     end

    % Changed for AD because Toeplitz not implemented in ADiMat
    
    % Make v into a row vector
    v = v(:).';
    n = length(v);
    C = zeros(n, n);
    for k = 1 : n
        C(k, [k:n, 1:(k-1)]) = v;
    end
    
end

