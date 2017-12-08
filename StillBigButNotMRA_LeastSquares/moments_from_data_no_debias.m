function [M1, M2, M3] = moments_from_data_no_debias(y, W)
% Inputs:
%     y is the long observation
%     W is the window length (maximum shift used in computing correlations)

    % First order moment
    M1 = sum(y);

    % Second order moment
    M2 = zeros(W, 1);
    for k = 0:(W-1)
        M2(k+1) = sum(y.*circshift_ad(y, k));
    end

    % Third order moment
    M3 = zeros(W, W);
    for k = 0:(W-1)
        for j = 0:(W-1)
            M3(k+1, j+1) = sum(y.*circshift_ad(y, k).*circshift_ad(y, j));
        end
    end

end
