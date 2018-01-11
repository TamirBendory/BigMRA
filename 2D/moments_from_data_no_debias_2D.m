function [M1, M2, M3] = moments_from_data_no_debias_2D(Y, list2, list3)
% Y is an image (a matrix)
% list2 has size n2 x 2
% list3 has size n3 x 4
% Each row of list2 and list3 contains integers (shifts)

    % Initializations for AD
    k1 = 0;
    k2 = 0;
    l1 = 0;
    l2 = 0;
    shiftedY = zeros(size(Y));
    shiftedYY = zeros(size(Y));
    
    M1 = sum(Y(:));
    
    n2 = size(list2, 1);
    M2 = zeros(n2, 1);
    for k = 1 : n2
        k1 = list2(k, 1);
        k2 = list2(k, 2);
        shiftedY = circshift_ad_2D(Y, [k1, k2]);
        M2(k) = Y(:)'*shiftedY(:);
    end
    
    n3 = size(list3, 1);
    M3 = zeros(n3, 1);
    for k = 1 : n3
        k1 = list3(k, 1);
        k2 = list3(k, 2);
        l1 = list3(k, 3);
        l2 = list3(k, 4);
        shiftedYY = circshift_ad_2D(Y, [k1, k2]) .* circshift_ad_2D(Y, [l1, l2]);
        M3(k) = Y(:)'*shiftedYY(:);
    end
    
end
