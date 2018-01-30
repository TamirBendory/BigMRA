function [M1, M2, M3] = moments_from_data_no_debias_2D(Y, list2, list3)
% Y is an image (a matrix)
% list2 has size n2 x 2
% list3 has size n3 x 4
% Each row of list2 and list3 contains integers (shifts)

    [N1, N2] = size(Y);

    M1 = sum(Y(:));
    
    n2 = size(list2, 1);
    M2 = zeros(n2, 1);
    parfor k = 1 : n2
        
        shift1 = list2(k, :);
        
        vals1 = [0, shift1(1)];
        range1 = (1+max(vals1)) : (N1+min(vals1));
        vals2 = [0, shift1(2)];
        range2 = (1+max(vals2)) : (N2+min(vals2));
        
        X1 = Y(range1, range2); %#ok<PFBNS>
        X2 = Y(range1-shift1(1), range2-shift1(2));
        
        M2(k) = sum(X1(:) .* X2(:));
        
    end
    
    n3 = size(list3, 1);
    M3 = zeros(n3, 1);
    parfor k = 1 : n3
        
        shifts = list3(k, :);
        shift1 = shifts([1, 2]);
        shift2 = shifts([3, 4]);

        vals1 = [0, shift1(1), shift2(1)];
        range1 = (1+max(vals1)) : (N1+min(vals1));
        vals2 = [0, shift1(2), shift2(2)];
        range2 = (1+max(vals2)) : (N2+min(vals2));
        
        X1 = Y(range1, range2); %#ok<PFBNS>
        X2 = Y(range1-shift1(1), range2-shift1(2));
        X3 = Y(range1-shift2(1), range2-shift2(2));
        
        M3(k) = sum(X1(:) .* X2(:) .* X3(:));
        
    end
    
end
