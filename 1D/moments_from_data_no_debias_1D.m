function [M1, M2, M3] = moments_from_data_no_debias_1D(y, list2, list3)
% y is a signal (a vector)
% list2 has size n2 x 1
% list3 has size n3 x 2
% Each row of list2 and list3 contains integers (shifts) -- can be negative

    y = y(:);
    N = size(y, 1);

    M1 = sum(y);
    
    n2 = size(list2, 1);
    M2 = zeros(n2, 1);
    parfor k = 1 : n2
        
        shift1 = list2(k);
        
        vals1 = [0, shift1];
        range1 = (1+max(vals1)) : (N+min(vals1));
        
        y1 = y(range1); %#ok<PFBNS>
        y2 = y(range1-shift1);
        
        M2(k) = sum(y1 .* y2);
        
    end
    
    n3 = size(list3, 1);
    M3 = zeros(n3, 1);
    parfor k = 1 : n3
        
        shifts = list3(k, :);
        shift1 = shifts(1);
        shift2 = shifts(2);

        vals1 = [0, shift1, shift2];
        range1 = (1+max(vals1)) : (N+min(vals1));
        
        y1 = y(range1); %#ok<PFBNS>
        y2 = y(range1-shift1);
        y3 = y(range1-shift2);
        
        M3(k) = sum(y1 .* y2 .* y3);
        
    end
    
end
