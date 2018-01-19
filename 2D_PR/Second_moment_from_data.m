function M2 = Second_moment_from_data(Y,W,sigma)
% Y is an image (a matrix)

[N1, N2] = size(Y);
M2 = zeros(W);

for shift1 = -(W-1)/2:(W-1)/2
    for shift2 = 0:(W-1)/2        
        vals1 = [0, shift1];
        range1 = (1+max(vals1)) : (N1+min(vals1));
        vals2 = [0, shift2];
        range2 = (1+max(vals2)) : (N2+min(vals2));        
        X1 = Y(range1, range2);
        X2 = Y(range1-shift1, range2-shift2);
        M2(shift1+(W+1)/2,shift2+(W+1)/2) = sum(X1(:) .* X2(:));
    end
end

% construct the entire 2nd moment with symmetries

M2 = circshift(M2,[-(W-1)/2,-(W-1)/2]);
M2((W+1)/2+1:end,(W+1)/2+1:end) = rot90(M2(2:(W+1)/2,2:(W+1)/2),2);
M2(2:(W+1)/2,(W+1)/2+1:end) = rot90(M2((W+1)/2+1:end,2:(W+1)/2),2);
M2(1,(W+1)/2+1:end) = fliplr(M2(1,2:(W+1)/2));
% unbiasing
M2(1,1) = M2(1,1) - sigma^2*N1*N2;
end

