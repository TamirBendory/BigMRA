function [A2, A3] = AC_est(y,W)

A2 = zeros(W,1);
A3 = zeros(2*W+1);

% computing second moment
for i = 0:W-1
    A2(i+1) = sum(y.*circshift(y,i));
end

% computing third moment

for i = -(W-1):W-1
    for j = -(W-1):W-1
        A3(i+W,j+W) = sum(y.*circshift(y,i).*circshift(y,-j));
    end
end

