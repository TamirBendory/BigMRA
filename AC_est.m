function [A1, A2, A3] = AC_est(y,W,m,sigma)

A2 = zeros(W,1);
A3 = zeros(W);

A1 = sum(y)/m; 

% computing second moment
for i = 0:W-1
        A2(i+1) = sum(y.*circshift(y,i));
end

A2(1) = A2(1) - length(y)*sigma^2; %debiasing
A2 = A2/m; 


% computing third moment

for i = 0:W-1
    for j = 0:W-1
        A3(i+1,j+1) = sum(y.*circshift(y,i).*circshift(y,j));
    end
end

A3 = A3/m;

%debiasing
Bmat = eye(W);
Bmat(1,:) = Bmat(1,:)+1; Bmat(:,1) = Bmat(:,1)+1;
A3 = A3 - A1*(sigma^2)*Bmat; % using the sum of the true signal

