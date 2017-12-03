function [A1, A2, A3] = AC_est(y,W,L,m,sigma)

A2 = zeros(W,1);
A3 = zeros(2*W+1);


A1 = sum(y)/m; 

% computing second moment
for i = 0:W-1
        A2(i+1) = sum(y.*circshift(y,i));
end

A2(1) = A2(1) - length(y)*sigma^2; %debiasing
A2 = A2/m; 


% computing third moment

for i = -(W-1):W-1
    for j = -(W-1):W-1
        A3(i+W,j+W) = sum(y.*circshift(y,i).*circshift(y,j));
    end
end

A3 = A3(W+1-L:W+L-1,W+1-L:W+L-1)/m;
Bmat = eye(2*L-1);
Bmat(1,:) = Bmat(1,:)+1; Bmat(:,1) = Bmat(:,1)+1;
Bmat = circshift(Bmat,[L-1,L-1]); % assuming W is odd 
A3 = A3 - A1*(sigma^2)*Bmat; % using the sum of the true signal

