clear; close all; clc;

L = 4; % length of signal
x = randn(L,1); x = x - mean(x);

N = 100;
k = 1;
sigma = 0;

y = gen_data(x,N,k,L,sigma); 

%% rearraging the matrix data
y_matrix = zeros(N,L);
y_stretch = [y ; y(1:L-1)];
for i = 1:N
    y_matrix(i,:) = y_stretch(i:i+L-1);
end
clear y_stretch;

%%

%s_est = zeros(N,1); %sind = randi([L,N-L+1],k,1);
%s_est(sind) = 1; %1/k; %clear sind;
s_est = rand(N,1);
s_est = s_est/sum(s_est);
s_est(1:L-1) = 0; s_est(N-L+1:N) = 0;


num_iter = 50;
x_est = zeros(4,1); err = zeros(num_iter,1); 
tau = 0.5; 

for iter = 1:num_iter
    
    for i = 1:N
        x_est = x_est + y_matrix(i,:)'*circshift(s_est(i),2);
    end
    
    if norm(x_est)~=0
    x_est = x_est/norm(x_est)*norm(x);
    end
    
    %re-estimating s
    
    for i = 1:N
        s_est(i) = y_matrix(i,:)*x_est;
    end
   % s_est = circshift(s_est,1);
    
    for i = L:N-L-1
        if abs(s_est(i)) == max(s_est([i-L+1:i+L-1]))
            continue;
        else
            s_est(i) = s_est(i)*tau;
        end
        
    end
    s_est = s_est/sum(s_est);
    
    err(iter) = norm(x_est- x)/norm(x);
    
end

figure; plot(err);
%figure; hold on; stem(x); stem(x_est,'.r'); xlim([1,L])
%figure; hold on; stem(s_est); stem(y,'*r');
