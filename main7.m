clear; close all; clc;
dbstop if error
L = 4; % length of signal
x = randn(L,1); 

N = 50; 
k = 4;
sigma = 0;


[y ind] = gen_data(x,N,k,L,sigma); 


%% running over indices
max_iter = 1000;
x_est = randn(L,1);
y_s = [y ; y(1:L)];
batch_size = N/2; 
 
for iter = 1:max_iter

    rand_ind = randi(N,[N/2,1]);
   
    for i = 1:batch_size
      
    w(i,iter) = y_s( rand_ind(i):rand_ind(i)+L-1)'*x_est;         
    end
    
w(:,iter) = w(:,iter)/sum(w(:,iter)); %normalization
x_new = zeros(L,1);
for i = 1:batch_size
    x_new = x_new + y_s(rand_ind(i):rand_ind(i)+L-1)*w(i,iter);         
end
x_est = x_new/norm(x_new)*norm(x);
err(iter) = norm(x - x_est)/norm(x);
end

figure; plot(err)
%% 

figure; hold on;  
stem(x); stem(x_est,'*r'); 
legend('signal','recovered');