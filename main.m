clear; close all; clc;

%% Defining the problem

L = 11;      %length of signal
sigma = 1;   % noise level
W = 21;     %window length (should be an odd number)
Nfactor = 4; % Sparsity factor
K = 1;  %number signals

% generating non-centered signals
%x = zeros(L,K);
% for k=1:K
%     x(:,k) =  randn(L,1);
% end
x = randn(L,1); %x = x - mean(x) + 0.1;
Lm = 5;
m_vec = round(logspace(1,5,Lm));
num_rep = 100;

err1 = zeros(Lm,num_rep);
err2 = err1;
err3 = err1;

for i = 1:Lm
m = m_vec(i)
N = W*m*Nfactor; % # of measurements

    for j = 1:num_rep
j

%% Generating data

[y,yc, ind, class] = gen_data(x,N,m,sigma,W,K);

m_eff = zeros(K,1);
for k = 1:K
    m_eff(k) = sum(class==k);
end

if isempty(gcp('nocreate'))
    parpool(2, 'IdleTimeout', 240);
end

%%  Computing and comparing auto-correlations

% estimating the mean 
A1_est = sum(y)/m_eff; 
AC1x = sum(x);

[A2_est,A3_est] = AC_est(y,W);  
A2_est(1) = A2_est(1) - N*sigma^2; %unbiasing
A2_est = A2_est/m_eff; 
A3_est = A3_est(W+1-L:W+L-1,W+1-L:W+L-1)/m_eff;
Bmat = eye(W); Bmat(1,:) = Bmat(1,:)+1; Bmat(:,1) = Bmat(:,1)+1;
Bmat = circshift(Bmat,[(W-1)/2,(W-1)/2]); % assuming W is odd 
A3_est = A3_est - AC1x/L*(sigma^2)*Bmat; % using the sum of the true signal

AC2x = xcorr([x; zeros(W-L,1)]); AC2x  = AC2x(W:end); 
AC3x = comp_A3x([x; zeros(W-L,1)],W); 
AC3x = AC3x(W+1-L:W+L-1,W+1-L:W+L-1);

err1(i,j) = norm(AC1x - A1_est)/norm(AC1x);
err2(i,j) = norm(AC2x - A2_est)/norm(AC2x);
err3(i,j) = norm(AC3x - A3_est,'fro')/norm(AC3x,'fro');

save('err1','err1');
save('err2','err2');
save('err3','err3');

%figure; imagesc(A3_est - AC3x); colorbar

    end
end

%% Plotting

if 1
    
figure; 
subplot(311); plot(m_vec(1:7),mean(err1(1:7,:),2));
set(gca,'XScale','log'); set(gca,'YScale','log');
title('Mean estimation'); xlabel('m'); ylabel('error'); axis tight
subplot(312); plot(m_vec(1:7),mean(err2(1:7,:),2));
set(gca,'XScale','log'); set(gca,'YScale','log');
title('Second moment estimation'); xlabel('m'); ylabel('error');  axis tight
subplot(313); plot(m_vec(1:7),mean(err3(1:7,:),2));
set(gca,'XScale','log'); set(gca,'YScale','log');
title('Third moment estimation'); xlabel('m'); ylabel('error'); axis tight

end