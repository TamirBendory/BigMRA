clear; close all; 
%clc;

%% Defining the problem

L = 11;      %length of signal
sigma = 1;   % noise level
W = 2*L;     %window length
Nfactor = 4; % Sparsity factor
overlapping = 1; %maximal overlap
K = 2;  %number signals

% generating non-centered signals
x = zeros(L,K); 
for k=1:K
x(:,k) = randn(L,1) + randn(1);
end

m_vec = 1000; %round(logspace(2,4,10)); % # of signal's repetitions (maximal number)
num_rep = 1;

M1_err = zeros(length(m_vec),num_rep); M2_err = M1_err; M3_err = M1_err; M2nd_err = M1_err;

for mm = 1:length(m_vec)
    m = m_vec(mm);
    N = W*m*Nfactor; % # of measurements
%    display(strcat('m=',num2str(m)));
    for iter = 1:num_rep

        
        
        %% Generating data
        %tic
        [y,yc, ind, class] = gen_data(x,N,m,sigma,W,K);
        %snr = norm(yc)^2/norm(y-yc)^2; % The problem's SNR
        % we assume to know how many signals at each class 
        m_eff = zeros(K,1);
        for k = 1:K
        m_eff(k) = sum(class==k);
        end
        
        %fprintf('Measurement length  = %e \n',N);
        %fprintf('SNR = %.4f \n',snr);
        %fprintf('Generating data time = %.2f [sec] \n',toc);
        
        
        
        %% rearraging the matrix data
        %tic
        y_mat = gen_data_mtx(y,W,overlapping);
        %fprintf('Constructing data matrix time = %.2f [sec] \n',toc);
        
        if isempty(gcp('nocreate'))
            parpool(2, 'IdleTimeout', 240);
        end
        
        %%  Computing and comparing the mean
        
        % computing the empirical invariants of the data
        [M1, M2, M3] = invariants_from_data_no_debias(y_mat);
        
        M2nd = compute_2M(y_mat);
        
        % Expectations by analytical formulas 
        M1_est = sum(x)*m_eff/N;
        M2_est = psx(x,W,sigma,m_eff,N);
        M3_est = bsx(x,W,sigma,m_eff,N,M1);
        M2nd_est = M2ndx(x,W,sigma,m_eff,N);
        
        % computing estimation errors
        M1_err(mm,iter) = abs(M1 - M1_est)/abs(M1); 
        M2_err(mm,iter) = norm(M2 - M2_est)/norm(M2); 
        M3_err(mm,iter) = norm(M3 - M3_est,'fro')/norm(M3,'fro'); 
        M2nd_err(mm,iter) = norm(M2nd_est - M2nd,'fro')/norm(M2nd,'fro'); 
        
    end
end

%% plotting 

% Average errors
M1_err_mean = mean(M1_err,2);
M2_err_mean = mean(M2_err,2);
M3_err_mean = mean(M3_err,2);
M2nd_err_mean = mean(M2nd_err,2);

display([M1_err_mean,M2_err_mean,M3_err_mean,M2nd_err_mean]);

if 1
figure; 
subplot(221); plot(m_vec,M1_err_mean);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
title(strcat('Mean estimation, \sigma = ',num2str(sigma))); 
ylabel('error'); xlabel('m');
%ylim([10^(-3),10^(-1)])
subplot(222); plot(m_vec,M2_err_mean); 
title('Power spectrum estimation'); ylabel('error'); xlabel('m');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log'); 
%ylim([10^(-3),10^(-1)])
subplot(223); plot(m_vec,M3_err_mean); 
title('Bispectrum estimation'); ylabel('error'); xlabel('m');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
%ylim([10^(-3),10^(-1)])

subplot(224); plot(m_vec,M2nd_err_mean); 
title('2nd moment estimation'); ylabel('error'); xlabel('m');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

end