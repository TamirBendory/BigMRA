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
Lm = 4;
m_vec = round(logspace(1,4,Lm));
num_rep = 100;

err1 = zeros(Lm,num_rep); err2 = err1; err3 = err1;

for i = 1:Lm
    m = m_vec(i)
    N = W*m*Nfactor; % # of measurements
    
    for j = 1:num_rep
        
        %% Generating data
        
        [y,yc, ind, class] = gen_data(x,N,m,sigma,W);
        
        m_eff = zeros(K,1);
        for k = 1:K
            m_eff(k) = sum(class==k);
        end
        
        if isempty(gcp('nocreate'))
            parpool(2, 'IdleTimeout', 240);
        end
        
        %%  Computing and comparing auto-correlations
        
        %estimating the AC functions from the data
        [A1_est, A2_est,A3_est] = AC_est(y,W,m_eff,sigma);
        
        % computing the signal's AC functions
        AC1x = sum(x);
        AC2x = xcorr([x; zeros(W,1)]);
        AC2x  = AC2x(W+L:2*W+L-1);
        AC3x = comp_A3x([x; zeros(W,1)],W);        
        
        err1(i,j) = norm(AC1x - A1_est)/norm(AC1x);
        err2(i,j) = norm(AC2x - A2_est)/norm(AC2x);
        err3(i,j) = norm(AC3x - A3_est,'fro')/norm(AC3x,'fro');
        
        
        %
        % save('err1','err1');
        % save('err2','err2');
        % save('err3','err3');
        
        %figure; imagesc(A3_est - AC3x); colorbar
        
    end
end

%% Plotting

if 1
    
    figure;
    subplot(311); plot(m_vec,mean(err1,2));
    set(gca,'XScale','log'); set(gca,'YScale','log');
    title('Mean estimation'); xlabel('m'); ylabel('error'); axis tight
    subplot(312); plot(m_vec,mean(err2,2));
    set(gca,'XScale','log'); set(gca,'YScale','log');
    title('Second moment estimation'); xlabel('m'); ylabel('error');  axis tight
    subplot(313); plot(m_vec,mean(err3,2));
    set(gca,'XScale','log'); set(gca,'YScale','log');
    title('Third moment estimation'); xlabel('m'); ylabel('error'); axis tight
    
end