clear; close all;
%clc;
dbstop if error
%% Defining the problem

L = 11;      %length of signal
sigma = 0;   % noise level
W = 21;     %window length
Nfactor = 4; % Sparsity factor
overlapping = 1; %maximal overlap
K = 1;  %number signals

% generating non-centered signals
x = zeros(L,K);
for k=1:K
    x(:,k) = randn(L,1) + randn(1);
end

m_vec = 100; %round(logspace(2,4,10)); % # of signal's repetitions (maximal number)
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

        m_eff = zeros(K,1);
        for k = 1:K
            m_eff(k) = sum(class==k);
        end
        
        
        %% rearraging the matrix data
        %tic
        y_mat = gen_data_mtx(y,W,overlapping);
        
        if isempty(gcp('nocreate'))
            parpool(2, 'IdleTimeout', 240);
        end
        
        %%  Computing and comparing the mean
        
        %[M1, M2, M3] = invariants_from_data_no_debias(y_mat);
         ac_flag = 1;
         M2nd = compute_2M(y_mat,0);      
         A2_data = compute_A2(y_mat);
         Q = (linspace(ceil(L/2),L,L))';
         Q = [flipud(Q) ; ones(W-L,1)];
         A2_data = A2_data(1:W);
         A2_data = A2_data/N./Q/W*L;
 
         [A3_data, A3_comp] = compute_A3(y_mat,L);
        M3nd = compute_3M(y_mat,0);
        count = check_dis_number(M3nd);
        %S1 = squeeze(M3nd(11,:,:));
       
        %figure; subplot(211); imagesc(A3_data); colorbar
        %subplot(212); imagesc(S1); colorbar
                
        M2nd_est = M2ndx(x,W,sigma,m_eff,N,0);
        A2_est = A2x(x,W,sigma,m_eff,N);
        
               
    end
end

%% plotting

% Average errors
M1_err_mean = mean(M1_err,2);
M2_err_mean = mean(M2_err,2);
M3_err_mean = mean(M3_err,2);
M2nd_err_mean = mean(M2nd_err,2);

%display([M1_err_mean,M2_err_mean,M3_err_mean,M2nd_err_mean]);

if 0
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% computing estimation errors
      %  M1_err(mm,iter) = abs(M1 - M1_est)/abs(M1);
       % M2_err(mm,iter) = norm(M2 - M2_est)/norm(M2);
       % M3_err(mm,iter) = norm(M3 - M3_est,'fro')/norm(M3,'fro');
        %M2nd_err(mm,iter) = norm(M2nd_est - M2nd,'fro')/norm(M2nd,'fro');
     
        
            
        % Expectations by analytical formulas
        %M1_est = sum(x)*m_eff/N;
        %M2_est = psx(x,W,sigma,m_eff,N);
        %M3_est = bsx(x,W,sigma,m_eff,N,M1);
      %   M2nd_est = M2ndx(x,W,sigma,m_eff,N,0);
        %A2_est = A2x(x,W,sigma,m_eff,N);
     