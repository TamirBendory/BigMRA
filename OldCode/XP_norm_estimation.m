clear; close all; clc;

%% Defining the problem

L = 11; % length of signal
window_size = 4*L;
Nfactor = 6; % Sparsity factor, should be ~6
overlapping_factor = 1; % windows are overlapped by window_size/overlapping_factor

Nfactor_vector = 6; % Sparsity factor, should be ~6
k_vector = logspace(2,6,20);  % # of signal's repetitions (maximal number)
num_iter = 100;
sigma_vec = [0.5 , 1, 2];  % noise level2

Err_normX = zeros(length(k_vector),num_iter,3);



for s = 1:length(sigma_vec)
    
    sigma = sigma_vec(s)
    
    for i = 1:length(k_vector)
        k = k_vector(i)
        N = window_size*k*Nfactor; % # of measurements
        
        for j = 1:num_iter
            
            %% Generating data
            %tic
            x = randn(L,1);
            [y,yc, ind] = gen_data(x,N,k,sigma,window_size);
            snr = norm(yc)^2/norm(y-yc)^2; % The problem's SNR
            k_eff = length(ind); % The actual nunber of signal's repetitions
            
            %estimating the signal's norm from the data
            normX = sqrt((norm(y)^2 - sigma^2*N)/k_eff);
            Err_normX (i,j,s) = abs(norm(x) - normX)/norm(x);
            %fprintf('Error of norm''s estimation = %.4f  \n',Err_normX);
            
        end
    end
    save('Err_normX','Err_normX');
end

    %%
    ln = 1.5; 
    figure; hold on;
    loglog(k_vector,mean(squeeze(Err_normX(:,:,1)),2),'linewidth',ln)
    loglog(k_vector,mean(squeeze(Err_normX(:,:,2)),2),'linewidth',ln)
    loglog(k_vector,mean(squeeze(Err_normX(:,:,3)),2),'linewidth',ln)
    legend('\sigma=0.1','\sigma=1','\sigma=3');
    ylabel('error'); xlabel('K');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    axis tight
    pdf_print_code(gcf,'NormError', 11)