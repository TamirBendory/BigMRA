function [B, L1_list, L2_list, L3_list, blk_id] = precomp_bispectrum_coeffs_windows(maxL, s_list, L, W)

beta_PSWF = 1; Trunc = 1e-1;
[Wt, ang_freq] = precomp_pswf_t_windows(L-1, beta_PSWF, Trunc);
Wt = fft2(Wt, 3*L-2, 3*L-2);
Wt = permute(Wt, [1,2,4,3]);

[psi_Nn, n_list] = PSWF_2D_full_cart(maxL, floor(L/2), beta_PSWF, Trunc);
beta = sph_Bessel_to_2D_PSWF_factors(maxL, n_list, s_list(1), floor(L/2)); % l x N x s x n

psi_lNs = cell(maxL+1,1);
for l = 0:maxL
    psi_lNs{l+1} = cell(2*l+1, 1); % N
    for N = -l:l
        psi_lNs{l+1}{N+l+1} = (-1)^(N*(N<0))*reshape(psi_Nn{N+maxL+1}*beta{l+1}{abs(N)+1}(1:s_list(l+1),:).', L, L, []); % r x s
    end
end

% clear psi_Nn beta n_list
maxN = max(ang_freq);
q_list = zeros(maxN+1, 1);
blk_id = [];
for ii = 1:length(q_list)
    q_list(ii) = sum(ang_freq == ii-1);
    blk_id(end+1:end+q_list(ii)^2, 1) = ii;
end
q_cumsum = cumsum([0; q_list(:)]);
q_sq_cumsum = cumsum([0; q_list(:).^2]);
kq_len = length(blk_id);

L1_list = []; L2_list = []; L3_list = [];
for ii = 1:(maxL+1)^2
    [L1, L2] = ind2sub([maxL+1, maxL+1], ii);
    L1 = L1-1; L2 = L2-1;
    
    L3_vals = abs(L1-L2):min(L1+L2, maxL);
    for jj = 1:length(L3_vals)
        L1_list(end+1) = L1;
        L2_list(end+1) = L2;
        L3_list(end+1) = jj;
    end
end

vec = @(x) x(:);
B = cell(length(L1_list), 1);
num_freqs = length(ang_freq);

parfor ll = 1:length(L1_list)
    L1 = L1_list(ll);
    L2 = L2_list(ll);
    ii_3 = L3_list(ll); L3 = abs(L1-L2) + ii_3 - 1;
    
    s1_len = s_list(L1+1);
    s2_len = s_list(L2+1);
    s3_len = s_list(L3+1);
    
    acc_vec = zeros(kq_len*s3_len*s2_len, s1_len, 'double');
    
    for N1 = -L1:L1
        for N2 = max(-L2, -L3-N1):min(L2, L3-N1)
            N3 = N1+N2;
            norm_fact = W{L1+1, L2+1}(N2+L2+1, N1+L1+1, ii_3)/L^2;
            
            I1 = fft2(psi_lNs{L2+1}{N2+L2+1}, 3*L-2, 3*L-2);
            I2 = fft2(psi_lNs{L3+1}{N3+L3+1}, 3*L-2, 3*L-2);
            
            T = zeros(3*L-2, 3*L-2, s2_len + s3_len, num_freqs, 'double');
            T(:,:,1:s2_len,:) = bsxfun(@times, I1, Wt); 
            T(:, :, s2_len+1:end, :) = bsxfun(@times, I2, Wt);
            
            T = ifft2(T);
            T = T(L:end-L+1, L:end-L+1, :, :); % row x col x (s2 + s3) x (k,q)
            T = reshape(T, L^2, s2_len + s3_len, num_freqs);
            
            for col = 1:L
                T_acc = zeros(kq_len*s3_len*s2_len, L, 'double');
                P1 = reshape(psi_lNs{L1+1}{N1+L1+1}(:, col, :), L, s1_len);
                for row = 1:L
                    T_add = zeros(kq_len,s3_len*s2_len,'double');
                    idx = row + (col-1)*L;
                    for N = 0:maxN
                        T_N = T(idx, :, q_cumsum(N+1)+1: q_cumsum(N+2));
                        T_N = reshape(T_N, s2_len+s3_len ,q_list(N+1));
                        T_N = vec(T_N(1:s2_len,:))*vec(T_N(s2_len+1:end,:))'; % (s2, q1) x (s3, q2)
                        T_N = reshape(T_N, s2_len, q_list(N+1), s3_len, q_list(N+1));
                        T_N = permute(T_N, [2,4,3,1]); % (q1, q2, s3, s2)
                        T_add(q_sq_cumsum(N+1)+1: q_sq_cumsum(N+2), :) = ...
                                                        reshape(T_N, q_list(N+1)^2, s3_len*s2_len);
                    end
                    T_acc(:, row) = T_add(:);
                end
                acc_vec = acc_vec + norm_fact*T_acc*P1; % (q1, q2, s3, s2) x s1
            end
        end
    end
    B{ll} = reshape(acc_vec, kq_len, s3_len, s2_len, s1_len); % (q1, q2) x s3 x s2 x s1
end

