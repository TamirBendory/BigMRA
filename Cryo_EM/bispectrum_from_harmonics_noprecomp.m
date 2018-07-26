function M = bispectrum_from_harmonics_noprecomp(a_lms, L, W)

Lr = floor(L/2);
maxL = length(a_lms)-1;
maxS = size(a_lms{1},1);

[PSWF_Kq_p, PSWF_quad_int_Kq] = precomp_pswf_t(L, 1, 1e-1, 1);

n_list = precomp_Nn_list(Lr);
psi_Nn = PSWF_2D_full_cart(maxL, n_list, Lr);
beta = sph_Bessel_to_2D_PSWF_factors(maxL, n_list, maxS, Lr); % l x N x s x n

b_LmNn = cell(maxL+1);
for l = 0:maxL
    b_LmNn{l+1} = cell(2*l+1, 1);
    for N = 0:l
        b_LmNn{l+1}{N+l+1} = beta{l+1}{N+1}(1:size(a_lms{l+1},1),:).'*a_lms{l+1}; % n x m
    end
    for N = -l:-1
        b_LmNn{l+1}{N+l+1} = conj(b_LmNn{l+1}{abs(N)+l+1});
    end
end

psi_lNm = cell(size(psi_Nn));
for l = 0:maxL
    psi_lNm{l+1} = cell(2*l+1, 1); % N
    for N = -l:l
        psi_lNm{l+1}{N+l+1} = reshape(psi_Nn{N+maxL+1}*b_LmNn{l+1}{N+l+1}, 2*Lr, 2*Lr, []); % r x m
    end
end

blk_id = [];
for ii = 1:max(PSWF_Kq_p.ang_freq)+1
    num_freqs = length(find(PSWF_Kq_p.ang_freq == ii-1));
    blk_id = [blk_id; ii*ones(num_freqs^2, 1)];
end

M = 0;
parfor ii = 1:L^2*(maxL+1)^2
    [row, col, L1, L2] = ind2sub([L, L, maxL+1, maxL+1], ii);
    L1 = L1-1; L2 = L2-1;
    
    for L3 = abs(L1-L2):min(L1+L2, maxL)
        
        for jj = 1:(2*L1+1)^2*(2*L2+1)^2
            [m1, N1, m2, N2] = ind2sub([2*L1+1, 2*L1+1, 2*L2+1, 2*L2+1], jj);
            m1 = m1-L1-1; N1 = N1-L1-1;
            m2 = m2-L2-1; N2 = N2-L2-1;
            
            m3 = m1 + m2; N3 = N1 + N2;
            if abs(m3) > L3 || abs(N3) > L3
                continue;
            end
            
            norm_fact = W{L1+1,L2+1}(L3-abs(L1-L2)+1, m1+L1+1, m2+L2+1)*W{L1+1,L2+1}(L3-abs(L1-L2)+1, N1+L1+1, N2+L2+1);
            
            patch = extract_patch(psi_lNm{L2+1}{N2+L2+1}(:,:,m2+L2+1), row, col, L);
            patch(:,:,2) = extract_patch(psi_lNm{L3+1}{N3+L3+1}(:,:,m3+L3+1), row, col, L);
            
            T  = pswf_t_f_fast(patch, L, 1, 1e-1, 1, PSWF_Kq_p, PSWF_quad_int_Kq, '../nfft');
            T = PSWF_coeff_convert(T, PSWF_Kq_p.ang_freq);
            
            T = cellfun(@(x) psi_lNm{L1+1}{N1+L1+1}(row,col,m1+L1+1)*x(:,1)*x(:,2)', T, 'UniformOutput', 0);
            T = cellfun(@(x)x(:), T, 'UniformOutput', 0);
            
            M = M + norm_fact*vertcat(T{:})./L^2;
        end
    end
end

M = accumarray(blk_id, M, [max(PSWF_Kq_p.ang_freq)+1, 1], @(x) {reshape(x, sqrt(length(x)), [])});