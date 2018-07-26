function [B, L1_list, L2_list, L3_list, blk_id] = precomp_bispectrum_coeffs_LS(maxL, s_list, L, W, psi_lNs, Psi, ang_freq, pts_in_disc)

% clear psi_Nn beta n_list

q_list = zeros(max(ang_freq)+1, 1);
blk_id = [];
for ii = 1:length(q_list)
    q_list(ii) = sum(ang_freq == ii-1);
    blk_id(end+1:end+q_list(ii)^2, 1) = ii;
end
kq_len = length(blk_id);
q_cell = num2cell(q_list);

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

for ll = 1:length(L1_list)
    L1 = L1_list(ll);
    L2 = L2_list(ll);
    ii_3 = L3_list(ll); L3 = abs(L1-L2) + ii_3 - 1;
    
    s1_len = s_list(L1+1);
    s2_len = s_list(L2+1);
    s3_len = s_list(L3+1);
    
    acc_vec = zeros(kq_len*s3_len*s2_len, s1_len,'double');
    psi_lNs1 = cellfun(@(x) reshape(x, L, L, []), psi_lNs{L1+1}, 'UniformOutput', 0);
    psi_lNs2 = cellfun(@(x) reshape(x, L, L, []), psi_lNs{L2+1}, 'UniformOutput', 0);
    psi_lNs3 = cellfun(@(x) reshape(x, L, L, []), psi_lNs{L3+1}, 'UniformOutput', 0);
    
    parfor ii = 1:L^2*(2*L1+1)*(2*L2+1)
        [row, col, N1, N2] = ind2sub([L, L, 2*L1+1, 2*L2+1], ii);
        N1 = N1-L1-1; N2 = N2-L2-1;
        N3 = N1 + N2;
        if abs(N3) > L3, continue; end
        
        norm_fact = W{L1+1, L2+1}(N2+L2+1, N1+L1+1, ii_3)/L^2;
        
        patch = extract_patch(psi_lNs2{N2+L2+1}, row, col, L-1);
        patch(:,:,end+1:end+s3_len) = extract_patch(psi_lNs3{N3+L3+1}, row, col, L-1);
        patch = reshape(patch(1:end-1,1:end-1,:), 4*(L-1)^2, size(patch,3));
        
        x = Psi\patch(pts_in_disc, :);
        T = x(1:length(ang_freq),:);
        T(ang_freq>0, :) = T(ang_freq>0, :) + 1i*x(length(ang_freq)+1:end, :);
        T = mat2cell(T, q_list, s2_len + s3_len);
        
        T = cellfun(@(x) vec(x(:,1:s2_len))*vec(x(:,s2_len+1:end))', T, 'UniformOutput', 0); % (q1, s2) x (q2, s3)
        
        T = cellfun(@(x,y) reshape(permute(reshape(x, y, s2_len, y, s3_len), [1,3,4,2]), y^2, []), T, q_cell, 'UniformOutput', 0); %(q1, q2, s3, s2, s1)
        
        acc_vec = acc_vec + norm_fact*vec(vertcat(T{:}))*reshape(psi_lNs1{N1+L1+1}(row, col, :), 1, s1_len); % (q1, q2, s3, s2) x s1
    end
    B{ll} = reshape(acc_vec, kq_len, s3_len, s2_len, s1_len); % (q1, q2) x s3 x s2 x s1
end

