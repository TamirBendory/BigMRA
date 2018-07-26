function M = bispectrum_from_harmonics_par(a_lms, L, B, W, L1_list, L2_list, L3_list, s_lens)

maxL = length(a_lms)-1;
maxS = size(a_lms{1},1);

if ~exist('W', 'var') || isempty(W)
    W = precomp_wigner_weights(maxL);
end

if ~exist('B','var') || isempty(B)
    B = precomp_bispectrum_coeffs(maxL, maxS, L, W);
    [B, L1_list, L2_list, L3_list] = form_B_cell(B);
end

maxK = length(B)-1;
M = cell(maxK+1, 1);
M = cellfun(@(x) 0, M, 'UniformOutput', 0);
for k = 0:maxK
    M_acc = 0;
    B_k = B{k+1};
    q_k = size(B_k{1},1);
    parfor ii = 1:length(B_k)
        L1 = L1_list(ii); s1_curr = s_lens(L1+1); 
        L2 = L2_list(ii); s2_curr = s_lens(L2+1);
        jj_3 = L3_list(ii); L3 = abs(L1-L2) + jj_3 - 1; s3_curr = s_lens(L3+1);
        
        % Perform summation over s1, s2, s3:
        tmp = reshape(B_k{ii}(:,:,1:s1_curr,1:s2_curr,1:s3_curr), q_k^2*s1_curr*s2_curr, s3_curr)*conj(a_lms{L3+1}); % (q1, q2, s1, s2) x m3
        tmp = reshape(tmp.', [2*L3+1, q_k, q_k, s1_curr, s2_curr]); % m3 x q1 x q2 x s1 x s2
        
        tmp = reshape(tmp, (2*L3+1)*q_k^2*s1_curr, s2_curr)*a_lms{L2+1}; % m3 x q1 x q2 x s1 x m2
        tmp = reshape(tmp.', [2*L2+1, 2*L3+1, q_k, q_k, s1_curr]); % m2 x m3 x q1 x q2 x s1
        
        tmp = reshape(tmp, (2*L2+1)*(2*L3+1)*q_k^2, s1_curr)*a_lms{L1+1}; % (m2, m3, q1, q2) x m1
        tmp = reshape(tmp.', [2*L1+1, 2*L2+1, 2*L3+1, q_k, q_k]); % m1 x m2 x m3 x q1 x q2
        
        % Perform summation over m2, m3:
        for m2 = -L2:L2
            for m1 = max(-L1, -L3-m2):min(L1, L3-m2)
                m3 = m1+m2;
                M_acc = M_acc + W{L1+1, L2+1}(jj_3, m1+L1+1, m2+L2+1)*reshape(tmp(m1+L1+1, m2+L2+1, m3+L3+1,:,:), q_k, q_k); % q1 x q2
            end
        end
    end
    M{k+1} = real(M_acc + M_acc')/2;
end
