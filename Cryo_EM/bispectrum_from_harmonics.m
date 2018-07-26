function M = bispectrum_from_harmonics(a_lms, L, B, W)

maxL = length(a_lms)-1;
maxS = size(a_lms{1},1);

if ~exist('W', 'var') || isempty(W)
    W = precomp_wigner_weights(maxL);
end

if ~exist('B','var') || isempty(B)
    B = precomp_bispectrum_coeffs(maxL, maxS, L, W);
end

M = cell(length(B{1}{1}), 1);
M = cellfun(@(x) 0, M, 'UniformOutput', 0);
for ii = 1:(maxL+1)^2
    [L1, L2] = ind2sub([maxL+1, maxL+1], ii);
    L1 = L1-1; L2 = L2-1;
    
    L3_vals = abs(L1-L2):min(L1+L2, maxL);
    for jj = 1:length(L3_vals)
        L3 = L3_vals(jj);
        
        for k = 0:length(B{ii}{jj})-1
            sz_B = size(B{ii}{jj}{k+1}); % q1 x q2 x s1 x s2 x s3
            % Perform summation over s1, s2, s3:
            sz = size(a_lms{L3+1});
            tmp = reshape(B{ii}{jj}{k+1}(:,:,:,:,1:sz(1)), [], sz(1))*conj(a_lms{L3+1}); % (q1, q2, s1, s2) x m3
            tmp = reshape(tmp.', [sz(2), sz_B(1:4)]); % m3 x q1 x q2 x s1 x s2
            
            sz_B = size(tmp);
            sz = size(a_lms{L2+1});
            tmp = reshape(tmp(:,:,:,:,1:sz(1)), [], sz(1))*a_lms{L2+1}; % m3 x q1 x q2 x s1 x m2
            tmp = reshape(tmp.', [sz(2), sz_B(1:4)]); % m2 x m3 x q1 x q2 x s1
            
            sz_B = size(tmp);
            sz = size(a_lms{L1+1});
            tmp = reshape(tmp(:,:,:,:,1:sz(1)), [], sz(1))*a_lms{L1+1}; % (m2, m3, q1, q2) x m1
            tmp = reshape(tmp.', [sz(2), sz_B(1:4)]); % m1 x m2 x m3 x q1 x q2
            
            % Perform summation over m2, m3:
            for m1 = -L1:L1
                for m2 = -L2:L2
                    m3 = m1+m2;
                    if abs(m3) > L3, continue; end
                    M{k+1} = M{k+1} + W{ii}(jj, m1+L1+1, m2+L2+1)*squeeze(tmp(m1+L1+1, m2+L2+1, m3+L3+1,:,:)); % q1 x q2
                end
            end
        end
    end
end

M = cellfun(@(x) real(x + x')/2, M, 'UniformOutput', 0); % enforce symmetry and realness 
