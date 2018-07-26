function [M1, M2, M3] = moments_from_harmonics(a_lmi, L0)
% Evaluate the three first moments of the dataset on a L x L cartesian grid

% load('triple_factors.mat', 'C')
maxL = length(a_lmi) - 1;
[x,y] = meshgrid(-(L0-1):L0-1,-(L0-1):L0-1);
x_sum = bsxfun(@plus, x(:), x(:).'); y_sum = bsxfun(@plus, y(:), y(:).');
[phi_grid,k_grid] = cart2pol(x, y);
[phi_grid_sum,k_grid_sum] = cart2pol(x_sum, y_sum);

Y_l = YN2YL(getSH(maxL, [phi_grid(:), (pi/2)*ones(numel(phi_grid), 1)], 'real'));
Y_l = cellfun(@(x) x(1:2:end,:).', Y_l, 'UniformOutput', 0);
Y_l_sum = YN2YL(getSH(maxL, [phi_grid_sum(:), (pi/2)*ones(numel(phi_grid_sum), 1)], 'real'));
Y_l_sum = cellfun(@(x) x(1:2:end,:).', Y_l_sum, 'UniformOutput', 0);

j_l = generate_spherical_bessel_basis(0, size(a_lmi{1},1), 1/2, 0);
M1 = j_l{1}*a_lmi{1}/(sqrt(4*pi)*L0^2);

j_l = generate_spherical_bessel_basis(length(a_lmi)-1, cellfun(@(x)size(x,1),a_lmi), 1/2, k_grid(:)./(2*(L0-1)));
j_l_sum = generate_spherical_bessel_basis(length(a_lmi)-1, cellfun(@(x)size(x,1),a_lmi), 1/2, k_grid_sum(:)./(2*(L0-1)));
A_lm = cellfun(@mtimes, j_l, a_lmi, 'UniformOutput', 0);
A_lm_sum = cellfun(@mtimes, j_l_sum, a_lmi, 'UniformOutput', 0);

clear j_l j_l_sum phi_grid_sum k_grid_sum

M2 = cellfun(@(x) sum(abs(x).^2, 2), A_lm, 'UniformOutput', 0);
M2 = sum(cat(2, M2{:}), 2);
M2 = reshape(M2, size(k_grid))./(4*pi*L0^2);

M3 = zeros((2*L0-1)^2, (2*L0-1)^2);
for L1 = 0:maxL
    for L2 = 0:maxL
        for L3 = abs(L1-L2):min(L1+L2, maxL)
            list1 = -L1:2:L1; list2 = -L2:2:L2;
            parfor ii = 1:(2*L1+1)*(L1+1)*(2*L2+1)*(L2+1)
                [m1_ind, m2_ind, m1p_ind, m2p_ind] = ind2sub([2*L1+1, 2*L2+1, L1+1, L2+1], ii);
                m1 = m1_ind - L1-1; m2 = m2_ind - L2-1; 
                m1p = list1(m1p_ind); m2p = list2(m2p_ind); m3p_ind = find(-L3:2:L3 == m1p+m2p);
                if abs(m1+m2) <= L3 && abs(m1p + m2p) <= L3 && ~isempty(m3p_ind)
                    m3_ind = m1 + m2 + L3+1;
                    prod1 = (A_lm{L1+1}(:,m1_ind).*Y_l{L1+1}(:,m1p_ind)) * (A_lm{L2+1}(:,m2_ind).*Y_l{L2+1}(:,m2p_ind)).';
                    prod2 = conj(A_lm_sum{L3+1}(:,m3_ind).*Y_l_sum{L3+1}(:,m3p_ind));
                    prod3 = prod1.*reshape(prod2, size(prod1));
                    const = (-1)^(m1+m2+m1p+m2p)*double(w3j(L1, L2, L3, m1, m2, -m1-m2)*w3j(L1, L2, L3, m1p, m2p, -m1p-m2p));
                    M3 = M3 + prod3*const;
                end
            end
        end
    end
end
M3 = reshape(M3, [2*L0-1,2*L0-1,2*L0-1,2*L0-1])./(L0^2);