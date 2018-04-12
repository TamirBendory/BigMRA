% See notes April 12, 2018
clear;
% close all;
clc;

figure(1);
clf;

L = 11;
range = -(L-1) : (L-1);
[ell2, ell1] = ndgrid(range);

% Start by keeping only shifts such that all three signals overlap
M = (abs(ell2 - ell1) <= (L-1));

% Reject biased terms
M = M & (ell1 ~= 0);
M = M & (ell2 ~= 0);
M = M & (ell1 ~= ell2);

% Reject symmetry ell1 == ell2
M = M & (ell1 >= ell2);

% Identify further symmetries (seems there are none: see printed list at command prompt)
for a = 1:(2*L-1)
    for b = (2*L-1)
        l1 = ell1(a, b);
        l2 = ell2(a, b);
        ind = find( ell1 == -l1 & ell2 == l2 - l1 );
        if ~isempty(ind)
            [c, d] = ind2sub(size(M), ind);
            if M(a, b) && M(c, d)
                fprintf('(%d, %d) ~ (%d, %d)\n', l1, l2, -l1, l2-l1);
            end
        end
    end
end

imagesc(range, range, M);
xlabel('ell_1');
ylabel('ell_2');
axis equal; axis tight;
colormap gray;
title('Bright pixels are good shifts');

% Are these unique?

z = randn(L, 1);
list2 = (1:(L-1))';
list3 = zeros(nnz(M), 2);
k = 0;
for k1 = 1 : (2*L-1)
    for k2 = 1 : (2*L-1)
        if M(k1, k2)
            k = k + 1;
            list3(k, :) = [ell1(k1, k2), ell2(k1, k2)];
        end
    end
end
[M1, M2, M3] = moments_from_data_no_debias_1D(z, list2, list3);

numel(uniquetol(M2, 1e-8))
numel(M2)

numel(uniquetol(M3, 1e-8))
numel(M3)

% Write in the picture which moments are the same (they get the same number)
hold all;
[A, B, C] = uniquetol(M3, 1e-8);
for k = 1 : length(C)
    ll = list3(k, :);
    text(ll(1), ll(2), num2str(C(k)));
end

set(gca, 'YDir', 'normal')