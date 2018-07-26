function [B_n, L1_list, L2_list, L3_list] = form_B_cell(B)

maxL = size(B,1)-1;
maxK = length(B{1}{1})-1;

B_n = cell(maxK+1, 1);
L1_list = []; L2_list = []; L3_list = [];
for ii = 1:(maxL+1)^2
    [L1, L2] = ind2sub([maxL+1, maxL+1], ii);
    L1 = L1-1; L2 = L2-1;
    
    L3_vals = abs(L1-L2):min(L1+L2, maxL);
    for jj = 1:length(L3_vals)
        L1_list = [L1_list;  L1];
        L2_list = [L2_list;  L2];
        L3_list = [L3_list;  jj];
    end
end

for k = 0:maxK
    B_n{k+1} = cell(length(L3_list), 1);
    for ii = 1:length(L3_list)
        B_n{k+1}{ii} = B{L1_list(ii)+1, L2_list(ii)+1}{L3_list(ii)}{k+1};
    end
end