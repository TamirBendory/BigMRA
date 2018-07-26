function B = form_B_cell_inverse(B_n, L1_list, L2_list, L3_list)

maxL = max(L1_list);
maxK = size(B_n,1)-1;

B = cell(maxL+1, maxL+1);

for ii = 1:length(L1_list)
    for k = 0:maxK
        B{L1_list(ii)+1, L2_list(ii)+1}{L3_list(ii)}{k+1} = B_n{k+1}{ii};
    end
end
