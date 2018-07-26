function s_list = gen_s_list(maxL, r_cut, r_select_ratio, N)

fname='Besselj_L200_S500.mat';
load(fname)
 
B = B( B(:, 3)< 2*pi*N*r_cut * r_select_ratio & B(:,1) <= maxL, :); %Nyquist criterion
l_grid = B(:, 1);

s_list = zeros(max(l_grid)+1,1);
for l = 0:max(l_grid)
    s_list(l+1) = length(find(l_grid == l));
end