function [Wt, ang_freq] = precomp_pswf_t_windows(L, beta, T)

[Mt, ang_freq] = precomp_pswf_t_mat(L, beta, T);
[x,y] = meshgrid(-L:L, -L:L); pts_notin_disc = sqrt(x.^2 + y.^2) > L;
Mt(:, pts_notin_disc) = 0;
Wt = reshape(Mt, size(Mt,1), 2*L+1, 2*L+1);
Wt = permute(Wt, [2,3,1]);
Wt = Wt(end:-1:1, end:-1:1, :);