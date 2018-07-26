function [R_sum, P_sum, R, P, w_q] = gen_shift_grid(quad_pts)

% Generate grid for r integral:
[xx, w_x] = lgwt(quad_pts, -1, 1);
[X, Y] = meshgrid(xx, xx);
[R, P] = cart2pol(X, Y);
w_q = kron(w_x, w_x);

X_sum = bsxfun(@plus, X(:), 2*X(:).'); % (x) x (dx)
Y_sum = bsxfun(@plus, Y(:), 2*Y(:).'); % (y) x (dy)
[P_sum, R_sum] = cart2pol(X_sum, Y_sum);