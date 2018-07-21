function [m2, m3] = debias_moments_steerable(m1, m2, m3, sigma, L)

[R_0q, alpha_0q] = PSWF_radial_2D(0, size(m3{1},1)-1, pi*(L-1), 0);
R = (R_0q(:)/sqrt(2*pi))*(4*R_0q./(sqrt(2*pi)*alpha_0q.'));
R = R + R.';

m3 = cellfun(@(x) x - sigma^2*m1*eye(size(x)), m3, 'UniformOutput', 0);
m3{1} = m3{1} - sigma^2*m1*R;

m2 =  m2 - sigma^2*R_0q(:)/sqrt(2*pi);